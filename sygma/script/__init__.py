"""SyGMa: Systematically Generating potential Metabolites"""

import argparse
import sygma
import sys, logging, shutil
from rdkit import Chem
from rdkit.Chem import AllChem


def run_sygma(args):
    scenario = sygma.Scenario([
        [sygma.ruleset['phase1'], args.phase1],
        [sygma.ruleset['phase2'], args.phase2]
    ])

    parent = Chem.MolFromSmiles(args.parentmol)
    Chem.AllChem.Compute2DCoords(parent)
    metabolites_network = scenario.run(parent)
    metabolites_network.calc_scores()
    metabolites_network.add_coordinates()
    if args.outputtype == "sdf":
        metabolites_network.write_sdf()
    elif args.outputtype == "smiles":
        print metabolites_network.to_smiles()[:-1]

    return None

def main():
    "Entry point for magma script"
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--version', action='version', version='%(prog)s ' + sygma.version())
    ap.add_argument('-o', '--outputtype', help="Molecule output type (default: %(default)s)", default="sdf", type=str)
    ap.add_argument('-1', '--phase1', help="Number of phase 1 cycles", type=int)
    ap.add_argument('-2', '--phase2', help="Number of phase 2 cycles", type = int)
    ap.add_argument('parentmol', help="Smiles string of parent molecule structure", type=str)
    """Parse arguments and run subcommand"""
    args = ap.parse_args(sys.argv[1:])
    return run_sygma(args)

if __name__ == "__main__":
    main()
