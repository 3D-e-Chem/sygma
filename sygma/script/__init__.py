"""SyGMa: Systematically Generating potential Metabolites"""

from builtins import str
import argparse
import sygma
import sys
from rdkit import Chem, RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)
import logging
logging.basicConfig()
logger = logging.getLogger('sygma')

def run_sygma(args, file=sys.stdout):
    logger.setLevel(args.loglevel.upper())
    scenario = sygma.Scenario([
        [sygma.ruleset['phase1'], args.phase1],
        [sygma.ruleset['phase2'], args.phase2]
    ])

    parent = Chem.MolFromSmiles(args.parentmol)
    metabolic_tree = scenario.run(parent)
    metabolic_tree.calc_scores()
    if args.outputtype == "sdf":
        metabolic_tree.write_sdf(file)
    elif args.outputtype == "smiles":
        file.write("\n".join([m+" "+str(s) for m,s in metabolic_tree.to_smiles()])+'\n')
    return None

def get_sygma_parser():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--version', action='version', version='%(prog)s ' + sygma.__version__)
    ap.add_argument('-o', '--outputtype', help="Molecule output type (default: %(default)s)", default="sdf", type=str)
    ap.add_argument('-1', '--phase1', help="Number of phase 1 cycles (default: %(default)s)", default=1, type=int)
    ap.add_argument('-2', '--phase2', help="Number of phase 2 cycles (default: %(default)s)", default=1, type=int)
    ap.add_argument('-l', '--loglevel', help="Set logging level (default: %(default)s)", default='info',
                    choices=['debug', 'info', 'warn',' error'])
    ap.add_argument('parentmol', help="Smiles string of parent molecule structure", type=str)
    return ap

def main():
    """Entry point for magma script"""

    # Parse arguments and run subcommand
    ap = get_sygma_parser()
    args = ap.parse_args(sys.argv[1:])
    return run_sygma(args)

if __name__ == "__main__":
    main()
