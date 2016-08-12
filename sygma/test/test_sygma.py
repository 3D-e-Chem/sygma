import sygma
from rdkit import Chem
import argparse
from StringIO import StringIO

def test_predict_phenol_metabolites():
    """Test prediction of phenol metabolites by sygma module"""

    # Each step in a scenario lists the ruleset and the number of reaction cycles to be applied
    scenario = sygma.Scenario([
        [sygma.ruleset['phase1'], 1],
        [sygma.ruleset['phase2'], 1]])

    # An rdkit molecule, optionally with 2D coordinates, is required as parent molecule
    parent = Chem.MolFromSmiles("c1ccccc1O")

    metabolic_tree = scenario.run(parent)
    metabolic_tree.calc_scores()

    metabolite_list = metabolic_tree.to_list()
    assert len(metabolite_list) == 12
    assert metabolite_list[0]['SyGMa_score'] == 1
    assert metabolite_list[1]['SyGMa_pathway'] == 'O-glucuronidation_(aromatic_hydroxyl); \n'

def test_sygma_command_line1():
    """Test prediction of phenol metabolites by sygma command, output smiles"""

    args = argparse.Namespace()
    args.phase1 = 1
    args.phase2 = 2
    args.parentmol = 'c1ccccc1O'
    args.outputtype = 'smiles'
    args.loglevel = 'INFO'
    out = StringIO()
    sygma.script.run_sygma(args, out)
    assert len(out.getvalue().split('\n')) == 28

def test_sygma_command_line2():
    """Test prediction of phenol metabolites by sygma command, output sdf"""

    args = argparse.Namespace()
    args.phase1 = 1
    args.phase2 = 0
    args.parentmol = 'c1ccccc1O'
    args.outputtype = 'sdf'
    args.loglevel = 'INFO'
    out = StringIO()
    sygma.script.run_sygma(args, out)
    assert len(out.getvalue().split('$$$$\n')) == 4