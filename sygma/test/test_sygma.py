import sygma
from rdkit import Chem
from rdkit.Chem import AllChem


def test_sygma():
    scenario = sygma.Scenario([
        [sygma.ruleset['phase1'], 1],
        [sygma.ruleset['phase2'], 1]
    ])

    parent = Chem.MolFromSmiles("c1ccccc1O")
    Chem.AllChem.Compute2DCoords(parent)
    metabolites_network = scenario.run(parent)
    metabolites_network.calc_scores()
    metabolites_network.add_coordinates()

    metabolite_list = metabolites_network.to_list()
    assert len(metabolite_list) == 12
