import sygma
from rdkit import Chem

def test_prediction_phenol():
    """Integration test of sygma module"""
    scenario = sygma.Scenario([
        [sygma.ruleset['phase1'], 1],
        [sygma.ruleset['phase2'], 1]
    ])

    parent = Chem.MolFromSmiles("c1ccccc1O")

    metabolites_network = scenario.run(parent)
    metabolites_network.calc_scores()
    metabolites_network.add_coordinates()

    metabolite_list = metabolites_network.to_list()
    assert len(metabolite_list) == 12
    assert metabolite_list[0]['SyGMa_score'] == 1
    assert metabolite_list[1]['SyGMa_path'] == 'O-glucuronidation_(aromatic_hydroxyl); \n'
