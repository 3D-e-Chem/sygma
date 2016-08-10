from rdkit.Chem import AllChem
from sygma.network import Network


class Scenario:
    """
    Class to read and process metabolic scenario

    :param scenario:
        A list of lists, each representing a metabolic phase as
        [name_of_file_containing_rules, number_of_cycles_to_apply]
    """

    def __init__(self, scenario):
        self.rules = {}
        for step in scenario:
            name, cycles = step
            step.append(self.read_reaction_rules(name))
        self.scenario = scenario

    def read_reaction_rules(self, filename):
        rules = []
        for l in open(filename, "r"):
            if l != "\n" and l[0] != "#":
                smarts, probability, name, rest = l.split("\t")
                rules.append(Rule(name, probability, smarts))
        return rules

    def run(self, parentmol):
        """
        :param parentmol:
            An RDKit molecule
        :return:
            A sygma.Network object
        """
        if parentmol.GetNumConformers() == 0:
            # make sure the parentmolecule has coordinates
            AllChem.Compute2DCoords(parentmol)
        network = Network(parentmol)
        for name, cycles, rules in self.scenario:
            network.metabolize_all_nodes(rules, cycles)
        return network


class Rule:
    """
    Class to contain a metabolic rule

    :param rulename:
        A string containing a unique name of the rule
    :param probability:
        A probability value between 0 and 1 indicating the empirical success rate of the rule
    :param smarts:
        A reaction smarts describing the chemical transformation of the rule
    """

    def __init__(self, rulename, probability, smarts):
        self.rulename = rulename
        self.probability = probability
        self.reaction = AllChem.ReactionFromSmarts(smarts)
