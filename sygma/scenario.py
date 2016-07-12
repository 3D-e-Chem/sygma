from rdkit.Chem import AllChem
from sygma.network import Network


class Scenario:
    """class to read and process metabolic scenario"""

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
        network = Network(parentmol)
        for name, cycles, rules in self.scenario:
            network.metabolize_all_nodes(rules, cycles)
        return network


class Rule:
    """class describing a metabolic rule"""

    def __init__(self, rulename, probability, smarts):
        self.rulename = rulename
        self.probability = probability
        self.reaction = AllChem.ReactionFromSmarts(smarts)
