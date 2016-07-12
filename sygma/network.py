from rdkit import Chem
from rdkit.Chem import AllChem
import copy
from sygma.treenode import TreeNode


class Network:
    def __init__(self, parentmol=None):
        self.nodes = {}
        if parentmol:
            parentnode = TreeNode(parentmol, parent=None, rule=None, score=1, path="")
            self.nodes[parentnode.ikey] = parentnode
            self.parentkey = parentnode.ikey

    def react(self, reactant, reaction):
        """ Apply reaction to reactant and return products """
        ps = reaction.RunReactants([reactant])
        products = []
        for product in ps:
            frags = (Chem.GetMolFrags(product[0], asMols=True))
            for p in frags:
                q = copy.copy(p)
                Chem.SanitizeMol(q)
                products.append(q)
        return products

    def metabolize_node(self, node, rules):
        for rule in rules:
            products = self.react(node.mol, rule.reaction)
            for x in products:
                ikey = AllChem.InchiToInchiKey(AllChem.MolToInchi(x))[:14]
                node.children.append(ikey)
                if ikey in self.nodes:
                    if node.ikey not in self.nodes[ikey].parents or \
                                    self.nodes[ikey].parents[node.ikey].probability < rule.probability:
                        self.nodes[ikey].parents[node.ikey] = rule
                else:
                    self.nodes[ikey] = TreeNode(x, parent=node.ikey, rule=rule)

    def metabolize_all_nodes(self, rules, cycles=1):
        for i in range(cycles):
            ikeys = self.nodes.keys()
            for ikey in ikeys:
                self.metabolize_node(self.nodes[ikey], rules)

    def add_coordinates(self):
        for node in self.nodes.itervalues():
            node.gen_coords()

    def calc_scores(self):
        for key in self.nodes:
            self.calc_score(self.nodes[key])

    def calc_score(self, node):
        if node.score > 0:
            return
        else:
            for pkey in node.parents:
                if self.nodes[pkey].score is None:
                    self.calc_score(self.nodes[pkey])
                newscore = self.nodes[pkey].score * float(node.parents[pkey].probability)
                if node.score < newscore:
                    node.score = newscore
                    node.path = self.nodes[pkey].path + node.parents[pkey].rulename + "; \n"

    def to_list(self, parent_key_name=None):
        def sortkey(rowdict):
            return rowdict['SyGMa_score']

        list = []
        for key in self.nodes:
            path = "parent" if key == self.parentkey else self.nodes[key].path
            list.append({"parent": self.nodes[self.parentkey].mol,
                         "SyGMa_reaction": path,
                         "SyGMa_metabolite": self.nodes[key].mol,
                         "SyGMa_score": self.nodes[key].score})
        list.sort(key=sortkey, reverse=True)
        return list
