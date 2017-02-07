from rdkit import Chem
from rdkit.Chem import AllChem
import copy
import sys
from sygma.treenode import TreeNode
import logging
logger = logging.getLogger('sygma')


class Tree(object):
    """
    Class to build and analyse a metabolic tree

    :param parentmol:
        An RDKit molecule
    """
    def __init__(self, parentmol=None):
        self.nodes = {}
        if parentmol:
            parentnode = TreeNode(parentmol, parent=None, rule=None, score=1, pathway="")
            self.nodes[parentnode.ikey] = parentnode
            self.parentkey = parentnode.ikey

    def _react(self, reactant, reaction):
        """Apply reaction to reactant and return products"""
        ps = reaction.RunReactants([reactant])
        products = []
        for product in ps:
            frags = (Chem.GetMolFrags(product[0], asMols=True, sanitizeFrags=False))
            for p in frags:
                q = copy.copy(p)
                try:
                    Chem.SanitizeMol(q)
                    products.append(q)
                except:
                    pass # Ignore fragments that cannot be sanitized
        return products

    def metabolize_node(self, node, rules):
        for rule in rules:
            products = self._react(node.mol, rule.reaction)
            for x in products:
                ikey = AllChem.InchiToInchiKey(AllChem.MolToInchi(x))[:14]
                x.SetProp("_Name", ikey)
                node.children.append(ikey)
                if ikey in self.nodes:
                    if node.ikey not in self.nodes[ikey].parents or \
                                    self.nodes[ikey].parents[node.ikey].probability < rule.probability:
                        self.nodes[ikey].parents[node.ikey] = rule
                else:
                    self.nodes[ikey] = TreeNode(x, parent=node.ikey, rule=rule)

    def metabolize_all_nodes(self, rules, cycles=1):
        """
        Metabolize all nodes according to [rules], for [cycles] number of cycles

        :param rules:
            List of rules
        :param cycles:
            Integer indicating the number of subsequent steps to apply the rules
        """
        for i in range(cycles):
            logger.info('Cycle ' + str(i + 1))
            ikeys = list(self.nodes.keys())
            for ikey in ikeys:
                self.metabolize_node(self.nodes[ikey], rules)

    def add_coordinates(self):
        """
        Add missing atomic coordinates to all metabolites
        """
        for node in self.nodes.values():
            node.gen_coords()

    def calc_scores(self):
        """
        Calculate probability scores for all metabolites
        """
        for key in self.nodes:
            self.calc_score(self.nodes[key])

    def calc_score(self, node):
        if node.score != None:
            return
        else:
            node.score = -1 # Indicating the score is "waiting" for a result
            for pkey in node.parents:
                if self.nodes[pkey].score is None:      # No calculation is requested for parents with score -1,
                    self.calc_score(self.nodes[pkey])   # to avoid closed loops ....
                newscore = self.nodes[pkey].score * float(node.parents[pkey].probability)
                if node.score < newscore:   # This implies that parents with a newscore of -1 are ignored,
                    node.score = newscore   # to avoid closed loops ...
                    node.pathway = self.nodes[pkey].pathway + node.parents[pkey].rulename + "; \n"


    def to_list(self, filter_small_fragments = True, parent_column = 'parent'):
        """
        Generate a list of metabolites

        :param filter_small_fragments:
            Boolean to activate filtering all metabolites with less then 15% of original atoms (of the parent)
        :param parent_column:
            String containing the name for the column with the parent molecule
        :return:
            A list of dictionaries for each metabolites, containing the SyGMa_metabolite (an RDKit Molecule),
            SyGMa_pathway and SyGMa_score, sorted by decreasing probability.
        """
        def sortkey(rowdict):
            return rowdict['SyGMa_score']
        output_list = []
        n_parent_atoms = self.nodes[self.parentkey].mol.GetNumAtoms()
        for key in self.nodes:
            if filter_small_fragments and float(self.nodes[key].n_original_atoms) <= 0.15 * n_parent_atoms:
                continue
            pathway = "parent;" if key == self.parentkey else self.nodes[key].pathway
            output_list.append({parent_column: self.nodes[self.parentkey].mol,
                         "SyGMa_pathway": pathway,
                         "SyGMa_metabolite": self.nodes[key].mol,
                         "SyGMa_score": self.nodes[key].score})
        output_list.sort(key=sortkey, reverse=True)
        return output_list

    def to_smiles(self, filter_small_fragments = True):
        """
        Generate a smiles list of metabolites

        :param filter_small_fragments:
            Boolean to activate filtering all metabolites with less then 15% of original atoms (of the parent)
        :return:
            A list of metabolites as list ``[[SyGMa_metabolite as smiles, SyGMa_score]]``
            sorted by decreasing probability score.
        """
        output_list = self.to_list(filter_small_fragments=filter_small_fragments)
        smiles_list = []
        for entry in output_list:
            smiles_list.append([Chem.MolToSmiles(entry['SyGMa_metabolite']),entry['SyGMa_score']])
        return smiles_list

    def write_sdf(self, file=sys.stdout, filter_small_fragments = True):
        """
        Generate an SDFile with metabolites including the SyGMa_pathway and the SyGMa score as properties

        :param file:
            The SDF file to write to
        :param filter_small_fragments:
            Boolean to activate filtering all metabolites with less then 15% of original atoms (of the parent)
        """
        output_list = self.to_list(filter_small_fragments=filter_small_fragments)
        sdf = Chem.SDWriter(file)
        for entry in output_list:
            mol = entry['SyGMa_metabolite']
            mol.SetProp("Pathway", entry['SyGMa_pathway'][:-1])
            mol.SetProp("Score", str(entry['SyGMa_score']))
            sdf.write(mol)
