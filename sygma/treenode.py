from rdkit import Geometry
from rdkit.Chem import AllChem, rdMolTransforms


class TreeNode:
    def __init__(self, mol, parent="", rule=None, score=None, path=""):
        self.mol = mol
        self.parents = {parent: rule}
        self.children = []
        self.ikey = AllChem.InchiToInchiKey(AllChem.MolToInchi(mol))[:14]
        self.score = score
        self.path = path
        self.n_original_atoms = None

    def gen_coords(self):
        """ Calculate 2D positions for atoms without coordinates """
        conf = self.mol.GetConformer(0)
        coord_dict = {}
        # Put known coordinates in coordDict
        for i in range(self.mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            if pos.x != 0.0 or pos.y != 0.0:
                coord_dict[i] = Geometry.Point2D(pos.x, pos.y)
        self.n_original_atoms = len(coord_dict)
        if self.n_original_atoms > 1:
            # calculate average length of all bonds with coordinates
            total = 0
            n = 0
            for bond in self.mol.GetBonds():
                b = bond.GetBeginAtomIdx()
                e = bond.GetEndAtomIdx()
                if b in coord_dict and e in coord_dict:
                    n += 1
                    total += rdMolTransforms.GetBondLength(conf, b, e)
            av = total / n
            # compute coordinates for new atoms, keeping known coordinates
            AllChem.Compute2DCoords(self.mol, coordMap=coord_dict, bondLength=av)
        else:
            AllChem.Compute2DCoords(self.mol)
