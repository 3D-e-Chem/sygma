import sygma
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem


def test_treenode_init():
    """Test init of TreeNode class"""

    # generate phenol molecule
    mol = Chem.MolFromSmiles('c1ccccc1O')
    node = sygma.TreeNode(mol)

    # test if inchikey is generated correctly
    assert node.ikey == "ISWSIDIOOBJBQZ"

def test_treenode_gen_coords():
    """Test gen_coords method of TreeNode class"""

    # generate phenol molecule with coordinates
    mol = Chem.MolFromSmiles('c1ccccc1O')
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer(0)
    pos2x = conf.GetAtomPosition(2).x
    # remove coordinates of oxygen atom
    conf.SetAtomPosition(6, Geometry.Point3D(0.0, 0.0, 0.0))
    # create TreeNode
    node = sygma.TreeNode(mol)
    # generate coordinates for atoms without atom positions (i.e. all coordinates are zero)
    node.gen_coords()

    # test if existing coordinates are unchanged in gen_coords method
    assert conf.GetAtomPosition(2).x == pos2x

    # test if additional coordinates are generated in gen_coords method
    conf = node.mol.GetConformer(0)
    assert conf.GetAtomPosition(6).x != 0.0
