import sygma
import argparse
from StringIO import StringIO


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