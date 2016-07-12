# SyGMa
SyGMa is a python library for the  **Sy**stematic **G**eneration of potential **M**et**a**bolites.

## Requirements
SyGMa requires RDKit with INCHI support

## Installation
* See http://www.rdkit.org/docs/Install.html for RDKit installation instructions.
* python setup.py install

## Example: generating metabolites of phenol
```
import sygma
from rdkit import Chem

# each step in a scenario lists the ruleset and the number of reaction cycles to be applied 
scenario = sygma.Scenario([
    [sygma.ruleset['phase1'], 2],
    [sygma.ruleset['phase2'], 1]
    ])

# An rdkit molecule, optionally with 2D coordinates, is required as parent molecule
parent = Chem.MolFromSmiles("c1ccccc1O")
Chem.AllChem.Compute2DCoords(parent)

metabolites_network = scenario.run(parent)
metabolites_network.calc_scores()
metabolites_network.add_coordinates() # Only if a parent molecule was provided with coordinates

print metabolites_network.to_list()
```



