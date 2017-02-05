SyGMa
=====
SyGMa is a python library for the **Sy**\ stematic **G**\ eneration of potential **M**\ et\ **a**\ bolites.
It is a reimplementation of the metabolic rules outlined in
`Ridder, L., & Wagener, M. (2008)
SyGMa: combining expert knowledge and empirical scoring in the prediction of metabolites.
ChemMedChem, 3(5), 821-832
<http://onlinelibrary.wiley.com/doi/10.1002/cmdc.200700312/full>`_.

.. image:: https://travis-ci.org/3D-e-Chem/sygma.svg?branch=master
    :target: https://travis-ci.org/3D-e-Chem/sygma
.. image:: https://api.codacy.com/project/badge/Grade/7f18ab1d1a80437f8e28ac1676c70519
    :target: https://www.codacy.com/app/3D-e-Chem/sygma?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=3D-e-Chem/sygma&amp;utm_campaign=Badge_Grade
.. image:: https://api.codacy.com/project/badge/Coverage/7f18ab1d1a80437f8e28ac1676c70519
    :target: https://www.codacy.com/app/3D-e-Chem/sygma?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=3D-e-Chem/sygma&amp;utm_campaign=Badge_Coverage
.. image:: https://img.shields.io/badge/docker-ready-blue.svg
    :target: https://hub.docker.com/r/3dechem/sygma
.. image:: https://anaconda.org/3d-e-chem/sygma/badges/installer/conda.svg
    :target: https://conda.anaconda.org/3d-e-chem

Requirements
------------
SyGMa requires RDKit with INCHI support

Installation
------------
* Install with Anaconda: ``conda install -c 3d-e-Chem -c rdkit sygma``

OR

* Install RDKit following the instructions in http://www.rdkit.org/docs/Install.html

AND

* ``pip install sygma`` OR, after downloading sygma, ``python setup.py install``

Example: generating metabolites of phenol
-----------------------------------------
.. code-block:: python

    import sygma
    from rdkit import Chem

    # Each step in a scenario lists the ruleset and the number of reaction cycles to be applied
    scenario = sygma.Scenario([
        [sygma.ruleset['phase1'], 1],
        [sygma.ruleset['phase2'], 1]])

    # An rdkit molecule, optionally with 2D coordinates, is required as parent molecule
    parent = Chem.MolFromSmiles("c1ccccc1O")

    metabolic_tree = scenario.run(parent)
    metabolic_tree.calc_scores()

    print metabolic_tree.to_smiles()


Docker
------
SyGMa can be executed in a Docker (https://www.docker.com/) container as follows:

.. code-block:: bash

    docker run 3dechem/sygma c1ccccc1O
