.. sygma documentation master file, created by
   sphinx-quickstart on Tue Aug  9 13:41:23 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SyGMa's documentation!
=================================

Contents:

.. toctree::
   :maxdepth: 2

   Introduction <self>
   api
   command

.. include:: introduction.rst

Example
-------
.. include:: ../sygma/test/test_sygma.py
    :code: python

Docker
------
SyGMa can be executed in a `Docker <https://www.docker.com/>`_ container as follows:

.. code-block:: bash

    docker run 3dechem/sygma c1ccccc1O

Rulesets
--------
.. automodule:: ruleset

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

