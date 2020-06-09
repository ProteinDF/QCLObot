Setting up calculation model
============================

Prepare Oxytocin model
----------------------

You can get Oxytocin model (1NPO_model.pdb) from here.
Or you save the following text as 1NPO_model.pdb.

.. literalinclude:: 1NPO_model.pdb
   :language: none
   :encoding: utf-8



Check the model structure
----------------------------

You should check the model structure by using viewing software such as PyMol and Chimera.
The current model has a dislufide bond between CYC1 and CYC6.


              
Convert PDB file to Bridge format file
--------------------------------------

.. code-block:: bash

  $ ${PDF_HOME}/bin/pdb2brd.py 1NPO_model.pdb 1NPO_model.brd

   
.. note::

  You can check the bridge file contents using the following utility.

.. code-block:: bash

  $ ${PDF_HOME}/bin/brd2txt.py 1NPO_model.brd
   

  
Create input file for QCLObot
-----------------------------

You save the following text as 1NPO_model.QCLO.yaml.

.. literalinclude:: 1NPO_model.QCLO.yaml
   :language: yaml
   :encoding: utf-8


              
Run QCLObot
-----------

.. code-block:: bash

   $ ${PDF_HOME}/bin/QCLObot.py 1NPO_mode.QCLO.yaml



