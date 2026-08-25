Setting up modeling proteins
============================

Prepare Applications
--------------------

Before running QCLObot-modeler, AmberTools is required.
Please set the environment variable, `AMBERHOME`. 


Prepare Oxytocin model
----------------------

You can get Oxytocin model (2MGO.pdb) from PDB or :download:`here <_static/2MGO.pdb>`.

The file (2,875 lines) begins as follows; download the full file via the link above.

.. literalinclude:: _static/2MGO.pdb
   :language: none
   :encoding: utf-8
   :lines: 1-20

Create input file for QCLObot-modeler
-------------------------------------

You save the following text as :download:`2MGO_modeling.yaml <_static/2MGO_modeling.yaml>`.

.. literalinclude:: _static/2MGO_modeling.yaml
   :language: yaml
   :encoding: utf-8

Run QCLObot-modeler
-------------------

.. code-block:: bash

   $ ${PDF_HOME}/bin/QCLObot-modeler.py 2MGO_modeling.yaml




