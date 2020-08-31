****************************************
explanation of QCLObot-modeler YAML File
****************************************

QCLObot-modeler playbook is separated by following two parts:

- vars section

  definition of variables
  
- tasks section

  definition of task for general and for frame molecule

vars section
============

Declare the variable. 
The declared variables can be used with "with_items" and so on.


tasks section
=============

The tasks section is consist of three parts:

- protonate


protonate
^^^^^^^^^

Protonate the molecule.

- src

Specify the PDB file.


.. code-block:: yaml

  - name: 2MGO
    protonate:
      src: 2mgo.pdb


neutralize
^^^^^^^^^^

Add ions to neutralize.
Ignore ions on the ion-pairs.

.. code-block:: yaml

  - name: 2MGO_optx
    neutralize:
      reference: 2MGO_opt
      dest: 2MGO_optx.pdb


opt
^^^^^^^^^

Do energy minimization simulation.

- reference

Specify the task to be the input structure.

- belly_mask:

(option)

Use the belly method.
The force of a specified molecule is considered to be zero.


.. code-block:: yaml

  - name: 2MGO_optx_mdwatx_optwatx
    opt:
      reference: 2MGO_optx_mdwatx
      belly_mask:
        - water
        - ions
      dest: 2MGO_optx_mdwatx.pdb


md
^^^^^

Run molecular dynamics simulation.

- name

(mandatory)
Give this task a name. 
A directory of this name is created. The task is executed in this directory.

- reference

Specify the task to be the input structure.

- solvation:

(option)
Specify the solvation model.

- belly_mask:

(option)

Use the belly method.
The force of a specified molecule is considered to be zero.


.. code-block:: yaml
                
  - name: 2MGO_optx_mdwatx
    md:
      reference: 2MGO_optx
      solvation:
        method: cap
        model: TIP3PBOX
      belly_mask:
        - water
        - ions