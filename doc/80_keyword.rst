********************************
explanation of QCLObot YAML File
********************************

QCLObot playbook is separated by following two parts:
- vars section
  definition of variables
  
- tasks section
  definition of task for general and for frame molecule

vars section
============


  
tasks section
=============

The tasks section is consist of three parts:
- condition
- task
- frame definition


condition
---------

with_items
^^^^^^^^^^

repeat the task.

  
when
^^^^

do the task if the condition is satisfied.
  
  
task
----

debug
^^^^^


mail
^^^^


frame
-----

frame definition
^^^^^^^^^^^^^^^^

name(mandatory)
"""""""""""""""

All frame require the name.
Based of this name value, the working directory is created on the current directory. 


fragments
"""""""""

A frame molecule consists of fragment(s).
Fragment is defined by following keywords.
All fragment requires the "name" attribute, which is used as name of the fragment.


* atomlist

'atomlist' directive makes fragment from atom list

.. code-block:: yaml
                
   atomlist:
     - "N  0.000000   0.000000   0.000000"
     - "N  1.000000   0.000000   0.000000"

        
The atomlist is an array object.
Each atom is defined by string separated white space, 
or by array object as following:


.. code-block:: yaml

   atomlist:
     - [O, -7.328, -30.909,  17.923]
     - [H, -6.026, -31.757,  17.909]                

       
* add_CH3

If the keyword is defined as "yes",
a methyl group is add as fragment.

  * displacement

    This atom is substitute with methyl carbon.
    This value is specified by string as Bridge path.

  * root

    This atom is indicated to the next atom of the displacement atom.
    This value is specified by string as Bridge path.

  
* add_ACE

  
* add_NME


* reference

The fragment is created by using the previous calculation result.

  * frame

    This value indicates the name of the frame molecule.

  * fragment

    The name of the fragment in the frame molecule.


* brd_select

The group which is indicated by the value of "brd_select" keyword is add as fragment.
This value is specified by string as Bridge path.


frame action
^^^^^^^^^^^^

The following keyword indicates for the frame object to do. 

* pre_scf

If "pre_scf" is defined as "yes",
then the processing calculation before SCF loop is carried out in the frame molecule.


* guess

Creation of the initial guess is executed.
How to create guess depends on the value of "guess" keyword.

  * harris

    The initial guess is created by using Harris functional method.
    This is default.

  * QCLO

    The inigial guess is made of the QCLOs of corresponding fragments by using QCLO method.
    If the QCLO of the child fragment has not been created,
    it is computed automatically.

    
* sp

If the "sp" is defined as "yes",
the single-point calculation of the frame molecule is carried out.
If "pre_scf" and "guess" keywords are not indicated,
these operations are automatically performed.


* gradient

If the "gradient" is defined as "yes",
the energy gradient is gained in the frame molecule.


calculation configuration
^^^^^^^^^^^^^^^^^^^


default frame
^^^^^^^^^^^^^

if the name section is 'default', this frame parameters are used as default values.
