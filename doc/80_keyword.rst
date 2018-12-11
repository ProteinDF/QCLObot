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

Declare the variable. 
The declared variables can be used with "with_items" and so on.


.. code-block:: yaml
  
  vars:
    residues: [27, 28, 48, 97, 101, 108, 248, 547, 558, 559]
  
  tasks:
    - name: res_{{ item }}
      sp: true
      fragments:
        - name: res_{{ item }}
          brd_select: /model_1/A/{{ item }}/
        - name: ACE_{{ item }}
          add_ACE: true
          brd_select: /model_1/A/{{ item -1 }}/
        - name: NME_{{ item }}
          add_NME: true
          brd_select: /model_1/A/{{ item +1 }}/
      with_items: residues


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

mail
^^^^

Send mail.

.. code-block:: yaml

  - name: send_mail  
    mail:
      smtp_server: xxx.xxx.xxx.xxx
      smtp_port: xx
      use_SSL: yes
      smtp_account: xxxx
      smtp_password: xxxx
      from_address: xxxx@xxxx.xxxx
      to_address: xxxx@xxxx.xxxx
      subject: 'test mail'
      msg: 'job finished.'


debug
^^^^^

do nothing.


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
The fragment is defined by following keywords.
All fragment requires the "name" attribute, which is used as name of the fragment.


* atomlist

'atomlist' directive makes fragment from atom list

.. code-block:: yaml
                
  - name: N2
    sp: true
    fragments:
      - name: N2
        atomlist:
          - "N  0.000000   0.000000   0.000000"
          - "N  1.000000   0.000000   0.000000"

        
The atomlist is an array object.
Each atom is defined by string separated white space, 
or by array object as following:


.. code-block:: yaml

  - name: OH
    sp: true
    fragments:
      - name: OH
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

.. code-block:: yaml

  - name: small_mol
    sp: true
    fragments:
      - name: frag1
        add_CH3: true
        displacement: "/model_1/A/100/100_C1"
        root: "/model_1/A/100/100_C2"


* add_ACE

Place the acetyl group in the specified place and add it as a fragment.

.. code-block:: yaml

    - name: res_3
      sp: true
      fragments:
        - name: res_3
          brd_select: /model_1/A/3/
        - name: ACE_3
          add_ACE: true
          brd_select: /model_1/A/2/


* add_NME

Place the N-methyl group in the specified place and add it as a fragment.

.. code-block:: yaml

    - name: res_3
      sp: true
      fragments:
        - name: res_3
          brd_select: /model_1/A/3/
        - name: NME_3
          add_NME: true
          brd_select: /model_1/A/4/


* reference

The fragment is created by using the previous calculation result.

  * frame

    This value indicates the name of the frame molecule.

  * fragment

    The name of the fragment in the frame molecule.

.. code-block:: yaml

  - name: res_3
    sp: true
    fragments:
      - name: res_3
        brd_select: /model_1/A/3/
      - name: ACE_3
        add_ACE: true
        brd_select: /model_1/A/3/
      - name: NME_3
        add_NME: true
        brd_select: /model_1/A/3/
  
  - name: res_3-7
    sp: true
    guess: QCLO
    fragments:
      - name: referenced_res_3
        reference:
          frame: res_3
          fragment: res_3


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


.. code-block:: yaml

  - name: res_3-7
    pre_scf: yes
    guess: QCLO
    sp: yes
    gradient: yes


default frame
^^^^^^^^^^^^^

The "default frame" is a special frame.

If the name section is 'default', this frame parameters are used as default values.

In the following example, the frame is calculated as DZVP2 as the basisset and the exchange correlation functional is used by B3LYP.


.. code-block:: yaml

  tasks:
    - name: default
      brd_file: sample.brd
      basis_set: DZVP2
      XC_functional: B3LYP

    - name: N2
      atomlist:
        - "N  0.000000   0.000000   0.000000"
        - "N  1.000000   0.000000   0.000000"

