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
    - name: frame_{{ item }}
      sp: true
      fragments:
        - name: AA_{{ item }}
          brd_select: /model_1/A/{{ item }}/
        - name: ACE_{{ item }}
          add_ACE: true
          brd_select: /model_1/A/{{ item -1 }}/
        - name: NME_{{ item }}
          add_NME: true
          brd_select: /model_1/A/{{ item +1 }}/
      with_items: residues


.. warning::

  Note that when using templates, especially when operating within a template, there is a distinction between strings and numbers.


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
  

include
^^^^^^^

Insert an external file. Fill in the external file with the "tasks" statement.

.. code-block:: yaml
  
  tasks:
    - name: include_step1
      include_tasks: "step1.yaml"


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


cmd_alias
"""""""""

Accept values of the dictionary type. Replaces the default external commands. The following commands are currently supported:

.. It can be found in _get_default_cmds() of qcframe.py.

- archive
- mat-extend
- mat-mul
- mat-select
- mat-symetrize
- mat-transpose
- mat-diagonal


.. db_filename
.. """""""""""

.. Specify the name of the DB file to be used.


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

  - name: frame_3
    sp: true
    fragments:
      - name: AA_3
        brd_select: /model_1/A/3/
      - name: ACE_3
        add_ACE: true
        brd_select: /model_1/A/3/
      - name: NME_3
        add_NME: true
        brd_select: /model_1/A/3/
  
  - name: frame_3-7
    sp: true
    guess: QCLO
    fragments:
      - name: AA_3
        reference:
          frame: frame_3
          fragment: AA_3


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
  Use mapping format in "guess" section.

  * method

    Specified method for guess. Possible values are as follows.

    * harris

      The initial guess is created by using Harris functional method.
      This is default.

    * QCLO

      The inigial guess is made of the QCLOs of corresponding fragments by using QCLO method.
      If the QCLO of the child fragment has not been created,
      it is computed automatically.

  * force

    boolean. The default value is False. If you want to force execution even if it has already been done, specify True.
    
* sp

  If the "sp" is defined as "yes",
  the single-point calculation of the frame molecule is carried out.
  If "pre_scf" and "guess" keywords are not indicated,
  these operations are automatically performed.


* force

  If the "force" is defined as "yes",
  the energy force is gained in the frame molecule.


  .. code-block:: yaml

    - name: res_3-7
      pre_scf: yes
      guess: QCLO
      sp: yes
      force: yes


* summary

  Displays a summary of the calculation. 
  There are three different methods depending on the data format.

  * boolean

    Outputs a standard summary (True).

    .. code-block:: yaml

      - name: res_3-7
        pre_scf: yes
        guess: QCLO
        sp: yes
        summary: yes


  * string

    Output according to the given string.
    Specific strings are replaced by the corresponding content.


    ============== ===================================
    keyword        content
    ============== ===================================
    {NUM_OF_ATOMS} number of atoms
    {NUM_OF_AO}    number of AOs
    {NUM_OF_MO}    number of MOs
    {METHOD}       method
    {IS_CONVERGED} Whether the SCF is converged or not
    {ITERATION}    iteration
    {TOTAL_ENERGY} total energy
    {GRADIENT_RMS} gradient RMS
    ============== ===================================


    .. code-block:: yaml

      - name: res_3-7
      pre_scf: yes
      guess: QCLO
      sp: yes
      summary: "atoms: {NUM_OF_ATOMS} iterations: {ITERATION}"



  * dict

    If you want to export to a file, you can use this format. The output file is written in appendix mode.

    * format

      Output according to the format string.

    * filepath

      Specify the file path to be output.

      .. code-block:: yaml

        - name: res_3-7
          pre_scf: yes
          guess: QCLO
          sp: yes
          summary:
            format: "atoms: {NUM_OF_ATOMS} iterations: {ITERATION}"
            filepath: "summary.txt"


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
