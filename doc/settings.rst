.. role:: xml(code)
   :language: xml

.. role:: python(code)
   :language: python

###################
Configuration file
###################

The description of backmapping configuration is stored in XML file format. The format of file is very similarly
to what is used by VOTCA_. In next two examples we will explain step-by-step how to create that file.


Simple mapping - dodecane
==========================

The original VOTCA_ mapping file is shown below. It represents mapping where two united-atoms are replaced by the one coarse-grained bead.

.. code-block:: xml

  <cg_molecule>
  <name>DOD</name> <!-- the name of the molecule after mapping -->
  <ident>DOD</ident> <!-- the name of the molecule in atomistic topology --> 
  <topology>
    <cg_beads>
      <cg_bead>
        <name>A1</name>
        <type>A</type>
        <mapping>A</mapping>
        <beads>
          1:DOD:C1 1:DOD:C2
        </beads>
      </cg_bead>
      <cg_bead>
        <name>B1</name>
        <type>B</type>
        <mapping>B</mapping>
        <beads>
          1:DOD:C3 1:DOD:C4
        </beads>
      </cg_bead>
      <cg_bead>
        <name>B2</name>
        <type>B</type>
        <mapping>B</mapping>
        <beads>
          1:DOD:C5 1:DOD:C6
        </beads>
      </cg_bead>
      <cg_bead>
        <name>B3</name>
        <type>B</type>
        <mapping>B</mapping>
        <beads>
          1:DOD:C7 1:DOD:C8
        </beads>
      </cg_bead>
      <cg_bead>
        <name>B4</name>
        <type>B</type>
        <mapping>B</mapping>
        <beads>
          1:DOD:C9 1:DOD:C10
        </beads>
      </cg_bead>
      <cg_bead>
        <name>A2</name>
        <type>A</type>
        <mapping>A</mapping>
        <beads>
          1:DOD:C11 1:DOD:C12
        </beads>
      </cg_bead>
    </cg_beads>
    <cg_bonded>
      <bond>
        <name>bond</name>
        <beads>
          A1 B1
          A2 B4
          B1 B2
          B2 B3
          B3 B4
        </beads>
      </bond>
      <angle>
        <name>angle_ch3</name>
        <beads>
          A1 B1 B2
          B3 B4 A2
        </beads>
      </angle>
      <angle>
        <name>angle_ch2</name>
        <beads>
          B1 B2 B3
          B2 B3 B4
        </beads>
      </angle>
      <dihedral>
        <name>dihedral</name>
        <beads>
          A1 B1 B2 B3
          B1 B2 B3 B4
          B2 B3 B4 A2
        </beads>
      </dihedral>
    </cg_bonded>
  </topology>
  <maps>
    <map>
      <name>A</name>
      <weights>15.035 14.027</weights>
    </map>
    <map>
      <name>B</name>
      <weights>14.027 14.027</weights>
    </map>
  </maps>
  </cg_molecule>

`<settings>`
+++++++++++++++

First we have to surround whole file with :xml:`<settings>` tag:

.. code-block:: xml

  <settings>
    <cg_molecule>
    ...
    </cg_molecule>
  </settings>

`<cg_molecule>`
++++++++++++++++++

Next, we have to set up source coordinate file for each of molecule. We add tag :xml:`<source_file>` and :xml:`<source_topology>`
which defines topology of atomistic fragment.

.. code-block:: xml

    <cg_molecule>
      ...
      <source_coordinate><file>dodecane_single.gro</file></file></source_coordinate>
      <source_topology><file>topol.top</file></source_topology>
      ...
    </cg_molecule>

There could be multiple :xml:`<cg_molecule>` sections that will describe different types of atomistic molecule.

`<cg_configuration>`
+++++++++++++++++++++++++

After the list of CG molecules, we have to define
input format and a name of file of CG coordinate file that is reverse mapped.

.. code-block:: xml

  <cg_molecule>
  ...
  </cg_molecule>

  <cg_configuration>
    <coordinate>cg_conf.gro</coordinate>
    <topology>cg_topol.top</topology>
  </cg_configuration>
  ...

`<hybrid_configuration>`
++++++++++++++++++++++++++

This section describes output hybrid coordinate file. Currently we support only GROMACS output file.

.. code-block:: xml
  
  <hybrid_configuration>
    <file>hyb_conf.gro</file>
  </hybrid_configuration>

In this block above we define that the output hybrid coordinate file will be stored in `hyb_conf.gro`.


`<hybrid_topology>`
+++++++++++++++++++++++++++

The options for output hybrid topology file are defined in tag :xml:`<hybrid_topology>`. Very simple example is shown below.

.. code-block:: xml

  <hybrid_topology>
    <file>hyb_topol.top</file>
    <include>
      #include "/usr/share/gromacs/top/oplsaa.ff/forcefield.itp"
    </include>
    <molecule_type>
        <name>PE</name>
        <exclusion>3</exclusion>
    </molecule_type>
    <system>PE</system>
  </hybrid_topology>

Here we define that output file is `hyb_topol.top`. We define include section in :xml:`<include></include>`. :xml:`<molecule_type>` tag defines the name of molecule in `<name>` tag and the exclusion rule. In this case, every particle that is separated by 1, 2, 3 bonds will be excluded from non-bonded interactions. This follows the GROMACS_ rules (see `GROMACS manual`_ for further information). 

.. Complex mapping - EPON-828 network
   -------------------------------

.. _xml-structure:


**References**

.. target-notes::

.. _VOTCA: http://votca.org/
.. _GROMACS: http://www.gromacs.org/
.. _LAMMPS: http://lammps.sandia.gov/
.. _data file: http://lammps.sandia.gov/doc/2001/data_format.html
.. _input files: http://lammps.sandia.gov/doc/99/input_commands.html
.. _GROMACS manual: http://www.gromacs.org/Documentation/Manual
