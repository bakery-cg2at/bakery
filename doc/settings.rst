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
      <source_file>dodecane_single.gro</source_file>
      <source_topology>topol.top</source_topology>
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
    <format>GROMACS</format>
    <file>cg_conf.gro</file>
  </cg_configuration>
  ...

In this example our input CG coordinate file is stored in `cg_conf.gro` file which is in GROMACS_ format. `Baker` supports also LAMMPS_ files. The definition is a bit different:

.. code-block:: xml
  
  <cg_molecule>
  ...
  </cg_molecule>
  <cg_configuration>
    <data_files>data.1000000</data_files>
    <input_files>in.epoxy</input_files>
    <format>LAMMPS</format>
    <type2chain>
      1:EPO
      2:EPO
      3:EPO
      4:EPO
      5:DET
      6:DET
    </type2chain>
    <name_seq chain_name="EPO">
      A1 B1 C1 D1 C2 B2 A2
    </name_seq>
    <name_seq chain_name="DET">
      E1 F1 E2
    </name_seq>
  </cg_configuration>

This is an example of input CG configuration from LAMMPS_ simulation. Here we define input `data file`_ in :xml:`<data_files>` tag (there can be several files, it will be read one by one) and `input files`_ which can also hold several files, that are read one by one.
The :xml:`<type2chain>` tag describes mapping between type of particle to molecule name, e.g. molecular chain and :xml:`<name_seq>` defines names of coarse-grained particle names. LAMMPS_ data files does not store such information (in contrast to input files of GROMACS).

`<hybrid_configuration>`
++++++++++++++++++++++++++

This section describes output hybrid coordinate file. Currently we support only GROMACS output file.

.. code-block:: xml
  
  <hybrid_configuration>
    <file>hyb_conf.gro</file>
    <format>GRO</format>
  </hybrid_configuration>

In this block above we define that the output hybrid coordinate file will be stored in `hyb_conf.gro` and will be in format of GROMACS.


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

XML structure
======================
Here we describe XML structure. Every **dot** describes the tree structure, e.g. **root.a** represents :xml:`<root><a></a></root>`. 
The notation *1* or *1..N* express the number of elements in the tree.


settings
  *1*, root

cg_molecule
  *1..N*, different types of atomistic molecules

cg_molecule.name
  name of the molecule after mapping

cg_molecule.ident
  name of the molecule in atomistic topology

cg_molecule.source_file
  name of coordinate file that stores atomistic fragments which represents molecule

cg_molecule.source_topology
  name of topology file (GROMACS format) that describes molecule

cg_molecule.topology.cg_beads
  *1*, definition of the coarse-grained beads and corresponding atomistic particles

cg_molecule.topology.cg_beads.cg_bead
  *1..N*, definition of coarse-grained bead

cg_molecule.topology.cg_beads.cg_bead.name
  the name of CG bead

cg_molecule.topology.cg_beads.cg_bead.type
  the type of CG bead

cg_molecule.topology.cg_beads.cg_bead.mapping
  the corresponding mass mapping scheme

cg_molecule.topology.cg_beads.cg_bead.beads
  the list of atomistic particles from which the CG bead will be built

cg_molecule.topology.cg_bonded
  *1*, the bonded terms of CG molecule

Here we can have three types of bonded terms: `bond`, `angle`, `dihedral`. In the next lines we use term `TYPE`.

cg_molecule.topology.cg_bonded.`TYPE`
  *1..N*, the definition of the coarse-grained bonded term

cg_molecule.topology.cg_bonded.`TYPE`.name
  the name of the bonded term

cg_molecule.topology.cg_bonded.`TYPE`.params
  the parameters of the CG bonded term (that has to be in the GROMACS format)

cg_molecule.topology.cg_bonded.`TYPE`.beads
  the list of pairs, triplets or quadruplets that will form this bonded
  term


cg_molecule.maps
  *1*, definition of mass mapping terms

cg_molecule.maps.map
  *1..N*, particular mass definition

cg_molecule.maps.map.name
  the name of map that will be used in `<mapping>` tag in the definition of CG beads

cg_molecule.maps.map.weights
  the list of numbers that defines mass of atomistic particles


cg_configuration
 *1*, definition of the CG configuration

cg_configuration.format
  the format of the coarse-grained configurato, currently it can be GROMACS or LAMMPS

cg_configuration.format.file
  the name of input file with the coarse-grained coordinate file

hybrid_configuration.file
  the output hybrid coordinate file

hybrid_configuration.format
  the format of output hybrid coordinate file, currently only GRO is accepted

hybrid_topology
  *1*, the definition of output hybrid topology. The format is similar to the one
  used by GROMACS

hybrid_topology.file
  the name of output file

hybrid_topology.include
  the include section of hybrid top

hybrid_topology.molecule_type
  *1*, the molecule_type section

hybrid_topologoy.molecule_type.name
  the name of molecule

hybrid_topology.molecule_type.exclusion
  the exclusion

hybrid_topology.system
  the name of system


**References**

.. target-notes::

.. _VOTCA: http://votca.org/
.. _GROMACS: http://www.gromacs.org/
.. _LAMMPS: http://lammps.sandia.gov/
.. _data file: http://lammps.sandia.gov/doc/2001/data_format.html
.. _input files: http://lammps.sandia.gov/doc/99/input_commands.html
.. _GROMACS manual: http://www.gromacs.org/Documentation/Manual
