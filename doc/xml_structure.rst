XML structure
======================
Here we describe structure of XML settings file. Every **dot** describes the tree structure, e.g. **root.a** represents :xml:`<root><a></a></root>`.
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