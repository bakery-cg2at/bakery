XML structure
======================
Here we describe structure of XML settings file. Every **dot** describes the tree structure, e.g. **root.a** represents :xml:`<root><a></a></root>`.
The notation *1* or *1..N* express the number of elements in the tree.

settings
  *1*, root

cg_molecule
  *1..N*, different types of atomistic molecules
  defines the list of fragments in given CG molecule.

  **attributes**
     - equilibrate_charges: if set to "1" then the partial charges will be adjusted to make the fragment neutral.

cg_molecule.name
  coarse-grained name of the molecule

cg_molecule.ident
  name of the molecule after backmapping

cg_molecule.source_coordinate.file
  name of coordinate file that stores atomistic fragments which represents molecule

  **attributes**
     - molecule_degree: declares that given coordinate file will be used when the molecule is in *degree*

cg_molecule.source_topology.file
  name of topology file (GROMACS format) that describes molecule

  **attributes**
     - molecule_degree: declares that given coordinate file will be used when the molecule is in *degree*

cg_molecule.cg_bead
  *1..N*, definition of coarse-grained bead

cg_bead.name
  the name of CG bead

cg_bead.type
  the type of CG bead

cg_bead.beads
  the list of atomistic particles from which the CG bead will be built

  **attributes**
     - degree: the degree of particular CG bead
     - molecule_degree: the degree of the whole CG molecule
     - active_site: the atom that be used to create atomistic bond (format: `<chain_name>:<atom_name>`

cg_bead.beads.charge_map
  the list of partial charges that will be assigned to the atoms.
  e.g. 0.0 0.0 * * 1.0; by using *** user can declare that the partial charge for that atom
  will be taken from topology file.

cg_configuration
 *1*, definition of the CG configuration

cg_configuration.coordinate
  the name of input file with the coarse-grained coordinate file
cg_configuration.topology
  the name of input file with the coarse-grained coordinate file

hybrid_configuration.file
  the output hybrid coordinate file

hybrid_topology
  *1*, the definition of output hybrid topology. The format is similar to the one
  used by GROMACS

hybrid_topology.file
  the name of output file

hybrid_topology.(bonds|angles|dihedrals)
  *1..N* defines the atomistic cross-bonds
  The content contains a list fo pairs in the format <molecule name>:<atom name>

  **attributes**
     - params: parameters of the bonds

hybrid_topology.include
  the include section of hybrid top  (do not add #include part)

hybrid_topology.molecule_type
  *1*, the molecule_type section

hybrid_topologoy.molecule_type.name
  the name of molecule

hybrid_topology.molecule_type.exclusion
  the exclusion

hybrid_topology.system
  the name of system