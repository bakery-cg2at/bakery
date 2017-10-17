XML settings file
=================

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
      coarse-grained name of the residue

    cg_molecule.ident
      name of the residue after backmapping

cg_molecule.source_coordinate
  name of the coordinate file if there is no need to define it per molecule degree.

    cg_molecule.source_coordinate.file
      name of the coordinate file that stores atomistic fragment

      **attributes**
         - molecule_degree: declares that given coordinate file will be used when the molecule is in *degree*

      **remark**
         For every *molecule_degree* a separate coordinate file can be defined

cg_molecule.source_topology
  name of the topology file if there is no need to define it per molecule degree.

    cg_molecule.source_topology.file
      name of topology file (GROMACS format) that stores atomistic fragment

      **attributes**
         - molecule_degree: declares that given coordinate file will be used when the molecule is in *degree*

      **remark**
         For every *molecule_degree* a separate topology file can be defined

cg_molecule.cg_bead
  *1..N*, definition of coarse-grained beads

    cg_bead.name
      the name of CG bead

    cg_bead.type
      the type of CG bead

    cg_bead.beads
      the list of atomistic particles from which the CG bead will be built

      **attributes**
         - degree: the degree of particular CG bead
         - molecule_degree: the degree of the whole CG molecule
         - active_site: the atom that be used to create atomistic bond (format: `<chain_name>:<atom_name>:<max_bonds>`

      **remark**
         *max_bonds* in *active_site* definition defines the maximum number of bonds that given active site can have.
         This take into account existing atomistic bonds. Keep in mind that the package does not know anything about
         double bonds. Therefore, if your active site is already connected withing the fragment with other atoms
         by double bonds, you should reduce the maximum number of bonds. For instance, if your carbon has already double bond
         with other atom, then *max_bonds* will be 3 instead of 4.

    cg_bead.beads.charge_map
      the list of partial charges that will be assigned to the atoms.
      e.g. 0.0 0.0 * * 1.0; by using *** user can declare that the partial charge for that atom
      will be taken from topology file.

    cg_bead.beads.remove
      the list of atoms to remove whenever a given active site is used.

      **attributes**
         - active_site: the name of active site, it should be the same as defined in *active_site* attribute of *cg_beads.beads*

cg_configuration
 *1*, definition of the CG configuration

cg_configuration.coordinate
  the name of input file with the coarse-grained coordinate file
cg_configuration.topology
  the name of input file with the coarse-grained coordinate file

hybrid_configuration.file
  the output hybrid coordinate file

hybrid_configuration.format
  the output of hybrid coordinate file, currently only GROMACS is supported (.gro file)

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
         - typeid: optional parameter with the numeric label of the bond, useful for converting output topology file to LAMMPS format.

      **remark**
         the parameters in *params* are directly injected into output topology file. Therefore you can put comments
         after *;* character as in the normal GROMACS topology file.

    hybrid_topology.include
      the include section of hybrid top  (do not add #include part) Put every entry in the new line.

    hybrid_topology.molecule_type
      *1*, the molecule_type section

        hybrid_topologoy.molecule_type.name
          the name of molecule

        hybrid_topology.molecule_type.exclusion
          the global exclusion value, based on this the exclusion for non-bonded interactions will be calculated.

        hybrid_topology.molecule_type.exclusion_cg
          the exclusion only for CG particles, this value will be used to generate exclusion list only for CG particles

        hybrid_topology.molecule_type.exclusion_at
          the exclusion only for atomistic particles, this value will be used to generate exclusion list only for atomistic particles

hybrid_topology.system
  the name of system