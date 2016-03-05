#  Copyright (C) 2015-2016
#      Jakub Krajniak (jkrajniak at gmail.com)
#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# The file was modyfied to support non-standard entries that appears in gromacs topology
# file.

from collections import namedtuple
import espressopp
from topology_helper import *

__doc__ = """This Python module allows one to use GROMACS data files as the
input to an ESPResSo++ simulation, set interactions for given
particle types and convert GROMACS potential tables into
ESPResSo++ tables.
It containts functions: read(), setInteractions(), convertTable()
"""

GromacsSystem = namedtuple(
    'GromacsSystem', [
        'defaults',
        'types',
        'masses',
        'charges',
        'res_ids',
        'atomtypeparams',
        'bondtypes',
        'bondtypeparams',
        'angletypes',
        'angletypeparams',
        'dihedraltypes',
        'dihedraltypeparams',
        'pairtypes',
        'pairtypeparams',
        'nonbond_params',
        'exclusions',
        'x', 'y', 'z',
        'vx', 'vy', 'vz',
        'Lx', 'Ly', 'Lz'
    ])


def read(gro_file, top_file="", doRegularExcl=True, defines=None):
    """ Read GROMACS data files.

    Keyword arguments:
    gro_file -- contains coordinates of all particles, the number of particles, velocities and box size.
    top_file -- contains topology information. Included topology files (.itp) are also read
    doRegularExcl -- if True, exclusions are generated automatically based on the nregxcl parameter (see gromacs manual)
    """

    if defines is None:
        defines = {}

    # read gro file
    if gro_file != "":
        f = open(gro_file)
        f.readline() # skip comment line
        total_num_particles = int(f.readline())

        # store coordinates and velocities
        x, y, z = [], [], []
        vx, vy, vz = [], [], []
        for i in range(total_num_particles):
            s = f.readline()[20:69]
            # coordinates
            x.append(float(s[0:8]))
            y.append(float(s[8:16]))
            z.append(float(s[16:24]))

            if len(s.split()) > 3:
                # velocities
                vx.append(float(s[24:32]))
                vy.append(float(s[32:40]))
                vz.append(float(s[40:49]))

        # store box size
        Lx, Ly, Lz = map(float, f.readline().split()) # read last line, convert to float
        f.close()

    # read top and itp files
    masses, charges = [], [] # mases and charges of the whole configuration
    types=[] # tuple: atomindex(int) to atomtypeid(int)
    bonds={} # dict: key bondtypeid value: tuple of bond pairs
    angles={} # dict: key angletype value: tuple of triples
    dihedrals = {} #same...
    pairs_1_4 = {}  # dict: key pairtype value: tuple of pairs
    exclusions=set() #list of atom pairs no considered in non-bonded interactions

    defaults={} # gromacs default values
    atomtypeparams={} # a dict: key atomtypeid , value : class storing actual parameters of each type e.g. c6, c12, etc..
    use_atomtypeparams = {}  # dict with the atomtypes that are use in the topology
    nonbond_params = {}
    use_nonbond_params = {}
    bondtypeparams={} # same for bonds
    angletypeparams={} # same for angles
    dihedraltypeparams={} # same for dihedrals
    pairtypeparams={}
    use_pairtypeparams={}

    if top_file != "":
        #f = open(top_file)
        # FileBuffer: a class which behaves like a file, but all lines are in memory
        # we use this for emulating a 'preprocessor' which handles the #include
        # statements in the .top and .itp files
        fb=FileBuffer()

        defines = {}
        FillFileBuffer(top_file, fb, defines=defines)
        f = PostProcessFileBuffer(fb, defines)

        print "Reading top file: "+top_file
        line = ''
        a, p = 0, 0
        bondtypecount, angletypecount, dihedraltypecount=0,0,0
        readdefaults, readattypes, readnonbondtypes, readpairtypes, readbdtypes, readantypes, readdhtypes  = (
                False, False, False, False, False, False, False)
        defaults={} # gromacs default values
        atomtypes={} # a dict: key atomtypename(str) value: atomtypeid(int)
        bondtypes={} # a dict: key atomindex(int),atomindex(int)  value: bondtypeid(int)
        angletypes={} # a dict: key atomindex(int), atomindex(int),atomindex(int) value: angletypeid(int)
        dihedraltypes={} # a dict: key atomtindex(int), atomindex(int), atomindex(int),atomindex(int) value: dihedraltypeid(int)
        nonbonds = {}

        atnum_attype = {}
        wildcard_type = None  # Use for dihedrals. Special name 'X'.

        # it was moved out of "if" statement
        #atomtypeparams={} # a dict: key atomtypeid , value : class storing actual parameters of each type e.g. c6, c12, etc..
        #bondtypeparams={} # same for bonds
        #angletypeparams={} # same for angles
        #dihedraltypeparams={} # same for dihedrals

        #atomparams={} # key: atomindex(int) value: per atom parameters e.g. q, mass
        molecules=[]
        #molecules = {} # key: moleculeid value: name (string)
        readmolecules = False

        skip_section = False

        for line in f.lines:
            if line[0] == ";":  # skip comment line
                continue

            if skip_section and line.startswith('#end'):
                skip_section = False
                continue

            if skip_section:
                continue

            if line.startswith('#ifdef'):
                define_tmp = line.split()
                if len(define_tmp) > 1:
                    skip_section = defines.get(define_tmp[1], False)
                else:
                    skip_section = True
                continue

            if line.startswith('#else'):
                skip_section = True
                continue

            if line.startswith("#define"):
                define_tmp = line.split()
                defines[define_tmp[1]] = True
                continue

            if 'defaults' in line: # store some gromacs default values
                readdefaults =True
                continue

            if readdefaults:
                if line.strip() == "" or '[' in line: # end of defaults section
                    readdefaults=False
                    print "Defaults: ", defaults
                else:
                    fields=line.split()
                    if len(fields)==5:
                        defaults={"nbtype":fields[0], "combinationrule":int(fields[1]),
                        "genpairs":fields[2], "fudgeLJ":float(fields[3]), "fudgeQQ":float(fields[4])}
                    else:
                        defaults={"nbtype":fields[0], "combinationrule":fields[1]}

            if 'atomtypes' in line: # map atom types (espressopp++ uses ints)
                readattypes = True
                print "Reading atomtypes (GROMACS: ESPResSo++): "
                continue

            if readattypes:
                if line.strip() == "" or '[' in line: # end of atomtypes section
                    readattypes = False
                    atomtypes.update({'X': a})
                    wildcard_type = a
                    atnum_attype['X'] = 'X'
                else:
                    fields=line.split(';')[0].split()
                    attypename = fields[0]

                    #make a map containing the properties
                    # sig, eps may be c6 and c12: this is specified in the defaults
                    # and converted later
                    if fields[0].startswith('opls'):
                        tmpprop = {
                            'atnum': fields[1],
                            'mass': float(fields[3]),
                            'charge': float(fields[4]),
                            'particletype': fields[5],
                            'sig': float(fields[6]),
                            'eps': float(fields[7])
                            }
                        atnum_attype[attypename] = fields[1]
                    elif len(fields)==7:
                        tmpprop={
                                "atnum":int(fields[1]),
                                "atname": fields[0],
                                "mass": float(fields[2]),
                                "charge":float(fields[3]),
                                "particletype":fields[4],
                                "sig":float(fields[5]),
                                "eps":float(fields[6])}
                    elif len(fields) == 8:
                        tmpprop = {
                            'atnum': fields[1],
                            'mass': float(fields[3]),
                            'charge': float(fields[4]),
                            'sig': float(fields[6]),
                            'eps': float(fields[7])
                            }
                    else:
                        print('AA other: {}'.format(fields))
                        tmpprop={
                            "atnum": fields[0],
                            "mass": float(fields[1]),
                            "charge": float(fields[2]),
                            "particletype": fields[3],
                            "sig": float(fields[4]),
                            "eps":float(fields[5])
                        }

                if attypename not in atomtypes:
                    atomtypes.update({attypename:a}) # atomtypes is used when reading the "atoms" section
                    atomtypeparams.update({a:tmpprop})
                    a += 1

            if 'nonbond_params' in line:
                print 'Reading [ nonbond_params ]'
                readnonbondtypes = True
                continue
            if readnonbondtypes:
                l = line.strip()
                if l == '' or '[' in l:
                    readnonbondtypes = False
                else:
                    fields = l.split(';')[0].split()
                    if len(fields) > 0:
                        a1, a2, fn, c6, c12 = fields[:5]
                        if int(fn) != 1:
                            continue
                        at1, at2 = sorted([atomtypes.get(a1), atomtypes.get(a2)])
                        if (at1, at2) not in nonbond_params:
                            nonbond_params[(at1, at2)] = {
                                'sig': float(c6),
                                'eps': float(c12)
                                }

            if 'pairtypes' in line:
                print 'Reading [ pairtypes ]'
                readpairtypes = True
                continue
            if readpairtypes:
                l = line.strip()
                if l == '' or '[' in l:
                    readpairtypes = False
                else:
                    fields = l.split(';')[0].split()
                    if len(fields) > 0:
                        a1, a2, fn, c6, c12 = fields[:5]
                        if int(fn) != 1:
                            continue
                        at1, at2 = sorted([atomtypes.get(a1), atomtypes.get(a2)])
                        if at1 and at2:
                            if (at1, at2) not in pairtypeparams:
                                pairtypeparams[(at1, at2)] = {
                                    'sig': float(c6),
                                    'eps': float(c12)
                                    }

            if 'bondtypes' in line:
                print 'Reading [ bondtypes ]'
                readbdtypes = True
                continue

            if readbdtypes:
                if line.strip() == "" or '[' in line: # end of bondtypes section
                    readbdtypes = False
                else:
                    tmp = line.split()
                    i, j = tmp[:2]
                    p=ParseBondTypeParam(line)
                    #check if this type has been defined before
                    bdtypeid=FindType(p, bondtypeparams)
                    if bdtypeid==None:
                        bdtypeid=len(bondtypeparams)
                        bondtypeparams.update({bdtypeid:p})

                    if i in bondtypes:
                        bondtypes[i].update({j:bdtypeid})
                    else:
                        bondtypes.update({i:{j:bdtypeid}})

            if 'angletypes' in line:
                print 'Reading [ angletypes ]'
                readantypes = True
                continue

            if readantypes:
                if line.strip() == "" or '[' in line: # end of angletypes section
                    readantypes = False
                else:
                    tmp = line.split()
                    i, j, k = tmp[:3]
                    p=ParseAngleTypeParam(line)

                    atypeid=FindType(p, angletypeparams)
                    if atypeid==None:
                        atypeid=len(angletypeparams)
                        angletypeparams.update({atypeid:p})

                    if i in angletypes:
                        if j in angletypes[i]:
                            angletypes[i][j].update({k:atypeid})
                        else:
                            angletypes[i].update({j:{k:atypeid}})
                    else:
                        angletypes.update({i:{j:{k:atypeid}}})
                    #print "FOUND angletype: ", angletypecount, " : ", p.parameters
                    #angletypeparams.update({angletypecount:p})
                    #angletypecount+=1

            if 'dihedraltypes' in line:
                print 'Reading [ dihedraltypes ]'
                readdhtypes = True
                continue

            if readdhtypes:
                if line.strip() == "" or '[' in line: # end of angletypes section
                    readdhtypes = False
                else:
                    tmp = line.split()
                    try:
                        int(tmp[4])
                    except ValueError:
                        continue
                    i, j, k, l = tmp[:4]
                    p=ParseDihedralTypeParam(line)
                    if p is False:
                        print('Skip dihedral line: {}'.format(line))
                        continue

                    dtypeid=FindType(p, dihedraltypeparams)
                    if dtypeid==None:
                        dtypeid=len(dihedraltypeparams)
                        dihedraltypeparams.update({dtypeid:p})
                    if i in dihedraltypes:
                        if j in dihedraltypes[i]:
                            if k in dihedraltypes[i][j]:
                                dihedraltypes[i][j][k].update({l:dtypeid})
                            else:
                                dihedraltypes[i][j].update({k:{l:dtypeid}})
                        else:
                            dihedraltypes[i].update({j:{k:{l:dtypeid}}})
                    else:
                        dihedraltypes.update({i:{j:{k:{l:dtypeid}}}})

            if 'molecules' in line: # store number of chains
                readmolecules = True
                continue

            if readmolecules:
                if line.strip() == "" or '[' in line: # end of molecules section
                    readmolecules = False
                else:
                    mol, nrmol = line.strip().split()
                    #we have to check if the same molecules comes multiple times in the molecules section
                    if len(molecules) == 0:
                      molecules.append({'name':mol, 'count':int(nrmol)})
                    elif molecules[-1]['name'] == mol: #check if mol was added earlier already
                        molecules[-1]['count'] = molecules[-1]['count'] + int(nrmol) #update count
                    else: molecules.append({'name':mol, 'count':int(nrmol)}) #if mol newly added


        molstartindex=0 #this is the index of the first atom in the molecule being parsed
        res_idx = 0  # index of molecule like single polymer chain.


        f.seek(0) # Now we search for bonds, angles definitions and start from the beginning of the file buffer

        for mol in molecules:
            print "Preparing molecule %s... " % mol['name']
            print "-----------------------------"

            # find and store number of molecules
            num_molecule_copies=mol['count']
            # this does not what the name suggests....
            nrexcl = storeMolecules(f, molecules, mol)
            # find and store atom types
            types, masses, charges, num_atoms_molecule, res_ids = \
                storeAtoms(f, defaults, types, atomtypes, atomtypeparams, use_atomtypeparams,
                           nonbond_params, use_nonbond_params, masses, charges,
                           num_molecule_copies,
                           res_idx)

            # find and store bonds
            bonds = storeBonds(f, types, bondtypes, bondtypeparams, bonds,
                               num_atoms_molecule, num_molecule_copies,
                               molstartindex, exclusions, nrexcl, doRegularExcl)

            # find and store angles
            angles = storeAngles(f, types, angletypes, angletypeparams, angles,
                                 num_atoms_molecule, num_molecule_copies, molstartindex)

            # find and store dihedrals
            dihedrals = storeDihedrals(f, types, dihedraltypes, dihedraltypeparams, dihedrals,
                                       num_atoms_molecule, num_molecule_copies,
                                       molstartindex, atomtypeparams,
                                       wildcard_type)

            pairs_1_4 = storePairs(f, defaults, types, pairtypeparams, use_pairtypeparams,
                                   atomtypeparams, pairs_1_4,
                                   num_atoms_molecule, num_molecule_copies,
                                   molstartindex)
            if doRegularExcl:
                storeExclusions(exclusions, nrexcl, bonds, angles, dihedrals)

            molstartindex+=num_molecule_copies*num_atoms_molecule
            res_idx += num_molecule_copies


    # Update typeparams
    use_keys = [s[0] for s in bonds]
    bondtypeparams = {k: v for k, v in bondtypeparams.iteritems() if k in use_keys}
    use_keys = [s[0] for s in angles]
    angletypeparams = {k: v for k, v in angletypeparams.iteritems() if k in use_keys}
    use_keys = [s[0] for s in dihedrals]
    dihedraltypeparams = {k: v for k, v in dihedraltypeparams.iteritems() if k in use_keys}

    params = []

    unpackvars=[]

    # The data is packed into a touple, unpackvars contains a string which
    # tells the user which kind of data was read.

    print 'Found default values', defaults
    print 'Found {} types'.format(len(types))
    print 'Found {} nonbonded_pairs'.format(len(use_nonbond_params))
    print 'Found {} masses'.format(len(masses))
    print 'Found {} charges'.format(len(charges))
    print 'Found {} atom type parameters'.format(len(use_atomtypeparams))
    print 'Found {} bonds'.format(len(bonds))
    print 'Found {} bond type parameters'.format(len(bondtypeparams))
    print 'Found {} angles'.format(len(angles))
    print 'Found {} angle type parameters'.format(len(angletypeparams))
    print 'Found {} dihedrals'.format(len(dihedrals))
    print 'Found {} dihedral type parameters'.format(len(dihedraltypeparams))
    print 'Found {} 1-4 pair type parameters'.format(len(use_pairtypeparams))
    print 'Found {} 1-4 pairs'.format(len(pairs_1_4))
    print 'Found {} bond exclusions'.format(len(exclusions))
    print 'Found box: {}'.format([Lx, Ly, Lz])
    print 'Found {} particles'.format(len(x))

    gromacs_system = GromacsSystem(
        defaults, types, masses, charges, res_ids, use_atomtypeparams,
        bonds, bondtypeparams, angles, angletypeparams,
        dihedrals, dihedraltypeparams, pairs_1_4, use_pairtypeparams,
        use_nonbond_params, exclusions,
        x, y, z, vx, vy, vz, Lx, Ly, Lz)

    return gromacs_system


def storeMolecules(f, molecules, mol=""):
    nrexcl=0
    line = ''
    line=f.readlastline()
    while not 'moleculetype' in line:
        line = f.readline()
        if not line: break # break out of while if EOF
    line = f.readline()
    while(len(line) > 1 and not '[' in line):
        if line[0] == ";":   # skip comment lines
            #print "skipping line: "+line.strip("\n")
            line = f.readline()
            continue
        fields=line.split()
        #mol = fields[0]
        nrexcl=int(fields[1])
        line = f.readline()
    return nrexcl


def storeAtoms(f, defaults, types, atomtypes,
               atomtypeparams,
               use_atomtypeparams,
               nonbondedparams,
               use_nonbond_params,
               masses,
               charges,
               num_molecule_copies,
               res_idx):
    line = ''
    types_tmp = []
    charge_tmp =[]
    mass_tmp=[]
    molecule_index = []
    pos = f.tell()

    combinationrule = defaults['combinationrule']

    line=f.readlastline()
    while not 'atoms' in line:
        line = f.readline()
        if not line: break # break out of while if EOF
    line = f.readline()
    while(len(line) > 1 and not '[' in line):
        if line[0] == ";":   # skip comment lines
            line = f.readline()
            continue
        fields=line.split()
        attypeid=atomtypes[fields[1]] # map str type to int type
        types_tmp.append(attypeid)
        if len(fields) > 6:
            # this atom has a charge different from its atomtype
            charge_tmp.append(float(fields[6]))
        else:
            #look up default values for this atom type
            charge_tmp.append(atomtypeparams[attypeid]['charge'])
        if len(fields) > 7:
            # also has a special mass
            mass_tmp.append(float(fields[7]))
        else:
            mass_tmp.append(atomtypeparams[attypeid]['mass'])

        use_atomtypeparams.update({attypeid: atomtypeparams[attypeid]})

        line = f.readline()

    # Convert to sigma/epsilon
    if combinationrule == 1:
        for k, v in use_atomtypeparams.iteritems():
            c6, c12 = float(v['sig']), float(v['eps'])
            sig, eps = convertc6c12(c6, c12)
            print '{}, Convert C6({}), C12({}) to sig({}), eps({})'.format(
                k, c6, c12, sig, eps)
            use_atomtypeparams[k]['sig'] = sig
            use_atomtypeparams[k]['eps'] = eps

    # Prepare nonbonded params to contains only that store in atomtypes
    for k, v in nonbondedparams.iteritems():
        if k[0] in use_atomtypeparams and k[1] in use_atomtypeparams:
            use_nonbond_params.update({k: v})
            if combinationrule == 1:
                c6, c12 = float(v['sig']), float(v['eps'])
                sig, eps = convertc6c12(c6, c12)
                print '{}, Convert C6({}), C12({}) to sig({}), eps({})'.format(
                    k, c6, c12, sig, eps)
                use_nonbond_params[k]['sig'] = sig
                use_nonbond_params[k]['eps'] = eps

    # Reset position
    f.seek(pos)

    # extend copies of this molecule
    num_atoms_molecule = len(types_tmp)
    for i in range(num_molecule_copies):
        types.extend(types_tmp)
        charges.extend(charge_tmp)
        masses.extend(mass_tmp)
        molecule_index.extend([res_idx+i]*num_atoms_molecule)

    return types, masses, charges, num_atoms_molecule, molecule_index

def storePairs(f, defaults, types, pairtypeparams,
               use_pairtypeparams,
               atomtypeparams, pairs, num_atoms_molecule, num_molecule_copies, molstartindex):
    line = ''
    pairs_tmp = []
    top = False
    file_pos = f.tell()
    line = f.readlastline()
    fudgeLJ = float(defaults.get('fudgeLJ', 1.0))
    print('Using fudgeLJ: {}'.format(fudgeLJ))
    combinationrule = defaults['combinationrule']
    types_pairtypeid = {}

    line = f.readline().strip()
    in_section = False
    cross_pairs = False
    while line:
        if line.startswith('['):
            if 'pairs' in line or 'cross_pairs' in line:
                in_section = True
                cross_pairs =  'cross_pairs' in line
            else:
                in_section = False
            line = f.readline().strip()
            continue
        elif line.startswith(';'):
            line = f.readline().strip()
            continue
        else:
            if in_section:
                tmp = line.split(';')[0].split()
                lookup = len(tmp) <= 3
                pid1, pid2 = sorted(map(int, tmp[0:2]))
                t1, t2 = sorted([types[pid1-1], types[pid2-1]])
                pairtypeid = max(use_pairtypeparams) + 1 if use_pairtypeparams else 0
                if lookup:  # Look for parameters
                    at1 = atomtypeparams[t1]
                    at2 = atomtypeparams[t2]
                    if (t1, t2) in pairtypeparams:
                        if types_pairtypeid:
                            pairtypeid = types_pairtypeid.setdefault((t1, t2), max(types_pairtypeid.values())+1)
                        else:
                            pairtypeid = 0
                            types_pairtypeid[(t1, t2)] = 0
                        pair_params = pairtypeparams[(t1, t2)]
			use_pairtypeparams[pairtypeid] = pairtypeparams[(t1, t2)]
                    else:
                        sig_1, eps_1 = at1['sig'], at1['eps']
                        sig_2, eps_2 = at2['sig'], at2['eps']
                        eps = fudgeLJ*(eps_1*eps_2)**(1.0/2.0)
                        if combinationrule == 2:
                            sig = 0.5*(sig_1 + sig_2)
                        else:
                            sig = (sig_1*sig_2)**(1.0/2.0)
                        pairtypeid = max(use_pairtypeparams) + 1 if use_pairtypeparams else 0
                        use_pairtypeparams[pairtypeid] = {'sig': sig, 'eps': eps}
                        pairtypeparams[(t1, t2)] = use_pairtypeparams[pairtypeid]
                        types_pairtypeid[(t1, t2)] = pairtypeid
                    pairs_tmp.append((pid1, pid2, pairtypeid, cross_pairs))
                else:  # Params provided
                    if int(tmp[2]) != 1:
                        print('Warning! Supported only pair with type 1, given: {}'.format(
                            tmp[2]))
                        line = f.readline()
                        continue
                    sig = float(tmp[3])
                    eps = float(tmp[4])
                    if combinationrule == 1:
                        c6, c12 = sig, eps
                        sig, eps = convertc6c12(c6, c12)
		    use_pairtypeparams.update({pairtypeid: {
		        'sig': sig,
			'eps': eps
			}})
		    pairs_tmp.append((pid1, pid2, pairtypeid, cross_pairs))
		pairtypeid += 1
            line = f.readline()

    # Extend pairs to copies of molecule
    pairs_per_molecule = len(pairs_tmp)
    for i in range(num_molecule_copies):
        for j in range(pairs_per_molecule):
            pid1, pid2, pairtypeid, cross_pairs = pairs_tmp[j]
            ia = molstartindex + pid1 + (i * num_atoms_molecule)
            ib = molstartindex + pid2 + (i * num_atoms_molecule)
            if (pairtypeid, cross_pairs) in pairs:
                pairs[(pairtypeid, cross_pairs)].append((ia, ib))
            else:
                pairs.update({(pairtypeid, cross_pairs): [(ia, ib)]})
    return pairs

def storeBonds(f, types, bondtypes, bondtypeparams, bonds, num_atoms_molecule,\
    num_molecule_copies, molstartindex, exclusions, nregxcl, doRegularExcl=True):
    line = ''
    bonds_tmp = []
    top = False
    pos = f.tell()
    line=f.readlastline()
    local_exclusions=[] # excluded pairs of atoms within this mol (local ids)

    line = f.readline().strip()
    in_section = False
    cross_bond = False
    while line:
        if line.startswith('['):
            if 'bonds' in line or 'cross_bonds' in line:
                in_section = True
                cross_bond = 'cross_bond' in line
            else:
                in_section = False
            line = f.readline().strip()
            continue
        elif line.startswith(';'):
            line = f.readline().strip()
            continue
        else:
            if in_section:
                tmp = line.split(';')[0].split()
                lookup = len(tmp) <= 3
                pid1, pid2 = map(int, tmp[0:2])
                if lookup:
                    t1, t2 = types[pid1-1], types[pid2-1]
                    if t1 > t2:
                        t1, t2 = t2, t1
                    bdtypeid = bondtypes[t1][t2]
                else:
                    temptype = ParseBondTypeParam(line)
                    bdtypeid = FindType(temptype, bondtypeparams)
                    if bdtypeid == None:
                        bdtypeid = len(bondtypeparams)
                        bondtypeparams.update({bdtypeid: temptype})
                bonds_tmp.append((pid1, pid2, bdtypeid, cross_bond))
                if bondtypeparams[bdtypeid].automaticExclusion():
                    local_exclusions.append((pid1, pid2))
            line = f.readline().strip()
    f.seek(pos)

    # extend bonds to copies of this molecule
    bonds_per_mol = len(bonds_tmp)
    for i in range(num_molecule_copies):
        for j in range(bonds_per_mol):
            pid1, pid2, bdtypeid, cross_bond = bonds_tmp[j][0:4]
            ia=molstartindex+pid1 + (i * num_atoms_molecule) # index of copy atom i
            ib=molstartindex+pid2 + (i * num_atoms_molecule) # index of copy atom j

            if (bdtypeid, cross_bond) in bonds:
                bonds[(bdtypeid, cross_bond)].append((ia, ib))
            else:
                bonds.update({(bdtypeid, cross_bond): [(ia, ib)]})

    return bonds

def storeExclusions(exclusions, nrexcl, bonds, angles, dihedrals):
    print('Processing exclusion lists for nrexcl={}'.format(nrexcl))
    if nrexcl > 3:
        raise RuntimeError('Currently nrexcl > 3 is not supported')
    excl2list = {
        1: [bonds],
        2: [bonds, angles],
        3: [bonds, angles, dihedrals]
    }

    for l in excl2list[nrexcl]:
        for a in l.values():
            for t in a:
                exclusions.add(tuple(sorted([t[0], t[-1]])))

def storeAngles(f, types, angletypes, angletypeparams, angles, num_atoms_molecule, num_molecule_copies, molstartindex):
    line = ''
    angles_tmp = []
    pos = f.tell()
    line = f.readlastline()

    line = f.readline().strip()
    in_section = False
    cross_angle = False
    while line:
        if line.startswith('['):
            if 'angles' in line or 'cross_angles' in line:
                in_section = True
                cross_angle = 'cross_angles' in line
            else:
                in_section = False
            line = f.readline().strip()
            continue
        elif line.startswith(';'):
            line = f.readline().strip()
            continue
        else:
            if in_section:
                tmp = line.split(';')[0].split()
                lookup = len(tmp) <= 4
                pid1, pid2, pid3 = map(int, tmp[0:3])
                if lookup:
                    t1, t2, t3 = types[pid1-1], types[pid2-1], types[pid3-1]
                    try:
                        typeid = angletypes[t1][t2][t3]
                    except KeyError:
                        t1, t3 = t3, t1
                        typeid = angletypes[t1][t2][t3]
                else:
                    # Checks if we need to make new type.
                    temptype = ParseAngleTypeParam(line)
                    typeid = FindType(temptype, angletypeparams)
                    if typeid == None:
                        typeid = len(angletypeparams)
                        angletypeparams.update({typeid:temptype})
                angles_tmp.append((pid1, pid2, pid3, typeid, cross_angle))
        line = f.readline()

    f.seek(pos)

    # extend angles to copies of this molecule
    angles_per_mol = len(angles_tmp)
    for i in range(num_molecule_copies):
        for j in range(angles_per_mol):
            pid1, pid2, pid3, antypeid, cross_angle = angles_tmp[j][0:5]
            ia=molstartindex+pid1 + (i * num_atoms_molecule) # index of copy atom i
            ib=molstartindex+pid2 + (i * num_atoms_molecule) # index of copy atom j
            ic=molstartindex+pid3 + (i * num_atoms_molecule) # index of copy atom k
            if (antypeid, cross_angle) in angles:
                angles[(antypeid, cross_angle)].append((ia, ib, ic))
            else:
                angles.update({(antypeid, cross_angle): [(ia, ib, ic)]})
    return angles

def storeDihedrals(f, types, dihedraltypes, dihedraltypeparams, dihedrals,
                   num_atoms_molecule, num_molecule_copies, molstartindex,
                   atomtypeparams, wildcard_type):
    line = ''
    dihedrals_tmp = []
    pos = f.tell()
    line=f.readlastline()

    line = f.readline().strip()
    in_section = False
    cross_dih = False

    def check_type(t1, t2, t3, t4):
        wt = 'X'
        combinations = [
            (t1, t2, t3, t4),
            (wt, t2, t3, t4),
            (t1, t2, t3, wt),
            (wt, t2, t3, wt),
            (t4, t3, t2, t1),
            (wt, t3, t2, t1),
            (t4, t3, t2, wt),
            (wt, t3, t2, wt)
        ]
        for n1, n2, n3, n4 in combinations:
            try:
                return dihedraltypes[n1][n2][n3][n4]
            except KeyError:
                continue
        return dihedraltypes[n1][n2][n3][n4]

    while line:
        if line.startswith('['):
            if 'dihedrals' in line or 'cross_dihedrals' in line:
                in_section = True
                cross_dih = 'cross_dihedrals' in line
            else:
                in_section = False
            line = f.readline().strip()
            continue
        elif line.startswith(';'):
            line = f.readline().strip()
            continue
        else:
            if in_section:
                # Skip improper dihedrals, not supported yet
                if 'improper' in line:
                    line = f.readline().strip()
                    continue
                tmp = line.split(';')[0].split()
                lookup = len(tmp) <= 5
                pid1, pid2, pid3, pid4 = map(int, tmp[0:4])
                if lookup:
                    t1, t2, t3, t4 = (types[x-1] for x in map(int, tmp[0:4]))
                    try:
                        dihtypeid = check_type(t1, t2, t3, t4)
                    except KeyError:
                        t1, t2, t3, t4 = (
                            atomtypeparams[t1]['atnum'],
                            atomtypeparams[t2]['atnum'],
                            atomtypeparams[t3]['atnum'],
                            atomtypeparams[t4]['atnum']
                            )
                        try:
                            dihtypeid = check_type(t1, t2, t3, t4)
                        except KeyError as ex:
                            print(('Dihedral\n\t- {}\nnot found.'
                                   'Please define parameters in topology file').format(line))
                            print('{} {} {} {}'.format(t1, t2, t3, t4))
                            raise ex
                else:
                    #check if we need to make new type
                    temptype=ParseDihedralTypeParam(line)
                    dihtypeid=FindType(temptype, dihedraltypeparams)
                    if dihtypeid == None:
                        dihtypeid = len(dihedraltypeparams)
                        dihedraltypeparams.update({dihtypeid:temptype})
                dihedrals_tmp.append((pid1, pid2, pid3,pid4, dihtypeid, cross_dih))
        line = f.readline()

    f.seek(pos)

    # extend angles to copies of this molecule
    dihedrals_per_mol = len(dihedrals_tmp)
    for i in range(num_molecule_copies):
        for j in range(dihedrals_per_mol):
            pid1, pid2, pid3, pid4, dihtypeid, cross_dih = dihedrals_tmp[j][0:6]
            ia=molstartindex+pid1 + (i * num_atoms_molecule) # index of copy atom i
            ib=molstartindex+pid2 + (i * num_atoms_molecule) # index of copy atom j
            ic=molstartindex+pid3 + (i * num_atoms_molecule) # index of copy atom k
            id=molstartindex+pid4 + (i * num_atoms_molecule) # index of copy atom l
            if (dihtypeid, cross_dih) in dihedrals:
                dihedrals[(dihtypeid, cross_dih)].append((ia, ib, ic, id))
            else:
                dihedrals.update({(dihtypeid, cross_dih): [(ia, ib, ic, id)]})
    return dihedrals


def genParticleList(input_conf, use_velocity=False, use_charge=False, use_adress=False):
    """Generates particle list
    Args:
        input_conf: The tuple generate by read method.
        use_velocity: If set to true then velocity will be read.
        use_charge: If set to true then charge will be read.
    Returns:
        List of property names and particle list.
    """
    props = ['id', 'type', 'pos', 'res_id']
    use_mass = bool(input_conf.masses)
    use_velocity = use_velocity and bool(input_conf.vx)
    use_charge = use_charge and bool(input_conf.charges)
    if use_mass:
        props.append('mass')
    if use_velocity:
        props.append('v')
    if use_charge:
        props.append('q')

    Particle = namedtuple('Particle', props)

    particle_list = []
    num_particles = len(input_conf.types)
    for pid in range(num_particles):
        tmp = [pid+1,
               input_conf.types[pid],
               espressopp.Real3D(input_conf.x[pid], input_conf.y[pid], input_conf.z[pid]),
               input_conf.res_ids[pid]
              ]
        if use_mass:
            tmp.append(input_conf.masses[pid])
        if use_velocity:
            tmp.append(espressopp.Real3D(input_conf.vx[pid], input_conf.vy[pid], input_conf.vz[pid]))
        if use_charge:
            tmp.append(input_conf.charges[pid])
        particle_list.append(Particle(*tmp))
    return props, particle_list


def setBondedInteractions(system, bonds, bondtypeparams, ftpl=None):
    ret_list={}
    for (bid, _), bondlist in bonds.iteritems():
        if ftpl:
            fpl = espressopp.FixedPairListAdress(system.storage, ftpl)
        else:
            fpl = espressopp.FixedPairList(system.storage)
        fpl.addBonds(bondlist)
        bdinteraction=bondtypeparams[bid].createEspressoInteraction(system, fpl)
        if bdinteraction:
            system.addInteraction(bdinteraction)
            ret_list.update({bid: bdinteraction})
    return ret_list

def setPairInteractions(system, pairs, pairtypeparams, cutoff, ftpl=None):
    ret_list = {}
    for (pid, _), pair_list in pairs.iteritems():
        if ftpl:
            fpl = espressopp.FixedPairListAdress(system.storage, ftpl)
        else:
            fpl = espressopp.FixedPairList(system.storage)
        params = pairtypeparams[pid]
        if params['sig'] > 0.0 and params['eps'] > 0.0:
            fpl.addBonds(pair_list)
            print 'Pair interaction', params, ' num pairs:', len(pair_list)
            interaction = espressopp.interaction.FixedPairListLennardJones(
                system,
                fpl,
                espressopp.interaction.LennardJones(
                    sigma=params['sig'],
                    epsilon=params['eps'],
                    shift='auto',
                    cutoff=cutoff))
            system.addInteraction(interaction)
            ret_list[pid] = interaction
    return ret_list


def setAngleInteractions(system, angles, angletypeparams, ftpl=None):
    ret_list={}

    for (aid, _), anglelist in angles.iteritems():
        if ftpl:
            fpl = espressopp.FixedTripleListAdress(system.storage, ftpl)
        else:
            fpl = espressopp.FixedTripleList(system.storage)
        fpl.addTriples(anglelist)
        angleinteraction=angletypeparams[aid].createEspressoInteraction(system, fpl)
        if angleinteraction:
            system.addInteraction(angleinteraction)
            ret_list.update({aid: angleinteraction})
    return ret_list

def setDihedralInteractions(system, dihedrals, dihedraltypeparams, ftpl=None):
    ret_list={}

    for (did, _), dihedrallist in dihedrals.iteritems():
        if ftpl:
            fpl = espressopp.FixedQuadrupleListAdress(system.storage, ftpl)
        else:
            fpl = espressopp.FixedQuadrupleList(system.storage)
        fpl.addQuadruples(dihedrallist)
        dihedralinteraction=dihedraltypeparams[did].createEspressoInteraction(system, fpl)
        if dihedralinteraction:
            system.addInteraction(dihedralinteraction)
            ret_list.update({did: dihedralinteraction})
    return ret_list

def setLennardJonesInteractions(system, defaults, atomtypeparams, verletlist, cutoff, nonbonded_params=None, hadress=False, ftpl=None, table_groups=None):
    """ Set lennard jones interactions which were read from gromacs based on the atomypes"""
    if table_groups is None:
        table_groups = []
    if ftpl:
        if hadress:
            interaction = espressopp.interaction.VerletListHadressLennardJones(verletlist, ftpl)
        else:
            interaction = espressopp.interaction.VerletListAdressLennardJones(verletlist, ftpl)
    else:
        interaction = espressopp.interaction.VerletListLennardJones(verletlist)

    if nonbonded_params is None:
        nonbonded_params = {}

    combinationrule = int(defaults['combinationrule'])
    print "Setting up Lennard-Jones interactions"

    type_pairs = sorted({
        tuple(sorted([type_1, type_2]))
            for type_1, pi in atomtypeparams.iteritems()
            for type_2, pj in atomtypeparams.iteritems()
            if ((pi['atnum'] not in table_groups and pj['atnum'] not in table_groups) and \
            (pi.get('atname') not in table_groups and pj.get('atname') not in table_groups))
        })
    print('Number of pairs: {}'.format(len(type_pairs)))
    for type_1, type_2 in type_pairs:
        pi = atomtypeparams[type_1]
        pj = atomtypeparams[type_2]
        if pi['particletype'] == 'V' or pj['particletype'] == 'V':
            print('Skip {}-{}'.format(type_1, type_2))
            continue
        param = nonbonded_params.get((type_1, type_2))
        if param:
            print 'Using defined non-bonded cross params', param
            sig, eps = param['sig'], param['eps']
        else:
            sig_1, eps_1 = float(pi['sig']), float(pi['eps'])
            sig_2, eps_2 = float(pj['sig']), float(pj['eps'])
            if combinationrule == 2:
                sig = 0.5*(sig_1 + sig_2)
                eps = (eps_1*eps_2)**(1.0/2.0)
            else:
                sig = (sig_1*sig_2)**(1.0/2.0)
                eps = (eps_1*eps_2)**(1.0/2.0)
        if sig > 0.0 and eps > 0.0:
            print "Setting LJ interaction for", type_1, type_2, "to sig ", sig, "eps", eps, "cutoff", cutoff
            ljpot = espressopp.interaction.LennardJones(epsilon=eps, sigma=sig, shift='auto', cutoff=cutoff)
            if ftpl:
                interaction.setPotentialAT(type1=type_1, type2=type_2, potential=ljpot)
            else:
                interaction.setPotential(type1=type_1, type2=type_2, potential=ljpot)

    system.addInteraction(interaction)
    return interaction


def setCoulombInteractions(system, verletlist, rc, atomtypeparams,
                           epsilon1, epsilon2, kappa, hadress=False, adress=False, ftpl=None,
                           pot=None, interaction=None
                           ):

    pref=138.935485 # we want gromacs units, so this is 1/(4 pi eps_0) ins units of kJ mol^-1 e^-2

    type_pairs = sorted({
        tuple(sorted([type_1, type_2]))
            for type_1, pi in atomtypeparams.iteritems()
            for type_2, pj in atomtypeparams.iteritems()
            if ((pi.get('charge', 0.0) != 0.0 and pj.get('charge', 0.0) != 0.0) and \
            (pi['particletype'] != 'V' and pj['particletype'] != 'V'))
        })

    print('Number of coulombic pairs: {}'.format(len(type_pairs)))
    if len(type_pairs) > 0:
        if pot is None:
            pot = espressopp.interaction.ReactionFieldGeneralized(
                prefactor=pref, kappa=kappa, epsilon1=epsilon1, epsilon2=epsilon2, cutoff=rc)
        if interaction is None:
            if hadress and adress:
                raise RuntimeError('Ambiguous option, it is only possible to use Adress or HAdress.')

            if hadress:
                interaction = espressopp.interaction.VerletListHadressReactionFieldGeneralized(
                    verletlist, ftpl)
            elif adress:
                interaction = espressopp.interaction.VerletListAdressReactionFieldGeneralized(
                    verletlist, ftpl)
            else:
                interaction = espressopp.interaction.VerletListReactionFieldGeneralized(verletlist)

        if hadress or adress:
            setPotential_fn = interaction.setPotentialAT
        else:
            setPotential_fn = interaction.setPotential

        for type_1, type_2 in type_pairs:
            print('Set coulomb interaction: {}-{}'.format(type_1, type_2))
            setPotential_fn(type1=type_1, type2=type_2, potential=pot)
        system.addInteraction(interaction)
        return interaction
    else:
        return None


def setCoulombInteractionsProtein(system, verletlist, rc, types, epsilon1, epsilonprot,epsilonwat,kappa,otype,htype, hadress=False, adress=False, ftpl=None):

    print "# Setting up Coulomb reaction field interactions"
    print "# Using ",epsilonwat," for water and wat-prot and ",epsilonprot," for protein"

    pref=138.935485 # we want gromacs units, so this is 1/(4 pi eps_0) ins units of kJ mol^-1 e^-2

    potwat = espressopp.interaction.ReactionFieldGeneralized(prefactor=pref, kappa=kappa, epsilon1=epsilon1, epsilon2=epsilonwat, cutoff=rc)
    potprot = espressopp.interaction.ReactionFieldGeneralized(prefactor=pref, kappa=kappa, epsilon1=epsilon1, epsilon2=epsilonprot, cutoff=rc)

    if hadress and adress:
      print "Error! In gromacs.setCoulombInteractions, you cannot use adress and hadress at the same time"
      return
    if hadress:
        interaction=espressopp.interaction.VerletListHadressReactionFieldGeneralized(verletlist, ftpl)
    elif adress:
        interaction=espressopp.interaction.VerletListAdressReactionFieldGeneralized(verletlist, ftpl)
    else:
	interaction=espressopp.interaction.VerletListReactionFieldGeneralized(verletlist)

    for i in range(max(types)+1):
        for k in range(i, max(types)+1):
            if i==otype or i==htype or k==otype or k==htype:
	        if hadress or adress:
	            interaction.setPotentialAT(type1=i, type2=k, potential=potwat)
	        else:
	            interaction.setPotential(type1=i, type2=k, potential=potwat)
            else:
	        if hadress or adress:
	            interaction.setPotentialAT(type1=i, type2=k, potential=potprot)
	        else:
	            interaction.setPotential(type1=i, type2=k, potential=potprot)

    system.addInteraction(interaction)
    return interaction


def setCoulomb14Interactions(system, defaults, onefourlist, rc, types):

    #in Gromas, 1-4 interactions don't have reaction field correction
    print "# Setting up 1-4 Coulomb interactions"

    if defaults:
        fudge=float(defaults['fudgeQQ'])
        print "# Using electrostatics 1-4 fudge factor ",fudge

    pref=138.935485*fudge # we want gromacs units, so this is 1/(4 pi eps_0) ins units of kJ mol^-1 e^-2, scaled by fudge factor

    #pot = espressopp.interaction.CoulombRSpace(prefactor=pref, alpha=0.0, cutoff=rc)
    pot = espressopp.interaction.CoulombTruncated(prefactor=pref, cutoff=rc)

    #interaction=espressopp.interaction.FixedPairListTypesCoulombRSpace(system,onefourlist)
    interaction=espressopp.interaction.FixedPairListTypesCoulombTruncated(system,onefourlist)

    for i in range(max(types)+1):
        for k in range(i, max(types)+1):
            interaction.setPotential(type1=i, type2=k, potential=pot)

    system.addInteraction(interaction)
    return interaction



def setTabulatedInteractions(system, atomtypeparams, vl, cutoff, interaction=None, ftpl=None, table_groups=None):
    """Sets tabulated potential for types that has particletype set to 'V'."""
    spline_type = 2
    if table_groups is None:
        table_groups = []

    type_pairs = {
            tuple(sorted([type_1, type_2]))
            for type_1, v1 in atomtypeparams.iteritems()
            for type_2, v2 in atomtypeparams.iteritems()
            if (v1['atnum'] in table_groups and v2['atnum'] in table_groups)
            }
    if len(type_pairs) > 0:
        if interaction is None:
            if ftpl:
                interaction = espressopp.interaction.VerletListAdressTabulated(vl, ftpl)
            else:
                interaction = espressopp.interaction.VerletListTabulated(vl)
        else:
            if not ftpl:
                interaction = espressopp.interaction.VerletListTabulated(vl)

        for type_1, type_2 in type_pairs:
            print('Set tabulated potential {}-{}'.format(type_1, type_2))
            name_1 = atomtypeparams[type_1]['atnum']
            name_2 = atomtypeparams[type_2]['atnum']
            name_1, name_2 = sorted([name_1, name_2])
            table_name = '{}-{}.tab'.format(name_1, name_2)
            orig_table_name = 'table_{}_{}.xvg'.format(name_1, name_2)
            print('Converting table_{name1}_{name2}.xvg to {name1}-{name2}.tab'.format(
                name1=name_1, name2=name_2))
            convertTable(orig_table_name, table_name)
            if ftpl:
                interaction.setPotentialCG(
                    type1=type_1,
                    type2=type_2,
                    potential=espressopp.interaction.Tabulated(
                        itype=spline_type,
                        filename=table_name,
                        cutoff=cutoff))
            else:
                interaction.setPotential(
                    type1=type_1,
                    type2=type_2,
                    potential=espressopp.interaction.Tabulated(
                        itype=spline_type,
                        filename=table_name,
                        cutoff=cutoff))
        if interaction and not ftpl:
            system.addInteraction(interaction)

        return interaction
    else:
        return None



