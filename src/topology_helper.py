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


# Some helper classes usefull when parsing the gromacs topology

import espressopp
import math
import os


def convertTable(gro_in_file, esp_out_file, sigma=1.0, epsilon=1.0, c6=1.0, c12=1.0):
    """Convert GROMACS tabulated file into ESPResSo++ tabulated file (new file
    is created). First column of input file can be either distance or angle.
    For non-bonded files, c6 and c12 can be provided. Default value for sigma, epsilon,
    c6 and c12 is 1.0. Electrostatics are not taken into account (f and fd columns).

    Args:
        gro_in_file: the GROMACS tabulated file name (bonded, nonbonded, angle
            or dihedral).
        esp_out_file: filename of the ESPResSo++ tabulated file to be written.
        sigma: optional, depending on whether you want to convert units or not.
        epsilon: optional, depending on whether you want to convert units or not.
        c6: optional
        c12: optional
    """
    # determine file type
    bonded, angle, dihedral = False, False, False
    if gro_in_file[6] == "b":
        bonded = True
    if gro_in_file[6] == "a":
        angle = True
        bonded = True
    if gro_in_file[6] == "d":
        dihedral = True
        bonded = True

    fin = open(gro_in_file, 'r')
    fout = open(esp_out_file, 'w')

    if bonded:  # bonded has 3 columns
        for line in fin:
            if line.startswith('#'):
                continue

            columns = line.split()
            r = float(columns[0])
            f = float(columns[1])  # energy
            fd = float(columns[2])  # force

            # convert units
            if angle or dihedral:  # degrees to radians
                r = math.radians(r)
                fd = fd*180/math.pi
            else:
                r = r / sigma
            e = f / epsilon
            f = fd*sigma / epsilon

            if (not angle and not dihedral and r != 0) or \
                (angle and r <= math.pi and r >= 0) or \
                    (dihedral and r >= -math.pi and r <= math.pi):
                fout.write("%15.8g %15.8g %15.8g\n" % (r, e, f))

    else:  # non-bonded has 7 columns
        for line in fin:
            if line.startswith('#'):  # skip comment lines
                continue
            # TODO: Skiped columns 1 and 2, electrostatics is not implemented yet.
            columns = line.split()
            r = float(columns[0])
            g = float(columns[3])  # dispersion
            gd = float(columns[4])
            h = float(columns[5])  # repulsion
            hd = float(columns[6])

            e = c6*g + c12*h
            f = c6*gd + c12*hd

            # convert units
            r = r / sigma
            e = e / epsilon
            f = f*sigma / epsilon

            if r != 0:  # skip 0
                fout.write("%15.8g %15.8g %15.8g\n" % (r, e, f))

    fin.close()
    fout.close()


class FileBuffer():
    def __init__(self):
        self.linecount = 0
        self.lines = []
        self.pos = 0

    def appendline(self, line):
        self.lines.append(line)
        self.linecount = len(self.lines)

    def readline(self):
        try:
            line = self.lines[self.pos]
        except:
            return ''
        self.pos += 1
        return line

    def readlastline(self):
        try:
            line = self.lines[self.pos-1]
        except:
            return ''
        return line

    def seek(self, p):
        self.pos = p

    def tell(self):
        return self.pos

    def eof(self):
        return self.pos == len(self.lines)


def FillFileBuffer(fname, filebuffer, cwd=None, defines=None):
    if cwd is None:
        cwd = '.'
    if defines is None:
        defines = {}
    f = open(os.path.join(cwd, fname), 'r')
    for line in f:
        if line.startswith(';'):
            continue
        if "include" in line:
            name = line.split()[1].strip('\"')
            cwd_name = os.path.dirname(name)
            if cwd_name != '':
                cwd = cwd_name
            FillFileBuffer(name, filebuffer, cwd, defines)
        elif 'define' in line:
            t = line.strip().split()
            if len(t) > 2:
                defines[t[1]] = ' '.join(t[2:])
        else:
            l = line.rstrip('\n')
            if l:
                filebuffer.appendline(l)

    f.close()
    return


def PostProcessFileBuffer(filebuffer, defines):
    """Replace all defines with the value from the dictionary."""
    ret_fb = FileBuffer()
    define_keys = set(defines)
    for line in filebuffer.lines:
        line = line.strip()
        if line:
            if not (line.startswith(';') or line.startswith('#define')
                    or line.startswith('#include') or line.startswith('#ifdef')
                    or line.startswith('#ifndef')):
                def_key = set.intersection(set(map(str.strip, line.split())), define_keys)
                if def_key:
                    def_key = def_key.pop()
                    ret_fb.appendline(
                        line.replace(def_key, defines[def_key]))
                else:
                    ret_fb.appendline(line)
            else:
                ret_fb.appendline(line)
    return ret_fb


def FindType(proposedtype, typelist):
    list=[typeid for (typeid,atype) in typelist.iteritems() if atype==proposedtype ]
    if len(list)>1:
        print "Error: duplicate type definitons", proposedtype.parameters
        exit()
    elif len(list)==0:
        return None
    return list[0]


def convertc6c12(c6, c12):
    if c12 == 0.0:
        return 1.0, 0.0
    sig = pow(c12/c6, 1.0/6.)
    if sig > 0.0:
        eps = 0.25*c6*pow(sig, -6.0)
    else:
        eps = 0.0
    return sig, eps


class InteractionType:
    def __init__(self, parameters):
        self.parameters=parameters
    def __eq__(self,other):
        # interaction types are defined to be equal if all parameters are equal
        for k, v in self.parameters.iteritems():
            if k not in other.parameters: return False
            if other.parameters[k]!=v: return False
        return True
    def createEspressoInteraction(self, system, fpl, is_cg=None):
        print "WARNING: could not set up interaction for", self.parameters, ": Espresso potential not implemented"
        return None
    def automaticExclusion(self):
        #overwrite in derrived class if the particular interaction is automatically excluded
        return False

class HarmonicBondedInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl, is_cg=None):
        # interaction specific stuff here
        # spring constant kb is half the gromacs spring constant
        pot = espressopp.interaction.Harmonic(self.parameters['kb']/2.0, self.parameters['b0'])
        print 'setting harmonic bond k=', self.parameters['kb']/2.0, 'b0=', self.parameters['b0']
        if is_cg is not None:
            return espressopp.interaction.FixedPairListAdressHarmonic(system, fpl, pot, is_cg)
        else:
            return espressopp.interaction.FixedPairListHarmonic(system, fpl, pot)

    def automaticExclusion(self):
        return True

class MorseBondedInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl, is_cg=None):
        # interaction specific stuff here
        if is_cg is not None:
            raise RuntimeError('Morse potential is not implemented to support Adress')
        pot = espressopp.interaction.Morse(self.parameters['D'], self.parameters['beta'], self.parameters['rmin'])
        interb = espressopp.interaction.FixedPairListMorse(system, fpl, pot)
        return interb
    def automaticExclusion(self):
        return True

class FENEBondedInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl, is_cg=None):
        # interaction specific stuff here
        # spring constant kb is half the gromacs spring constant
        if is_cg is not None:
            raise RuntimeError('FENE potential is not implemented to support Adress')
        pot = espressopp.interaction.Fene(self.parameters['kb']/2.0, self.parameters['b0'])
        interb = espressopp.interaction.FixedPairListFene(system, fpl, pot)
        return interb
    def automaticExclusion(self):
        return True

class HarmonicAngleInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl, is_cg=None):
        # interaction specific stuff here
        # spring constant kb is half the gromacs spring constant. Also convert deg to rad
        K = self.parameters['k'] / 2.0
        theta0 = self.parameters['theta']*math.pi/180.0
        print 'setting angular harmonic k=', K, 'theta=', theta0
        pot = espressopp.interaction.AngularHarmonic(K=K, theta0=theta0)
        if is_cg is not None:
            return espressopp.interaction.FixedTripleListAdressAngularHarmonic(system, fpl, pot, is_cg)
        else:
            return espressopp.interaction.FixedTripleListAngularHarmonic(system, fpl, pot)


class TabulatedBondInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl, is_cg=None):
        spline=1
        fg = "table_b"+str(self.parameters['tablenr'])+".xvg"
        fe = fg.split(".")[0]+".pot" # name of espressopp file
        if not os.path.exists(fe):
            convertTable(fg, fe)
        print('Tabulated bond: {}'.format(fe))
        potTab = espressopp.interaction.Tabulated(itype=spline, filename=fe)
        if is_cg is not None:
            return espressopp.interaction.FixedPairListAdressTabulated(system, fpl, potTab, is_cg)
        else:
            return espressopp.interaction.FixedPairListTabulated(system, fpl, potTab)
    def automaticExclusion(self):
        return self.parameters['excl']

class TabulatedAngleInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl, is_cg=None):
        spline=1
        fg = "table_a"+str(self.parameters['tablenr'])+".xvg"
        fe = fg.split(".")[0]+".pot" # name of espressopp file
        if not os.path.exists(fe):
            convertTable(fg, fe)
        print('Tabulated angular: {}'.format(fe))
        potTab = espressopp.interaction.TabulatedAngular(itype=spline, filename=fe)
        if is_cg is not None:
            return espressopp.interaction.FixedTripleListAdressTabulatedAngular(system, fpl, potTab, is_cg)
        else:
            return espressopp.interaction.FixedTripleListTabulatedAngular(system, fpl, potTab)

class TabulatedDihedralInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl, is_cg=None):
        spline = 1
        fg = "table_d"+str(self.parameters['tablenr'])+".xvg"
        fe = fg.split(".")[0]+".pot" # name of espressopp file
        if not os.path.exists(fe):
            convertTable(fg, fe)
        print('Tabulated dihedral: {}'.format(fe))
        potTab = espressopp.interaction.TabulatedDihedral(itype=spline, filename=fe)
        if is_cg is not None:
            return espressopp.interaction.FixedQuadrupleListAdressTabulatedDihedral(system, fpl, potTab, is_cg)
        else:
            return espressopp.interaction.FixedQuadrupleListTabulatedDihedral(system, fpl, potTab)

class RyckaertBellemansDihedralInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl, is_cg=None):
        print('RyckaertBellemans: {}'.format(self.parameters))
        pot = espressopp.interaction.DihedralRB(**self.parameters)
        if is_cg is not None:
            return espressopp.interaction.FixedQuadrupleListAdressDihedralRB(system, fpl, pot, is_cg)
        else:
            return espressopp.interaction.FixedQuadrupleListDihedralRB(system, fpl, pot)


class HarmonicNCosDihedralInteractionType(InteractionType):
    def createEspressoInteraction(self, system, fpl, is_cg=None):
        theta0 = self.parameters['theta0'] * math.pi/180.0
        print('HarmonicNCosDihedral, theta0: {}, k:{}, n: {}'.format(
            theta0, self.parameters['k'], self.parameters['n']))
        pot = espressopp.interaction.DihedralHarmonicNCos(
            K=self.parameters['k'], phi0=theta0, multiplicity=self.parameters['n'])
        if is_cg is not None:
            return espressopp.interaction.FixedQuadrupleListAdressDihedralHarmonicNCos(system, fpl, pot, is_cg)
        else:
            return espressopp.interaction.FixedQuadrupleListDihedralHarmonicNCos(system, fpl, pot)


def ParseBondTypeParam(line):
    tmp = line.split()
    btype= tmp[2]
    # TODO: handle exclusions automatically
    if btype == "8":
        p=TabulatedBondInteractionType({"tablenr":int(float(tmp[3])),"k":float(tmp[4]), 'excl':True})
    elif btype == "9":
        p=TabulatedBondInteractionType({"tablenr":int(float(tmp[3])), "k":float(tmp[4]), 'excl':False})
    elif btype == "1":
        p=HarmonicBondedInteractionType({"b0":float(tmp[3]), "kb":float(tmp[4])})
    elif btype == "2":
        # Converts from GROMOS Harmonic to Harmonic
        # k_harm = 2*kb*b0**2
        kb = 0.5*2*float(tmp[4])*(float(tmp[3])**2)
        p=HarmonicBondedInteractionType({"b0":float(tmp[3]), "kb":kb})
    elif btype == "3":
        p=MorseBondedInteractionType({"b0":float(tmp[3]), "D":float(tmp[4]), "beta":float(tmp[5])})
    elif btype == "7":
        p=FENEBondedInteractionType({"b0":float(tmp[3]), "kb":float(tmp[4])})
    elif btype == "9":
        p=TabulatedBondInteractionType({"tablenr":int(float(tmp[3])), "k":float(tmp[4])})
    else:
        print "Unsupported bond type", tmp[2], "in line:"
        print line
        exit()
    return p

def ParseAngleTypeParam(line):
    tmp = line.split()
    type= float(tmp[3])
    if type == 1:
        p=HarmonicAngleInteractionType({"theta":float(tmp[4]), "k":float(tmp[5])})
    elif type == 2:
        theta = float(tmp[4])*math.pi/180.0
        k = 0.5*float(tmp[5])*(math.sin(theta)**2)
        p=HarmonicAngleInteractionType({"theta":float(tmp[4]), "k":k})
    elif type == 8:
        p=TabulatedAngleInteractionType({"tablenr":int(float(tmp[4])),"k":float(tmp[5])})
    else:
        print "Unsupported angle type", type, "in line:"
        print line
        exit()
    return p

def ParseDihedralTypeParam(line):
    #
    tmp = line.split(';')[0].split()
    type = int(tmp[4])
    if type == 8:
        p=TabulatedDihedralInteractionType({"tablenr":int(float(tmp[5])), "k":float(tmp[6])})
    elif type == 3:
        tmp[5:11] = map(float, tmp[5:11])
        p = RyckaertBellemansDihedralInteractionType({'K0': tmp[5], 'K1': tmp[6], 'K2': tmp[7], 'K3': tmp[8], 'K4': tmp[9], 'K5': tmp[10]})
    elif type == 1:
        p = HarmonicNCosDihedralInteractionType({'theta0': float(tmp[5]), 'k': float(tmp[6]), 'n': int(tmp[7])})
    else:
        print "Unsupported dihedral type", type, "in line:"
        print line
        return False
    return p



# Usefull code for generating the regular exclusions

class Node():
    def __init__(self, id):
	self.id=id
	self.neighbours=[]
    def addNeighbour(self, nb):
	self.neighbours.append(nb)

def FindNodeById(id, nodes):
    list=[n for n in nodes if n.id==id ]
    if len(list)>1:
        print "Error: duplicate nodes", id
        exit()
    elif len(list)==0:
        return None
    return list[0]

def FindNNextNeighbours(startnode, numberNeighbours, neighbours, forbiddenNodes):
    if numberNeighbours==0:
	return neighbours

    #avoid going back the same path
    forbiddenNodes.append(startnode)

    # Loop over next neighbours and add them to the neighbours list
    # Recursively call the function with numberNeighbours-1
    for n in startnode.neighbours:
	if not n in forbiddenNodes:
	    if n not in neighbours: neighbours.append(n) # avoid double counting in rings
	    FindNNextNeighbours(n, numberNeighbours-1, neighbours, forbiddenNodes)


def GenerateRegularExclusions(bonds, nrexcl, exclusions):
    nodes=[]
    # make a Node object for each atom involved in bonds
    for b in bonds:
        bids=b[0:2]
        for i in bids:
            if FindNodeById(i, nodes)==None:
               n=Node(i)
               nodes.append(n)

    # find the next neighbours for each node and append them
    for b in bonds:
        permutations=[(b[0], b[1]), (b[1], b[0])]
        for p in permutations:
            n=FindNodeById(p[0], nodes)
            nn=FindNodeById(p[1], nodes)
            n.addNeighbour(nn)

    # for each atom, call the FindNNextNeighbours function, which recursively
    # seraches for nrexcl next neighbours
    for n in nodes:
        neighbours=[]
        FindNNextNeighbours(n, nrexcl, neighbours, forbiddenNodes=[])
        for nb in neighbours:
            # check if the permutation is already in the exclusion list
            # this may be slow, but to do it in every MD step is even slower...
            # TODO: find a clever algorithm which does avoid permuations from the start
            if not (n.id, nb.id) in exclusions:
                if not (nb.id, n.id) in exclusions:
                    exclusions.append((n.id, nb.id))

    return exclusions
