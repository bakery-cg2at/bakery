import files_io
import sys

top1 = files_io.GROMACSTopologyFile(sys.argv[1])
top1.read()
top2 = files_io.GROMACSTopologyFile(sys.argv[2])
top2.read()

what_to_compare = ['atoms', 'dihedrals', 'angles', 'bonds',
                   'cross_bonds', 'cross_angles', 'cross_dihedrals',
                   'cross_pairs', 'pairs']

for s in what_to_compare:
    s1 = getattr(top1, s)
    s2 = getattr(top2, s)
    print((s, set(s1) == set(s2)))
