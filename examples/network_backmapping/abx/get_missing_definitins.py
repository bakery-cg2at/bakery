import sys

input_f = sys.argv[1]

missing_angles = set()
missing_dihedrals = set()

for l in open(input_f, 'r'):
    t = l.split()
    if len(t) == 6:
        missing_angles.add(tuple([x.split(':')[1] for x in t[:3]]))
    elif len(t) == 8:
        missing_dihedrals.add(tuple([x.split(':')[1] for x in t[:4]]))

print('Missing angles: ')
with open('missing_angles.txt', 'w') as oa:
    for a in missing_angles:
        print(' '.join(a))
        oa.write(' '.join(a))
        oa.write('\n')

print('Missing dihedrals: ')
with open('missing_dihedrals.txt', 'w') as od:
    for d in missing_dihedrals:
        print(' '.join(d))
        od.write(' '.join(d))
        od.write('\n')
