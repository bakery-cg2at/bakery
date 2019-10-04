from md_libs import files_io
import sys

top_f = sys.argv[1]
missing_def = sys.argv[2]

top = files_io.GROMACSTopologyFile(top_f)
top.read()

angles_params = {tuple([top.atoms[x].name for x in k]): v for k, v in top.angles.items()}
dihedrals_params = {tuple([top.atoms[x].name for x in k]): v for k, v in top.dihedrals.items()}

for m in open(missing_def, 'r'):
    t = m.split()
    if len(t) == 3:
        param = angles_params.get(tuple(t), angles_params.get(tuple(reversed(t))))
        print('{}: {}'.format(t,param))
    elif len(t) == 4:
        param = dihedrals_params.get(tuple(t), dihedrals_params.get(tuple(reversed(t))))
        print('{}: {}'.format(t,param))
