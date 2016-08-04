#!/usr/bin/env python
"""
Copyright (C) 2015 Jakub Krajniak <jkrajniak@gmail.com>

This file is distributed under free software licence:
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Scan for .conv files and print convergence values for each 
of non-bonded potentials.

"""

import os

steps = [x for x in os.listdir('.') if x.startswith('step')]
steps.sort()
data = []

header = None
for step in steps:
    conv = [x for x in os.listdir(step) if x.endswith('dist.conv')]
    tmp = {}
    for k in conv:
        tmp[k.split('.')[0]] = float(open(os.path.join(step, k)).read())
    data.append(tmp)
    if not header and conv:
        header = [x.split('.')[0] for x in conv]

output = open('convergence.csv', 'w')

header_format = '{:^15}' + '{:<15}' * (len(header))
row_format = '{:^15}' + '{:<15.4}' * (len(header))

print(header_format.format('step', *header))

output.write(header_format.format('step', *header))
output.write('\n')

for step_id, conv in enumerate(data):
    print(row_format.format(step_id, *map(conv.get, header)))
    output.write('{}\n'.format(row_format.format(step_id, *map(conv.get, header))))

output.close()
print('Saved to convergence.csv ...')
