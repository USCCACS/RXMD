import os, sys

def get_mdframe(filename):

    mdframe = dict()
    mdframe['FILENAME'] = filename
    with open(filename, 'r') as fin:
       for line in fin:
           data = line.split()

           if 'total energy' in line:
               mdframe['ENERGY'] = float(data[-1])

           if 'PRIMVEC' in line:
               lattice = list()
               lattice.append(float(fin.readline().split()[0]))
               lattice.append(float(fin.readline().split()[1]))
               lattice.append(float(fin.readline().split()[2]))
               mdframe['LATTICE'] = lattice

           if 'PRESSURE' in line:
               mdframe['PRESSURE'] = [float(p) for p in fin.readline().split()]

           if 'PRIMCOORD' in line:
               num_atoms = int(fin.readline().split()[0])

               elems = list()
               coordinates = list()
               forces = list()
               for n in range(num_atoms):
                   data = fin.readline().split()
                   elems.append(data[0])
                   coordinates.append([float(f) for f in data[1:4]])
                   forces.append([float(f) for f in data[4:7]])
               mdframe['ELEMENTS'] = elems
               mdframe['COORDINATES'] = coordinates
               mdframe['FORCES'] = forces

    return mdframe

def save_mdframe(mdframe, filename=""):

    if filename == "": filename = mdframe['FILENAME']+'.xyz'

    with open(filename,'w') as fout:
        fout.write(f'{len(mdframe["ELEMENTS"])}\n')
        string = [str(f) for f in mdframe['LATTICE']]
        if len(string) == 3: string+= ["90"]*3
        fout.write(' '.join(string) + '\n')

        elems = mdframe['ELEMENTS']
        coords = mdframe['COORDINATES']
        forces = mdframe['FORCES']

        for n in range(len(elems)):
            string = [elems[n]]
            string += [str(c) for c in coords[n]]
            string += [str(c) for c in forces[n]]
            fout.write(' '.join(string) + '\n')

filename = sys.argv[1]

mdframe = get_mdframe(filename)

for k,v in mdframe.items():
    print(k,': ',v)

save_mdframe(mdframe)
