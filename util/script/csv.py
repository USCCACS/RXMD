import sys

filename = sys.argv[1]

data = []
with open(filename,'r') as f:
    for line in f:
        print(line)
        data.append(line.split())

with open(filename+'.csv','w') as f:
    for d in data:
        line = ''
        for d0 in d:
            line = line + (d0+',') 
        f.write(line+'\n')
