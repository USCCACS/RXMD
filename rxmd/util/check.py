import subprocess
import os
import difflib
import shutil

curDir=os.getcwd()

rxmdCom=[curDir+"/rxmd"]
geninitCom=[curDir+"/init/geninit",'input.xyz','ffield']
cleanCom=['make','-f',curDir+"/init/Makefile",'clean']

refNames=["rdx","sic"]
logFile="logfile.txt"
refFile="reffile.txt"

for d in refNames:
    print '==================================='
    refDir = curDir + "/refs/" + d
    cwd = os.chdir(refDir)
    print 'changed dir to {}'.format(refDir)

# create a new system
    subprocess.call(geninitCom)
    shutil.move('./rxff.bin','DAT/rxff.bin')

# run a short run and compare the output against reference output from previous run
    logPath = refDir + "/" + logFile
    refPath = refDir + "/" + refFile

    with open(logPath, 'w') as log:
        p = subprocess.Popen(rxmdCom, shell=True, universal_newlines=True, stdout=log)
        p.wait()
    subprocess.call(['diff',refPath,logPath])

# cleanup 
    subprocess.call(cleanCom)

