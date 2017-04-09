import subprocess
import os
import shutil

curDir=os.getcwd()

rxmdCom=[curDir+"/rxmd"]
geninitCom=[curDir+"/init/geninit",'-inputxyz','input.xyz','-ffield','ffield']
cleanCom=['make','-f',curDir+"/init/Makefile",'clean']

# reference tests 
refNames=['rdx','sic','FeS','rdx-exL','sic-exL']
extraArgs={key:[] for key in refNames}

# create 2x2x2 unit cells for FeS test
extraArgs['FeS']+=['-mc','2','2','2']

logFile="logfile.txt"
refFile="reffile.txt"

for ref in refNames:
    print '==================================='
    print ref
    print '==================================='
    refDir = curDir + "/refs/" + ref
    cwd = os.chdir(refDir)
    print 'changed dir to {}'.format(refDir)

# create a new system
    subprocess.call(geninitCom+extraArgs[ref])
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

