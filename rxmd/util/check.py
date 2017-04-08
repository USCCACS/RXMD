import subprocess
import os
import difflib

curDir=os.getcwd()
rxmdExec=curDir+"/"+"rxmd"

refNames=["rdx","sic"]
logFile="logfile.txt"
refFile="reffile.txt"

def getInitData(workDir):
    os.chdir("init")


for d in refNames:
    refDir = curDir + "/refs/" + d
    cwd = os.chdir(refDir)
    print 'changed dir to {}'.format(refDir)

    logPath = refDir + "/" + logFile
    refPath = refDir + "/" + refFile

    with open(logPath, 'w') as log:
        p = subprocess.Popen(rxmdExec, shell=True, universal_newlines=True, stdout=log)
        p.wait()
    subprocess.call(['diff',refPath,logPath])

