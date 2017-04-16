import subprocess
import os, sys
import shutil

def FileExist(filePath):
    if(os.path.exists(filePath)):
        return True
    else:
        print 'file {} does not exist'.format(filePath)
        sys.exit(1)

curDir=os.getcwd()

rxmdPath=os.path.join(curDir,"rxmd")
FileExist(rxmdPath)

geninitPath=os.path.join(curDir,"init/geninit")
FileExist(geninitPath)

# command to run 'rxmd' with run-profile data
rxmdCom=[rxmdPath,"--profile"]

# command for system initialization and cleanup
geninitCom=[geninitPath,'-inputxyz','input.xyz','-ffield','ffield']

# command to cleanup working dir
cleanCom=['make','-f',os.path.join(curDir,"init/Makefile"),'clean']

# reference tests 
refNames=['rdx','sic','FeS','rdx-exL','sic-exL']
extraArgs={key:[] for key in refNames}

# create 2x2x2 unit cells for FeS test
extraArgs['FeS']+=['-mc','2','2','2']

runProfileFile="profile.dat"
refFile="ref.dat"

for ref in refNames:
    print '==================================='
    print ref
    print '==================================='
    refDir = os.path.join(curDir,"refs",ref)
    cwd = os.chdir(refDir)
    #print 'changed dir to {}'.format(refDir)

# run short simulation, compare run profile of this run and previous reference value (profile.dat vs ref.dat)
# remove profile.dat at the end.  
    runProfilePath = os.path.join(refDir,runProfileFile)
    refPath = os.path.join(refDir,refFile)

# All stdout is stored in 'DAT/log'
    with open('DAT/log', 'w') as log:

# create new system
        subprocess.call(geninitCom+extraArgs[ref], stdout=log)
        shutil.move('./rxff.bin','DAT/rxff.bin')

        p = subprocess.call(rxmdCom, stdout=log)
        
        proc = subprocess.Popen(['diff',refPath,runProfilePath], stdout=subprocess.PIPE)
        stdOut = proc.stdout.read()

        if(len(stdOut)==0):
            print "Pass!"
        else:
            print "Failed!\n",stdOut

        os.remove(runProfilePath)

# cleanup 
        subprocess.call(cleanCom, stdout=log)
