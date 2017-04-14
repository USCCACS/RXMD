import subprocess
import os
import shutil

curDir=os.getcwd()

# command to run 'rxmd' with summary data
rxmdCom=[curDir+"/rxmd","--summary"]

# command for system initialization and cleanup
geninitCom=[curDir+"/init/geninit",'-inputxyz','input.xyz','-ffield','ffield']
cleanCom=['make','-f',curDir+"/init/Makefile",'clean']

# reference tests 
refNames=['rdx','sic','FeS','rdx-exL','sic-exL']
extraArgs={key:[] for key in refNames}

# create 2x2x2 unit cells for FeS test
extraArgs['FeS']+=['-mc','2','2','2']

summaryFile="summary.dat"
refFile="ref.dat"

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

# run a short run, compare the output against reference output from previous run, 
# then cleanup temporary summary file. 
    summaryPath = refDir + "/" + summaryFile
    refPath = refDir + "/" + refFile

    with open('DAT/log', 'w') as log:
        p = subprocess.call(rxmdCom, stdout=log)
        
        proc = subprocess.Popen(['diff',refPath,summaryPath], stdout=subprocess.PIPE)
        stdOut = proc.stdout.read()

        print '-----------------------------------'
        if(len(stdOut)==0):
            print "Pass!"
        else:
            print "Error!\n",stdOut
        print '-----------------------------------'

        os.remove(summaryPath)

# cleanup 
    subprocess.call(cleanCom)
