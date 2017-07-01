import os, glob
import subprocess
import filecmp
import argparse

def CleanBuildDir(path):
    subprocess.call(['make','-C',path,'clean'])

def BuildExec(path):
    subprocess.call(['make','-C',path,'-j','12'])

def BuildGetinit(rootPath):

    initPath=os.path.join(rootPath,'init')
    print "initPath = %s"%(initPath)

    print "\n\nbuilding getinit..\n"
    CleanBuildDir(initPath)
    BuildExec(initPath)

def Clean(cwdPath,rootPath,unitTests):

    srcPath=os.path.join(rootPath,'src')
    initPath=os.path.join(rootPath,'init')
    CleanBuildDir(srcPath)
    CleanBuildDir(initPath)

    for t in unitTests:
        print "================================="
        print t
        print "================================="

        unitTestPath = os.path.join(cwdPath,t)
        runDir = os.path.join(unitTestPath,'run')
        CleanBuildDir(runDir)

def RunUnitTests(cwdPath,unitTests):

    for t in unitTests:
        print "================================="
        print t
        print "================================="

        unitTestPath = os.path.join(cwdPath,t)
        runDir = os.path.join(unitTestPath,'run')
        print "  Entering %s"%(runDir)
        os.chdir(runDir)
    
        logFileName = t+'.log'
        print "  Running test"
        with open (logFileName, 'w') as output:
            subprocess.call(['make','run'], stdout=output, stderr=output)
    
        refDir = os.path.join(unitTestPath,'refs')
        refFiles = glob.glob(refDir+'/*.txt')
        print "  refDir = %s "%(refDir)
    
        for refPath in refFiles:
            ref = os.path.basename(refPath) 
            currentRefPath = os.path.join(runDir,ref)
            if(not filecmp.cmp(currentRefPath, refPath)):
                print "    ---------------------------------"
                print "    %s is different: "%(ref)
                print "    ---------------------------------"
                subprocess.call(['diff',currentRefPath,refPath])
            else:
                print "    %s Pass"%(ref)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('-r','--run', help="run all unit tests.", action="store_true", default=False)
    parser.add_argument('-c','--clean',help="clean up unit tests directories.", action="store_true", default=False)
    parser.add_argument('-g','--getinit',help="build getinit that is used in each test.", action="store_true", default=False)

    cwdPath=os.getcwd()
    rootPath=os.path.dirname(os.path.join('..',cwdPath))

    unitTests = filter(os.path.isdir, os.listdir(cwdPath))

    args = parser.parse_args()

    if args.run:
        RunUnitTests(cwdPath,unitTests)

    if args.getinit:
        BuildGetinit(rootPath)

    if args.clean:
        Clean(cwdPath,rootPath,unitTests)

