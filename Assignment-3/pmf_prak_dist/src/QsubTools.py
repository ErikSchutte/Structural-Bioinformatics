"""
Sample SARA/LISA job scripts:

##############################################################
simple job, requiring one node:
##############################################################
#PBS -lwalltime=2:00:00
                         # 2 hours wall-clock
                         # time allowed for this job
#PBS -lnodes=1:ppn=1
                         # 1 node for this job
#PBS -S /bin/bash
n=`wc -l < $PBS_NODEFILE`
echo start of job in directory $PBS_O_WORKDIR
echo number of nodes is $n
echo the allocated nodes are:
cat $PBS_NODEFILE
##############################################################


##############################################################
Example job to run gromacs on two nodes, 8 processes per node:
##############################################################

#PBS -lnodes=2:ppn=8 -lwalltime=1:00:00
echo start of gromacs job
module load gromacs openmpi fortran
cd my-gromacs-workdir             # put here your actual work directory
mpiexec mdrun -------             # put after mdrun the relevant parameters
echo end of gromacs job

"""

import os
import sys
import datetime as dt

import Env as env
import OsTools as ot

def generate_jobscript(commandlist,argumentlist,jobname,ppn,walltime,ClusterName):
    if(ClusterName=='lisa'):
        header    = '################################\n'
        header   += '# Job: '+jobname+'\n'
        header   += '# Generated automatically by: '+os.path.abspath(__file__)+'\n'
        header   += '# User: '+os.getlogin()+'\n'
        now = dt.datetime.now()
        header   += '# Date: '+str(now.date())+'\n'
        header   += '# Time: '+str(now.time())+'\n'
        header   += '################################'
        walltime  = '#PBS -lwalltime='+walltime
        nodes_ppn = '#PBS -lnodes=1:ppn='+str(ppn)
        gmxlib    = 'export GMXLIB=/home/repool/gromacs/top'
        modules   = 'module load gromacs openmpi fortran'
        cdpwd     = 'cd '+os.getcwd()
    #   write to file
        fw = open(jobname,'w')
        fw.write(header)
        fw.write('\n\n')
        fw.write(walltime)
        fw.write('\n\n')
        fw.write(nodes_ppn)
        fw.write('\n\n')
        fw.write(gmxlib)
        fw.write('\n\n')
        fw.write(modules)
        fw.write('\n\n')
        fw.write(cdpwd)
        fw.write('\n\n')

        if (len(commandlist)==len(argumentlist)):
            for i in range(len(commandlist)):
                command   = commandlist[i]+' \\\n'
                arguments = argumentlist[i]
                for j in range(len(arguments)-1):
                    line         = '\t'
                    line        += str(arguments[j])
                    line        += ' \\\n'
                    arguments[j] = line
                line         = '\t'
                line        += str(arguments[-1])
                arguments[-1] = line
    #           write to file
                fw.write(command)
                for line in arguments:
                    fw.write(line)
                fw.write('\\\n')
                outfile = jobname+'.'+str(i)+'.out'
                fw.write('\t> '+outfile+' 2>&1 &')
                fw.write('\n\n')
            fw.write('\n')
            fw.write('wait\n')
        fw.close()
    elif(ClusterName=='ibis'):
        content  = '#!/bin/bash\n'
        content += '#\n'
        content += '# ###########################\n'
        content += '# Job: '+jobname+'.sh\n'
        content += '# Generated automatically by: '+os.path.abspath(__file__)+'\n'
        content += '# User: '+os.getlogin()+'\n'
        now      = dt.datetime.now()
        content += '# Date: '+str(now.date())+'\n'
        content += '# Time: '+str(now.time())+'\n'
        content += '# ###########################\n'
        content += '#\n'
        content += '# ---------------------------\n'
        content += '# set the name of the job\n'
        content += '#$ -N '+jobname+'.sh\n'
        content += '# ---------------------------\n'
        content += '#\n'
        content += '# ---------------------------\n'
        content += '# set up the parameters for qsub\n'
        content += '# ---------------------------\n'
#        content += '#'
#        content += '#  Mail to user at beginning/end/abort/on suspension'
#        content += '#$ -m beas'
#        content += '#  By default, mail is sent to the submitting user'
#        content += '#  Use  $ -M username    to direct mail to another userid'
        content += '#\n'
        content += '# Execute the job from the current working directory\n'
        content += '# Job output will appear in this directory\n'
        content += '#$ -cwd\n'
        content += '#   can use -o dirname to redirect stdout\n'
        content += '#   can use -e dirname to redirect stderr\n'
        content += '#\n'
        content += '# to request resources at job submission time\n'
        content += '# use #-l resource=value\n'
        content += '# For instance, the commented out\n'
        content += '# lines below request a resource of \'express\'\n'
        content += '# and a hard CPU time of 10 minutes\n'
        content += '####$ -l express\n'
        content += '####$ =l h_cpu='+walltime+'\n'
        content += '#\n'
        content += '#  Export these environment variables\n'
        content += '#$ -v PATH\n'
        content += '#\n'
        content += 'export PATH=$TMPDIR:$PATH\n'
        content += 'export GMXLIB='+env.home+'/packages/gromacs-4.0.5/share/gromacs/top\n'
        content += '# MARTINI FORCE FIELD\n'
        content += 'export MFFHOME=\"$HOME/MARTINI/MFF\"\n'
        content += '# PYTHONPATH\n'
        content += '# BioPython:\n'
        content += 'export PYTHONPATH=\"'+env.home+'/packages/python/biopython-1.53/lib/python2.5/site-packages\"\n'
        content += '# modeller:\n'
        content += 'export PYTHONPATH=$PYTHONPATH:\"'+env.home+'/bin/modeller9v7/modlib\"\n'
        content += 'export PYTHONPATH=$PYTHONPATH:\"'+env.home+'/bin/modeller9v7/lib/x86_64-intel8\"\n'
        content += 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:\"'+env.home+'/bin/modeller9v7/lib/x86_64-intel8\"\n'
        content += '# Pull module:\n'
        content += 'export PYTHONPATH=$PYTHONPATH:\"'+env.home+'/PyWork/trunk/Pull/src\"\n'
        content += '#\n'
        content += '# ---------------------------\n'
        content += '# run the job\n'
        content += '# ---------------------------\n'
        content += '#\n'
#       write content to file
        fw = open(jobname+'.sh','w')
        fw.write(content)
#       write actual commands and arguments
        if (len(commandlist)==len(argumentlist)):
            for i in range(len(commandlist)):
                command   = commandlist[i]+' \\\n'
                arguments = argumentlist[i]
                for j in range(len(arguments)-1):
                    line         = '\t'
                    line        += str(arguments[j])
                    line        += ' \\\n'
                    arguments[j] = line
                line         = '\t'
                line        += str(arguments[-1])
                arguments[-1] = line
#               write to file
                fw.write(command)
                for line in arguments:
                    fw.write(line)
                fw.write('\\\n')
                outfile = jobname+'.sh.'+str(i)+'.out'
                fw.write('\t> '+outfile+' 2>&1')
                fw.write('\n\n')
            fw.write('\n')
        fw.close()
    else:
        print 'ERROR (generate_jobscript): len(commandlist)!=len(argumentlist)'
        sys.exit(1)


def qsub(jobname):
    if ot.system('qsub', options=jobname):
        print "qsub failed for", jobname
        exit(-1)
