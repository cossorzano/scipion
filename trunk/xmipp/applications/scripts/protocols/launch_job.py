# The job should be launched from the working directory!
def launch_job(programname,
               params,
               log,
               DoParallel,
               NumberOfMpiProcesses,
               NumberOfThreads,
               SystemFlavour):

    import os,sys

    if not DoParallel:
        command = programname + ' ' + params

    else:
        mpiprogramname=programname.replace('xmipp','xmipp_mpi')

        if (SystemFlavour=='SLURM-BSC'):
            mpicommand = 'srun '
        elif (SystemFlavour=='PBS-FT'):
            mpicommand = 'mpirun -np ' + str(NumberOfMpiProcesses) + ' -machinefile ' + os.environ.get('PBS_NODEFILE')
        elif (SystemFlavour=='XMIPP_MACHINEFILE'):
            mpicommand = 'mpirun -np ' + str(NumberOfMpiProcesses) + ' -machinefile ' + os.environ.get('XMIPP_MACHINEFILE')
        elif (SystemFlavour==''):
            mpicommand = 'mpirun -np ' + str(NumberOfMpiProcesses)
        else:
            message= "Error: unrecognized SystemFlavour: ", SystemFlavour
            print '* ',message
            print '*********************************************************************'
            log.info(message)
            sys.exit(1)

        command = mpicommand + ' `which '+ str(mpiprogramname) +'` ' + params

    print '* ',command,'\n'
    log.info(command)
    os.system(command)
