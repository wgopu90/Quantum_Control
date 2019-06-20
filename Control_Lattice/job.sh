      #PBS -N Quantum Control
      #PBS -S /bin/bash
      ## #PBS -e NapaValley:/home/gopi/Desktop/Matlab_codes/matlab/Control_of_lattice/benchmark.err
      ## #PBS -o NapaValley:/home/gopi/Desktop/Matlab_codes/matlab/Control_of_lattice/benchmark.log
      #PBS -l walltime=300:00:00
      #PBS -l nodes=1:ppn=14,mem=10gb
      #PBS -m abe -M wariergopu@gmail.com
      #PBS -V
      echo Working directory is $PBS_O_WORKDIR
      cd $PBS_O_WORKDIR
      
      ### Run some informational commands.
      echo Running on host `hostname`
      echo Time is `date`
      echo Directory is `pwd`
      echo This jobs runs on the following processors:
      echo `cat $PBS_NODEFILE`
      
      ### Define number of processors
      NPROCS=`wc -l < $PBS_NODEFILE`
      echo This job has allocated $NPROCS cpus
      
      ### Run Matlab script for benchmarking.
      matlab -nodisplay -nodesktop < trial_opt.m > outputfile_w.o
      
