#!/bin/bash

# Grid Engine options
#$ -cwd

#$ -l h_rt=24:00:00
#$ -l h_vmem=10G
#$ -N hriV06_3rd





#id=`printf %04d $SGE_TASK_ID`
echo "==========================================================="
echo Running task 0193 on $HOSTNAME
echo "==========================================================="

~/HighlanderLab/share/slim  -d "REPID='0193'" selAndNeut_allGTsV06.in
