#!/bin/bash

# Grid Engine options
#$ -cwd
#$ -t 1-200
#$ -l h_rt=12:00:00
#$ -l h_vmem=10G
#$ -N hriV02_200


id=`printf %04d $SGE_TASK_ID`
echo "==========================================================="
echo Running task $id on $HOSTNAME
echo "==========================================================="

~/HighlanderLab/share/slim  -d "REPID='$id'" selAndNeut_allGTs.in
