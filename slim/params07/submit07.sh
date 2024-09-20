#!/bin/bash

# Grid Engine options
#$ -cwd
#$ -l h_rt=6:00:00
#$ -t 1-200
#$ -l h_vmem=10G
#$ -N hriV07_200


id=`printf %04d $SGE_TASK_ID`
echo "==========================================================="
echo Running task $id on $HOSTNAME
echo "==========================================================="

~/HighlanderLab/share/slim  -d "REPID='$id'" selAndNeut_allGTsV07.in
