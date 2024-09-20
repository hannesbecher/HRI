#!/bin/bash

# Grid Engine options
#$ -cwd
#$ -t 1-8
#$ -l h_rt=24:00:00
#$ -l h_vmem=10G
#$ -N hriV06_8



rr=(0 79 82 88 184 186 187 188 193)

id=`printf %04d $SGE_TASK_ID`
echo "==========================================================="
echo Running task $id on $HOSTNAME
echo "==========================================================="

~/HighlanderLab/share/slim  -d "REPID='${rr[$id]}'" selAndNeut_allGTsV06.in
