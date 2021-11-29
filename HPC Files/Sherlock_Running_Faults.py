#!/usr/bin/env python

import os

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

job_directory = os.getcwd()

# Make top level directories
mkdir_p(job_directory)

#indexes = [2,4,5,6,10,11,12,13,14,15,19,20,21,22,33,34,40,49,50,51,52,54,55,56,57,58,60,61,64,65,66,67,72,73,74,75,76,78,81,85,87,88,96,97,99,101,102,105,106,107,108,109,110,114,115,116,118,119,126,127,192,219,224,229,232,233,234,235,237,241,254,259,273,274,325,632,636,638,639,640,648,664,681,687,688,689,692,693,696,698,701,704,706,710,715,859,860,864,865,867,1001,	1090,1099,1108,1175,1176,1179,1184,1185,1187,1256,1364,1387,1445,1450,1464,1473,1478,1487,1488,1489,1543,1553,1557,1558,1559,1563]

indexes = range(0,1992)

i=1
for set in range(91,101):

    first = 20*(set-1)
    last = 20*set+1
    if last > len(indexes):
        continue
    extra_time = (i-1)*600
    runs = indexes[first:last]
    i +=1
    for ind_run in runs:

        job_file = os.path.join(job_directory, "%s.job" % str(ind_run))
        with open(job_file,'wb') as fh:
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name=%s.2hrs\n" % str(ind_run))
            fh.writelines("#SBATCH --output=ModuleT.log\n")
            fh.writelines("#SBATCH --time=24:00:00\n")
            fh.writelines("#SBATCH --mem=32000\n")
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --ntasks-per-node=1\n")
            fh.writelines("#SBATCH --cpus-per-task=16\n")
            fh.writelines("#SBATCH --qos=normal\n")
            fh.writelines("#SBATCH --mail-type=FAIL\n")
            fh.writelines("#SBATCH --mail-user=rsilval@stanford.edu\n")
            fh.writelines("#SBATCH --partition=cee\n")
            fh.writelines("#SBATCH --begin=now+" + str(extra_time)+"\n")
            fh.writelines("ml python/2.7.13\n")
            fh.writelines("python ModuleT_Single_Scenario.py %s\n" % ind_run)

        os.system("sbatch %s" % job_file)