#! /bin/bash
#$ -S /bin/bash
#$ -N DIDA-qmake
#$ -V
#$ -wd /genesis/scratch/hmohamadi/DIDA/chr14/
#$ -j y
#$ -o $JOB_ID.o

export PATH=/genesis/scratch/hmohamadi/DIDA/bin/:$PATH
export reads=/scratch/hmohamadi/partition/chr14/treads.fastq


echo "hostname: $(hostname)"
echo "temp dir: $TMPDIR"
echo "pwd:      $(pwd)"
echo "job id:   $JOB_ID"

echo "Job started at: $(date)"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

cd /genesis/scratch/hmohamadi/DIDA/chr14/
touch ABySSTest-3-{1..12}.tmp
ssh login4 "cd /genesis/scratch/hmohamadi/DIDA/chr14/ && export reads=/scratch/hmohamadi/partition/chr14/treads.fastq && qmake -v PATH -inherit -cwd -C /genesis/scratch/hmohamadi/DIDA/chr14/ -- -j 13 -f dida_makefile all CTIGS=ABySSTest-3 PARTS=12 DMERS=28"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "Job ended at:   $(date)"

