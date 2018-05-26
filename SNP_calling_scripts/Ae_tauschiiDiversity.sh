#!/bin/bash -l

#SBATCH --job-name=Ae_tauschiiDiversity
#SBATCH -p killable.q

#SBATCH --mem-per-cpu=10G   					# Memory per core, use --mem= for memory per node
#SBATCH --nodes=1	--ntasks-per-node=10
#SBATCH --time=02-00:00:00   					# Use the form DD-HH:MM:SS
#SBATCH --mail-type=ALL --output="%x.%j.o" --error="%x.%j.e"

# set memory and cores to use in downstream analysis (generally no need to change)
mem=95; cores="$SLURM_TASKS_PER_NODE"

## remove comment (#) from which dbPath you want to use (first is tauschii then wheat)
	dbPath="/bulk/user/genomes/Ae_tauschii/Aet_v4.0/Aet_v4.0.fasta"

# exit script if dbPath is not selected
	[[ -z "$dbPath" ]] && { echo "ERROR: dbPath not selected" ; exit 1; }

#############################################
## NO NEED TO CHANGE ANYTHING FROM HERE ON ##
#############################################
name="$SLURM_JOB_NAME"
user="$(whoami)"
keyFile="/bulk/${user}/gbs/jobs/${name}.txt"
seqDir="/bulk/user/sequence"
tasselPath="/homes/user/softwares/tassel5/run_pipeline.pl"

mkdir "/bulk/${user}/gbs/projects/${name}"
mkdir "/bulk/${user}/gbs/projects/${name}/keyFileSh"
cd "/bulk/${user}/gbs/projects/${name}"

# Path for required software
export PATH=$PATH:/homes/user/usr/bin:/homes/user/usr/bin/bin

# Set JAVA VM Version
module load Java

echo " Exit code 0 is good and 1 is bad " >> z_exitcode.log
echo "==================================" >> z_exitcode.log

## GBSSeqToTagDBPlugin - RUN Tags to DB
$tasselPath -Xms${mem}G -Xmx${mem}G -fork1 -GBSSeqToTagDBPlugin -e PstI-MspI \
    -i $seqDir -db ${name}.db -k $keyFile -kmerLength 64 -minKmerL 64 -mnQS 20 -mxKmerNum 250000000 -endPlugin -runfork1 >> z_pipeline.log
echo "1. GBSSeqToTagDBPlugin: The exit code is $?" >> z_exitcode.log

## TagExportToFastqPlugin - export Tags
$tasselPath -Xms${mem}G -Xmx${mem}G -fork1 -TagExportToFastqPlugin \
    -db ${name}.db -o ${name}_tagsForAlign.fa.gz -c 10 -endPlugin -runfork1 >> z_pipeline.log
echo "2. TagExportToFastqPlugin: The exit code is $?" >> z_exitcode.log

## RUN BOWTIE
bowtie2 -p ${cores} --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.25 \
    -x $dbPath -U ${name}_tagsForAlign.fa.gz -S ${name}.sam >> z_pipeline.log
echo "3. bowtie2: The exit code is $?" >> z_exitcode.log

## SAMToGBSdbPlugin - SAM to DB
$tasselPath -Xms${mem}G -Xmx${mem}G -fork1 -SAMToGBSdbPlugin \
    -i ${name}.sam -db ${name}.db -aProp 0.0 -aLen 0 -endPlugin -runfork1 >> z_pipeline.log
echo "4. SAMToGBSdbPlugin: The exit code is $?" >> z_exitcode.log

## DiscoverySNPCallerPluginV2 - RUN DISCOVERY SNP CALLER
$tasselPath -Xms${mem}G -Xmx${mem}G -fork1 -DiscoverySNPCallerPluginV2 \
    -db ${name}.db -mnLCov 0.1 -mnMAF 0.01 -deleteOldData true -endPlugin -runfork1 >> z_pipeline.log
echo "5. DiscoverySNPCallerPlugin: The exit code is $?" >> z_exitcode.log

## SNPQualityProfilerPlugin - RUN QUALITY PROFILER
$tasselPath -Xms${mem}G -Xmx${mem}G -fork1 -SNPQualityProfilerPlugin \
    -db ${name}.db -statFile ${name}_SNPqual_stats.txt -endPlugin -runfork1 >> z_pipeline.log
echo "6. SNPQualityProfilerPlugin all: The exit code is $?" >> z_exitcode.log

## UpdateSNPPositionQualityPlugin - UPDATE DATABASE WITH QUALITY SCORE
$tasselPath -Xms${mem}G -Xmx${mem}G -fork1 -UpdateSNPPositionQualityPlugin \
    -db ${name}.db -qsFile ${name}_SNPqual_stats.txt -endPlugin -runfork1 >> z_pipeline.log
echo "7. UpdateSNPPositionQualityPlugin: The exit code is $?" >> z_exitcode.log

## ProductionSNPCallerPluginV2 - RUN PRODUCTION PIPELINE - output .h5
$tasselPath -Xms${mem}G -Xmx${mem}G -fork1 -ProductionSNPCallerPluginV2 \
    -db ${name}.db -i $seqDir -k $keyFile -o ${name}.vcf -e PstI-MspI -kmerLength 64 -endPlugin -runfork1 >> z_pipeline.log
echo "8. ProductionSNPCallerPlugin: The exit code is $?" >> z_exitcode.log

## Convert to Hapmap format
$tasselPath -Xms${mem}G -Xmx${mem}G -fork1 -vcf ${name}.vcf -export ${name} -exportType Hapmap
echo "9. HapmapConversion: The exit code is $?" >> z_exitcode.log

## Get errors in one file
grep "ERROR" z_pipeline.log >> z_ERROR.log

## Move all the files to the project directory
mv /bulk/${user}/gbs/jobs/${name}.* keyFileSh/
