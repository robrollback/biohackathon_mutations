#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# DnaSeq BatchScheduler Job Submission Bash script
# Version: 2.2.1-beta
# Created on: 2017-02-28T11:25:59
# Steps:
#   trimmomatic: 1 job
#   merge_trimmomatic_stats: 1 job
#   bwa_mem_picard_sort_sam: 2 jobs
#   picard_merge_sam_files: 1 job
#   gatk_indel_realigner: 1 job
#   merge_realigned: 1 job
#   fix_mate_by_coordinate: 2 jobs
#   picard_mark_duplicates: 2 jobs
#   recalibration: 3 jobs
#   gatk_haplotype_caller: 1 job
#   merge_and_call_individual_gvcf: 1 job
#   TOTAL: 16 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=~/
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/DnaSeq_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

SEPARATOR_LINE=`seq -s - 80 | sed 's/[0-9]//g'`

#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic.gcat_set_025_chr22
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.gcat_set_025_chr22
JOB_DONE=job_output/trimmomatic/trimmomatic.gcat_set_025_chr22.4a37a871f139f32ca3b5a06d7196c0ba.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/gcat_set_025_chr22 && \
`cat > trim/gcat_set_025_chr22/gcat_set_025_chr22.trim.adapters.fa << END
>Prefix/1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>Prefix/2
TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 1 \
  -phred33 \
  /gs/project/mugqic/projects/GCAT_benchmark/raw_reads/gcat_set_025/gcat_set_025.pair1.chr22.fastq.gz \
  /gs/project/mugqic/projects/GCAT_benchmark/raw_reads/gcat_set_025/gcat_set_025.pair2.chr22.fastq.gz \
  trim/gcat_set_025_chr22/gcat_set_025_chr22.trim.pair1.fastq.gz \
  trim/gcat_set_025_chr22/gcat_set_025_chr22.trim.single1.fastq.gz \
  trim/gcat_set_025_chr22/gcat_set_025_chr22.trim.pair2.fastq.gz \
  trim/gcat_set_025_chr22/gcat_set_025_chr22.trim.single2.fastq.gz \
  ILLUMINACLIP:trim/gcat_set_025_chr22/gcat_set_025_chr22.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/gcat_set_025_chr22/gcat_set_025_chr22.trim.log
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.0c7a52fdaf47b84011634a399f55c417.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Paired Reads #	Surviving Paired Reads #	Surviving Paired Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/gcat_set_025_chr22/gcat_set_025_chr22.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/gcat_set_025_chr22	gcat_set_025_chr22	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
cut -f1,3- metrics/trimReadsetTable.tsv | awk -F"	" '{OFS="	"; if (NR==1) {if ($2=="Raw Paired Reads #") {paired=1};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"} else {if (paired) {$2=$2*2; $3=$3*2}; raw[$1]+=$2; surviving[$1]+=$3}}END{for (sample in raw){print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}' \
  > metrics/trimSampleTable.tsv && \
mkdir -p report && \
cp metrics/trimReadsetTable.tsv metrics/trimSampleTable.tsv report/ && \
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"} else {print $1, $2, sprintf("%\47d", $3), sprintf("%\47d", $4), sprintf("%.1f", $5)}}' metrics/trimReadsetTable.tsv` && \
pandoc \
  /home/reveleigh/github/tumor_fresh_new/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --template /home/reveleigh/github/tumor_fresh_new/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --variable trailing_min_quality=30 \
  --variable min_length=50 \
  --variable read_type=Paired \
  --variable trim_readset_table="$trim_readset_table_md" \
  --to markdown \
  > report/Illumina.merge_trimmomatic_stats.md
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: bwa_mem_picard_sort_sam
#-------------------------------------------------------------------------------
STEP=bwa_mem_picard_sort_sam
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam.gcat_set_025_chr22
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.gcat_set_025_chr22
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.gcat_set_025_chr22.d5ecc6c25f7eb4cbaaff9a6d94631943.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/gcat_set_025_chr22/gcat_set_025_chr22 && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:gcat_set_025_chr22	SM:gcat_set_025_chr22	LB:run01	PU:run1934_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/bwa_index/Homo_sapiens.GRCh37.fa \
  trim/gcat_set_025_chr22/gcat_set_025_chr22.trim.pair1.fastq.gz \
  trim/gcat_set_025_chr22/gcat_set_025_chr22.trim.pair2.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/localscratch/ \
 INPUT=/dev/stdin \
 OUTPUT=alignment/gcat_set_025_chr22/gcat_set_025_chr22/gcat_set_025_chr22.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=3750000
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_report
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam_report
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam_report.f46a10282ed6eec2998316b25218fada.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /home/reveleigh/github/tumor_fresh_new/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  --variable scientific_name="Homo_sapiens" \
  --variable assembly="GRCh37" \
  /home/reveleigh/github/tumor_fresh_new/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  > report/DnaSeq.bwa_mem_picard_sort_sam.md
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: symlink_readset_sample_bam.gcat_set_025_chr22
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.gcat_set_025_chr22
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.gcat_set_025_chr22.b76fffd05b872d4ca6c3e0108aad3cc2.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
mkdir -p alignment/gcat_set_025_chr22 && \
ln -s -f gcat_set_025_chr22/gcat_set_025_chr22.sorted.bam alignment/gcat_set_025_chr22/gcat_set_025_chr22.sorted.bam && \
ln -s -f gcat_set_025_chr22/gcat_set_025_chr22.sorted.bai alignment/gcat_set_025_chr22/gcat_set_025_chr22.sorted.bai
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: gatk_indel_realigner
#-------------------------------------------------------------------------------
STEP=gatk_indel_realigner
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner.gcat_set_025_chr22
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.gcat_set_025_chr22
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.gcat_set_025_chr22.11a0d0fa7da5fc9d95dec25ff553bcf2.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.7 && \
mkdir -p alignment/gcat_set_025_chr22/realign && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx3200M -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator  \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/gcat_set_025_chr22/gcat_set_025_chr22.sorted.bam \
  --intervals 22 \
  --known /gs/project/mugqic/analyste_dev/phase2/genomes/Homo_sapiens/hg1k_v37/annotations/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
  --out alignment/gcat_set_025_chr22/realign/all.intervals && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx3200M -jar $GATK_JAR \
  --analysis_type IndelRealigner  \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --intervals 22 \
  --input_file alignment/gcat_set_025_chr22/gcat_set_025_chr22.sorted.bam \
   \
  --targetIntervals alignment/gcat_set_025_chr22/realign/all.intervals \
  --knownAlleles /gs/project/mugqic/analyste_dev/phase2/genomes/Homo_sapiens/hg1k_v37/annotations/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
   \
  --out alignment/gcat_set_025_chr22/realign/all.bam \
  --maxReadsInMemory 500000 && \
ln -s -f realign/all.bam alignment/gcat_set_025_chr22/gcat_set_025_chr22.realigned.qsorted.bam
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: merge_realigned
#-------------------------------------------------------------------------------
STEP=merge_realigned
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_realigned_report
#-------------------------------------------------------------------------------
JOB_NAME=merge_realigned_report
JOB_DONE=job_output/merge_realigned/merge_realigned_report.042f4615b2029cb7107cf37efebedc9d.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
mkdir -p report && \
cp \
  /home/reveleigh/github/tumor_fresh_new/bfx/report/DnaSeq.gatk_indel_realigner.md \
  report/DnaSeq.gatk_indel_realigner.md
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: fix_mate_by_coordinate
#-------------------------------------------------------------------------------
STEP=fix_mate_by_coordinate
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: fix_mate_by_coordinate.gcat_set_025_chr22
#-------------------------------------------------------------------------------
JOB_NAME=fix_mate_by_coordinate.gcat_set_025_chr22
JOB_DONE=job_output/fix_mate_by_coordinate/fix_mate_by_coordinate.gcat_set_025_chr22.5caa7d3739760830ae0c49b02be72d68.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/picard/1.123 && \
java -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $BVATOOLS_JAR \
  groupfixmate \
  --level 1 \
  --bam alignment/gcat_set_025_chr22/gcat_set_025_chr22.realigned.qsorted.bam \
  --out alignment/gcat_set_025_chr22/gcat_set_025_chr22.matefixed.sorted.tmp.bam && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/localscratch/ \
 INPUT=alignment/gcat_set_025_chr22/gcat_set_025_chr22.matefixed.sorted.tmp.bam \
 OUTPUT=alignment/gcat_set_025_chr22/gcat_set_025_chr22.matefixed.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=3750000
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: fix_mate_by_coordinate_report
#-------------------------------------------------------------------------------
JOB_NAME=fix_mate_by_coordinate_report
JOB_DONE=job_output/fix_mate_by_coordinate/fix_mate_by_coordinate_report.a7ccc99a84daab81a3eadea52795b99f.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
mkdir -p report && \
cp \
  /home/reveleigh/github/tumor_fresh_new/bfx/report/DnaSeq.fix_mate_by_coordinate.md \
  report/DnaSeq.fix_mate_by_coordinate.md
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates.gcat_set_025_chr22
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.gcat_set_025_chr22
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.gcat_set_025_chr22.c00fa75615586d4d4ee8653d11206b5f.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/localscratch/ \
 INPUT=alignment/gcat_set_025_chr22/gcat_set_025_chr22.matefixed.sorted.bam \
 OUTPUT=alignment/gcat_set_025_chr22/gcat_set_025_chr22.sorted.dup.bam \
 METRICS_FILE=alignment/gcat_set_025_chr22/gcat_set_025_chr22.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates_report
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates_report.6498b27e93615e5d986e8cb8d37fb841.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
mkdir -p report && \
cp \
  /home/reveleigh/github/tumor_fresh_new/bfx/report/DnaSeq.picard_mark_duplicates.md \
  report/DnaSeq.picard_mark_duplicates.md
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: recalibration
#-------------------------------------------------------------------------------
STEP=recalibration
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: interval_list.gcat_set_025.b37.chr22.bed
#-------------------------------------------------------------------------------
JOB_NAME=interval_list.gcat_set_025.b37.chr22.bed
JOB_DONE=job_output/recalibration/interval_list.gcat_set_025.b37.chr22.bed.8f088815f42143bb1ce0f443307add01.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
module load mugqic/mugqic_tools/2.1.6 mugqic/perl/5.22.1 && \
bed2IntervalList.pl \
  --dict /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.dict \
  --bed /gs/project/mugqic/projects/GCAT_benchmark/raw_reads/gcat_set_025/gcat_set_025.b37.chr22.bed \
  > /gs/project/mugqic/projects/GCAT_benchmark/raw_reads/gcat_set_025/gcat_set_025.b37.chr22.interval_list
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: recalibration.gcat_set_025_chr22
#-------------------------------------------------------------------------------
JOB_NAME=recalibration.gcat_set_025_chr22
JOB_DONE=job_output/recalibration/recalibration.gcat_set_025_chr22.e87cd7a3253f919accb134efb951e4ea.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.7 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx11G -jar $GATK_JAR \
  --analysis_type BaseRecalibrator \
  --num_cpu_threads_per_data_thread 6 \
  --input_file alignment/gcat_set_025_chr22/gcat_set_025_chr22.sorted.dup.bam \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa  \
  --intervals /gs/project/mugqic/projects/GCAT_benchmark/raw_reads/gcat_set_025/gcat_set_025.b37.chr22.interval_list \
  --knownSites /gs/project/mugqic/analyste_dev/phase2/genomes/Homo_sapiens/hg1k_v37/annotations/dbSnp-138.vcf.gz \
  --knownSites /gs/project/mugqic/analyste_dev/phase2/genomes/Homo_sapiens/hg1k_v37/annotations/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
  --out alignment/gcat_set_025_chr22/gcat_set_025_chr22.sorted.dup.recalibration_report.grp && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx11G -jar $GATK_JAR \
  --analysis_type PrintReads \
  --num_cpu_threads_per_data_thread 6 \
  --input_file alignment/gcat_set_025_chr22/gcat_set_025_chr22.sorted.dup.bam \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --BQSR alignment/gcat_set_025_chr22/gcat_set_025_chr22.sorted.dup.recalibration_report.grp \
  --out alignment/gcat_set_025_chr22/gcat_set_025_chr22.sorted.dup.recal.bam && \
md5sum alignment/gcat_set_025_chr22/gcat_set_025_chr22.sorted.dup.recal.bam > alignment/gcat_set_025_chr22/gcat_set_025_chr22.sorted.dup.recal.bam.md5
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: recalibration_report
#-------------------------------------------------------------------------------
JOB_NAME=recalibration_report
JOB_DONE=job_output/recalibration/recalibration_report.37566d0ed0df60bff452f87f3f80c877.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
mkdir -p report && \
cp \
  /home/reveleigh/github/tumor_fresh_new/bfx/report/DnaSeq.recalibration.md \
  report/DnaSeq.recalibration.md
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: gatk_haplotype_caller
#-------------------------------------------------------------------------------
STEP=gatk_haplotype_caller
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller.gcat_set_025_chr22
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.gcat_set_025_chr22
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.gcat_set_025_chr22.d8e26b63eba2178925a2df5a0fbbcd88.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.7 && \
mkdir -p alignment/gcat_set_025_chr22/rawHaplotypeCaller && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -dt none -nct 1 \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --intervals 22 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/gcat_set_025_chr22/gcat_set_025_chr22.sorted.dup.recal.bam \
  --out alignment/gcat_set_025_chr22/rawHaplotypeCaller/gcat_set_025_chr22.hc.g.vcf.gz
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: merge_and_call_individual_gvcf
#-------------------------------------------------------------------------------
STEP=merge_and_call_individual_gvcf
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_and_call_individual_gvcf.gcat_set_025_chr22
#-------------------------------------------------------------------------------
JOB_NAME=merge_and_call_individual_gvcf.gcat_set_025_chr22
JOB_DONE=job_output/merge_and_call_individual_gvcf/merge_and_call_individual_gvcf.gcat_set_025_chr22.77cb39d2b1a36f43170c651dd22736da.mugqic.done
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.7 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -cp $GATK_JAR \
  org.broadinstitute.gatk.tools.CatVariants  \
  --reference /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --variant alignment/gcat_set_025_chr22/rawHaplotypeCaller/gcat_set_025_chr22.hc.g.vcf.gz \
  --outputFile alignment/gcat_set_025_chr22/gcat_set_025_chr22.hc.g.vcf.gz && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx32G -jar $GATK_JAR \
  --analysis_type GenotypeGVCFs -nt 3 \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --intervals 22 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --variant alignment/gcat_set_025_chr22/gcat_set_025_chr22.hc.g.vcf.gz \
  --out alignment/gcat_set_025_chr22/gcat_set_025_chr22.hc.vcf.gz
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r17-n02&ip=10.241.129.12&pipeline=DnaSeq&steps=trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,picard_merge_sam_files,gatk_indel_realigner,merge_realigned,fix_mate_by_coordinate,picard_mark_duplicates,recalibration,gatk_haplotype_caller,merge_and_call_individual_gvcf&samples=1&AnonymizedList=7992c18ba68c4a2705a534d5d4981de9" --quiet --output-document=/dev/null

