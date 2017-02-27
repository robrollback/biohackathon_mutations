# Biohackathon: Rapid indentification of cancer mutations
***By Robert Eveleigh, MSc., Mathieu Bourgey, PhD., Guillaume Bourque, PhD.***

In this biohackathon module, you are working at the Genome Quebec sequencing patient's tumors to guide decisions in medical treatment.  Your challenge is to successfully analyze the given tumor sample and identify the mutations in the tumor as quickly as possible.  The new Illumina seqeuencer, the Illumina SuperDuper, provides files containing all the reads (fastq or bam) and your goal is to generate the file with the mutations (vcf).

We will be working with one sample NA12878 from the CEPH 1463 pedigree sequenced by the [Genome in a Bottle](http://jimb.stanford.edu/giab/) (GIAB) consortium.

Due to time and resource restraints, we will provide reads from chromosome 22 only.

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.

## Original Setup

The initial structure of your folders should look like this:
```
<ROOT>
|--raw_reads/                
    `--gcat_set_25/          # raw reads of chr22 from the SuperDuper sequencer. Contains fastqs, initial bam, and bed files.
|-- bin/                     # bin directory. Contains bash scripts to execute standard pipeline, and benchmarking.
|-- benchmarking/            # benchmarking results. Contains high quality variants from the GIAB
    
```

### Environment setup
The team will utilizing Mcgill's HPC cluster: Guillimin.  Accounts have been setup and credientials to login in will be provided.

Once the team has access to the HPC, the team will need to setup the mugqic pipeline environment to allow the load of modules for the programs you will be employing.

To setup enviroment please add the following lines to your bash_profile.
See [Quick setup for abacus, guillimin and mammouth users](https://bitbucket.org/mugqic/mugqic_pipelines#markdown-header-quick-setup-for-abacus-guillimin-and-mammouth-users)

#Data processing and variant discovery

![Data processing diagram](img/dnaseq_pipeline.jpg)
    
The NA12878 samples from the CEPH pedigree has been sequenced. Given the raw data, C3G has processed the data using the GATK's 'best practices' as ground truth. This includes trimming, quality reads mapped to the reference genome and the resulting bam files have been further processed using indel realigner, fixmates, and base recalibration to reduce known systematic biases. Variants were called with GATK haplotype caller in GVCF mode.

###Goal

The team's goal is to achieve similar F1 scores while decreasing computational time.

Processing time: ~30 mins

| True-pos | False-pos | False-neg | Precision | Sensitivity | F-measure |
|:-------- |:--------- |:--------- |:--------- |:----------- |:--------- |
|553       |36         |167        |0.9389     |0.7681       |0.8449     |
