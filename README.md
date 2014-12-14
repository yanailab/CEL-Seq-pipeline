Overview
==========
This is a pipeline for the CEL-Seq method.
CEL-Seq is a HTS method published hashimshony et. al(2012), based on RNAseq. 
It typcially reads only the 3'UTR end of mRNA transcripts and is stranded.

The pipeline inpuuts read files in FASTQ format, alongside a reference genome and annotation and outputs a read count.

Dependencies
==============
General: bowtie2  
Python: python-dev, htseq, argparse

Orientation (local copy-YanaiLab)
===============
You are in one of three:  
* The full-repositories/work directory /data/tools/pipeline_dev  
      Houses dev and master and sub-branches 
* The dev directory /data/tools/pipeline_dev:
      Only has the stable development branch
* The stable directory /data/tools/pipeline_stable  
      Only houses the stable master branch

Synching across (local copy-YanaiLab)
===================
After done making changes at pipeline the stable/master branch;  
* Sync to github:  
        `cd /data/tools/pipeline`
        `git push origin master`  
* Sync to stable pipeline:  
        `cd /data/tools/pipeline_stable`  
        `git pull origin master`  
For updating changes at the dev branch;
* Sync to dev pipeline:
        `cd /data/tools/pipeline_dev`
        `git pull origin dev`

