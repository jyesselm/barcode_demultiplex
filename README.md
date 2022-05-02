barcode_demultiplex
-------------------
Barcode_demultiplex is a lightweight program for demultiplexing RNA sequencing reads with internal barcodes in the form of helices or hairpins.

### Install 
Requires `seqkit` to do the sequence matching. If you have conda you can install it with 

```shell
conda install -c bioconda seqkit 
```
If you don't have conda or are on the new ARM macs. you can download it here: https://bioinf.shenwei.me/seqkit/download/

There are also a few packages that need to be installed from GitHub repos.

```shell
# dreem which handles the pipeline for doing the downstream alignment and mutation counting for DMS-MaPseq data. 
# In theory, could probably make this a
# conditional requirement but haven't done that yet
https://github.com/jyesselm/dreem

# rna library for manipulation of rna secondary structure 
https://github.com/YesselmanLab/rna_library
# to get all the requirements
pip install -r requirements.txt 

# other packages required at the moment written in the lab 
pip install vienna
pip install rna_seq_tools
```

### How to run 
```shell
python barcode_demultiplex/run.py -csv examples/ires_minittrs_2/C009C.csv -fq1 examples/ires_minittrs_2/test_R2_001.fastq -fq2 examples/ires_minittrs_2/test_R1_001.fastq -helix 2 0 8
```

