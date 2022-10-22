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
pip install git+https://https://github.com/jyesselm/dreem.git@jyesselm/dev

# rna library for manipulation of rna secondary structure 
pip install git+https://github.com/YesselmanLab/rna_library
```

### How to run 
Here is a test run using example data. Will demultiplex the first barcode which uses the second helix in the RNA construct design. `-m 1` will only run the first barcode. 
```shell
barcode_demultiplex -csv examples/ires_minittrs_2/C009C.csv -fq1 examples/ires_minittrs_2/test_R2_001.fastq -fq2 examples/ires_minittrs_2/test_R1_001.fastq -helix 2 0 8 -m 1 
```

If you would like to also run dreem to get the mutation fraction of each residue you can add `--run-dreem` to the command. For each construct dreem will be run and will generate the summaried output in the form of a pickle file. Adding `--include-all-dreem-outputs` will generate all the standard outputs of dreem. This is not recommended for large libraries as it will yield thousands of output files.

```shell
barcode_demultiplex -csv examples/ires_minittrs_2/C009C.csv -fq1 examples/ires_minittrs_2/test_R2_001.fastq -fq2 examples/ires_minittrs_2/test_R1_001.fastq -helix 2 0 8 -m 1 --run-dreem  
```

