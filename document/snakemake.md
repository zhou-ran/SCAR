
## Initiate snakemke pipeline using SCAR

### Prepare configure information for snakemake

```
library(SCAR)

config_lst <- list(
'output' = 'data', # the file name of output foder
'samples' = c('LL'), # the sample label
'FASTQ' = 'fq', # the file location of fastq files
'STAR' = list('mix_index' = 'index', 
             'cores' = as.integer(8)), # the parameter for STAR
'soft' = list(
 'STAR' = 'STAR', # the soft location of STAR
 'seqkit' = 'seqkit', # the soft location of seqkit (https://github.com/shenwei356/seqkit)
 'split_script' = 'bin/split_read_id.py', # don't modify this
 'seqkt' = 'seqkt', # the soft location of seqkt (https://github.com/lh3/seqtk)
 'python' = 'python' # the soft location of python)
 )

# generate snakemake pipeline
init_smk_script(dir = 'smk_test/',
             config_list = config_lst,
             overwrite = T)



```


### Run the pipeline using snakemake 

```
pip install snakemake 

cd smk_test/

# dry-run
snakemake -s Split_PDX.py -np 

# Job stats:
# job                  count    min threads    max threads
# -----------------  -------  -------------  -------------
# align_to_mix_ref         4              1              1
# all                      1              1              1
# get_read_hg38            4              1              1
# get_read_mm10            4              1              1
# load_share_memo          1              1              1
# remove_share_memo        1              1              1
# split_fq                 1              1              1
# split_reads              4              1              1
# total                   20              1              1
# 
# This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.

```

