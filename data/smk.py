import os
import sys
import re


split_inds = ['001','002','003','004']

## Configuration file
if len(config) == 0:
	if os.path.isfile('./config.yaml'):
		configfile: './config.yaml'
	else:
		sys.exit('Make sure there is a config.yaml file in ' + os.getcwd() + ' or specify one with the --configfile commandline parameter.')

samples = config['samples']

def getpath(str):
	if str in ['', '.', './']:
		return ''
	if str.startswith('./'):
		regex = re.compile('^\./?')
		str = regex.sub('', str)
	if not str.endswith('/'):
		str += '/'
	return str

output_dir = getpath(config['output'])
fastq_dir = getpath(config['FASTQ'])


rule all:
	input:
		output_dir + 'STAR/remove_share_memo/ShareExit.log',
		expand(output_dir + 'split_reads/{sample}.part_{split_ind}/hg38.tsv.gz', sample = samples, split_ind = split_inds),
		expand(output_dir + 'hg38/{sample}/fq/{sample}_S1_L{split_ind}_R2_001.fastq.gz', sample = samples, split_ind=split_inds),
		expand(output_dir + 'mm10/{sample}/fq/{sample}_S1_L{split_ind}_R2_001.fastq.gz', sample = samples, split_ind=split_inds),

rule split_fq:
	input:
		fastq_dir + '{sample}_S1_L001_R2_001.fastq.gz'

	output:
		temp(expand(output_dir + '{{sample}}/{{sample}}_S1_L001_R2_001.part_{split_ind}.fastq.gz', split_ind = split_inds))

	params:
		prefix = 'data/{sample}/',
		soft = config['soft']['seqkit']

	shell:
		'''

		{params.soft} split2 \
		-1 {input} \
		-p 4 \
		-O {params.prefix}

		'''	
rule load_share_memo:
	output:
		output_dir + 'STAR/load_share_memo/ShareLoadInfo.log'
	params:
		STAR = config['soft']['STAR'],
		index = config['STAR']['mix_index'],
		outFileNamePrefix = output_dir + 'STAR/load_share_memo/'
	shell:
		'''

		{params.STAR} --genomeLoad LoadAndExit \
		--genomeDir {params.index} \
		--outFileNamePrefix {params.outFileNamePrefix} > {output}

		'''

rule align_to_mix_ref:
	input:
		fq = output_dir + '{sample}/{sample}_S1_L001_R2_001.part_{split_ind}.fastq.gz',
		log = rules.load_share_memo.output
	output:
		output_dir + '{sample}/{split_ind}.Aligned.out.bam'
	params:
		STAR = config['soft']['STAR'],
		cores = config['STAR']['cores'],
		prefix = output_dir + '{sample}/{split_ind}.',
		index = config['STAR']['mix_index']
	shell:
		'''

		{params.STAR} \
		 --runThreadN {params.cores} \
		 --genomeDir {params.index} \
		 --outSAMtype BAM Unsorted \
		 --outBAMcompression 9 \
		 --limitBAMsortRAM 60000000000 \
		 --readFilesCommand zcat \
		 --readFilesIn {input.fq} \
		 --outFileNamePrefix {params.prefix}

		'''

rule split_reads:
	input:
		output_dir + '{sample}/{split_ind}.Aligned.out.bam'
	output:
		mm = output_dir + 'split_reads/{sample}.part_{split_ind}/mm10.tsv.gz',
		hg = output_dir + 'split_reads/{sample}.part_{split_ind}/hg38.tsv.gz'

	params:
		prefix = output_dir + 'split_reads/{sample}.{split_ind}',
		script = config['soft']['split_script'],
		python = config['soft']['python']

	shell:
		'''

		{params.python} {params.script} \
		{input} {params.prefix}

		'''
    
rule get_read_mm10:
	input:
		fq_id = output_dir + 'split_reads/{sample}.part_{split_ind}/mm10.tsv.gz',
		r1 = fastq_dir + '{sample}_S1_L001_R1_001.fastq.gz',
		r2 = fastq_dir + '{sample}_S1_L001_R2_001.fastq.gz'
	output:
		r1 = output_dir + 'mm10/{sample}/fq/{sample}_S1_L{split_ind}_R1_001.fastq.gz',
		r2 = output_dir + 'mm10/{sample}/fq/{sample}_S1_L{split_ind}_R2_001.fastq.gz'

	params:
		seqtk = config['soft']['seqkt']

	shell:
		'''

		{params.seqtk} subseq {input.r1} {input.fq_id} | gzip > {output.r1};
		{params.seqtk} subseq {input.r2} {input.fq_id} | gzip > {output.r2}

		'''


rule get_read_hg38:
	input:
		fq_id = output_dir + 'split_reads/{sample}.part_{split_ind}/hg38.tsv.gz',
		r1 = fastq_dir + '{sample}_S1_L001_R1_001.fastq.gz',
		r2 = fastq_dir + '{sample}_S1_L001_R2_001.fastq.gz'
	output:
		r1 = output_dir + 'hg38/{sample}/fq/{sample}_S1_L{split_ind}_R1_001.fastq.gz',
		r2 = output_dir + 'hg38/{sample}/fq/{sample}_S1_L{split_ind}_R2_001.fastq.gz'

	params:
		seqtk = config['soft']['seqkt']

	shell:
		'''

		{params.seqtk} subseq {input.r1} {input.fq_id} | gzip > {output.r1};
		{params.seqtk} subseq {input.r2} {input.fq_id} | gzip > {output.r2}

		'''
rule remove_share_memo:
	input:
		expand(output_dir + '{sample}/{split_ind}.Aligned.out.bam', sample = samples, split_ind = split_inds)
	output:
		output_dir + 'STAR/remove_share_memo/ShareExit.log'
	params:
		STAR = config['soft']['STAR'],
		index = config['STAR']['mix_index'],
		outFileNamePrefix = output_dir + 'STAR/remove_share_memo/'

	shell:
		'''

		{params.STAR} --genomeLoad Remove \
		--genomeDir {params.index} \
		--outFileNamePrefix {params.outFileNamePrefix} > {output}

		'''




