import pysam
import sys
import fire
import os
import gzip

def main(bamfile, out_prefix):
	try:
		os.makedirs(out_prefix)
	except OSError:
		pass

	bam_fh = pysam.AlignmentFile(bamfile, 'rb')
	mm_out = gzip.open(f'{out_prefix}/mm10.tsv.gz', 'w')
	hg_out = gzip.open(f'{out_prefix}/hg38.tsv.gz', 'w')
	mix_out = gzip.open(f'{out_prefix}/mix.tsv.gz', 'w')

	c_id = ''
	chrom_set = set()

	for line in bam_fh.fetch(until_eof=True):
		if not c_id:
			c_id = line.query_name
			if line.reference_name.startswith('hg38'):
				chrom_set.add('hg38')
			else:
				chrom_set.add('mm10')
		else:
			if c_id != line.query_name:
				if len(chrom_set) == 2:
					mix_out.write(f'{c_id}\n'.encode())
				elif len(chrom_set) == 1 and 'mm10' in chrom_set:
					mm_out.write(f'{c_id}\n'.encode())
				elif len(chrom_set) == 1 and 'hg38' in chrom_set:
					hg_out.write(f'{c_id}\n'.encode())
				else:
					sys.exit(f'Error found in processing {c_id}')

				c_id = line.query_name
				if line.reference_name.startswith('hg38'):
					chrom_set = set(['hg38'])
				else:
					chrom_set = set(['mm10'])

			else:
				if line.reference_name.startswith('hg38'):
					chrom_set.add('hg38')
				else:
					chrom_set.add('mm10')

	if len(chrom_set) == 2:
		mix_out.write(f'{c_id}\n'.encode())
	elif len(chrom_set) == 1 and 'mm10' in chrom_set:
		mm_out.write(f'{c_id}\n'.encode())
	elif len(chrom_set) == 1 and 'hg38' in chrom_set:
		hg_out.write(f'{c_id}\n'.encode())
	else:
		sys.exit(f'Error found in processing {c_id}')

	mm_out.close()
	hg_out.close()
	mix_out.close()

if __name__ == '__main__':
	fire.Fire(main)
