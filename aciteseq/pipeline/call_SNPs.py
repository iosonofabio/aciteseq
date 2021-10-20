'''
Author:	        Fabio Zanini
Date:	        2021/02/03
Description:	Call SNPs for the demux
'''
import os
import sys
import argparse
import pandas as pd
import subprocess as sp


if __name__ == '__main__':

    pa = argparse.ArgumentParser()
    pa.add_argument('--samples', default='iMSC,AdMSC')
    pa.add_argument('--steps', required=True)
    args = pa.parse_args()

    samples = args.samples.split(',')
    steps = args.steps.split(',')
    
    fdn_data = os.path.abspath(
            '../../data/sequencing/experiments/202012_pilot1'
            )
    fdn_transcriptome = '/home/fabio/software/cellranger_reference_data/refdata-gex-GRCh38-2020-A'
 
    for sample in samples:
        print(sample)
        fdn_in = f'{fdn_data}/gene_expression_raw/{sample}'
        fdn_out = f'{fdn_data}/gene_expression_bam/{sample}'
        os.makedirs(fdn_out, exist_ok=True)

        if 'map' in steps:
            cwd = os.path.abspath(os.curdir)
            try:
                call = (
                    'cellranger-4.0.0 count '
                    #'--dry ' #FIXME
                    '--nosecondary '
                    f'--id {sample} '
                    f'--transcriptome {fdn_transcriptome} '
                    f'--fastqs {fdn_in} '
                    f'--sample {sample} '
                    '--force-cells 8000 '
                    '--localcores 32 '
                    '--localmem 128 '
                    )
                print(call)
                os.chdir(fdn_out)
                sp.run(call, shell=True, check=True)
            finally:
                os.chdir(cwd)

        if 'NtoD' in steps:
            print('Convert N CIGARs into D for freebayes - sigh!')
            import pysam

            fn_bam = f'{fdn_data}/gene_expression_bam/{sample}/{sample}/outs/possorted_genome_bam.bam'
            fn_newbam = f'{fdn_data}/gene_expression_vcf/{sample}/possorted_genome_bam_NtoD.bam'
            os.makedirs(os.path.dirname(fn_newbam), exist_ok=True)
            with pysam.AlignmentFile(fn_bam, 'r') as f:
                with pysam.AlignmentFile(fn_newbam, 'wb', template=f) as fout:
                    for ir, read in enumerate(f):
                        if ir % 100000 == 0:
                            print(ir, end='\r')
                        if read.cigarstring is None:
                            continue
                        if 'N' in read.cigarstring:
                            newcigar = []
                            for (code, length) in read.cigar:
                                # 3 is N (intron), 2 is D (deletion)
                                if code == 3:
                                    code = 2
                                newcigar.append((code, length))
                            read.cigar = newcigar
                        fout.write(read)
                    print()

        if 'indexBAM' in steps:
            print('Make BAM index for freebayes')
            fn_newbam = f'{fdn_data}/gene_expression_vcf/{sample}/possorted_genome_bam_NtoD.bam'
            call = f'samtools index {fn_newbam}'
            sp.run(call, shell=True, check=True)

        if 'freebayes' in steps:
            print('Call SNV via freebayes')
            fn_newbam = f'{fdn_data}/gene_expression_vcf/{sample}/possorted_genome_bam_NtoD.bam'
            fn_snv = f'{fdn_data}/gene_expression_vcf/{sample}/freebayes_snv.vcf'
            fn_freebayes_genome = f'{fdn_transcriptome}/fasta/genome.fa'
            fn_freebayes_regions = f'{fdn_data}/gene_expression_vcf/genome_regions_freebayes.fa.fai'
            call = ('freebayes-parallel '
                    f'{fn_freebayes_regions} '
                    '32 ' # n CPUS
                    f'-f {fn_freebayes_genome} '
                    '-iXu ' # ignore indels
                    '-C 2 -q 1 ' # 2 reads, quality 1
                    f'{fn_newbam} '
                    f'> {fn_snv}'
                    )
            print(call)
            sp.run(call, shell=True, check=True)

        if 'count_alleles' in steps:
            fn_bam = f'{fdn_data}/gene_expression_bam/{sample}/{sample}/outs/possorted_genome_bam.bam'
            fn_barcodes = f'{fdn_data}/gene_expression_bam/{sample}/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
            fn_snv = f'{fdn_data}/gene_expression_vcf/{sample}/freebayes_snv.vcf'
            fn_scSplit_bin = f'{fdn_data}/gene_expression_snv_counts/scSplit'
            fn_common_snps = f'{fdn_data}/gene_expression_snv_counts/scSplit_common_snvs_hg38.tar.gz'
            fdn_alleles = f'{fdn_data}/gene_expression_snv_counts/{sample}'
            os.makedirs(fdn_alleles, exist_ok=True)
            call = (f'{fn_scSplit_bin} count '
                    f'-v {fn_snv} '
                    f'-i {fn_bam} '
                    f'-c {fn_common_snps} '
                    '-r ref_filtered.csv '
                    '-a alt_filtered.csv '
                    f'-o {fdn_alleles} '
                    f'-b {fn_barcodes}'
                    )
            print(call)
            sp.run(call, shell=True, check=True)

            #import vcf
            #from collections import Counter
            #pos = Counter()
            #with open(fn_snv, 'r') as f:
            #    vcffile = vcf.Reader(f)
            #    for ir, record in enumerate(vcffile):
            #        pos[record.CHROM+':'+str(record.POS)] += record.samples[0].data.DP
            #        if ir == 10000:
            #            break

        if 'demux' in steps:
            fn_scSplit_bin = f'{fdn_data}/gene_expression_snv_counts/scSplit'
            fdn_alleles = f'{fdn_data}/gene_expression_snv_counts/{sample}'
            fn_ref = f'{fdn_alleles}/ref_filtered.csv'
            fn_alt = f'{fdn_alleles}/alt_filtered.csv'
            fdn_demux = f'{fdn_data}/gene_expression_snv_demux/{sample}'
            os.makedirs(fdn_demux, exist_ok=True)
            call = (f'{fn_scSplit_bin} run '
                    f'-r {fn_ref} '
                    f'-a {fn_alt} '
                    '-n 2 '
                    f'-o {fdn_demux} '
                    )
            print(call)
            sp.run(call, shell=True, check=True)


