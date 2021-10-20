'''
Author:	        Fabio Zanini
Date:	        2021/02/03
Description:	Run CITE-seq-Count on Avani's pilot CITE-Seq experiment.
'''
import os
import sys
import argparse
import pandas as pd
import subprocess as sp


if __name__ == '__main__':

    pa = argparse.ArgumentParser()
    pa.add_argument('--samples', default='iMSC,AdMSC')
    args = pa.parse_args()

    samples = args.samples.split(',')
    
    fdn_data = '../../data/sequencing/experiments/202012_pilot1'
    
    sample_suffixes = {'AdMSC': 'S1', 'iMSC': 'S2'}
    for sample in samples:
        print(sample)
        fdn_in = f'{fdn_data}/Ab_raw/201116_A00152_0332_AHNCV5DSXY/201116_A00152_0332_AHNCV5DSXY/fastq_path/CITE/{sample}'
        fdn_out = f'{fdn_data}/Ab_counts/{sample}'
        os.makedirs(fdn_out, exist_ok=True)
    
        suffix = sample_suffixes[sample]
        fns = {
            'R1': [f"{fdn_in}/{sample}_{suffix}_L00{i}_R1_001.fastq.gz" for i in range(1, 5)],
            'R2': [f"{fdn_in}/{sample}_{suffix}_L00{i}_R2_001.fastq.gz" for i in range(1, 5)],
        }
        
        fnr1 = ','.join(fns['R1'])
        fnr2 = ','.join(fns['R2'])
    
        fn_tags = f'{fdn_data}/Ab_params/Avanis_CITEseq_Antibodies1.csv'
        fn_tags_new = f'{fdn_data}/Ab_params/Avanis_CITEseq_Antibodies1_for_pipeline.csv'
        df = pd.read_csv(fn_tags, sep=',')
    
        # At least one Ab has two clones in there (sigh!)
        df['unique_id'] = df['Marker'] + '; ' + df['Clone']
        df[['Barcode', 'unique_id']].to_csv(fn_tags_new, sep=',', index=False, header=None)
    
        fn_whitelist = f'{fdn_data}/Ab_params/barcodes_{sample}.csv'
    
        call = (f'CITE-seq-Count -R1 "{fnr1}" -R2 "{fnr2}" -t {fn_tags_new} '
               #'-n 100 ' #FIXME
               '-cbf 1 '
               '-cbl 16 '
               '-umif 17 '
               '-umil 28 '
               '--max-error 3 '
               '-T 32 '
               '--dense '
               '--expected_cells 8000 '
               '--bc_collapsing_dist 1 '
               '--umi_collapsing_dist 2 '
               f'--whitelist {fn_whitelist} '
               f'-o "{fdn_out}"')
    
        print(call)
        sp.run(call, shell=True, check=True)
