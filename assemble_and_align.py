#!/usr/bin/env python3

import os
import sys
import pandas
import subprocess
import shutil
import argparse
import functools
from datetime import datetime



print = functools.partial(print, flush=True)

##########################################################################################################
#
#
def run_command(cmd: list):
    ext_cmd = [ '/usr/bin/time', '-v'] + cmd
    
    print(f"\n*************** SUBPROCESS STARTED: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ***************\n")
    
    p = subprocess.Popen(cmd)
    p.communicate()
    
    print(f"\n***************  SUBPROCESS FINISHED: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ***************\n")
    
    


##########################################################################################################
#
#
def splash(fastq_list: str, working_dir: str, anchors_tsv: str):
    
    cwd = os.getcwd() + '/'
    
    if not os.path.exists(anchors_tsv):
        print(f'Running SPLASH in {working_dir}...')
        
        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
            
            f_in = open(fastq_list)
            f_out = open(f'{working_dir}/samples.list', 'w')
            fastqs = f_in.read().splitlines()
            
            for entry in fastqs:
                head, fname = os.path.split(entry)
                f_out.write(f'{fname}\t{os.path.abspath(entry)}\n')
            
            f_in.close()
            f_out.close()
            
            
        os.chdir(working_dir)
           
        cmd = ['splash',
            '--outname_prefix', 'result',
            '--anchor_len', '27',
            '--gap_len', '0',
            '--target_len', '27',
            '--poly_ACGT_len', '8',
            '--dump_Cjs',
            '--max_pval_opt_for_Cjs', '0.05',
            '--with_effect_size_cts',
            '--with_pval_asymp_opt',
            '--n_threads_stage_1', '16',
            '--n_threads_stage_1_internal', '4',
            '--n_threads_stage_2', '32',
            '--n_bins', '128',
            '--n_most_freq_targets', '10',
            '--kmc_max_mem_GB', '12',
            '--dump_sample_anchor_target_count_txt',
            '--dump_sample_anchor_target_count_binary',
            '--keep_top_n_target_entropy', '100000',
            '--keep_top_n_effect_size_bin', '100000',
            '--keep_top_target_entropy_anchors_satc',
            '--keep_top_effect_size_bin_anchors_satc',
            '--without_compactors',
            '--exclude_postprocessing_item', 'postprocessing/blast.json',
            f'{working_dir}/samples.list']
        
        run_command(cmd)
        
        os.chdir(cwd)
        
        shutil.move(working_dir + 'result.after_correction.scores.tsv', anchors_tsv) 
                    
          
        
##########################################################################################################
#
#       
def lookup(
    in_anchors_tsv: str, 
    pos_index_list: str, 
    pos_index_k: int,
    neg_index_list: str,
    neg_index_k: int,
    working_dir: str,
    out_anchors_tsv: str):
    
    if not os.path.exists(working_dir):
            os.mkdir(working_dir)
            
    
    indices = {
        f'{working_dir}/positive.index': (pos_index_list, pos_index_k),
        f'{working_dir}/negative.index': (neg_index_list, neg_index_k)
    }
     

    for index, index_params in indices.items():
        
        if not os.path.exists(index):
            print(f'Building index {index} from {index_params[1]}\n')
            
            cmd = ['build_lookup_table.py',  
                '--poly_ACGT_len', '6', 
                '--kmer_len', f'{index_params[1]}',
                '--outname', index,
                index_params[0]]
                
            run_command(cmd)
        
        
        out = f'{index}.tsv' 
    
        if not os.path.exists(out):
            print(f'Querying {in_anchors_tsv} against {index}\n')
      
            cmd = ['lookup_table', 'query', 
                '--input_fmt', 'extendors', 
                '--stats_fmt', 'with_stats',
                '--report_fmt', 'empty',
                '--output_fmt', 'extendors',
                '--truncate_paths', 
                in_anchors_tsv,
                index, 
                out]
            
            run_command(cmd)
    
    
    pattern = '1: 0, 2: 0, 3: 0, 4: 0'
            
    if not os.path.exists(out_anchors_tsv):
        print(f'Generating {out_anchors_tsv}')
        
        f_in = open(in_anchors_tsv)
        f_out = open(out_anchors_tsv, 'w')
        
        f_neg = open(f'{working_dir}/negative.index.tsv')
        f_pos = open(f'{working_dir}/positive.index.tsv')
                
        n_neg_passed = 0        
        n_out = 0

        for i, (line_in, line_neg, line_pos) in enumerate(zip(f_in, f_neg, f_pos)):
            if i == 0:
                f_out.write('anchor\n') # copy header
                continue
            else:
                if pattern not in line_neg: 
                    continue
                
                n_neg_passed += 1                

                q = 0
                while True:
                    q = line_pos.find('1:', q)
                    if q == -1:
                        break

                    if line_pos[q: q + len(pattern)] != pattern:
                        n_out += 1
                        f_out.write(line_in[0: line_in.find('\t')] + '\n')
                        #print('HIT: ' + line_pos)
                        break
                    q += 1
            

            if i % 100000 == 0:
                print(f'{i}: {n_neg_passed} negative passed, {n_out} positive passed')
               #break
            
        f_in.close()
        f_pos.close()
        f_neg.close()
        f_out.close()


##########################################################################################################
#
#       
def compactors(fastq_list: str, anchors_tsv: str, working_dir: str, compactors_fasta: str):
    
    cwd = os.getcwd() + '/'
    
    if not os.path.exists(working_dir):
        print(f'Running compactors on independently on {fastq_list}')

        os.mkdir(working_dir)
      
        f_in = open(fastq_list)
        f_out = open(f'{working_dir}/samples.list', 'w')
        fastqs = f_in.read().splitlines()
       
        for entry in fastqs:
            f_out.write(f'{os.path.abspath(entry)}\n')
        
        f_in.close()
        f_out.close()
           
        os.chdir(working_dir)
                
        cmd = [
            'compactors',
            '--num_threads', '64',
            '--beta', '0.5',
            '--epsilon', '0.001',
            '--lower_bound', '2',
            '--min_extender_specificity', '0.8',
            '--num_extenders', '3',
            '--extenders_shift', '5',
            '--no_subcompactors',
            '--independent_outputs',
            'samples.list',
            anchors_tsv,
            'compactors.tsv',
           ]
        
        run_command(cmd)
        
        os.chdir(cwd)


    if not os.path.exists(compactors_fasta):

        compactors_fasta = open(compactors_fasta, 'w')
            
        print(f'Extracting compactors into FASTA')    

        for ie,entry in enumerate(os.listdir(f'{working_dir}')):
            if '.tsv' in entry:
                run_name = entry[0: entry.find('.fastq')] 
                
                print(f'{ie}: {run_name}...', end='')
                df = pandas.read_csv(f'{working_dir}/{entry}', sep='\t')
                col = df['compactor']

                for i,seq in enumerate(col):
                    compactors_fasta.write(f'>{run_name}-{i}\n')
                    compactors_fasta.write(f'{seq}\n')
                    
                print(f'{len(col)}')
                
        
        compactors_fasta.close()
        

##########################################################################################################
#
#    
def mmseqs(query_fasta: str, target_fasta: str, working_dir: str, out_tsv: str):
    
    query_db = f'{working_dir}/query-db'
    target_db = f'{working_dir}/target-db'
    out_db = f'{working_dir}/out.db'    
    
    
    if not os.path.exists(out_tsv):
    
        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
      
     
        if not os.path.exists(target_db):
            cmd = ['mmseqs',
            'createdb', 
            target_fasta,
            target_db]
            
            run_command(cmd)

        if not os.path.exists(query_db):
            cmd = ['mmseqs',
            'createdb', 
            query_fasta,
            query_db]
            
            run_command(cmd)
             
           
        if not os.path.exists(out_db):
            cmd = ['mmseqs',
            'search',
            '--alignment-mode', '3',
            query_db,
            target_db,
            out_db,
            'tmp'
            ] 

            run_command(cmd)        

        
         
        cmd = ['mmseqs',
        'convertalis', 
        query_db,
        target_db,
        out_db,
        out_tsv
        ] 

        run_command(cmd)           
    

##########################################################################################################
#
#
if __name__ == '__main__':
     
    parser = argparse.ArgumentParser(
        prog = "assemble_and_align.py",
        description=f'Generate compactors and align against database\n\n',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    
    parser.add_argument('fastq_list', help='list of FASTQ files to analyze')
    parser.add_argument('target_db_fasta', help='FASTA file with sequences to be used as a target database')
    
    parser.add_argument('--out_dir', default='out-align', type=str, help='output directory with all the results')
    
    parser.add_argument('--positive_index_list', default=None, type=str, help='list of FASTA files for positive index')
    parser.add_argument('--positive_index_k', default=18, type=int, help='k-mer length for positive index')
  
    parser.add_argument('--negative_index_list', default=None, type=str, help='list of FASTA files for negative index')
    parser.add_argument('--negative_index_k', default=18, type=int, help='k-mer length for negative index')
    
    
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)
        
    args=parser.parse_args()     
    
    if bool(args.positive_index_list) ^ bool(args.negative_index_list):
        print('Arguments --positive_index_list and --negative_index_list must be given together')
        sys.exit(1)
    
        
  
    # output data
    out_dir = os.path.abspath(args.out_dir)
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    anchors_tsv = f'{out_dir}/anchors.scores.tsv'
    filtered_anchors_tsv = anchors_tsv
    compactors_fasta = f'{out_dir}/compactors.fasta'
    alignments_tsv = f'{out_dir}/alignments.tsv'
   
    # steps
    print('Identifying anchors...')
    splash(args.fastq_list, f'{out_dir}/splash/', anchors_tsv)
    
    if args.positive_index_list is not None:
        filtered_anchors_tsv = f'{out_dir}/anchors.filtered.tsv'
        print('Filtering anchors using lookup tables...')
        lookup(anchors_tsv, args.positive_index_list, args.positive_index_k, args.negative_index_list, args.negative_index_k, f'{out_dir}/lookup/', filtered_anchors_tsv)
    
    print('Generating compactors...')
    compactors(args.fastq_list, filtered_anchors_tsv, f'{out_dir}/compactors/', compactors_fasta)
    
    print('Aligning compactors against database...')
    mmseqs(compactors_fasta, args.target_db_fasta, f'{out_dir}/mmseqs/', alignments_tsv)
         
         
