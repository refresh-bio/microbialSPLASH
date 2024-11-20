#!/usr/bin/env python3

import pandas
import numpy
import os
import sys
import pickle
import math
import argparse
import functools
from pathlib import Path

print = functools.partial(print, flush=True)

N_UNIPROT_SYMBOLS = 88263858984 # this is to calculate normalized e-value with respect to the entire UniProtKB
MAX_OVERLAP = 0.5
INVALID_ACCESSIONS = {'Q9FS29', 'Q1KLZ2', 'Q76MX2', 'Q1KLZ1', 'P18503'}


########################################################################################################################
def is_cas(acc, c):
    return (c[0:3] == 'CAS') and (c[3:4].isnumeric()) and (acc not in INVALID_ACCESSIONS)

########################################################################################################################
def load_compactors_to_multirows(df: pandas.DataFrame, subdir: str, accessions2proteins: dict):
    file = f'{subdir}/compactors_to_multi_rows.pickle'

    if not os.path.exists(file):
        runs = set()

        queries = df['query']
        q_start = df['query_start']
        q_end = df['query_end']
        accessions = df['target']
        evals = df['e_value']

        print(f'Searching best hits for each compactor rows')
        compactors_to_rows = {}
        for rid in range(len(df)):
            if (rid > 0) and (rid % 1000000 == 0):
                print(f'\r{rid}', end='')

            q = queries[rid]
            acc = accessions[rid]
            t = accessions2proteins[acc][0]

            if (is_cas(acc, t)):
                ev = evals[rid]
                if q_start[rid] < q_end[rid]:
                    qs, qe = q_start[rid], q_end[rid]
                else:
                    qs, qe = q_end[rid], q_start[rid]

                if q not in compactors_to_rows:
                    compactors_to_rows[q] = {t : (rid, ev, qs, qe, acc)}
                else:
                    best_hits = compactors_to_rows[q]
                    if (t not in best_hits) or (ev < best_hits[t][1]):
                        best_hits[t] = (rid, ev, qs, qe, acc)

        print('done')

        handle = open(file, 'wb')
        pickle.dump(compactors_to_rows, handle)
        handle.close()
    else:
        print(f'Loading best hits for each compactor from {file}')
        handle = open(file, 'rb')
        compactors_to_rows = pickle.load(handle)
        handle.close()

    print(f'Compactors count: {len(compactors_to_rows)}')

    return compactors_to_rows

########################################################################################################################
def filter_overlaps(compactors_to_best_multirows: dict):

    n_singletons = 0
    n_multiple = 0
    n_total_hits = 0
    

    for compactor, best_hits in compactors_to_best_rows.items():

        hits_list = [*best_hits.items()]
        hits_list = sorted(hits_list, key= lambda x:x[1][1]) # sort by e-value

        for i, current_hit  in enumerate(hits_list):
            for j in range(i):
                ref_hit = hits_list[j]
                if (ref_hit[1] is not None) and (current_hit[1][2] <= ref_hit[1][3]) and (current_hit[1][3] >= ref_hit[1][2]):
                   
                    overlap_len = min(current_hit[1][3] - ref_hit[1][2], ref_hit[1][3] - current_hit[1][2])

                    if (overlap_len > MAX_OVERLAP * (ref_hit[1][3] - ref_hit[1][2]) or overlap_len > MAX_OVERLAP * (current_hit[1][3] - current_hit[1][2])):
                        hits_list[i] = (hits_list[i][0], None)
                    else:
                        pass

        for i in range(len(hits_list)):
            if hits_list[i][1] == None:
                best_hits.pop(hits_list[i][0])

        n_total_hits += len(best_hits)
        if len(best_hits) == 1:
            n_singletons += 1
        else:
            n_multiple += 1

    print('Filtering overlapping alignments with greedy approach\n' +
            f'Number of aligning compactors: {len(compactors_to_best_multirows)},\n' +
            f'\twith single CAS targets: {n_singletons},\n' +
            f'\twith multiple CAS targets: {n_multiple},\n' +
            f'Total number of hits: {n_total_hits}\n')

    return n_total_hits




########################################################################################################################
def load_multitargets(df: pandas.DataFrame,  subdir: str, compactors_to_best_multirows: dict):
    file = f'{subdir}/multitargets.pickle'

    if not os.path.exists(file):
        print('Identifying targets')
        targets = set()
        for compactor, best_hits in compactors_to_best_rows.items():
            targets.update(best_hits.keys())

        targets = sorted([*targets])
        handle = open(file, 'wb')
        pickle.dump(targets, handle)
        handle.close()
    else:
        print(f'Loading targets from {file}')
        handle = open(file, 'rb')
        targets = pickle.load(handle)
        handle.close()

    print(f'Targets count: {len(targets)}')
    return targets

########################################################################################################################
def map_accessions(fasta_handle):
    
    map = {}
    lines = fasta_handle.readlines()
    n_symbols = 0

    for line in lines:
        if line[0] == '>':
            a = line.find('|')
            a += 1
            b = line.find('|', a)
            accession = line[a:b]

            b += 1
            c = line.find('_', b)

            protein = line[b:c]

            c += 1
            d = line.find(' ')
            organism = line[c:d]

            map[accession] = (protein, organism)
        else:
            n_symbols += len(line)

    return map, n_symbols



########################################################################################################################
if __name__ == '__main__':
     
    parser = argparse.ArgumentParser(
        prog = "extract_cas.py",
        description=f'Extract CAS-targeted alignments\n\n',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    
    parser.add_argument('fastq_list', help='list of FASTQ files to analyze')
    parser.add_argument('target_db_fasta', help='FASTA file with sequences to be used as a target database')
    parser.add_argument('alignments_tsv', help='input TSV table with alignments of compactors against the target database produced by assemble_and_align.py')
    parser.add_argument('compactors_fasta', help='input FASTA with compactors sequences produced by assemble_and_align.py')
    
    parser.add_argument('--out_dir', default='out-cas', type=str, help='output directory with results of analyzes')
   
    if len(sys.argv) < 5:
        parser.print_help()
        sys.exit(1)
             
    args=parser.parse_args() 

    # output data
    out_dir = os.path.abspath(args.out_dir)
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    cas_alignments_tsv = f'{out_dir}/cas-alignments.tsv'
    cas_compactors_fasta = f'{out_dir}/cas-compactors.fasta'     
    
    print('Mapping accessions to proteins')
    f = open(args.target_db_fasta)
    acc2prot, n_symbols = map_accessions(f)
    f.close()
 
    if not os.path.exists(cas_alignments_tsv):

        norm_factor = N_UNIPROT_SYMBOLS / n_symbols

        print('Loading runs...')
        f = open(args.fastq_list)
        lines = f.readlines()
        runs = sorted([Path(p).stem for p in lines])
        runs_to_rows = {run: rid for rid,run in enumerate(runs)}

        print('Loading data frame...')
        df = pandas.read_csv(args.alignments_tsv, sep='\t', header=None,
                             names=['query','target','identity','alignment_length','num_mismatches','num_gapopen',
                                    'query_start','query_end','target_start','target_end','e_value', 'bitscore'])

        print(f'{len(df)} rows loaded')

        compactors_to_best_rows = load_compactors_to_multirows(df, out_dir, acc2prot)

        n_total_hits = filter_overlaps(compactors_to_best_rows)
        targets = load_multitargets(df, out_dir, compactors_to_best_rows)

        targets_to_cols = {t: cid for cid, t in enumerate(targets)}

        hits_rows = []
        hits_target_protein = []
        hits_target_organism = []
        hits_norm_evalue = []


        print('Filling final tables by iterating over compactors')
        for cid, (compactor, best_hits) in enumerate(compactors_to_best_rows.items()):
            if (cid > 0) and (cid % 1000 == 0):
                print(f'\r{cid}', end='')

            run = compactor[:compactor.find('-')]
            run_id = runs_to_rows[run]

            for target, (rid, ev, qstart, qend, acc) in best_hits.items():
                col_id = targets_to_cols[target]
                hits_rows.append(rid)
                hits_target_protein.append(target)
                hits_target_organism.append(acc2prot[acc][1])
                hits_norm_evalue.append(ev * norm_factor)

               
        sub_df = df.iloc[hits_rows, :]
        sub_df = sub_df.assign(
            target_protein=hits_target_protein,
            target_organism=hits_target_organism,
            norm_e_value=hits_norm_evalue
        )

        sub_df.to_csv(cas_alignments_tsv, sep='\t', index=False)
        
     
    if not os.path.exists(cas_compactors_fasta):
   
        df = pandas.read_csv(cas_alignments_tsv, sep='\t')

        compactors_set = set(df['query'])

        in_fasta = open(args.compactors_fasta)
        out_fasta = open(cas_compactors_fasta, 'wt')

        print('Analyzing compactors FASTA...')

        for rid,(header, seq) in enumerate(zip(in_fasta, in_fasta)):
            if (rid > 0) and (rid % 10000 == 0):
                print(f'\r{rid}', end='')

            q = header[1:-1]
            if q in compactors_set:
                out_fasta.write(header)
                out_fasta.write(seq)

        out_fasta.close()
        in_fasta.close()

