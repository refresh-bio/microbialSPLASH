#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys
import numpy
import multiprocessing as mp


def download_srr(srr: str):
  
    print(srr)
    cmd = ['prefetch', srr]
    proc = subprocess.Popen(cmd) 
    ec = proc.wait()
    
    if ec != 0: return ec
    
    cmd = ['fastq-dump', '--split-3', '--skip-technical', srr]
    proc = subprocess.Popen(cmd) 
    ec = proc.wait()
    
    beg = end
  
    return ec

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        prog = "srr-download.py",
        description='Welcome to Parallel SRR downloader',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    
    parser.add_argument('srr_list', help='list of SRR accessions')
    parser.add_argument('out_dir', help='output directory where to save files') 
    
    parser.add_argument('--n_threads', default=1, type=int, help='number of working threads')
    parser.add_argument('--n_random', default=0, type=int, help='select given number of SRRs randomly, 0 = take all')
    parser.add_argument('--n_first', default=0, type=int, help='select given number of SRRs from the beginning; 0 = take all ')
  
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)
   
    args=parser.parse_args()
    
    file = open(args.srr_list)
    runs = file.read().splitlines()
   
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    
  
    if args.n_random > 0:
        print(f'Downloading {args.n_random} randomly selected SRRs from {len(runs)} accessions listed in {args.srr_list}')
        indices = numpy.random.choice(len(runs), size=args.n_random, replace=False)
    elif args.n_first > 0:
        print(f'Downloading {args.n_first} first SRRs from {len(runs)} accessions listed in {args.srr_list}')
        indices = range(0,args.n_first)
    else:
        print(f'Downloading all SRRs from {len(runs)} accessions listed in {args.srr_list}')
        indices = range(0,len(runs))
   
    cwd = os.getcwd()
    os.chdir(args.out_dir)
    
    with mp.Pool(processes=args.n_threads) as pool: # pool only visible in the context below
        results = pool.starmap(download_srr, ((runs[i],) for i in indices))
        
        
    os.chdir(cwd)       


        
    