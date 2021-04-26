#!/usr/bin/env python

import os
import sys
import subprocess

import re
import pandas as pd
import numpy as np

def run_mGEMS(themisto_align_exec, mSWEEP_exec, mGEMS_exec, t, min_abun, r1, r2, ref_d, out_d, out_d_bin, msweep_abun):
    sys.stderr.write("Running mSWEEP/mGEMS pipeline on {},{} against reference set {}...\n".format(r1, r2, ref_d))

    ref_idx="{}/ref_idx".format(ref_d)
    ref=os.path.basename(ref_d)
    out=os.path.basename(out_d)

    if not os.path.isdir(out_d):
        os.makedirs(out_d)
    if not os.path.isdir(out_d_bin):
        os.makedirs(out_d_bin)
    
    log=open("{}/run.log".format(out_d), 'w')

    ref_clu="{}/ref_clu.txt".format(ref_d)
    
    r_ali1="{}/ali_1.aln.gz".format(out_d)
    r_ali2="{}/ali_2.aln.gz".format(out_d)
    r_ali1_p=re.sub(r'.gz$', '', r_ali1)
    r_ali2_p=re.sub(r'.gz$', '', r_ali2)
    msweep_abun_p=re.sub(r'_abundances.txt$', '', msweep_abun)
    msweep_abun_prob="{}/msweep_probs.csv".format(out_d)
    
    sys.stderr.write("Pseudoaligning r1 reads with themisto...\n")
    themisto_cmd="{} --index-dir {} --query-file {} --outfile {} --rc --temp-dir {}/tmp --n-threads {} --sort-output --gzip-output".format(themisto_align_exec, ref_idx, r1, r_ali1_p, ref_idx, t)
    std_result=subprocess.run(themisto_cmd, shell=True, check=True, capture_output=True, text=True)
    log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
    
    sys.stderr.write("Pseudoaligning r2 reads with themisto...\n")
    themisto_cmd="{} --index-dir {} --query-file {} --outfile {} --rc --temp-dir {}/tmp --n-threads {} --sort-output --gzip-output".format(themisto_align_exec, ref_idx, r2, r_ali2_p, ref_idx, t)
    std_result=subprocess.run(themisto_cmd, shell=True, check=True, capture_output=True, text=True)
    log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
    
    sys.stderr.write("Estimating abundances with mSWEEP...\n")
    mSWEEP_cmd="{} -t 1 --themisto-1 {} --themisto-2 {} -o {} -i {} --write-probs".format(mSWEEP_exec, r_ali1, r_ali2, msweep_abun_p, ref_clu)
    std_result=subprocess.run(mSWEEP_cmd, shell=True, check=True, capture_output=True, text=True)
    log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
    
    abun_tmp=pd.read_csv(msweep_abun, sep="\t", comment="#", names=["cluster", "abundance"])
    abun_tmp=abun_tmp.query("abundance >= {}".format(min_abun))

    clus_to_bin=abun_tmp["cluster"].to_list()
    clus_to_ext=["{}/".format(out_d_bin) + x for x in clus_to_bin]
    clus_to_ext=[x + ".bin" for x in clus_to_ext]

    clus_to_bin=",".join(clus_to_bin)
    clus_to_ext=",".join(clus_to_ext)

    sys.stderr.write("Binning reads with mGEMS...\n")
    mGEMS_cmd="{} bin --groups {} --themisto-alns {},{} -o {} --probs {} -a {} --index {} --min-abundance {}".format(mGEMS_exec, clus_to_bin, r_ali1, r_ali2, out_d_bin, msweep_abun_prob, msweep_abun, ref_idx, min_abun)
    std_result=subprocess.run(mGEMS_cmd, shell=True, check=True, capture_output=True, text=True)
    log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
    
    sys.stderr.write("Extracting reads with mGEMS...\n")
    mGEMS_cmd="{} extract --bins {} -r {},{} -o {}".format(mGEMS_exec, clus_to_ext, r1, r2, out_d_bin)
    std_result=subprocess.run(mGEMS_cmd, shell=True, check=True, capture_output=True, text=True)
    log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
    
    sys.stderr.write("mSWEEP/mGEMS pipeline on {},{} completed\n".format(r1, r2))

    log.close()

