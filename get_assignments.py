#!/usr/bin/env python

import os
import sys
import subprocess
import shutil

import re
import pandas as pd
import numpy as np

def run_mGEMS(themisto_align_exec, mSWEEP_exec, mGEMS_exec, t, min_abun, r1, r2, ref_d, out_d, out_d_bin, msweep_abun, keep):
    sys.stderr.write("Running mSWEEP/mGEMS pipeline on {},{} against reference set {}...\n".format(r1, r2, ref_d))

    ref_idx_d="{}/ref_idx".format(ref_d)
    ref_idx="{}/ref_idx/ref_idx".format(ref_d)
    ref=os.path.basename(ref_d)
    out=os.path.basename(out_d)
    tmp_out_d="{}/tmp".format(out_d)

    if not os.path.isdir(out_d):
        os.makedirs(out_d)
    if not os.path.isdir(out_d_bin):
        os.makedirs(out_d_bin)
    if not os.path.isdir(tmp_out_d):
        os.makedirs(tmp_out_d)

    log=open("{}/run.log".format(out_d), 'w')

    ref_clu="{}/ref_clu.txt".format(ref_d)
    
    r_ali1="{}/ali_1.aln.gz".format(out_d)
    r_ali2="{}/ali_2.aln.gz".format(out_d)
    r_ali1_p=re.sub(r'.gz$', '', r_ali1)
    r_ali2_p=re.sub(r'.gz$', '', r_ali2)
    msweep_abun_p=re.sub(r'_abundances.txt$', '', msweep_abun)
    msweep_abun_prob_t="{}/msweep_probs.tsv.gz".format(out_d)
    msweep_abun_prob_c="{}/msweep_probs.csv.gz".format(out_d)

    sys.stderr.write("Pseudoaligning r1 reads with themisto...\n")
    themisto_cmd="{} --index-prefix {} --query-file {} --outfile {} --rc --temp-dir {} --n-threads {} --sort-output --gzip-output".format(themisto_align_exec, ref_idx, r1, r_ali1_p, tmp_out_d, t)
    std_result=subprocess.run(themisto_cmd, shell=True, check=True, capture_output=True, text=True)
    log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
    
    sys.stderr.write("Pseudoaligning r2 reads with themisto...\n")
    themisto_cmd="{} --index-prefix {} --query-file {} --outfile {} --rc --temp-dir {} --n-threads {} --sort-output --gzip-output".format(themisto_align_exec, ref_idx, r2, r_ali2_p, tmp_out_d, t)
    std_result=subprocess.run(themisto_cmd, shell=True, check=True, capture_output=True, text=True)
    log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
    
    sys.stderr.write("Estimating abundances with mSWEEP...\n")
    mSWEEP_cmd="{} -t {} --themisto-1 {} --themisto-2 {} -o {} -i {} --write-probs --compress z".format(mSWEEP_exec, t, r_ali1, r_ali2, msweep_abun_p, ref_clu)
    std_result=subprocess.run(mSWEEP_cmd, shell=True, check=True, capture_output=True, text=True)
    log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
    
    # convert tsv to csv

    msweep_abun_prob_tmp=pd.read_csv(msweep_abun_prob_t, sep="\t")
    msweep_abun_prob_tmp.to_csv(msweep_abun_prob_c, index=False)

    abun_tmp=pd.read_csv(msweep_abun, sep="\t", comment="#", names=["cluster", "abundance"], dtype={'cluster': 'str'})
    abun_tmp=abun_tmp.query("abundance >= {}".format(min_abun))

    clu_to_bin=abun_tmp["cluster"].to_list()
    clu_to_ext=["{}/".format(out_d_bin) + x for x in clu_to_bin]
    clu_to_ext=[x + ".bin" for x in clu_to_ext]

    clu_to_bin=",".join(clu_to_bin)

    sys.stderr.write("Binning reads with mGEMS...\n")
    mGEMS_cmd="{} bin --groups {} --themisto-alns {},{} -o {} --probs {} -a {} --index {} --min-abundance {} -i {}".format(mGEMS_exec, clu_to_bin, r_ali1, r_ali2, out_d_bin, msweep_abun_prob_c, msweep_abun, ref_idx_d, min_abun, ref_clu)
    std_result=subprocess.run(mGEMS_cmd, shell=True, check=True, capture_output=True, text=True)
    log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
    
    for clu_bin_ext in clu_to_ext:
        if os.path.isfile(clu_bin_ext):
            sys.stderr.write("Extracting reads with mGEMS...\n")
            mGEMS_cmd="{} extract --bins {} -r {},{} -o {}".format(mGEMS_exec, clu_bin_ext, r1, r2, out_d_bin)
            std_result=subprocess.run(mGEMS_cmd, shell=True, check=True, capture_output=True, text=True)
            log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
    
    if not keep:
        os.remove(msweep_abun_prob_c)
        os.remove(msweep_abun_prob_t)
        os.remove(r_ali1)
        os.remove(r_ali2)
        shutil.rmtree(tmp_out_d)
    
    sys.stderr.write("mSWEEP/mGEMS pipeline on {},{} completed\n".format(r1, r2))

    log.close()

