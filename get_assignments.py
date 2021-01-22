#!/usr/bin/env python

import os
import sys
import subprocess

import re
import pandas as pd
import numpy as np

def run_mGEMS(themisto_align_exec, mSWEEP_exec, mGEMS_exec, t, min_abun, r1, r2, ref_d, out_d, out_d_bin, msweep_abun):
    
    ref_idx="{}/ref_idx".format(ref_d)
    ref=os.path.basename(ref_d)
    out=os.path.basename(out_d)

    if not os.path.isdir(out_d):
        os.makedirs(out_d)
    if not os.path.isdir(out_d_bin):
        os.makedirs(out_d_bin)

    ref_clu="{}/ref_clu.txt".format(ref_d)
    
    r_ali1="{}/ali_1.aln.gz".format(out_d)
    r_ali2="{}/ali_2.aln.gz".format(out_d)
    r_ali1_p=re.sub(r'.gz$', '', r_ali1)
    r_ali2_p=re.sub(r'.gz$', '', r_ali2)
    msweep_abun_p=re.sub(r'_abundances.txt$', '', msweep_abun)
    msweep_abun_prob="{}/msweep_probs.csv".format(out_d)

    themisto_cmd="{} --index-dir {} --query-file {} --outfile {} --rc --temp-dir {}/tmp --n-threads {} --sort-output --gzip-output".format(themisto_align_exec, ref_idx, r1, r_ali1_p, ref_idx, t)
    subprocess.run(themisto_cmd, shell=True, check=True)

    themisto_cmd="{} --index-dir {} --query-file {} --outfile {} --rc --temp-dir {}/tmp --n-threads {} --sort-output --gzip-output".format(themisto_align_exec, ref_idx, r2, r_ali2_p, ref_idx, t)
    subprocess.run(themisto_cmd, shell=True, check=True)
    
    mSWEEP_cmd="{} -t 1 --themisto-1 {} --themisto-2 {} -o {} -i {} --write-probs".format(mSWEEP_exec, r_ali1, r_ali2, msweep_abun_p, ref_clu)
    subprocess.run(mSWEEP_cmd, shell=True, check=True)

    mGEMS_cmd="{} -r {},{} --themisto-alns {},{} -o {} --probs {} -a {} --index {} --min-abundance {}".format(mGEMS_exec, r1, r2, r_ali1, r_ali2, out_d_bin, msweep_abun_prob, msweep_abun, ref_idx, min_abun)
    subprocess.run(mGEMS_cmd, shell=True, check=True)

