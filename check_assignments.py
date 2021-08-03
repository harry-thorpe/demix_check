#!/usr/bin/env python

import os
import sys
import subprocess

import re
import pandas as pd
import numpy as np

from sketch import *

def check_mGEMS(mash_exec, t, ss, m, min_abun, ref_d, out_d, binned_reads_d, msweep_abun):
    sys.stderr.write("Checking analysis {} against reference set {}...\n".format(binned_reads_d, ref_d))

    ref_idx="{}/ref_idx".format(ref_d)
    ref_clu="{}/ref_clu.tsv".format(ref_d)
    ref=os.path.basename(ref_d)
    out=os.path.basename(out_d)
    out_d_sketch="{}/binned_reads_sketches".format(out_d)
    
    if not os.path.isdir(out_d):
        os.makedirs(out_d)
    if not os.path.isdir(out_d_sketch):
        os.makedirs(out_d_sketch)
    
    log=open("{}/check.log".format(out_d), 'w')

    ref_msh="{}/ref.msh".format(ref_d)
    ref_clu_thr="{}/ref_clu_thr.tsv".format(ref_d)
    
    abun_out="{}/msweep_abundances.tsv".format(out_d)
    abun_out_filt="{}/msweep_abundances_filt.tsv".format(out_d)
    
    clu_score_out="{}/clu_score.tsv".format(out_d)
    
    abun_tmp=pd.read_csv(msweep_abun, sep="\t", comment="#", names=["cluster", "abundance"])
    abun_tmp.to_csv(abun_out, sep="\t", index=False)

    abun_tmp_filt=abun_tmp.query("abundance >= {}".format(min_abun))
    abun_tmp_filt.to_csv(abun_out_filt, sep="\t", index=False)

    clu=pd.read_csv(abun_out_filt, sep="\t")
    clu_count=len(clu)
    clu_score=clu
    clu_score["score"]=0
    
    sys.stderr.write("Found {} cluster/s in {}\n".format(clu_count, binned_reads_d))
    for i in range(clu_count):
        cluster=clu["cluster"][i]

        r1="{}/{}_1.fastq.gz".format(binned_reads_d, cluster)
        r2="{}/{}_2.fastq.gz".format(binned_reads_d, cluster)
        r_f="{} {}".format(r1, r2)

        msh_out="{}/{}.msh".format(out_d_sketch, cluster)
        msh_dis_out="{}/{}_msh_dis.tsv.gz".format(out_d_sketch, cluster)
        msh_dis_clu_out="{}/{}_msh_dis_clu.tsv.gz".format(out_d_sketch, cluster)
        msh_scr_dis_out="{}/{}_msh_scr_dis.tsv.gz".format(out_d_sketch, cluster)
        msh_scr_dis_clu_out="{}/{}_msh_scr_dis_clu.tsv.gz".format(out_d_sketch, cluster)
        
        sys.stderr.write("Creating mash sketches...\n")
        std_result=run_mash_sketch(mash_exec, t, r_f, msh_out, ss, m, in_type="fq")
        log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
        
        sys.stderr.write("Calculating mash distances...\n")
        std_result=run_mash_dist(mash_exec, t, ref_msh, msh_out, msh_dis_out)
        log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
        add_clusters(ref_clu, msh_dis_out, msh_dis_clu_out, ref=True, met=False)
        
        #std_result=run_mash_screen(mash_exec, t, ref_msh, r_f, msh_scr_dis_out, cluster)
        #log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
        #add_clusters(ref_clu, msh_scr_dis_out, msh_scr_dis_clu_out, ref=True, met=False)
        
        dis=pd.read_csv(msh_dis_clu_out, sep="\t", dtype={'ref_id': 'str', 'met_id': 'str'})
        dis=dis[(dis.ref_cluster == cluster)]

        thr=pd.read_csv(ref_clu_thr, sep="\t")
        thr=thr[(thr.cluster == cluster)]
        
        # quick check that there is only one threshold line per cluster
        if thr.shape[0] != 1:
            sys.stderr.write('ERROR: threshold file is in the wrong format\n')
            sys.exit(1)
        
        sys.stderr.write("Calculating scores...\n")
        if min(dis["distance"]) < max(thr["dis_same_max"]):
            clu_score['score']=np.where(clu_score['cluster'] == cluster, 1, clu_score['score'])
        elif min(dis["distance"]) < max(thr["threshold"]):
            clu_score['score']=np.where(clu_score['cluster'] == cluster, 2, clu_score['score'])
        elif abs(np.median(dis["distance"]) - max(thr["dis_same_med_all"])) < abs(np.median(dis["distance"]) - max(thr["dis_diff_med_all"])):
            clu_score['score']=np.where(clu_score['cluster'] == cluster, 3, clu_score['score'])
        else:
            clu_score['score']=np.where(clu_score['cluster'] == cluster, 4, clu_score['score'])

    clu_score.to_csv(clu_score_out, sep="\t", index=False)
    
    sys.stderr.write("Analysis {} checking completed\n".format(binned_reads_d))

    log.close()

