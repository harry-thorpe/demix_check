#!/usr/bin/env python

import os
import sys
import subprocess

import re
import pandas as pd
import numpy as np

from sketch import *

def check_mGEMS(mash_exec, seqtk_exec, t, ss, min_abun, ref_d, out_d, binned_reads_d, msweep_abun):
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
    ref_len_in="{}/ref_clu_comp.tsv".format(ref_d)
    
    abun_out="{}/msweep_abundances.tsv".format(out_d)
    abun_out_filt="{}/msweep_abundances_filt.tsv".format(out_d)

    clu_score_out="{}/clu_score.tsv".format(out_d)
    
    all_clu_msh_in_tmp=list()
    all_clu_msh_out="{}/all_clu.msh".format(out_d)
    all_clu_msh_dis_out="{}/all_clu_msh_dis.tsv.gz".format(out_d)
    
    abun_tmp=pd.read_csv(msweep_abun, sep="\t", comment="#", names=["cluster", "abundance"])
    abun_tmp.to_csv(abun_out, sep="\t", index=False)

    abun_tmp_filt=abun_tmp.query("abundance >= {}".format(min_abun))
    abun_tmp_filt.to_csv(abun_out_filt, sep="\t", index=False)
    
    ref_len=pd.read_csv(ref_len_in, sep="\t")
    ref_len=ref_len[["cluster", "length_ave"]]
    
    thr=pd.read_csv(ref_clu_thr, sep="\t")

    clu=pd.read_csv(abun_out_filt, sep="\t")
    clu_count=len(clu)

    clu_score=clu.merge(ref_len, how="left", on="cluster")

    clu_score_all=pd.DataFrame(columns=["cluster", "abundance", "length_ave", "score", "read_count", "total_bases", "coverage", "subsampled", "coverage_final", "notes"])
    
    sys.stderr.write("Found {} cluster/s in {}\n".format(clu_count, out_d))
    for i in range(clu_count):
        cluster=clu_score["cluster"].iloc[i]
        
        clu_score_tmp=clu_score[(clu_score.cluster == cluster)]
        clu_score_tmp=clu_score_tmp.reset_index()
        
        thr_tmp=thr[(thr.cluster == cluster)]
        
        # quick check that there is only one line per cluster
        if thr_tmp.shape[0] != 1 or clu_score_tmp.shape[0] != 1:
            sys.stderr.write('ERROR: file is in the wrong format\n')
            sys.exit(1)

        msh_out="{}/{}.msh".format(out_d_sketch, cluster)
        msh_dis_out="{}/{}_msh_dis.tsv.gz".format(out_d_sketch, cluster)
        msh_dis_clu_out="{}/{}_msh_dis_clu.tsv.gz".format(out_d_sketch, cluster)
        msh_scr_dis_out="{}/{}_msh_scr_dis.tsv.gz".format(out_d_sketch, cluster)
        msh_scr_dis_clu_out="{}/{}_msh_scr_dis_clu.tsv.gz".format(out_d_sketch, cluster)

        r1="{}/{}_1.fastq.gz".format(binned_reads_d, cluster)
        r2="{}/{}_2.fastq.gz".format(binned_reads_d, cluster)

        if os.path.isfile(r1) and os.path.isfile(r2):

            sys.stderr.write("Estimating coverage...\n")
            seqtk_cmd="{} fqchk {}".format(seqtk_exec, r1)
            std_result=subprocess.run(seqtk_cmd, shell=True, check=True, capture_output=True, text=True)
            
            subsampled=0
            total_bases=int(re.search(r'ALL\t(\d+)\t', std_result.stdout).group(1))*2
            read_len=float(re.search(r'avg_len:\s+([^;]+);', std_result.stdout).group(1))
            read_count=(total_bases/read_len)/2
            coverage=total_bases/clu_score_tmp["length_ave"].iloc[0]
            coverage_final=coverage
            notes=""

            if coverage > 100:
                sub_read_count=int(read_count/(coverage/100))
            
                sys.stderr.write("Coverage estimated to be {}x, so downsampling to 100x...\n".format(coverage))
                s_r1="{}/{}_1.fastq.gz".format(out_d_sketch, cluster)
                s_r2="{}/{}_2.fastq.gz".format(out_d_sketch, cluster)
            
                seqtk_cmd="{} sample -s 11 {} {} | gzip > {}".format(seqtk_exec, r1, sub_read_count, s_r1)
                std_result=subprocess.run(seqtk_cmd, shell=True, check=True, capture_output=True, text=True)
                seqtk_cmd="{} sample -s 11 {} {} | gzip > {}".format(seqtk_exec, r2, sub_read_count, s_r2)
                std_result=subprocess.run(seqtk_cmd, shell=True, check=True, capture_output=True, text=True)
                
                subsampled=1
                coverage_final=100
                notes="Warning: high coverage ({}) so {} reads were subsampled to reduce the coverage to 100".format(coverage, sub_read_count)

                r1=s_r1
                r2=s_r2
            elif coverage < 10:
                notes="Warning: low coverage ({}) so distances may not be accurate".format(coverage)

            r_f="{} {}".format(r1, r2)

            m=int(coverage_final/10)+2

            sys.stderr.write("Creating mash sketches...\n")
            std_result=run_mash_sketch(mash_exec, t, r_f, msh_out, ss, m, in_type="fq")
            log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
            
            all_clu_msh_in_tmp.append(msh_out)

            sys.stderr.write("Calculating mash distances...\n")
            std_result=run_mash_dist(mash_exec, t, ref_msh, msh_out, msh_dis_out)
            log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
            add_clusters(ref_clu, msh_dis_out, msh_dis_clu_out, ref=True, met=False)
            
            dis=pd.read_csv(msh_dis_clu_out, sep="\t", dtype={'ref_id': 'str', 'met_id': 'str'})
            dis=dis[(dis.ref_cluster == cluster)]
            
            sys.stderr.write("Calculating scores...\n")
            if min(dis["distance"]) < thr_tmp["dis_same_max"].iloc[0]:
                clu_score_tmp['score']=1
            elif min(dis["distance"]) < thr_tmp["threshold"].iloc[0]:
                clu_score_tmp['score']=2
            elif abs(np.median(dis["distance"]) - thr_tmp["dis_same_med_all"].iloc[0]) < abs(np.median(dis["distance"]) - thr_tmp["dis_diff_med_all"].iloc[0]):
                clu_score_tmp['score']=3
            else:
                clu_score_tmp['score']=4
            
            clu_score_tmp['read_count']=read_count
            clu_score_tmp['total_bases']=total_bases
            clu_score_tmp['coverage']=coverage
            clu_score_tmp['subsampled']=subsampled
            clu_score_tmp['coverage_final']=coverage_final
            clu_score_tmp['notes']=notes

        else:
            clu_score_tmp['score']=4

            clu_score_tmp['read_count']=None
            clu_score_tmp['total_bases']=None
            clu_score_tmp['coverage']=None
            clu_score_tmp['subsampled']=None
            clu_score_tmp['coverage_final']=None
            clu_score_tmp['notes']="Warning: no binned reads were found - probably a bad cluster assignment."

        clu_score_all=clu_score_all.append(clu_score_tmp, ignore_index=True)
    
    clu_score_all.to_csv(clu_score_out, sep="\t", index=False)
    
    if len(all_clu_msh_in_tmp) > 0:
        
        all_clu_msh_in=" ".join(all_clu_msh_in_tmp)

        std_result=run_mash_paste(mash_exec, all_clu_msh_in, all_clu_msh_out)
        log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
        std_result=run_mash_dist(mash_exec, t, all_clu_msh_out, all_clu_msh_out, all_clu_msh_dis_out)
        log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))

    sys.stderr.write("Analysis {} checking completed\n".format(binned_reads_d))

    log.close()

