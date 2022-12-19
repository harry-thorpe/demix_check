#!/usr/bin/env python

import os
import sys
import subprocess
import gzip

import re
import pandas as pd
import numpy as np

from sketch import *

def setup_reference(mash_exec, themisto_build_exec, seqtk_exec, d, t, ss, thr_prop_min, thr_abs_min, thr_prop_exp, redo_thr, no_build_index, no_build_fasta):
    sys.stderr.write("Setting up reference set {}...\n".format(d))

    seq_info_f="{}/ref_info.tsv".format(d)
    seq_info=pd.read_csv(seq_info_f, sep="\t", dtype={'id': 'str'})
    seq_count=len(seq_info)
    
    fa_out="{}/ref.fasta.gz".format(d)
    msh_out="{}/ref.msh".format(d)
    msh_dis_out="{}/ref_msh_dis.tsv.gz".format(d)
    msh_dis_clu_out="{}/ref_msh_dis_clu.tsv.gz".format(d)
    clu_thr_out="{}/ref_clu_thr.tsv".format(d)
    clu_out_m="{}/ref_clu.txt".format(d)
    clu_out="{}/ref_clu.tsv".format(d)
    comp_out="{}/ref_comp.tsv".format(d)
    clu_comp_out="{}/ref_clu_comp.tsv".format(d)
    paths_out="{}/ref_paths.txt".format(d)

    if redo_thr == True:
        
        sys.stderr.write("Calculating thresholds...\n")
        get_thresholds(clu_out, msh_dis_clu_out, thr_prop_min, thr_abs_min, thr_prop_exp, clu_thr_out)
    
    else:
        log=open("{}/ref.log".format(d), 'w')
        
        if not no_build_fasta:
            sys.stderr.write("Creating multifasta...\n")
            fo=gzip.open(fa_out, 'wt')
            for i in range(seq_count):
                seq_id=seq_info["id"][i]
                cluster=seq_info["cluster"][i]
                assembly=seq_info["assembly"][i]

                fo.write(">{}\n".format(seq_id))

                if assembly.endswith('.gz'):
                    f_tmp=gzip.open(assembly, 'rt')
                else:
                    f_tmp=open(assembly, 'r')
                    for l in f_tmp:
                        if not re.match(r'^>', l):
                            fo.write(l)
                    f_tmp.close()
            fo.close()

        clu_m=seq_info[["cluster"]]
        clu_m.to_csv(clu_out_m, sep="\t", index=False, header=None)

        clu=seq_info[["id", "cluster"]]
        clu.to_csv(clu_out, sep="\t", index=False)

        paths=seq_info[["assembly"]]
        paths.to_csv(paths_out, sep="\t", index=False, header=False)

        get_comp(seqtk_exec, clu_out, fa_out, comp_out, clu_comp_out, seq_info, no_build_fasta)

        sys.stderr.write("Creating mash sketches...\n")
        if not no_build_fasta:
            std_result=run_mash_sketch(mash_exec, t, fa_out, msh_out, ss, m=None, in_type="fa")
            log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
        else:
            std_result=run_mash_sketch(mash_exec, t, paths_out, msh_out, ss, m=None, in_type="list")
            log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))

        sys.stderr.write("Calculating mash distances...\n")
        std_result=run_mash_dist(mash_exec, t, msh_out, msh_out, msh_dis_out)
        log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))

        add_clusters(clu_out, msh_dis_out, msh_dis_clu_out, ref=True, met=True)
        
        sys.stderr.write("Calculating thresholds...\n")
        get_thresholds(clu_out, msh_dis_clu_out, thr_prop_min, thr_abs_min, thr_prop_exp, clu_thr_out)

        if not no_build_index:
            idx_d="{}/ref_idx".format(d)
            idx_d_tmp="{}/ref_idx/tmp".format(d)

            if not os.path.isdir(idx_d):
                os.makedirs(idx_d)
    
            if not os.path.isdir(idx_d_tmp):
                os.makedirs(idx_d_tmp)
        
            sys.stderr.write("Indexing reference set...\n")
            themisto_cmd="{} --k 31 --n-threads {} --input-file {} --auto-colors --index-dir {} --temp-dir {}".format(themisto_build_exec, t, fa_out, idx_d, idx_d_tmp)
            std_result=subprocess.run(themisto_cmd, shell=True, check=True, capture_output=True, text=True)
            log.write("{}\n\n{}\n{}\n\n".format(std_result.args, std_result.stderr, std_result.stdout))
        
            sys.stderr.write("Reference set {} setup completed\n".format(d))

        log.close()

def get_thresholds(in_clu, in_dis, thr_prop_min, thr_abs_min, thr_prop_exp, out_file):
    
    clu=pd.read_csv(in_clu, sep="\t", dtype={'id': 'str'})
    clu=clu.groupby("cluster", as_index=False).agg(n=("id", "count"))

    dis=pd.read_csv(in_dis, sep="\t", dtype={'ref_id': 'str', 'met_id': 'str'})

    dis_diff=dis.query("ref_cluster != met_cluster")
    dis_diff=dis_diff.query("ref_id != met_id")
    dis_diff=dis_diff.groupby("ref_cluster", as_index=False).agg(
            dis_diff_med=("distance", "median"),
            dis_diff_min=("distance", "min"))
    dis_diff_med_med=np.median(dis_diff["dis_diff_med"])

    dis_same=dis.query("ref_cluster == met_cluster")
    dis_same=dis_same.query("ref_id != met_id")
    dis_same=dis_same.groupby("ref_cluster", as_index=False).agg(
            dis_same_med=("distance", "median"),
            dis_same_min=("distance", "min"),
            dis_same_max=("distance", "max"))
    dis_same_med_med=np.median(dis_same["dis_same_med"])
    
    dis=dis_same.merge(dis_diff, how="left", on="ref_cluster")

    dis["threshold"]=dis["dis_same_max"]*(1+thr_prop_exp)
    dis["threshold"]=np.where(dis["threshold"] > dis["dis_diff_min"], dis["dis_same_max"], dis["threshold"])
    dis=dis.rename(columns={"ref_cluster":"cluster"})
    dis=dis[["cluster", "threshold", "dis_same_max"]]
    
    t_min=np.median(dis_diff["dis_diff_med"])*thr_prop_min
    if thr_abs_min:
        t_min=thr_abs_min
    t_med=np.median(dis["threshold"])
    t_comp=np.where(t_med < t_min, t_min, t_med)

    t_summary=clu.merge(dis, how="left", on="cluster")
    t_summary["threshold"]=np.where(pd.isna(t_summary["threshold"]), 0, t_summary["threshold"])
    t_summary["threshold"]=np.where(t_summary["threshold"] < t_comp, t_comp, t_summary["threshold"])
    t_summary["dis_same_med_all"]=dis_same_med_med
    t_summary["dis_diff_med_all"]=dis_diff_med_med

    t_summary.to_csv(out_file, sep="\t", index=False)

def get_comp(seqtk_exec, in_clu, in_fa, out_file, out_file_summary, seq_info, no_build_fasta):
    
    clu=pd.read_csv(in_clu, sep="\t", dtype={'id': 'str'})
    

    if no_build_fasta:
        seq_count=len(seq_info)
        with open(out_file, "w") as outfile:
            setup_cmd="echo \"length\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts\""
            std_result=subprocess.run(setup_cmd, shell=True, stdout=outfile)
            for i in range(seq_count):
                assembly=seq_info["assembly"][i]
                seqtk_cmd="{} comp {} | cut -f2-13".format(seqtk_exec, assembly)
                std_result=subprocess.run(seqtk_cmd, shell=True, check=True, text=True, stdout=outfile)
    else:
        seqtk_cmd="{} comp {} >> {}".format(seqtk_exec, in_fa, out_file)
        std_result=subprocess.run(seqtk_cmd, shell=True, check=True, capture_output=True, text=True)

    lens=pd.read_csv(out_file, sep="\t", dtype={'id': 'str'})
    lens[["id"]] = seq_info[["id"]]
    lens.to_csv(out_file, sep='\t', index=False)
    lens=lens[["id", "length"]]

    clu_lens=clu.merge(lens, how="left", on="id")
    clu_lens=clu_lens.groupby("cluster", as_index=False).agg(
            n=("id", "count"),
            length_ave=("length", "mean"),
            length_min=("length", "min"),
            length_max=("length", "max"))

    clu_lens.to_csv(out_file_summary, sep="\t", index=False)

