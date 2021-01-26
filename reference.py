#!/usr/bin/env python

import os
import sys
import subprocess

import re
import pandas as pd
import numpy as np

from sketch import *

def setup_reference(mash_exec, themisto_build_exec, d, t, ss, thr_prop_min, thr_prop_exp, redo_thr):

    seq_info_f="{}/ref_info.tsv".format(d)
    seq_info=pd.read_csv(seq_info_f, sep="\t")
    seq_count=len(seq_info)
    
    fa_out="{}/ref.fasta".format(d)
    msh_out="{}/ref.msh".format(d)
    msh_dis_out="{}/ref_msh_dis.tsv".format(d)
    msh_dis_clu_out="{}/ref_msh_dis_clu.tsv".format(d)
    clu_thr_out="{}/ref_clu_thr.tsv".format(d)
    clu_out_m="{}/ref_clu.txt".format(d)
    clu_out="{}/ref_clu.tsv".format(d)

    if redo_thr == True:

        get_thresholds(clu_out, msh_dis_clu_out, thr_prop_min, thr_prop_exp, clu_thr_out)
    
    else:
        fo=open(fa_out, 'w')
        for i in range(seq_count):
            seq_id=seq_info["id"][i]
            cluster=seq_info["cluster"][i]
            assembly=seq_info["assembly"][i]

            fo.write(">{}\n".format(seq_id))
        
            f_tmp=open(assembly, 'r')
            for l in f_tmp:
                if not re.match(r'^>', l):
                    fo.write(l)
        fo.close()

        clu_m=seq_info[["cluster"]]
        clu_m.to_csv(clu_out_m, sep="\t", index=False, header=None)

        clu=seq_info[["id", "cluster"]]
        clu.to_csv(clu_out, sep="\t", index=False)

        run_mash_sketch(mash_exec, t, fa_out, msh_out, ss, m=None, in_type="fa")
        run_mash_dist(mash_exec, t, msh_out, msh_out, msh_dis_out)
        add_clusters(clu_out, msh_dis_out, msh_dis_clu_out, ref=True, met=True)
        get_thresholds(clu_out, msh_dis_clu_out, thr_prop_min, thr_prop_exp, clu_thr_out)
    
        idx_d="{}/ref_idx".format(d)
        idx_d_tmp="{}/ref_idx/tmp".format(d)

        if not os.path.isdir(idx_d):
            os.makedirs(idx_d)
    
        if not os.path.isdir(idx_d_tmp):
            os.makedirs(idx_d_tmp)
    
        themisto_cmd="{} --k 31 --n-threads {} --input-file {} --auto-colors --index-dir {} --temp-dir {}".format(themisto_build_exec, t, fa_out, idx_d, idx_d_tmp)
        subprocess.run(themisto_cmd, shell=True, check=True)

def get_thresholds(in_clu, in_dis, thr_prop_min, thr_prop_exp, out_file):
    
    clu=pd.read_csv(in_clu, sep="\t")
    clu=clu.groupby("cluster", as_index=False).agg(n=("id", "count"))

    dis=pd.read_csv(in_dis, sep="\t")

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
    t_med=np.median(dis["threshold"])
    t_comp=np.where(t_med < t_min, t_min, t_med)

    t_summary=clu.merge(dis, how="left", on="cluster")
    t_summary["threshold"]=np.where(pd.isna(t_summary["threshold"]), 0, t_summary["threshold"])
    t_summary["threshold"]=np.where(t_summary["threshold"] < t_comp, t_comp, t_summary["threshold"])
    t_summary["dis_same_med_all"]=dis_same_med_med
    t_summary["dis_diff_med_all"]=dis_diff_med_med

    t_summary.to_csv(out_file, sep="\t", index=False)

