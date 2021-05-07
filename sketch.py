#!/usr/bin/env python

import os
import sys
import subprocess
import gzip

import re
import pandas as pd
import numpy as np

def run_mash_sketch(mash_exec, t, in_file, out_file, ss, m, in_type):
    std_result=""

    out_file_p=re.sub(r'.msh$', '', out_file)
    seq_id=os.path.basename(out_file_p)

    if in_type == "fa":
        mash_cmd="{} sketch -p {} -s {} -i -o {} {}".format(mash_exec, t, ss, out_file_p, in_file)
        std_result=subprocess.run(mash_cmd, shell=True, check=True, capture_output=True, text=True)
    elif in_type == "fq":
        mash_cmd="{} sketch -p {} -s {} -r -m {} -I {} -C - -o {} {}".format(mash_exec, t, ss, m, seq_id, out_file_p, in_file)
        std_result=subprocess.run(mash_cmd, shell=True, check=True, capture_output=True, text=True)
    
    return(std_result)

def run_mash_dist(mash_exec, t, ref_file, query_file, out_file):
    std_result=""

    out_file_tmp="{}.tmp.gz".format(out_file)

    mash_cmd="{} dist -p {} {} {} | gzip > {}".format(mash_exec, t, ref_file, query_file, out_file_tmp)
    std_result=subprocess.run(mash_cmd, shell=True, check=True, capture_output=True, text=True)

    data=pd.read_csv(out_file_tmp, sep="\t", names=["ref_id", "met_id", "distance", "p", "hashes"])
    data[["hashes", "ss"]]=data.hashes.str.split('/', expand=True)
    data=data[["ref_id", "met_id", "distance", "hashes", "ss", "p"]]
    data.to_csv(out_file, sep="\t", index=False)

    os.remove(out_file_tmp)
    
    return(std_result)

def run_mash_screen(mash_exec, t, ref_file, query_file, out_file, met_id):
    std_result=""

    out_file_tmp="{}.tmp.gz".format(out_file)

    mash_cmd="{} screen -p {} {} {} | gzip > {}".format(mash_exec, t, ref_file, query_file, out_file_tmp)
    std_result=subprocess.run(mash_cmd, shell=True, check=True, capture_output=True, text=True)

    data=pd.read_csv(out_file_tmp, sep="\t", names=["distance", "hashes", "coverage", "p", "ref_id", "TMP"])
    data[["hashes", "ss"]]=data.hashes.str.split('/', expand=True)
    data["distance"]=1-data["distance"]
    data["met_id"]=met_id
    data=data[["ref_id", "met_id", "distance", "hashes", "ss", "p", "coverage"]]
    data.to_csv(out_file, sep="\t", index=False)

    os.remove(out_file_tmp)
    
    return(std_result)

def add_clusters(in_clu, in_dis, out_file, ref=False, met=False):
    clu=pd.read_csv(in_clu, sep="\t")
    
    data=pd.read_csv(in_dis, sep="\t")
    if ref == True:
        data=data.merge(clu, how="left", left_on="ref_id", right_on="id")
        data=data.rename(columns={"cluster":"ref_cluster"})
        data=data.drop("id", axis=1)
    if met == True:
        data=data.merge(clu, how="left", left_on="met_id", right_on="id")
        data=data.rename(columns={"cluster":"met_cluster"})
        data=data.drop("id", axis=1)
        data["category"]=np.where(data["ref_cluster"] == data["met_cluster"], 'same cluster', 'different cluster')
    data.to_csv(out_file, sep="\t", index=False)

