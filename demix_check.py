#!/usr/bin/env python

import os
import io
import sys
import subprocess
import argparse
import shutil

import re
import pandas as pd
import numpy as np

# import other python functions
from get_assignments import *
from check_assignments import *
from sketch import *
from reference import *

# dependency commands
themisto_build_exec='themisto build'
themisto_align_exec='themisto pseudoalign'
mSWEEP_exec='mSWEEP'
mGEMS_exec='mGEMS'
mash_exec='mash'
seqtk_exec='seqtk'

dir_path=os.path.dirname(__file__)

# check dependencies
def cmd_exists(cmd):
    return shutil.which(cmd) is not None

if not cmd_exists('themisto'):
    sys.stderr.write("WARNING: can't find command: {}\n".format('themisto'))
if not cmd_exists(mSWEEP_exec):
    sys.stderr.write("WARNING: can't find command: {}\n".format(mSWEEP_exec))
if not cmd_exists(mGEMS_exec):
    sys.stderr.write("WARNING: can't find command: {}\n".format(mGEMS_exec))
if not cmd_exists(mash_exec):
    sys.stderr.write("WARNING: can't find command: {}\n".format(mash_exec))
if not cmd_exists(seqtk_exec):
    sys.stderr.write("WARNING: can't find command: {}\n".format(seqtk_exec))

# set up args ##########################
parser=argparse.ArgumentParser(description="Pipeline for assessing the cluster assignments from mGEMS", add_help=False)

# args to specify run mode
parser._optionals.title="Mode arguments"
parser._optionals.description="Arguments to select the run mode"
mode_group=parser.add_mutually_exclusive_group()
mode_group.add_argument('--mode_setup', action="store_true", default=False, help='set up one or more reference set/s')
mode_group.add_argument('--mode_check', action="store_true", default=False, help='check the results of an existing mGEMS analysis')
mode_group.add_argument('--mode_run', action="store_true", default=False, help='run mGEMS and then check the results')

# args specific to setup mode
setup_group=parser.add_argument_group('Setup arguments', 'Arguments for setup mode')
setup_group.add_argument('--redo_thr', action="store_true", default=False, help='quickly recalculate thresholds only (based on --thr_prop_exp and/or --thr_prop_min/--thr_abs_min). The reference set must have been previously set up before running with this option [default = off]')
setup_group.add_argument('--thr_prop_exp', type=float, default=0.2, help='proportion of maximum divergence within a cluster to expand the threshold by [default = %(default)s]')
setup_group.add_argument('--thr_prop_min', type=float, default=0.2, help='proportion of median divergence between clusters to set minimum threshold to [default = %(default)s]')
setup_group.add_argument('--thr_abs_min', type=float, help='absolute minimum threshold [default = not set] [overrides --thr_prop_min]')

# args specific to check mode
check_group=parser.add_argument_group('Check arguments', 'Arguments for check mode')
check_group.add_argument('--binned_reads_dir', type=str, required='--mode_check' in sys.argv, help='directory which contains the binned reads from an existing mGEMS analysis')
check_group.add_argument('--msweep_abun', type=str, required='--mode_check' in sys.argv, help='file which contains the abundance estimates from an existing mSWEEP analysis')

# args specific to run mode
run_group=parser.add_argument_group('Run arguments', 'Arguments for run mode')
run_group.add_argument('--r1', type=str, required='--mode_run' in sys.argv, help='r1 file')
run_group.add_argument('--r2', type=str, required='--mode_run' in sys.argv, help='r2 file')

# general args
general_group=parser.add_argument_group('General arguments')
general_group.add_argument('--out_dir', type=str, required='--mode_check' in sys.argv or '--mode_run' in sys.argv, help='output directory')
general_group.add_argument('--ref', type=str, required='--mode_setup' in sys.argv or '--mode_check' in sys.argv or '--mode_run' in sys.argv, help='reference set/s to use [either a string specifying the path to the reference directory or a file containing paths to the reference directories]')
general_group.add_argument('--min_abun', type=float, default=0.01, help='mSWEEP/mGEMS - only accept clusters with this abundance or greater [default = %(default)s]')
general_group.add_argument('--sketch_size', type=int, default=10000, help='mash - sketch size to use [default = %(default)d]')
general_group.add_argument('--plots', action="store_true", default=False, help='plot results (requires R and the tidyverse and cowplot packages) [default = off]')
general_group.add_argument('--threads', type=int, default=1, help='number of threads to use [default = %(default)d]')
general_group.add_argument('--keep', action="store_true", default=False, help='keep large intermediate files [default = off]')
general_group.add_argument('-h', '--help', action='help', help="show this help message and exit")

# set up arg vars
args=parser.parse_args()

mode_setup=args.mode_setup
mode_check=args.mode_check
mode_run=args.mode_run

redo_thr=args.redo_thr
thr_prop_exp=args.thr_prop_exp
thr_prop_min=args.thr_prop_min
thr_abs_min=args.thr_abs_min

binned_reads_d=args.binned_reads_dir
msweep_abun=args.msweep_abun

r1=args.r1
r2=args.r2

out_d=args.out_dir
ref_in=args.ref
min_abun=args.min_abun
ss=args.sketch_size
t=args.threads
plots=args.plots
keep=args.keep

# trim trailing slashes
if binned_reads_d:
    binned_reads_d=binned_reads_d.rstrip('/')
if ref_in:
    ref_in=ref_in.rstrip('/')
if out_d:
    out_d=out_d.rstrip('/')
########################################


# collect reference information ########
if ref_in:
    sys.stderr.write("Finding reference set/s...\n")
    if os.path.isfile(ref_in):
        ref_str=pd.read_csv(ref_in, sep="\t", header=None)
    elif os.path.isdir(ref_in):
        ref_str=pd.read_csv(io.StringIO(ref_in), sep="\t", header=None)
    else:
        sys.stderr.write("ERROR: can't find reference {}\n".format(ref_in))
        sys.exit(1)

    ref_str_c=ref_str.shape[1]
    ref_str_r=ref_str.shape[0]

    ref_str_long=ref_str.melt(value_name="value")
    ref_ds=list(set(ref_str_long.value.values))
    ref_ds_count=len(ref_ds)
    for tmp_ref_d in ref_ds:
        if not os.path.isdir(tmp_ref_d):
            sys.stderr.write("ERROR: can't find reference {}\n".format(tmp_ref_d))
            sys.exit(1)

    ref_str_d={}
    for c in range(1, ref_str_c):
        for r in range(0, ref_str_r):
            ref_str_d[ref_str[c][r]]=ref_str[(c-1)][r]

    sys.stderr.write("Found {} reference set/s:\n{}\n".format(ref_ds_count, ref_ds))
########################################


# run setup mode #######################
if mode_setup:
    
    sys.stderr.write("Running in setup mode...\n")
    for ref_d in ref_ds:
        if os.path.isdir(ref_d) and os.path.isfile("{}/ref_info.tsv".format(ref_d)):
            
            setup_reference(mash_exec, themisto_build_exec, seqtk_exec, ref_d, t, ss, thr_prop_min, thr_abs_min, thr_prop_exp, redo_thr)
            
            if plots:
                sys.stderr.write("Plotting output...\n")
                plot_cmd="Rscript {}/plot_reference.R {}".format(dir_path, ref_d)
                subprocess.run(plot_cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    sys.stderr.write("Setup mode completed\n")
########################################


# run check mode #######################
if mode_check:
    
    sys.stderr.write("Running in check mode...\n")
    ref_d=ref_ds[0]
    check_mGEMS(mash_exec, seqtk_exec, t, ss, min_abun, ref_d, out_d, binned_reads_d, msweep_abun)

    if plots:
        sys.stderr.write("Plotting output...\n")
        plot_cmd="Rscript {}/plot_sample_single.R {} {}".format(dir_path, out_d, ref_d)
        subprocess.run(plot_cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    sys.stderr.write("Check mode completed\n")
########################################


# run run mode #########################
if mode_run:
    
    sys.stderr.write("Running in run mode...\n")
    if not os.path.isdir(out_d):
        os.makedirs(out_d)

    out_f="{}/clu_out_summary.tsv".format(out_d)
    out_data=pd.DataFrame()

    # loop through reference directories
    for c in range(ref_str_c):
        l=ref_str.iloc[:, c].values
        l=list(set(l))

        for ref_d in l:
            ref=os.path.basename(ref_d)

            # use either starting reads or binned reads from mGEMS
            if c == 0:
                rr1=r1
                rr2=r2
                ref_p="input_reads"

                if not os.path.isfile(rr1):
                    sys.stderr.write("ERROR: can't find read file {}\n".format(rr1))
                    sys.exit(1)
                if not os.path.isfile(rr2):
                    sys.stderr.write("ERROR: can't find read file {}\n".format(rr2))
                    sys.exit(1)
            else:
                ref_d_p=ref_str_d[ref_d]
                ref_p=os.path.basename(ref_d_p)

                rr1="{}/{}/binned_reads/{}_1.fastq.gz".format(out_d, ref_p, ref)
                rr2="{}/{}/binned_reads/{}_2.fastq.gz".format(out_d, ref_p, ref)

            if os.path.isdir(ref_d) and os.path.isfile(rr1) and os.path.isfile(rr2):
                
                out_dr="{}/{}".format(out_d, ref)
                binned_reads_d="{}/binned_reads".format(out_dr)
                msweep_abun="{}/msweep_abundances.txt".format(out_dr)
                
                if not os.path.isdir(out_dr):
                    os.makedirs(out_dr)

                ref_clu="{}/ref_clu.tsv".format(ref_d)
                ref_msh="{}/ref.msh".format(ref_d)
                rr12="{} {}".format(rr1, rr2)
                msh_scr_dis_out="{}/msh_scr_dis.tsv.gz".format(out_dr)
                msh_scr_dis_clu_out="{}/msh_scr_dis_clu.tsv.gz".format(out_dr)

                # run mash screen
                sys.stderr.write("Running mash screen...\n")
                std_result=run_mash_screen(mash_exec, t, ref_msh, rr12, msh_scr_dis_out, ref_p)
                add_clusters(ref_clu, msh_scr_dis_out, msh_scr_dis_clu_out, ref=True, met=False)
                
                if plots:
                    sys.stderr.write("Plotting output...\n")
                    plot_cmd="Rscript {}/plot_mash_screen.R {} {}".format(dir_path, out_dr, ref_d)
                    subprocess.run(plot_cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

                # run the mSWEEP/mGEMS pipeline
                run_mGEMS(themisto_align_exec, mSWEEP_exec, mGEMS_exec, t, min_abun, rr1, rr2, ref_d, out_dr, binned_reads_d, msweep_abun, keep)

                # check the mGEMS bins
                check_mGEMS(mash_exec, seqtk_exec, t, ss, min_abun, ref_d, out_dr, binned_reads_d, msweep_abun)
                
                if plots:
                    sys.stderr.write("Plotting output...\n")
                    plot_cmd="Rscript {}/plot_sample_single.R {} {}".format(dir_path, out_dr, ref_d)
                    subprocess.run(plot_cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                
                out_tmp_f="{}/clu_score.tsv".format(out_dr)
                if os.path.isfile(out_tmp_f):
                    out_data_tmp=pd.read_csv(out_tmp_f, sep="\t", dtype={'cluster': 'str'})
                    out_data_tmp["ref"]=ref
                    out_data_tmp["level"]=c
                    out_data_tmp["out_d"]=out_dr
                    out_data_tmp["ref_d"]=ref_d

                    out_data=pd.concat([out_data, out_data_tmp], sort=False, ignore_index=True)

    out_data.to_csv(out_f, sep="\t", index=False)
    
    if plots and ref_str_c > 1:
        sys.stderr.write("Plotting output...\n")
        plot_cmd="Rscript {}/plot_sample_multi.R {}".format(dir_path, out_d)
        subprocess.run(plot_cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    sys.stderr.write("Run mode completed\n")
########################################

