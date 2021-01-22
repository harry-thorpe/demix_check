# demix_check

This pipeline checks the demixed binned reads from an mGEMS analysis. 

## Dependencies

* python3 with numpy and pandas
* R with tidyverse and cowplot (for plotting)
* build_index and pseudoalign
* mSWEEP
* mGEMS

## Running the pipeline

Clone the repository with:

```git clone https://github.com/harry-thorpe/demix_check.git```

The main script is:

```demix_check/demix_check.py```

and calling it with:

```python demix_check/demix_check.py --help```

prints the following help menu:

```
Pipeline for assessing the cluster assignments from mGEMS

Mode arguments:
  Arguments to select the run mode

  --mode_setup          Set up one or more references
  --mode_check          Check the results of an existing mGEMS analysis
  --mode_run            Run mGEMS and then check the results

Setup arguments:
  Arguments for setup mode

  --thr_prop_exp THR_PROP_EXP
                        proportion of maximum divergence within a cluster to expand the threshold by [default = 0.5]
  --thr_prop_min THR_PROP_MIN
                        proportion of median divergence between clusters to set minimum threshold to [default = 0.3]

Check arguments:
  Arguments for check mode

  --binned_reads_dir BINNED_READS_DIR
                        directory which contains the binned reads from an existing mGEMS analysis
  --msweep_abun MSWEEP_ABUN
                        file which contains the abundance estimates from an existing mSWEEP analysis

Run arguments:
  Arguments for run mode

  --r1 R1               r1 file
  --r2 R2               r2 file

General arguments:
  --out_dir OUT_DIR     output directory
  --ref REF             reference/s to use [either a string specifying the path to the reference directory or file containing paths to the reference directories]
  --min_abun MIN_ABUN   mSWEEP/mGEMS - only accept clusters with this abundance or greater [default = 0.01]
  --kmer_min_freq KMER_MIN_FREQ
                        mash - only use kmers with this frequency or greater [default = 3]
  --sketch_size SKETCH_SIZE
                        mash - sketch size to use [default = 10000]
  --plots               plot results (requires R and tidyverse ggplot cowplot) [default = off]
  --threads THREADS     number of threads to use [default = 1]
  -h, --help            show this help message and exit
```

## modes

There are 3 ways to run the pipeline.

## setup

The setup mode sets up the reference set/s for use with the pipeline. A reference set 
