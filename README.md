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

The setup mode sets up the reference set/s for use with the pipeline. A reference set consists of a set of assemblies with labels, as would normally be used with mSWEEP/mGEMS. For the demix_check pipeline, each reference set is stored in a separate folder, and the name of the folder defines the name of the reference set. This folder must contain a file named ref_info.tsv, which is a tab separated file with information about the isolate names, clusters, and the location of the reference assemblies. For example, a reference for Klebsiella pneumoniae may be stored as follows:

```
Kpne
Kpne/ref_info.tsv
```

where the contents of ```Kpne/ref_info.tsv``` have the header names [id, cluster, assembly]:
```
id  cluster assembly
isolate_1 SC1 path/to/isolate_1.fasta
isolate_2 SC1 path/to/isolate_2.fasta
isolate_3 SC2 path/to/isolate_3.fasta
```

The ```Kpne``` folder can then be given as input to demix_check in setup mode, which will do the following:

* Create a multifasta of the reference genomes
* Index this multifasta for use with themisto
* Generate mash sketches from the multifasta
* Calculate pairwise mash distances between the reference genomes
* Calculate appropriate thresholds for each cluster within the reference set
* Plot the within and between cluster distances along with the thresholds

This reference set is then ready for use.

## check

The check mode takes an existing mGEMS analysis, along with a reference set which has been set up with ```setup``` (this must be the same reference set as was used for the mGEMS analysis), and runs the checking part of the demix_check pipeline. This mode is useful if mGEMS has already been run. The folder containing the binned reads from mGEMS and the file containing the mSWEEP abundance estimations must both be specified as inputs, along with the reference set. The pipeline will then:

* Generate a mash sketch from each set of binned reads
* Calculate mash distances between the binned reads and the reference genomes
* Compare the query-ref mash distances to the thresholds (derived from the ref-ref comparisons) to decide whether the assignments are likely to be correct or not
* plot the ref-ref distances against the query-ref distances

## run

The run mode takes a set of mixed reads along with a reference set (which has been set up with ```setup```), and runs the mGEMS pipeline before then checking the results as described above. An output folder is specified, and then the results are placed in a folder which is specific to the reference set, so if the ```Kpne``` reference set is used, the results are placed in ```out_dir/Kpne```.

This mode can also be used to run mGEMS in a hierachical manner with several reference sets. This may be desirable, for example if we suspect that the mixed reads may contain multiple sequence clusters from several different species, we may first want to demix the reads into species bins, and then into sequence cluster bins from each species bin. To do this, we need one reference set which has clusters labelled by species, and then a separate reference set for each species that we want to separate into sequence clusters. These reference sets can be composed of the same isolates, but with different labels, for example:

```
ref_dir/species_ref
ref_dir/Ecoli
ref_dir/Kpne

ref_dir/species_ref/ref_info.tsv:
id  cluster assembly
isolate_1 Ecoli path/to/isolate_1.fasta
isolate_2 Ecoli path/to/isolate_2.fasta
isolate_3 Kpne path/to/isolate_3.fasta
isolate_4 Kpne path/to/isolate_4.fasta
isolate_5 Koxy path/to/isolate_5.fasta
isolate_6 Kgri path/to/isolate_6.fasta
isolate_7 Kqps path/to/isolate_7.fasta

ref_dir/Ecoli/ref_info.tsv:
id  cluster assembly
isolate_1 SC1 path/to/isolate_1.fasta
isolate_2 SC2 path/to/isolate_2.fasta

ref_dir/Kpne/ref_info.tsv:
id  cluster assembly
isolate_3 SC1 path/to/isolate_3.fasta
isolate_4 SC2 path/to/isolate_4.fasta

