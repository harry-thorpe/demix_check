# demix_check

This pipeline assess the demixed binned reads from an mGEMS analysis to help with interpreting the assigned clusters. This is an important step when running mGEMS on complex mixtures, when there is possible contamination, or when the species being analysed doesn't have a comprehensive reference set. In these cases the assigned clusters may be the closest available sequences from the reference set, but the reads may actually be from an unknown cluster which is not present in the reference set. To address this, demix_check calculates the genetic distances between reference isolates, and these are used to build distributions of within and between cluster distances. Genetic distances are then calculated between the demixed binned reads and the reference isolates. The query-ref distances are then compared to the ref-ref distances to determine whether the binned reads are from the assigned cluster or not. Mash is used to calculate distances as it has been shown to be very accurate and adds little computational burden to the whole pipeline.

## Dependencies

* python3 with numpy and pandas
* R with tidyverse and cowplot (for plotting)
* build_index and pseudoalign (from themisto)
* mSWEEP
* mGEMS
* mash

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

  --mode_setup          set up one or more reference set/s
  --mode_check          check the results of an existing mGEMS analysis
  --mode_run            run mGEMS and then check the results

Setup arguments:
  Arguments for setup mode
  
  --redo_thr            quickly recalculate thresholds only (based on --thr_prop_exp and/or --thr_prop_min/--thr_abs_min). The reference set/s must have been previously set up before running with this option [default = off]
  --thr_prop_exp THR_PROP_EXP
                        proportion of maximum divergence within a cluster to expand the threshold by [default = 0.2]
  --thr_prop_min THR_PROP_MIN
                        proportion of median divergence between clusters to set minimum threshold to [default = 0.2]
  --thr_abs_min  THR_ABS_MIN
                        absolute minimum threshold [default = not set] [overrides --thr_prop_min]

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
  --ref REF             reference set/s to use [either a string specifying the path to the reference directory or a file containing paths to the reference directories]
  --min_abun MIN_ABUN   mSWEEP/mGEMS - only accept clusters with this abundance or greater [default = 0.01]
  --kmer_min_freq KMER_MIN_FREQ
                        mash - only use kmers with this frequency or greater [default = 3]
  --sketch_size SKETCH_SIZE
                        mash - sketch size to use [default = 10000]
  --plots               plot results (requires R and tidyverse ggplot cowplot) [default = off]
  --threads THREADS     number of threads to use [default = 1]
  -h, --help            show this help message and exit
```

## Modes

There are 3 ways to run the pipeline - these are specified with the ```--mode_setup```, ```--mode_check```, ```--mode_run``` options.

### setup

The setup mode sets up the reference set/s for use with the pipeline. A reference set consists of a set of assemblies with labels, as would normally be used with mSWEEP/mGEMS. For the demix_check pipeline, each reference set is stored in a separate folder, and the name of the folder defines the name of the reference set. This folder must contain a file named ref_info.tsv, which is a tab separated file with information about the isolate names, clusters, and the location of the reference assemblies. For example, a reference set for Klebsiella pneumoniae may be stored as follows:

```
ref_dir/Kpne
ref_dir/Kpne/ref_info.tsv
```

where the contents of ```ref_dir/Kpne/ref_info.tsv``` has the header names [id, cluster, assembly]:

```
id  cluster assembly
isolate_1 SC1 path/to/isolate_1.fasta
isolate_2 SC1 path/to/isolate_2.fasta
isolate_3 SC2 path/to/isolate_3.fasta
```

The ```ref_dir/Kpne``` folder can then be given as input to demix_check in setup mode.

Example command:

```python demix_check/demix_check.py --mode_setup --ref ref_dir/Kpne```

This will do the following:

* Create a multifasta of the reference genomes
* Index this multifasta for use with themisto
* Generate mash sketches from the multifasta
* Calculate pairwise mash distances between the reference genomes
* Calculate appropriate thresholds for each cluster within the reference set
* Plot the within and between cluster distances along with the thresholds

This reference set is then ready for use. If the thresholds aren't appropriate, they can be adjusted with the ```--thr_prop_exp``` and ```--thr_prop_min``` options, and specifying ```--redo_thr``` will quickly recalculate the thresholds without doing the rest of the reference set up (creating and indexing the fasta file etc).

### check

The check mode takes an existing mGEMS analysis, along with a reference set which has been set up with ```setup``` (this must be the same reference set as was used for the mGEMS analysis), and runs the checking part of the demix_check pipeline. This mode is useful if mGEMS has already been run. The folder containing the binned reads from mGEMS and the file containing the mSWEEP abundance estimates must both be specified as inputs, along with the reference set.

Example command:

```python demix_check/demix_check.py --mode_check --binned_reads_dir mGEMS_analysis/binned_reads --msweep_abun mGEMS_analysis/msweep_abundances.txt --out_dir output_dir --ref ref_dir/Kpne```

This will do the following:

* Generate a mash sketch from each set of binned reads
* Calculate mash distances between the binned reads and the reference genomes
* Compare the query-ref distances to the ref-ref distances and thresholds (derived from the ref-ref comparisons) to calculate a score for each cluster
* plot the ref-ref distances against the query-ref distances

### run

The run mode takes a set of mixed reads along with a reference set (which has been set up with ```setup```), and runs the whole mSWEEP/mGEMS pipeline before then checking the results as described above. An output folder is specified, and then the results are placed in a folder which is specific to the reference set, so if the ```Kpne``` reference set is used, the results are placed in ```out_dir/Kpne```.

Example command:

```python demix_check/demix_check.py --mode_run --r1 reads_1.fastq.gz --r2 reads_2.fastq.gz --out_dir out_dir --ref ref_dir/Kpne```

This will do the following:

* Pseudoalign the reads to the reference set with themisto
* Run the mSWEEP abundance estimation
* Run mGEMS to bin the reads
* Run the checking part of demix_check

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
```

To tell the pipeline to run in a hierachical manner, a file containing the paths to the reference folders must be given to ```--ref```. This file must be tab separated (no header), with each column corresponding to a different level:

```
ref_file.tsv:
ref_dir/species_ref ref_dir/Ecoli
ref_dir/species_ref ref_dir/Kpne
```

For this to work, each reference set in the second level must be present as a cluster label in the first level e.g. ```Ecoli``` must be a cluster label in ```species_ref```.

Example command:

```python demix_check/demix_check.py --mode_run --r1 reads_1.fastq.gz --r2 reads_2.fastq.gz --out_dir out_dir --ref ref_file.tsv```

The pipeline will start by running the whole demix_check pipeline on each unique reference set in the first level, using the input fastq files as the mixed reads. It will then proceed to the second level, and for each reference set in that level, it will look for the binned reads for that cluster from the first level e.g. for ```Ecoli```, it will look for ```Ecoli_1.fastq.gz``` and ```Ecoli_2.fastq.gz``` in the mGEMS output from the first level. This continues until all levels have been completed.

## Output and interpretation

A score is calculated (1 - highest confidence, 4 - lowest confidence) for each cluster based on comparisons of the genetic distances between demixed binned reads (query) and known reference sequences (ref):

1. There is overlap between the query-ref distances and the ref-ref distances from the same cluster.
2. The query-ref distances are greater than the ref-ref distances from the same cluster, but are within the threshold.
3. The query-ref distances are greater than the threshold, but are closer to the median within-cluster distance than the median between-cluster distance.
4. The query-ref distances are greater than the threshold, and are closer to the median between-cluster distance than the median within-cluster distance.

The ref-ref distances are also plotted against the query-ref distances, along with the thresholds.

For each checking analysis the scores are saved in ```clu_score.tsv```, and the plots in ```sample_plot.pdf```.

If the pipeline has been run in a hierarchical manner, there are two additional output files, ```clu_out_summary.tsv``` and ```summary_plot.pdf```, which show the results for the first two levels in the hierarchy.
