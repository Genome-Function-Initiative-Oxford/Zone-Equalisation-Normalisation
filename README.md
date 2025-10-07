# Zone Equalisation Normalisation

ZEN-norm is a Python package for normalising bigWigs of genomic signal, such as ATAC-seq, ChIP-seq and TT-seq by Zone Equilisation Normalisation. Other features include: creating plots to compare normalisation method performance and reversing prior bigWig normalisation.

<p><img src="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/blob/assets/Images/ZEN_Overview_Figure.png" width="100%"></p>

---

<details open="open">
  <summary><b>Contents</b></summary>
  <ol>
    <li><a href="#installation">Installation</a></li>
    <li><a href="#zen_norm">ZEN bigWig Normalisation</a></li>
    <li><a href="#norm_compare">Plots for Evaluating Normalisation Method Performance</a></li>
    <li><a href="#reverse_norm">Reversing Prior bigWig Normalisation</a></li>
    <li><a href="#visualising_bigwigs">Visualising bigWigs</a></li>
  </ol>
</details>

<br>

---

<a id="installation"></a>
## Installation
ZEN-norm is designed to run on Python 3.10 and above and is installable from either PyPI or Conda.

<a id=""></a>
<details open="open">
  <summary><b>PyPI Installation</b></summary>
  To install the ZEN-norm package from <a href="https://test.pypi.org/project/ZEN-norm-test/">PyPI</a>, run the command below:

```
python -m pip install  --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple ZEN-norm-test
```
</details>

<a id=""></a>
<details open="open">
  <summary><b>Conda Installation</b></summary>
  To install the ZEN-norm package from <a href="">Conda</a>, active the conda environment you'd like to install the package into (<code>conda activate ...</code>) and run the command below:

```
conda install zen-norm
```

  Alternatively, if there are issues installing ZEN-norm, a conda environment with the required packages can be created using the <code>zen_environment.yml</code> file.

  ```
  conda env create --name zen_env --file=environment/zen_environment.yml
  conda activate zen_env
  ```
</details>

<br>

---

<a id="example_data"></a>
## Example Data
For the following tutorial, we demonstrate the features of ZEN-norm using the following publically avaliable data:

<a id=""></a>
<details open="open">
  <summary><b>Mouse Embryonic Stem Cell ATAC-seq</b></summary>
  This dataset contains ATAC-seq of E14 mouse embryonic stem cells (mESCs) treated with leukemia inhibitory factor (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3399494">LIF</a>) and retinoic acid (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3399495">RA</a>) (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120376">GSE120376</a>). Two bigWigs for these mapped to the mm10 genome can be downloaded from NCBI GEO using the following commands:

```
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3399495&format=file&file=GSM3399495%5FE14%5FATAC%5FRA%2Ebw" -O E14_ATAC_RA.bw
```
```
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3399494&format=file&file=GSM3399494%5FE14%5FATAC%5FLIF%2Ebw" -O E14_ATAC_LIF.bw
```

</details>

<br>

---

<a id="zen_norm"></a>
## ZEN bigWig Normalisation
The module `ZoneNorm` is for normalising genomic coverage with ZEN. Steps include: BAM to bigWig mapping, creating smoothed signals, distribution fitting, signal zone prediction and creating normalised bigWigs.

<br>

<a id=""></a>
<details open="open">
  <summary><b>Specifying Input BAMs or bigWigs</b></summary>
  ZEN-norm is designed to run on genomic signals from either BAM or bigWig files. Input files must be specified as either a list of file paths or a directory. If your input files are bigWigs that have been pre-normalised, then it is advisable to first remap them without normalisation, or to use ZEN-norm's <a href="#reverse_norm">reverse bigWig normalisation</a> feature.

  <br>

  <details>
    <summary><b><sub>Click for BAM Example</sub></b></summary>
    Say you have a folder containing two BAMs <code>E14_ATAC_RA.bam</code> and <code>E14_ATAC_LIF.bam</code> in a directory <code>mESC/mm10/BAMs</code> (see <a href="#example_data">Example Data</a>), you can specify the input as either:

```python
# Specify specific BAMs to process
bam_paths = ["mESC/mm10/BAMs/E14_ATAC_RA.bam", "mESC/mm10/BAMs/E14_ATAC_LIF.bam"]
# Or set directory containing BAMs
bam_directory = "mESC/mm10/BAMs"
```
    
  </details>

  <details>
    <summary><b><sub>Click for bigWig Example</sub></b></summary>
    Say you have a folder containing two bigWigs <code>E14_ATAC_RA.bw</code> and <code>E14_ATAC_LIF.bw</code> in a directory <code>mESC/mm10/bigWigs</code> (see <a href="#example_data">Example Data</a>), you can specify the input as either:
  
```python
# Specify specific bigWigs to process
bigwig_paths = ["mESC/mm10/bigWigs/E14_ATAC_RA.bw", "mESC/mm10/bigWigs/E14_ATAC_LIF.bw"]
# Or set directory containing bigWigs
bigwig_directory = "mESC/mm10/bigWigs"
```

  </details>

  <br>
</details>

<a id=""></a>
<details open="open">
  <summary><b>Initialisation of a ZoneNorm Object</b></summary>
  After deciding the input files, a <code>ZoneNorm</code> object must be created to specify the parameters for normalisation. For example: 
  
```python
from ZEN_norm.zone_norm import ZoneNorm

# Create ZoneNorm object
znorm = ZoneNorm(analysis_name = "mESC_Analysis",
                 bigwig_paths = ["mESC/mm10/bigWigs/E14_ATAC_RA.bw", "mESC/mm10/bigWigs/E14_ATAC_LIF.bw"],
                 n_cores = 4,
                 norm_method = "ZEN")
```

  <details>
    <summary open="open"><b><sub>Key Parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>analysis_name</code> | The name of the output folder to save results to. |
| <code>bam_paths</code> | List of paths of BAM files of interest. This takes priority over bam_directory, bigwig_paths and bigwig_directory. |
| <code>bam_directory</code> | Directory containing the BAM files of interest. This takes priority over bigwig_paths and bigwig_directory. |
| <code>bigwig_directory</code> | List of paths of bigWig and/or wig files of interest. This takes priority over bigwig_directory. |
| <code>bam_directory</code> | Directory containing the bigWig and/or wig files of interest. This is required if neither bam_paths, bam_directory or bigwig_paths are set. |
| <code>n_cores</code> | The number of cores / CPUs to use if using multiprocessing. |
| <code>norm_method</code> | Name of the normalisation method to apply. Options include: <code>ZEN</code> |

  </details>

  <details>
    <summary open="open"><b><sub>Additional Parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>chromosomes</code> | List of chromosomes to run analysis on. By default, all autosomal and sex chromosomes will be used. |
| <code>chrom_sizes_file</code> | Path to file of tab separated chromosomes and sizes for creating bigBed or converting wig to bigWig. |
| <code>interleave_sizes</code> | If <code>True</code> and using multiple cores, then process larger chromosomes alongside smaller ones to reduce memory usage. Otherwise process chromosomes from largest to smallest. |
| <code>sample_names</code> | Optionally set as a list of custom names corresponding to each file. e.g. 'cntrl_s1.bw' and 'flt3_inhibitor_s1.bw' could be set as ["Control Sample", "Treated Sample"]. This will be converted to a dictionary mapping original file names to the provided custom names. e.g. accessing sample_names would return {"cntrl_s1": "Control Sample", "flt3_inhibitor_s1": "Treated Sample"}. |
| <code>blacklist</code> | File path to blacklist file with chromosome coordinates to exclude. |
| <code>genome_size</code> | Required if using "RPGC" BAM normalisation. |
| <code>extend_reads</code> | Whether to enable read extension when creating bigWigs from BAMs with deepTools. |
| <code>filter_strand</code> | Specifies whether to separated reads by strand when creating bigWigs from BAMs with deepTools. Setting as True will create a forward and reverse strand bigWig per BAM (e.g. for nascent transcription). Leaving as False will disable strand filtering. Setting this as "forward" or "reverse" will extract the specified strand per BAM. |
| <code>exclude_zero</code> | Whether to ignore zeros when calculating statistics within signal zones. |
| <code>zone_remove_percent</code> | Set as a value greater than zero to ignore an upper and lower percentile when calculating statistics of signal in zones. e.g. zone_remove_percent = 100 for a signal of 1,000,000 and norm_method = "ZEN" allows the 100th upper and lower percentile to be removed when calculating mean and standard deviation. This prevents extreme values, such as those caused by technical artifacts in ATAC/ChIP-seq, from biasing the normalisation. |
| <code>norm_stats_type</code> | Type of statistics to use for ZEN normalisation. By default, this is "signal_padded_sample_zones", which will calculate averages and deviations from signal within sample-specific padded zones. Other options include: "signal_unpadded_sample_zones", which instead uses unpadded zones, "signal_padded_merged_zones" or "signal_unpadded_merged_zones", which use zone coordinates merged across all samples, or "signal", which considers all signal for statistics calculations. |
| <code>norm_power</code> | If norm_method is set as "Power", set this as a number to raise signal to the power of. |
| <code>deletion_size</code> | Threshold number of consecutive zeros in the signal for a region to be considered a potential deletion. |
| <code>downsample_size</code> | Number of positions of signal to select when fitting distributions. |
| <code>downsample_seed</code> | Integer seed for random sampling of signal when fitting distributions. |
| <code>kernel</code> | Smoothing kernel applied to predict signal scores. Custom kernel can be given as an array. |
| <code>test_distributions</code> | List of distributions to test for signal zone prediction. Supported distributions include: "norm", "laplace", "logistic", "cauchy", "gumbel_l", "gumbel_r". |
| <code>log_transform</code> | Whether to apply logarithm before transforming the smoothed signal during distribution testing. |
| <code>zone_distribution</code> | The distribution to fit for signal zone prediction. Can be set as either a distribution name, or as "infer" to automatically set the best tested distribution. |
| <code>zone_param_type</code> | The distribution parameter type to fit for signal zone distribution. Can be set as either "mean_fit", "median_fit", "scipy_fit" or "infer" to automatically set the best tested parameter type. |
| <code>zone_probability</code> | Probabilty threshold (between 0 to 1) from which the signal cut off is dervied from the distribution fitted for signal zone prediction. |
| <code>bin_size</code> | Number of base pairs to use in a full bin, e.g. 1,000. |
| <code>extend_n_bins</code> | An integer that alongside bin size, determines the amount of padding to add to predicted zones. E.g. if bin_size = 1,000 and extend_n_bins = 2, an unpadded zone coordinate of [5,678, 6,789] is first rounded to [5,000, 7,000], then extended either side by (2 * 1,000) as [3,000, 9,000]. |
| <code>merge_depth</code> | Minimum distance between two zones for them to be combined into one. |
| <code>min_region_bps</code> | Minimum size of a region initially predicted to have signal, i.e. the number of base pairs that must consecutively exceed the zone threshold. The zone is likely to be larger than this if rounding region coordinates to the nearest bins and adding one or more bins worth of padding. |
| <code>quality_filter</code> | Set a True to filter out regions of signal within zones with insufficient signal. |
| <code>min_different_bps</code> | If using the quality filter, set the minimum number of different base pairs that need to be found consecutively for a signal to be considered good quality. |
| <code>verbose</code> | Set as an integer greater than 0 to display progress messages. |

  </details>
  
  <br>
  
</details>

<a id=""></a>
<details open="open">
  <summary><b>Creating Convolved Signal</b></summary>
  After initialising the <code>ZoneNorm</code> object, the next step is to run <code>readBigWigSignals</code>. This method reads signal from bigWigs, calculates statistics about the signal and creates a convolved version of the signal to use for distribution fitting. 
  
```python
# Extract signal from each bigWig and created smoothed version
znorm.readBigWigSignals()
```

  <details>
    <summary><b><sub>Optional Parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>chromosomes</code> | List of chromosomes to read signal from per bigWig. |
| <code>replace_existing</code> | Set as <code>True</code> to overwrite previously created files. |

  </details>

  <br>

</details>

<a id=""></a>
<details open="open">
  <summary><b>Distribution Fitting</b></summary>
  After creating a smoothed version of the signal, distributions are fitted to it via <code>testDistributions</code>. This method further transforms the smoothed signal, and tests how well different distributions fit this transformed signal. The aim of this is for the selected distribution to be used after for signal zone prediction.
  
```python
# Fit distributions for signal zone prediction
znorm.testDistributions()
```

  <details>
    <summary><b><sub>Optional Parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>chromosomes</code> | List of chromosomes to read signal from per bigWig. |
| <code>replace_existing</code> | Set as <code>True</code> to overwrite previously created files. |

  </details>
  
  <br>

</details>

<a id=""></a>
<details open="open">
  <summary><b>Signal Zone Prediction</b></summary>
  Signal zones are predicted using one of the distributions fitted in the previous step.
  
```python
# Use distribution to predict signal zone coordinates
znorm.predictSignalZones()
```

  <details>
    <summary><b><sub>Optional Parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>chromosomes</code> | List of chromosomes to read signal from per bigWig. |
| <code>replace_existing</code> | Set as <code>True</code> to overwrite previously created files. |

  </details>

  <br>

</details>

<a id=""></a>
<details open="open">
  <summary><b>Creating Normalised bigWigs</b></summary>
  Finally, the signal is normalised and outputted as bigWigs. For ZEN normalisation, this defaults to dividing signal by standard deviation within signal zones.
  
```python
# Create normalised bigWigs
znorm.normaliseSignal()
```

  <details>
    <summary><b><sub>Optional Parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>chromosomes</code> | List of chromosomes to read signal from per bigWig. |
| <code>replace_existing</code> | Set as <code>True</code> to overwrite previously created files. |

  </details>

</details>

<br>

---

<a id="norm_compare"></a>
## Plots for Evaluating Normalisation Method Performance
To quantify genome-wide performance across normalisation methods, Wasserstein distribution plots or MA plots can be created.

<br>

<a id=""></a>
<details open="open">
  <summary><b>Initialisation</b></summary>

| Norm    | Sample | bigWig |
| -------- | ------- | ------- |
| RPKM | Cell_Type_A | path/to/bams/cell_type_A_RPKM.bw |
| RPKM | Cell_Type_B | path/to/bams/cell_type_B_RPKM.bw |
| ZEN | Cell_Type_A | path/to/bams/cell_type_A_ZEN.bw |
| ZEN | Cell_Type_B | path/to/bams/cell_type_B_ZEN.bw |
| ... | ... | ... |

```python
bw_df = pd.DataFrame({"norm": np.array([[m] * len(sample_names) for m in norm_methods]).flatten(), 
                      "sample": np.tile(sample_names, len(norm_methods)),
                      "bigwig": bigwig_files})

peak_df = pd.DataFrame({"sample": sample_names,
                        "peaks": peak_files})
```

```python
from ZEN_norm.norm_compare import NormCompare

norm_compare = NormCompare(bigwig_df = bw_df,
                           peaks_df = peak_df,
                           min_peak_score = 0.5,
                           min_consensus = 1,
                           n_cores = znorm.n_cores,
                           analysis_name = = "Analysis_Name",
                           verbose = 2)
```
  
  <br>
</details>

<a id=""></a>
<details open="open">
  <summary><b>Creating Wasserstein Distribution Plots</b></summary>

```python
norm_compare.plotWasserstein(reference_norm = "ZEN", pdf_name = "Wasserstein_Plot.pdf")
```

  <details>
    <summary open="open"><b><sub>Key Parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>reference_norm</code> |  |
| <code>pdf_name</code> |  |

  </details>

  <details>
    <summary open="open"><b><sub>Additional Parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code></code> |  |

  </details>
  
</details>

<a id=""></a>
<details open="open">
  <summary><b>Creating MA Plots</b></summary>

```python
norm_compare.MAPlot(norm_method, plot_samples = [], chromosomes = [], n_cols = 4,
               point_transparency = 0.3, point_colour = "grey", title = "", plot_width = 8, plot_height = 6,
               pdf_name = ""):
```

  <details>
    <summary open="open"><b><sub>Key Parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>norm_method</code> |  |
| <code>pdf_name</code> |  |

  </details>

  <details>
    <summary open="open"><b><sub>Additional Parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>plot_samples</code> |  |
| <code>chromosomes</code> |  |
| <code>n_cols</code> |  |
| <code>point_transparency</code> |  |
| <code>point_colour</code> |  |
| <code>title</code> |  |
| <code>plot_width</code> |  |
| <code>plot_height</code> |  |

  </details>
  
</details>

<br>

---

<a id="reverse_norm"></a>
## Reversing Prior bigWig Normalisation
The module `ReverseNorm` enables raw bigWigs to be created from bigWigs that have previously been normalised. This may be required to avoid double normalisation if renormalising bigWigs with ZEN. 

This works by estimating the coverage value that has been produced by a single read fragment, and dividing signal by this to obtain the coverage that would have been produced has linear normalisation not been applied. 

Sometimes it is not possible or convenient to obtain bigWigs without prior normalisation. For example, if using published data that only provides bigWigs after RPKM normalisation. However, double normalisation will occur if pre-normalised bigWigs are used directly in ZoneNorm. To avoid this, it is best to first reverse normalise pre-normalised bigWigs using ReverseNorm to obtain raw bigWigs.

<br>

```python
from ZEN_norm.reverse_norm import ReverseNorm

rev = ReverseNorm(analysis_name = "Analysis", # Set custom output folder name
                        bigwig_paths = inputs, # Specify a list of bigWig paths
                        n_cores = 1) # Set number of cores to use
```

  <details>
    <summary open="open"><b><sub>Key Parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>analysis_name</code> |  |
| <code>bigwig_paths</code> |  |
| <code>bigwig_directory</code> |  |
| <code>n_cores</code> |  |

  </details>

  <details>
    <summary open="open"><b><sub>Additional Parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>sample_names</code> |  |
| <code>blacklist</code> |  |
| <code>verbose</code> |  |

  </details>

```python
rev.reverseNorm(chunk_size = 10000000, chromosomes = ["chr22"])
```

  <details>
    <summary open="open"><b><sub>Key Parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>chromosomes</code> |  |

  </details>

  <details>
    <summary open="open"><b><sub>Additional Parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>chunk_size</code> |  |
| <code>replace_existing</code> |  |

  </details>

<a id="visualising_bigwigs"></a>
## Visualising bigWigs
To visualise signal both before and / or after normalisation, it is recommended that you view signal on a genome browser such as UCSC or MLV. In addition, regions of bigWig signal can be viewed using ZEN-norm track plots.

<a id=""></a>
<details open="open">
  <summary><b>UCSC Genome browser</b></summary>
  <a href="https://genome.ucsc.edu/index.html">UCSC</a>
</details>

<a id=""></a>
<details open="open">
  <summary><b>Multi Locus View (MLV)</b></summary>
  <a href="https://mlv.molbiol.ox.ac.uk/">MLV</a>
</details>

<a id=""></a>
<details open="open">
  <summary><b>ZEN-norm Signal Plots</b></summary>

```python
znorm.plotTracks(start_coord = 342825, end_coord = 357106, chromosome = "chr5", 
                 plot_samples = znorm.regexFindSamples(regex = "1_neg", ignore_case = True))
```
</details>
