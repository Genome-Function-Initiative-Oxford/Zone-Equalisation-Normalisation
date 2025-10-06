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
  This dataset contains ATAC-seq of E14 mouse embryonic stem cells (mESCs) treated with <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3399494">leukemia inhibitory factor (LIF)</a> and <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3399495">retinoic acid (RA)</a> (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120376">GSE120376</a>). Two bigWigs for these mapped to the mm10 genome can be downloaded from NCBI GEO using the following commands:

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
    <summary open="open"><b><sub>Key parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>analysis_name</code> | The name of the output folder to save results to. |
| <code>bam_paths</code> |  |
| <code>bam_directory</code> |  |
| <code>bigwig_directory</code> |  |
| <code>bam_directory</code> |  |
| <code>n_cores</code> | The number of  |
| <code>norm_method</code> | Name of the normalisation method to apply. Options include: <code>ZEN</code> |

  </details>

  <details>
    <summary open="open"><b><sub>Additional parameters</sub></b></summary>

| Parameter | Usage |
| -------- | ------- |
| <code>chromosomes</code> | List of chromosomes to run analysis on. By default, all autosomal and sex chromosomes will be used. |
| <code>chrom_sizes_file</code> | Path to file of tab separated chromosomes and sizes for creating bigBed or converting wig to bigWig. |
| <code>interleave_sizes</code> | If <code>True</code> and using multiple cores, then process larger chromosomes alongside smaller ones to reduce memory usage. Otherwise process chromosomes from largest to smallest. |
| <code>sample_names</code> |  |
| <code>blacklist</code> | File path to blacklist file with chromosome coordinates to exclude. |
| <code>genome_size</code> | Required if using "RPGC" BAM normalisation. |
| <code>extend_reads</code> | Whether to enable read extension when creating bigWigs from BAMs with deepTools. |
| <code>filter_strand</code> |  |
| <code>exclude_zero</code> |  |
| <code>zone_remove_percent</code> |  |
| <code>norm_stats_type</code> |  |
| <code>norm_power</code> |  |
| <code>deletion_size</code> |  |
| <code>downsample_size</code> |  |
| <code>downsample_seed</code> |  |
| <code>kernel</code> |  |
| <code>test_distributions</code> |  |
| <code>log_transform</code> |  |
| <code>zone_distribution</code> |  |
| <code>zone_param_type</code> |  |
| <code>zone_probability</code> |  |
| <code>bin_size</code> |  |
| <code>extend_n_bins</code> |  |
| <code>merge_depth</code> |  |
| <code>min_region_bps</code> |  |
| <code>quality_filter</code> |  |
| <code>min_different_bps</code> |  |
| <code>verbose</code> |  |

  </details>

</details>

<a id=""></a>
<details open="open">
  <summary><b>Creating Smoothed Signal</b></summary>
  
```python
# Open each bigWig and extract signal to arrays
znorm.readBigWigSignals()
```

</details>

<a id=""></a>
<details open="open">
  <summary><b>Distribution Fitting</b></summary>
  
```python
# Fit distributions for signal zone prediction
znorm.testDistributions()
```

</details>

<a id=""></a>
<details open="open">
  <summary><b>Signal Zone Prediction</b></summary>
  
```python
# Use distribution to predict signal zone coordinates
znorm.predictSignalZones()
```

</details>

<a id=""></a>
<details open="open">
  <summary><b>Creating Normalised bigWigs</b></summary>
  
```python
# Create normalised bigWigs
znorm.normaliseSignal()
```

</details>

<br>

---

<a id="norm_compare"></a>
## Plots for Evaluating Normalisation Method Performance

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

### Creating Wasserstein Plots
```python
norm_compare.plotWasserstein(reference_norm = "ZEN", pdf_name = "Wasserstein_Plot.pdf")
```

### Creating MA Plots
```python
norm_compare.MAPlot(norm_method, plot_samples = [], chromosomes = [], n_cols = 4,
               point_transparency = 0.3, point_colour = "grey", title = "", plot_width = 8, plot_height = 6,
               pdf_name = ""):
```

<br>

---

<a id="reverse_norm"></a>
## Reversing Prior bigWig Normalisation
Sometimes it is not possible or convenient to obtain bigWigs without prior normalisation. For example, if using published data that only provides bigWigs after RPKM normalisation. However, double normalisation will occur if pre-normalised bigWigs are used directly in ZoneNorm. To avoid this, it is best to first reverse normalise pre-normalised bigWigs using ReverseNorm to obtain raw bigWigs.

<br>

```python
from ZEN_norm.reverse_norm import ReverseNorm

rev = ReverseNorm(analysis_name = "Analysis", # Set custom output folder name
                        bigwig_paths = inputs, # Specify a list of bigWig paths
                        n_cores = 1) # Set number of cores to use
rev.reverseNorm(chunk_size = 10000000, chromosomes = ["chr22"])
```

<a id="visualising_bigwigs"></a>
## Visualising bigWigs
UCSC Genome browser or Multi Locus View

<a id=""></a>
<details open="open">
  <summary><b>ZEN-norm Signal Plots</b></summary>

```python
znorm.plotTracks(start_coord = 342825, end_coord = 357106, chromosome = "chr5", 
                 plot_samples = znorm.regexFindSamples(regex = "1_neg", ignore_case = True))
```
</details>

<a id=""></a>
<details open="open">
  <summary><b>UCSC Genome browser</b></summary>
  ...
</details>

<a id=""></a>
<details open="open">
  <summary><b>Multi Locus View (MLV)</b></summary>
  ...
</details>
