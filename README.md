# Zone Equalisation Normalisation

ZEN-norm is a Python package for normalising bigWigs by Zone Equilisation Normalisation. Additional features include modules for reversing prior bigWig normalisation and comparing normalisation methods using Wasserstein distance plots.

---

<br>

<details open="open">
  <summary><b>Contents</b></summary>
  <ol>
    <li><a href="#installation">Installation</a></li>
    <li><a href="#zen_norm">ZEN bigWig Normalisation</a></li>
    <li><a href="#norm_compare">Plots for Evaluating Normalisation Method Performance</a></li>
    <li><a href="#reverse_norm">Reversing Prior bigWig Normalisation</a></li>
  </ol>
</details>

---

<a id="installation"></a>
## Installation

<details open="open">
  <summary><b>PyPI Installation</b></summary>
  To install the ZEN-norm package from PyPI, run the command below:

```
python -m pip install  --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple ZEN-norm-test
```
</details>

<details open="open">
  <summary><b>Conda Installation</b></summary>
  To install the ZEN-norm package from Conda, run the command below:

```
TBC
```

  Alternatively, if there are issues installing ZEN-norm, a conda environment with the required packages can be created using the <code>zen_environment.yml</code> file.

  ```
  conda env create --name zen_env --file=environment/zen_environment.yml
  ```
</details>

<a id="zen_norm"></a>
## ZEN bigWig Normalisation
The module `ZoneNorm` is for normalising genomic coverage with ZEN. Steps include: BAM to bigWig mapping, creating smoothed signals, distribution fitting, signal zone prediction and creating normalised bigWigs.

<a id="zen_norm_inputs"></a>
### Specifying Input BAMs or bigWigs
ZEN-norm supports either BAM or bigWig files as an input. These contain the genomic signal of interest per sample.

```python
bam_paths = ["path/to/bams/cell_type_A.bam", "path/to/bams/cell_type_B.bam"]
```

If bigWigs have been pre-normalised, then it is advisable to remap them without normalisation, or to use ZEN-norm's <a href="#reverse_norm">reverse bigWig normalisation</a>.

```python
bigwig_paths = ["path/to/bams/cell_type_A.bw", "path/to/bams/cell_type_B.bw"]
```

### Initialisation of ZoneNorm Object


```python
from ZEN_norm.zone_norm import ZoneNorm

# Create ZoneNorm object
znorm = ZoneNorm(analysis_name = "Analysis_Name",
                 bigwig_paths = input_paths,
                 n_cores = n_cores,
                 norm_method = "ZEN")
```

```python
# Open each bigWig and extract signal to arrays
znorm.readBigWigSignals()
```

```python
# Fit distributions for signal zone prediction
znorm.testDistributions()
```

```python
# Use distribution to predict signal zone coordinates
znorm.predictSignalZones()
```

```python
# Create normalised bigWigs
znorm.normaliseSignal()
```

<a id="norm_compare"></a>
## Plots for Evaluating Normalisation Method Performance

### Initialisation

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

<a id="reverse_norm"></a>
## Reversing Prior bigWig Normalisation
Sometimes it is not possible or convenient to obtain bigWigs without prior normalisation. For example, if using published data that only provides bigWigs after RPKM normalisation. However, double normalisation will occur if pre-normalised bigWigs are used directly in ZoneNorm. To avoid this, it is best to first reverse normalise pre-normalised bigWigs using ReverseNorm to obtain raw bigWigs.

```python
from ZEN_norm.reverse_norm import ReverseNorm

rev = ReverseNorm(analysis_name = "Analysis", # Set custom output folder name
                        bigwig_paths = inputs, # Specify a list of bigWig paths
                        n_cores = 1) # Set number of cores to use
rev.reverseNorm(chunk_size = 10000000, chromosomes = ["chr22"])
```
