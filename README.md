# Zone Equalisation Normalisation

ZEN-norm is a Python package for normalising bigWigs by Zone Equilisation Normalisation. Additional features include modules for reversing prior bigWig normalisation and comparing normalisation methods using Wasserstein distance plots.

<br>

---

## Installation
The ZEN-norm package can either be installed from PyPI.

### PyPI
```
python -m pip install --no-cache-dir --prefer-binary --only-binary=pybigwig \
  --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple \
  ZEN-norm-test==0.0.5
```

### Conda
The conda environment to run ZEN-norm can be created by running the command:
```
conda env create --name zen_env --file=zen_environment.yml
```

## Normalisation With ZEN
The module `ZoneNorm` is designed for normalising genomic coverage with ZEN. Steps include: ...

### Specifying Input BAMs or BigWigs


```python
bam_paths = ["path/to/bams/cell_type_A.bam", "path/to/bams/cell_type_B.bam"]
```

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

## Evaluating Normalisation Methods

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
norm_compare.plotWasserstein(chromosomes = [], plot_samples = [], norm_methods = [], use_chrom_maxs = False,
                             reference_norm = "ZEN", log_scale = False, pdf_name = "TTseq_Gumbell")
```

### Creating MA Plots
```python
norm_compare.MAPlot(norm_method, plot_samples = [], chromosomes = [], n_cols = 4,
               point_transparency = 0.3, point_colour = "grey", title = "", plot_width = 8, plot_height = 6,
               pdf_name = ""):
```

## Reversing Prior Normalisation
Sometimes it is not possible or convenient to obtain bigWigs without prior normalisation. For example, if using published data that only provides bigWigs after RPKM normalisation. However, double normalisation will occur if pre-normalised bigWigs are used directly in ZoneNorm. To avoid this, it is best to first reverse normalise pre-normalised bigWigs using ReverseNorm to obtain raw bigWigs.

```python
from ZEN_norm.reverse_norm import ReverseNorm

rev = ReverseBigWigNorm(analysis_name = "Analysis", # Set custom output folder name
                        bigwig_paths = inputs, # Specify a list of bigWig paths
                        n_cores = 1) # Set number of cores to use
rev.reverseNorm(chunk_size = 10000000, chromosomes = ["chr22"])
```
