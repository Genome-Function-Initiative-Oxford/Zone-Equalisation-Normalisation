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

```python
from ZEN_norm.zone_norm import ZoneNorm

# Create ZoneNorm object
znorm = ZoneNorm(analysis_name = "Analysis_Name",
                 bigwig_paths = input_paths,
                 n_cores = n_cores,
                 norm_method = "ZEN")

# Open each bigWig and extract signal to arrays
znorm.readBigWigSignals()
# Fit distributions for signal zone prediction
znorm.testDistributions()
# Use distribution to predict signal zone coordinates
znorm.predictSignalZones()
# Create normalised bigWigs
znorm.normaliseSignal()
```

## Evaluating Normalisation Methods
```python
from ZEN_norm.norm_compare import NormCompare
```

## Reversing Prior Normalisation
Sometimes it is not possible or convenient to obtain bigWigs without prior normalisation. For example, if using published data that only provides bigWigs after RPKM normalisation. However, double normalisation will occur if pre-normalised bigWigs are used directly in ZoneNorm. To avoid this, it is best to first reverse normalise pre-normalised bigWigs using ReverseNorm to obtain raw bigWigs.

```python
from ZEN_norm.reverse_norm import ReverseNorm
```
