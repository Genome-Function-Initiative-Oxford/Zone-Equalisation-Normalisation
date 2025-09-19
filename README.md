# Zone Equalisation Normalisation

ZEN-norm is a Python package for normalising bigWigs and evaluating performance. Modules normalisation of bigWigs with ZEN (ZoneNorm), reverse normalisation (ReverseNorm) and normalisation comparision (NormCompare).

<br>

---

## Environment Set Up
The conda environment to run ZEN can be created by running the command:
```
conda env create --name zen_env --file=zen_environment.yml
```

## Normalisation With ZEN

```
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

## Reversing Prior Normalisation
Sometimes it is not possible or convenient to obtain bigWigs without prior normalisation. For example, if using published data that only provides bigWigs after RPKM normalisation. However, double normalisation will occur if pre-normalised bigWigs are used directly in ZoneNorm. To avoid this, it is best to first reverse normalise pre-normalised bigWigs using ReverseNorm to obtain raw bigWigs.
