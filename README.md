# Zone Equalisation Normalisation

TBC

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

## Evaluation Normalisation Methods

## Reversing Prior Normalisation
