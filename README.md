# Zone Equalisation Normalisation

ZEN-norm is a Python package for normalising bigWigs of genomic signal, such as ATAC-seq, ChIP-seq and TT-seq by Zone Equilisation Normalisation (ZEN). This also includes modules for reversing prior bigWig normalisation and creating plots to compare performance of normalisation methods genome-wide.

<p><img src="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/blob/assets/Images/ZEN_Overview_Figure.png" width="100%"></p>

---

<details open="open">
  <summary><b>Contents</b></summary>
  <ol>
    <li><a href="#installation">Installation</a></li>
    <li><a href="#reverse_norm">Reversing Prior bigWig Normalisation</a></li>
    <li><a href="#zen_norm">Normalising bigWigs With ZEN</a></li>
    <li><a href="#norm_compare">Evaluating Normalisation Method Performance</a></li>
  </ol>
</details>

<br>

---

<a id="installation"></a>
## Installation
ZEN-norm is designed to run on Python 3.10 and above. It is installable from either PyPI or Conda.

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

<a id="reverse_norm"></a>
## Reversing Prior bigWig Normalisation
The module `ReverseNorm` provides an optional step to enable non-normalised bigWigs to be created from bigWigs that have previously been normalised. This can be used to avoid double normalisation if renormalising bigWigs with ZEN. It works by estimating the coverage value that has been produced by a single read fragment, and dividing signal by this to obtain the coverage that would have been produced has linear normalisation not been applied. 

<br>


<a id=""></a>
<details open="open">
  <summary><b>Quick Example</b></summary>

```python
from ZEN_norm.reverse_norm import ReverseNorm

rev = ReverseNorm(analysis_name = "Example_Analysis", # Set custom output folder name
                  bigwig_paths = ["path/sample_A.bw", "path/sample_B.bw"], # Specify a list of bigWig paths
                  n_cores = 8) # Set number of cores to use
rev.reverseNorm(chromosomes = ["chr19"])
```

</details>

<br>

---

<br>

<a id="zen_norm"></a>
## Normalising bigWigs With ZEN
The module `ZoneNorm` is normalises genomic coverage with ZEN. Steps include: BAM to bigWig mapping, convolution to create smoothed signals, distribution fitting, signal zone prediction and generating bigWigs normalised with ZEN. It runs on genomic signals from either BAM or bigWig files. If your input files are bigWigs that have been pre-normalised, then it is advisable to first remap them without normalisation, or to use ZEN-norm's <a href="#reverse_norm">reverse bigWig normalisation</a> feature.

<br>

<a id=""></a>
<details open="open">
  <summary><b>Quick Example</b></summary>

```python

```

</details>

<br>

---

<a id="norm_compare"></a>
## Evaluating Normalisation Method Performance
To quantify genome-wide performance across normalisation methods, Wasserstein distribution plots or MA plots can be created.

<a id=""></a>
<details open="open">
  <summary><b>Quick Example</b></summary>

```python

```

</details>

<br>

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
