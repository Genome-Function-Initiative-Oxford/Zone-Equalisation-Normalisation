# Zone Equalisation Normalisation

ZEN-norm is a Python package for normalising bigWigs of genomic signal, such as ATAC-seq, ChIP-seq and TT-seq by Zone Equilisation Normalisation (ZEN). This also includes modules for reversing prior bigWig normalisation and creating plots to compare performance of normalisation methods genome-wide.

<p><img src="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/blob/assets/Images/ZEN_Overview_Figure.png" width="100%"></p>

**Citation:** T. Wilson, TA. Milne, SG. Riva and JR. Hughes, _Zone Equalisation Normalisation For Improved Alignment of Epigenetic Signal_, Unpublished, 2025

<br>

---

<details open="open">
  <summary><b>Contents</b></summary>
  <ol>
    <li><a href="#installation">Installation</a></li>
    <li><a href="#tutorials">Tutorials</a></li>
    <li><a href="#reverse_norm">Reversing Prior bigWig Normalisation</a></li>
    <li><a href="#zone_norm">Normalising bigWigs With ZEN</a></li>
    <li><a href="#compare_norm">Evaluating Normalisation Method Performance</a></li>
  </ol>
</details>

<br>

---

<a id="installation"></a>
## 1. Installation
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

<a id="tutorials"></a>
## 2. Tutorials

<a id=""></a>
<details open="open">
  <summary><b>Main Tutorial</b></summary>
  A detailed Jupyter notebook tutorial explaining how to use ZEN for reversing prior normalisation, normalising bigWigs with ZEN and evaluating normalisation methods via Wasserstein distance plots and MA plots is provided within the <a href="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/tree/main/tutorials/zen_tutorial">tutorial/zen_tutorials</a> folder of this repository. For a quick overview of the avaliable features, see <a href="#reverse_norm">sections 3 to 6</a> below.
</details>

<a id=""></a>
<details open="open">
  <summary><b>Publication Supplementary Figures</b></summary>
  A Jupyter notebook is provided in folder <a href="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/tree/main/tutorials/zen_tutorial">tutorial/supplementary_figures</a> to document how the figures were created in the Supplementary section of the ZEN publication.
</details>

<br>

---

<a id="reverse_norm"></a>
## 3. Reversing Prior bigWig Normalisation
Module `ReverseNorm` provides an optional step to enable non-normalised bigWigs to be created from pre-normalised bigWigs. This is not required if BAMs are available, but is designed to avoid double normalisation if renormalising bigWigs with ZEN. It works by estimating the coverage value that has been produced by a single read fragment and dividing signal by this to obtain the coverage that would have been produced had linear normalisation (e.g. RPKM, CPM) not been applied.

<br>


<a id=""></a>
<details open="open">
  <summary><b>Basic Example</b></summary>

```python
from ZEN_norm.reverse_norm import ReverseNorm

rev = ReverseNorm(analysis_name = "Example_Analysis", # Set custom output folder name
                  bigwig_paths = ["path/to/normalised_bigwigs/sample_A.bw", "path/to/normalised_bigwigs/sample_B.bw"], # Specify a list of bigWig paths
                  n_cores = 8) # Set number of cores to use
rev.reverseNorm(chromosomes = ["chr19"])
```

</details>

<br>

---

<a id="zone_norm"></a>
## 4. Normalising bigWigs With ZEN
Module `ZoneNorm` normalises genomic coverage with ZEN. Steps include: BAM to bigWig mapping, convolution to create smoothed signals, distribution fitting, signal zone prediction and generating bigWigs normalised with ZEN. It runs on genomic signals from either BAMs or bigWigs. If using bigWigs that have been pre-normalised, then it is advisable to first remap them without normalisation, or to use ZEN-norm's <a href="#reverse_norm">ReverseNorm</a> module.

<br>

<a id=""></a>
<details open="open">
  <summary><b>Basic Example</b></summary>

```python
# EITHER create bigWigs without normalisation from BAMs
znorm = ZoneNorm(analysis_name = "Example_Analysis", # Name of output folder
                 bam_paths = ["path/to/bams/sample_A.bam", "path/to/bams/sample_B.bam"], # List or directory of BAM files
                 n_cores = 8, # Number of processors
                 extend_reads = True, # Whether to extend reads during BAM to bigWig mapping (False recommended for transcriptional assays)
                 filter_strand = False) # Whether to separate by strand (True recommended for transcriptional assays)

# OR set bigWigs directly
znorm = ZoneNorm(analysis_name = "Example_Analysis", # Name of output folder
                 bigwig_paths = ["path/to/raw_bigwigs/sample_A.bw", "path/to/raw_bigwigs/sample_B.bw"], # List or directory of bigWig files
                 n_cores = 8) # Number of processors

# Create smoothed signal
znorm.convolveSignals()
# Test Laplace distribution
znorm.testDistributions()
# Use distribution to predict signal zone coordinates
znorm.predictSignalZones()
# Create normalised bigWigs
znorm.normaliseSignal()
```

</details>

<br>

---

<a id="compare_norm"></a>
## 5. Evaluating Normalisation Method Performance
To quantify genome-wide performance across normalisation methods, Wasserstein distribution plots or MA plots can be created.

<p><img src="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/blob/assets/Images/Erythroid_ATAC_Wasserstein_Plot.jpg" width="70%"></p>

<br>

<a id="visualising_bigwigs"></a>
## 6. Visualising bigWigs
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
