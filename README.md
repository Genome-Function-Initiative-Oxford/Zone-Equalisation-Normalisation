# Zone Equalisation Normalisation

ZEN-norm is a Python package for normalising bigWigs of genomic signal, such as ATAC-seq, ChIP-seq and TT-seq by Zone Equilisation Normalisation (ZEN). This also includes modules for reversing prior bigWig normalisation and creating plots to compare performance of normalisation methods genome-wide.

<p><img src="https://raw.githubusercontent.com/Genome-Function-Initiative-Oxford/ZEN-norm/refs/heads/assets/Images/ZEN_Overview_Figure.png?token=GHSAT0AAAAAADLQKKO3R5BOS4C5USDIGE2I2IPOK5Q" width="100%"></p>

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
python -m pip install ZEN-norm
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
  A detailed Jupyter notebook tutorial explaining how to use ZEN for reversing prior normalisation, normalising bigWigs with ZEN and evaluating normalisation methods via Wasserstein distance plots and MA plots is provided within the <a href="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/tree/main/tutorials/zen_tutorial">tutorial/zen_tutorials</a> folder of the GitHub repository. For a quick overview of the avaliable features, see <a href="#reverse_norm">sections 3 to 6</a> below.
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
Module `ZoneNorm` normalises genomic coverage with ZEN. Steps include: BAM to bigWig mapping, convolution to create smoothed signals, distribution fitting, signal zone prediction (coordinates of consistent regions of signal) and generating bigWigs normalised with ZEN. It runs on genomic signals from either BAMs or bigWigs. If using bigWigs that have been pre-normalised, then it is advisable to first remap them without normalisation, or to use ZEN-norm's <a href="#reverse_norm">ReverseNorm</a> module.

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

<a id=""></a>
<details open="open">
  <summary><b>Wasserstein Distance Plots</b></summary>
  Within a Wasserstein distance plot, min-max scaled pairwise sample Wasserstein distance (w) is measured over regions (e.g. peaks or zones) and plotted as violin and / or box plots per normalisation method. The normalisation method with the lowest average w therefore has the best alignment across the genome for regions of interest. For example in the plot below of erythroid ATAC-seq, ZEN has the lowest mean w and it is significantly lower than all other normalisation methods according to a t-test comparing the distributions.
  
<p><img src="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/blob/assets/Images/Erythroid_ATAC_Wasserstein_Plot.jpg" width="70%"></p>

</details>

<a id=""></a>
<details open="open">
  <summary><b>MA Plots</b></summary>
  MA plots compare differences in total signal (M) relative to signal intensity (A) between pairs of samples across regions (e.g. peaks or zones). This is useful to assess how effective a normalisation method is at reducing bias. For example in the HeLa TT-seq plot below, a point in each subplot represents the mean count of the signal over a region's coordinates for two samples after RPKM normalisation. The dotted red line is a reference, whereby the closer the points fall, the closer the average counts are for the samples.

<p><img src="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/blob/assets/Images/HeLa_TTseq_Reverse_MA_Plot.png" width="100%"></p>

</details>

<br>

<a id="visualising_signal"></a>
## 6. Visualising Signal and Zones
When normalising genomic signal with ZEN, using bigWigs as inputs and outputs allows some steps in the process to be visualised providing greater transparency than count based normalisation methods. For example, bigWigs can be saved after reverse normalisation, smoothing via convolution (useful to view thresholds against this signal) and after normalisation with ZEN. In additon, predicted signal zones (coordinates of consistent regions of signal) can be saved to BED files and visualised in the same way as peak calls. These signals can therefore be viewed using either track plots included in the package, or genome browsers.

<a id=""></a>
<details open="open">
  <summary><b>Track Plots</b></summary>
  When running ZEN-norm, regions of signal for one or more samples can be viewed as track plots. These are customisable as demonstrated in the examples below: 

<p></p>

  <ins>Viewing Non-Normalised Signals</ins>

  Signal from one or more samples can be overlaid to view alignment prior to normalisation.
  
<p><img src="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/blob/assets/Images/Erythroid_ATAC_Raw_Signal_Plot.png" width="100%"></p>
  
  <ins>Comparing Zone Thresholds Against Convoluted Signal</ins>

  After smoothing signal via convolution and distribution fitting, thresholds are derived from the distributions. These can be viewed against the convoluted signal to see how a threshold will separate signal from background noise during signal zone prediction.
  
<p><img src="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/blob/assets/Images/Erythroid_ATAC_Convoluted_Signal_Plot.png" width="100%"></p>

  <ins>Viewing Signal Zones</ins>

  After signal zone prediction, zones can be visualised as bars within the track plots.

<p><img src="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/blob/assets/Images/Erythroid_ATAC_Zone_Signal_Plot.png" width="100%"></p>

  <ins>Viewing Normalised Signals</ins>
  
  Signal from one or more samples can be overlaid to view alignment after normalisation.
  
<p><img src="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/blob/assets/Images/Erythroid_ATAC_ZEN_Signal_Plot.png" width="100%"></p>

</details>

<a id=""></a>
<details open="open">
  <summary><b>Genome Browsers</b></summary>
  After saving signals to bigWigs and zones as BED files, they can be uploaded to an interactive genome browser such as <a href="https://genome.ucsc.edu/index.html">UCSC Genome Browser</a> or <a href="https://mlv.molbiol.ox.ac.uk/">Multi Locus View (MLV)</a>.
</details>
