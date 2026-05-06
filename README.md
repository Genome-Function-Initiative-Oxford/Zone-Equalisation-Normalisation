# Zone Equalisation Normalisation
Zone Equilisation Normalisation (ZEN) is a method for normalising genomic signal bigWigs, such as ATAC-seq, ChIP-seq and TT-seq. ZEN is avaliable within the Python package ZEN-norm, which also includes modules for reversing prior bigWig normalisation and creating plots to compare performance of normalisation methods genome-wide.

<p><img src="https://github.com/Genome-Function-Initiative-Oxford/Zone-Equalisation-Normalisation/blob/assets/Images/ZEN_Overview_Figure.png?raw=True" width="100%"></p>

**Citation:** T. Wilson, TA. Milne, SG. Riva and JR. Hughes, <a href="https://www.biorxiv.org/content/10.64898/2025.12.10.693203v1">_Zone Equalisation Normalisation For Improved Alignment of Epigenetic Signal_</a>, bioRvix, 2025

<br>

---

<details open="open">
  <summary><b>Contents</b></summary>
  <ol>
    <li><a href="#installation">Installation</a></li>
    <li><a href="#tutorials">Tutorials</a></li>
  </ol>
</details>

<br>

---

<a id="installation"></a>
## 1. Installation
ZEN-norm is designed to run on Python 3.10 and above. It is installable from either Conda (recommended) or PyPI.

<a id=""></a>
<details open="open">
  <summary><b>Conda Installation</b></summary>
  To install the ZEN-norm package from <a href="https://anaconda.org/TomMakesThings/zen-norm">Conda</a>, active the conda environment you'd like to install the package into (<code>conda activate [environment_name]</code>) and run the command below:

```
conda install tommakesthings::zen-norm
```

  Or alternatively:

```
conda install -c https://conda.anaconda.org/TomMakesThings zen-norm
```
  
  If there are issues running either of the above, a conda environment with the required packages can be created with the <code>zen_environment.yml</code> file.

  ```
  conda env create --name zen_env --file=environment/zen_environment.yml
  conda activate zen_env
  ```
</details>

<a id=""></a>
<details open="open">
  <summary><b>PyPI Installation</b></summary>
  To install the ZEN-norm package from <a href="https://pypi.org/project/ZEN-norm/">PyPI</a>, run the command below:

```
pip install ZEN-norm
```

  Some ZEN-norm features require samtools, which is not distributed via PyPI and must be installed separately. 
</details>

<br>

---

<a id="tutorials"></a>
## 2. Tutorials
See the <a href="https://github.com/Genome-Function-Initiative-Oxford/Zone-Equalisation-Normalisation/tree/main/tutorials">Tutorials</a> page on GitHub for explanations of how to use ZEN-norm to:
* Normalise bigWigs with ZEN
* Reverse prior bigWig normalisation
* Compare and quanity normalisation method genome-wide alignment via Wasserstein distance
* Create plots included in the publication
