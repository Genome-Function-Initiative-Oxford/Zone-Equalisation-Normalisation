# Supplementary Figures

## About
<a href="">Supplementary_Figures.ipynb</a> is a Jupyter notebook with examples of how module `ZoneNorm` can be used for parameter tuning of the distribution fitting step. This includes creating the Kolmogorov–Smirnov test boxplots within the Supplementary section of the ZEN publication.


## Example Data
This notebook uses all the datasets detailed in the <a href="https://github.com/Genome-Function-Initiative-Oxford/Zone-Equalisation-Normalisation/tree/main/tutorials/zen_tutorial">main tutorial</a>, as well as additional datasets explained below. For these (except CATlas), BAMs can be created by downloaded FASTQ from NCBI GEO and aligning them to the hg38 genome using the <a href="https://github.com/Genome-Function-Initiative-Oxford/UpStreamPipeline/tree/main/genetics/CATCH-UP">CATCH-UP pipeline</a>.

<a id=""></a>
<details open="open">
  <summary><b>H1 Brain ATAC-seq</b></summary>
  This data by CE. Dundes et al. contains ATAC-seq of brain differentiation from the H1 cell line to forebrain, midbrain and hindbrain, with three replicates each (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286146">GSE286146</a>).

</details>

<a id=""></a>
<details open="open">
  <summary><b>CATlas Single-Cell ATAC-seq</b></summary>
  222 scATAC-seq bigWigs from human tissues can be downloaded from the <a href="https://decoder-genetics.wustl.edu/catlasv1/catlas_downloads/humantissues/Bigwig/">CATlas portal</a>. Note these need to be reverse normalised as they are pre-normalised by Signal Per Million Reads (SPMR).
  
</details>

<a id=""></a>
<details open="open">
  <summary><b>Erythroid CTCF ChIP-seq</b></summary>
  This data by E. Georgiades et al. includes CTCF ChIP-seq from three donors (donors 1, 2 and 30) without replicates (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244929">GSE244929</a>).

</details>

<a id=""></a>
<details open="open">
  <summary><b>HEL Pol II ChIP-seq</b></summary>
  This data by RG. Roeder et al. contains RNA Polymerase II Ser2P ChIP-seq from the HEL leukaemia cell line with three replicates (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84157">GSE84157</a>).

</details>

<a id=""></a>
<details open="open">
  <summary><b>A-375 TT-seq</b></summary>
  This data by ML. Insco et al. consists of TT-seq from the A-375 melanoma cell line of a CDK13 mutant clone R860Q and a control (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223888">GSE223888</a>).

</details>

<a id=""></a>
<details open="open">
  <summary><b>HEK293T TT-seq</b></summary>
  This data by CA. Mimoso et al. includes TT-seq of the HEK293T embryonic kidney cell line with four replicates (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218127">GSE218127</a>).

</details>






