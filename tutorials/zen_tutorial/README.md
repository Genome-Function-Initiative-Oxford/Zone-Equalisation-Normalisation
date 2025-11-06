# ZEN Tutorial

## About
The <a href="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/blob/main/tutorials/zen_tutorial/ZEN_Tutorial.ipynb">ZEN Tutorial</a> notebook contains detailed explainations for running bigWig normalisation with ZEN with module `ZoneNorm`, reversing prior bigWig normalisation with `ReverseNorm` and comparing performance of normalisation methods genome-wide with `CompareNorm`. Examples are given across a range of genomic assays using the datasets below.

<a id="example_data"></a>
## Example Data
The following datasets are all publically avaliable and files of URLs can be found in this repository under <a href="https://github.com/Genome-Function-Initiative-Oxford/ZEN-norm/tree/main/data">data</a>.

<a id=""></a>
<details open="open">
  <summary><b>Erythroid ATAC-seq</b></summary>
  This data consists of bulk ATAC-seq from two healthy day 13 erythroid donors (donors 2 and 3), each with seven technical replicates. Either BAMs (as used in the tutorial), or non-normalised bigWigs can be downloaded by command line:

  <br>
  

  <ins>BAMs</ins>
```
wget -i tests/data/Erythroid_ATAC_BAM_URLs.txt
```

  <ins>Non-Normalised bigWigs</ins>
```
wget -i tests/data/Erythroid_ATAC_bigWig_URLs.txt
```

<br>

</details>

<a id=""></a>
<details open="open">
  <summary><b>Erythroid RAD21 ChIP-seq</b></summary>
  This data by E. Georgiades et al (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244929">GSE244929</a>) consists of RAD21 ChIP-seq from three healthy day 13 erythroid donors (donors 1, 2, and 30), each with three replicates. Either BAMs, or non-normalised bigWigs (as used in the tutorial) can be downloaded by command line:

  <ins>BAMs</ins>
```
wget -i tests/data/Erythroid_RAD21_BAM_URLs.txt
```

  <ins>Non-Normalised bigWigs</ins>
```
wget -i tests/data/Erythroid_RAD21_bigWig_URLs.txt
```

<br>

</details>

<a id=""></a>
<details open="open">
  <summary><b>HeLa TT-seq</b></summary>
  This data by A. Fiszbein et al. (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284682">GSE284682</a>) contains TT-seq from HeLa cells treated with either U1-AMO or a control, each with two replicates. Either BAMs (as used in the tutorial), or non-normalised bigWigs split into forward (positive) and reverse (negative) strands can be downloaded by command line:

  <ins>BAMs</ins>
```
wget -i tests/data/HeLa_TTseq_BAM_URLs.txt
```

  <ins>Non-Normalised bigWigs</ins>
```
wget -i tests/data/HeLa_TTseq_bigWig_URLs.txt
```

<br>

</details>

<a id=""></a>
<details open="open">
  <summary><b>Mouse Embryonic Stem Cell ATAC-seq</b></summary>
  This ENCODE dataset (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120376">GSE120376</a>) is used to demonstrate how pre-normalised bigWigs can be reverse normalised prior to normalisation with ZEN. It contains ATAC-seq from E14 mouse embryonic stem cells (mESCs) treated with leukemia inhibitory factor (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3399494">LIF</a>) and retinoic acid (<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3399495">RA</a>). Two bigWigs mapped to the mm10 genome can be downloaded from NCBI GEO using the following commands.

  <br>
  
<ins>Pre-Normalised bigWigs</ins>

```
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3399495&format=file&file=GSM3399495%5FE14%5FATAC%5FRA%2Ebw" -O E14_ATAC_RA.bw
```
```
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3399494&format=file&file=GSM3399494%5FE14%5FATAC%5FLIF%2Ebw" -O E14_ATAC_LIF.bw
```

</details>





