---
title: "Analyzing Illumina RNA-seq Data with the CRI HPC Using The Dataset From The Project CRI-BIO-646-BMB-RKeenan"
author: '[Wen-Ching Chan](http://cri.uchicago.edu/people/#chan)'
date: "Nov 2018"
output:
  html_document:
    keep_md: yes
package: CRI-BIO-646-BMB-RKeenan
---




<!--
* Notes
    - add [MultiQC](http://multiqc.info/)
    - change grch38/hg38
    - incorporate other pipeline (ChIP-seq, WGBS, WGS/WXS)
    - Docker + WDL/CWL
-->


## <a name="Top"/>


## CRI RNAseq Pipelines
- [CRI-RNAseq-2016](https://github.com/riyuebao/CRI-Workshop-Nov2016-RNAseq/blob/master/Run_RNAseq.tutorial.ipynb) by [Riyue Bao](http://cri.uchicago.edu/people/#bao). Last modified on **November 12, 2016**.
- [CRI-RNAseq-2014](https://wiki.uchicago.edu/pages/viewpage.action?pageId=95855827) by [Lei Huang](http://cri.uchicago.edu/people/#huang). Last modified on **Apr 29, 2014**.


## **<span style="color:red">RUN AS PRACTICE</span>**

> 
> **Pipeline package is available on <span style="color:red">[CRI Bioinformatics FTP](ftp://logia.cri.uchicago.edu/bioinformatics/tutorials/Nov2018/CRI-BIO-646-BMB-RKeenan.tgz)</span>**
> 

The [Center for Research Informatics](http://cri.uchicago.edu/) (CRI) provides computational resources and expertise in biomedical informatics for researchers in the Biological Sciences Division (BSD) of the University of Chicago. 

<!--
This workshop is part of a series of [monthly training events](http://cri.uchicago.edu/seminar-series/) focusing on using the University of Chicago’s computational resources to analyze Next-Generation Sequencing and Microarray data.
-->

As a [bioinformatics core](http://cri.uchicago.edu/bioinformatics/), we are actively improving our pipelines and expanding pipeline functions. The tutorials will be updated in a timely manner but may not reflect the newest updates of the pipelines. Stay tuned with us for the latest pipeline release.

If you have any questions, comments, or suggestions, feel free to contact our core at bioinformatics@bsd.uchicago.edu or one of our bioinformaticians.


* [Quick Start](#Quick)
* [Introduction](#Intro)
* [Navigation Map](#Navi)
* [Work Flow](#Work)
* [Data Description](#Data)
* [Prerequisites](#Pre)
    + [Login and Setup Tutorial Working Directory](#Login)
* [Pipeline Steps](#Steps)
    + Step 1: [Quality Control](#ReadQC)
    + Step 2-1: [Read Alignment](#Aln)
    + Step 2-2: [Alignment QC](#AlnQC)
    + Step 3: [Expression Quantification](#Quant)
    + Step 4-1: [Identification of Differentially Expressed Genes (DEGs)](#DEG)
    + Step 4-2: [DEG Statistics](#LociStat)
    + Step 5: [Sample Correlation](#QuantQC)
    + Step 6: [Heat Map](#HeatMap)
    + Step 7: [Functional Enrichment Analysis](#GSEA)
* [BigDataScript Report](#BDS)
* [Reference](#Ref)



## <a name="Quick"/>Quick Start | [Top](#Top)

1. modify the generator script **Build_RNAseq.CRI-BIO-646.sh** accordingly
    1. project="PROJECT_AS_PREFIX" (e.g., **CRI-BIO-646** which is used as a prefix of metadata file **CRI-BIO-646**.metadata.txt and configuration file **CRI-BIO-646**.pipeline.yaml)
    2. padding="DIRECTORY_NAME_CONTAINING_PROJECT_DATA" (e.g., **CRI-BIO-646** which is the folder name to accommodate metadata file, configuration file, sequencing data folder, and references folder)
2. prepare metadata file as the example **[CRI-BIO-646.metadata.txt](CRI-BIO-646.metadata.txt)**
    1. **Single End (SE)** Library
        1. Set Flavor column as 1xReadLength (e.g., 1x50)
        2. Set Seqfile1 column as the file name of the repective sequencing file
    2. **Paired End (PE)** Library
        1. Set Flavor column as 2xReadLength (e.g., 2x50)
        2. Set Seqfile1 column as the file name of the repective read 1 (R1) sequencing file
        3. Set an additional column named '**Seqfile2**' as the file name of the repective read 2 (R2) sequencing file
    3. **Non strand-specific** Library
        1. Set LibType column to NS
    4. **Strand-specific** Library
        1. Inquire the library type from your seuqencing center and set LibType column to FR (the left-most end of the fragment (in transcript coordinates, or the first-strand synthesis) is the first sequenced) or RF (the right-most end of the fragment (in transcript coordinates) is the first sequenced, or the second-strand synthesis). You can read this [blog](http://onetipperday.sterding.com/2012/07/how-to-tell-which-library-type-to-use.html) for more details of strand-specific RNA-seq.
3. prepare reference files as the example reference hg38 under **[/gpfs/data/bioinformatics/cri_rnaseq_2018_ex/example/references/v28_92_GRCh38.p12)**
4. prepare pre-built STAR index files as the example reference hg38 under **[/gpfs/data/bioinformatics/cri_rnaseq_2018_ex/example/references/v28_92_GRCh38.p12/STAR)**



## <a name="Intro"/>Introduction | [Top](#Top)


RNA sequencing (RNA-seq) is a revolutionary approach that uses the next-generation sequencing technologies to detect and quantify expressed transcripts in biological samples. Compared to other methods such as microarrays, RNA-seq provides a more unbiased assessment of the full range of transcripts and their isoforms with a greater dynamic range in expression quantification.

In this tutorial, you will learn how to use the CRI's RNA-seq pipeline (available on both [CRI HPC cluster](http://cri.uchicago.edu/hpc/) and [GitHub](https://github.com/wenching/CRI-BIO-646-MED-YLi-CRIRnaSeqPipeline.git))) to analyze Illumina RNA sequencing data. The tutorial comprises the following Steps:

* Assess the sequencing data quality using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* Align the short reads to the target genome (e.g., Human) using [STAR](https://github.com/alexdobin/STAR)
* Measure alignment result using [Picard](https://broadinstitute.github.io/picard/) and [RSeQC](http://rseqc.sourceforge.net/)
* Quantify expression using [Subread](http://subread.sourceforge.net/)::[featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
* Identify differentially expressed genes using [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [limma](https://bioconductor.org/help/search/index.html?q=limma/)
* Heat Map using [pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html)
* Conduct functional enrichment analysis using [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)

By the end of this tutorial, you will:

* Understand the basics of RNA-seq experimental design
* Be familiar with the common data formats and standards in RNA-seq analysis
* Know computational tools for RNA-seq data processing
* Be able to align short sequence reads on a reference genome
* Be able to perform the basic analysis such as gene quantification and differential expression analysis

This tutorial is based on CRI's high-performance computing (HPC) cluster. If you are not familiar with this newly assembled cluster, a concise user's guide can be found [here](http://cri.uchicago.edu/wp-content/uploads/2017/04/Gardner-Part-1.pdf).



## <a name="Work"/>Work Flow | [Top](#Top)


The RNA-seq data used in this tutorial are from **CRI-BIO-646-BMB-RKeenan**.

In this tutorial, we use the sequencing reads in the project CRI-BIO-646 in mouse as example.
The sample information are saved in the file **<span style="color:red">`CRI-BIO-646.metadata.txt`</span>** (see [below](#Data)).

![Work Flow](IMG/RNAseq.workflow.png)



## <a name="Data"/>Data Description | [Top](#Top)

There are six (partial) single-end RNA-seq sequencing libraries will be used as the example dataset In this tutorial.
Their respective sample information is described in the metadata table CRI-BIO-646/**`CRI-BIO-646.metadata.txt`**.


<div style="border: 1px solid #ddd; padding: 5px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>Sample Description</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Sample </th>
   <th style="text-align:left;"> Library </th>
   <th style="text-align:left;"> ReadGroup </th>
   <th style="text-align:left;"> LibType </th>
   <th style="text-align:left;"> Platform </th>
   <th style="text-align:left;"> SequencingCenter </th>
   <th style="text-align:left;"> Date </th>
   <th style="text-align:right;"> Lane </th>
   <th style="text-align:left;"> Unit </th>
   <th style="text-align:left;"> Flavor </th>
   <th style="text-align:right;"> Encoding </th>
   <th style="text-align:right;"> Run </th>
   <th style="text-align:left;"> Genome </th>
   <th style="text-align:left;"> NucleicAcid </th>
   <th style="text-align:left;"> Group </th>
   <th style="text-align:left;"> Location </th>
   <th style="text-align:left;"> Seqfile1 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> WT01 </td>
   <td style="text-align:left;"> WT01 </td>
   <td style="text-align:left;"> BK-AA-1_S11_L007_R355 </td>
   <td style="text-align:left;"> RF </td>
   <td style="text-align:left;"> Illumina </td>
   <td style="text-align:left;"> GCF </td>
   <td style="text-align:left;"> 2018-11-07 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> K00242 </td>
   <td style="text-align:left;"> 1x50 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 355 </td>
   <td style="text-align:left;"> hg38 </td>
   <td style="text-align:left;"> rnaseq </td>
   <td style="text-align:left;"> WT </td>
   <td style="text-align:left;"> CRI-BIO-646/data/180216_K00242_0355_AHT7NTBBXX-RKeenan-AA-RS10/FastQ </td>
   <td style="text-align:left;"> BK-AA-1_S11_L007_R1_001.fastq.gz </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WT02 </td>
   <td style="text-align:left;"> WT02 </td>
   <td style="text-align:left;"> BK-AA-2_S12_L007_R355 </td>
   <td style="text-align:left;"> RF </td>
   <td style="text-align:left;"> Illumina </td>
   <td style="text-align:left;"> GCF </td>
   <td style="text-align:left;"> 2018-11-07 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> K00242 </td>
   <td style="text-align:left;"> 1x50 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 355 </td>
   <td style="text-align:left;"> hg38 </td>
   <td style="text-align:left;"> rnaseq </td>
   <td style="text-align:left;"> WT </td>
   <td style="text-align:left;"> CRI-BIO-646/data/180216_K00242_0355_AHT7NTBBXX-RKeenan-AA-RS10/FastQ </td>
   <td style="text-align:left;"> BK-AA-2_S12_L007_R1_001.fastq.gz </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WT03 </td>
   <td style="text-align:left;"> WT03 </td>
   <td style="text-align:left;"> BK-AA-3_S13_L007_R355 </td>
   <td style="text-align:left;"> RF </td>
   <td style="text-align:left;"> Illumina </td>
   <td style="text-align:left;"> GCF </td>
   <td style="text-align:left;"> 2018-11-07 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> K00242 </td>
   <td style="text-align:left;"> 1x50 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 355 </td>
   <td style="text-align:left;"> hg38 </td>
   <td style="text-align:left;"> rnaseq </td>
   <td style="text-align:left;"> WT </td>
   <td style="text-align:left;"> CRI-BIO-646/data/180216_K00242_0355_AHT7NTBBXX-RKeenan-AA-RS10/FastQ </td>
   <td style="text-align:left;"> BK-AA-3_S13_L007_R1_001.fastq.gz </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KO01 </td>
   <td style="text-align:left;"> KO01 </td>
   <td style="text-align:left;"> BK-AA-4_S14_L007_R355 </td>
   <td style="text-align:left;"> RF </td>
   <td style="text-align:left;"> Illumina </td>
   <td style="text-align:left;"> GCF </td>
   <td style="text-align:left;"> 2018-11-07 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> K00242 </td>
   <td style="text-align:left;"> 1x50 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 355 </td>
   <td style="text-align:left;"> hg38 </td>
   <td style="text-align:left;"> rnaseq </td>
   <td style="text-align:left;"> KO </td>
   <td style="text-align:left;"> CRI-BIO-646/data/180216_K00242_0355_AHT7NTBBXX-RKeenan-AA-RS10/FastQ </td>
   <td style="text-align:left;"> BK-AA-4_S14_L007_R1_001.fastq.gz </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KO02 </td>
   <td style="text-align:left;"> KO02 </td>
   <td style="text-align:left;"> BK-AA-5_S15_L007_R355 </td>
   <td style="text-align:left;"> RF </td>
   <td style="text-align:left;"> Illumina </td>
   <td style="text-align:left;"> GCF </td>
   <td style="text-align:left;"> 2018-11-07 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> K00242 </td>
   <td style="text-align:left;"> 1x50 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 355 </td>
   <td style="text-align:left;"> hg38 </td>
   <td style="text-align:left;"> rnaseq </td>
   <td style="text-align:left;"> KO </td>
   <td style="text-align:left;"> CRI-BIO-646/data/180216_K00242_0355_AHT7NTBBXX-RKeenan-AA-RS10/FastQ </td>
   <td style="text-align:left;"> BK-AA-5_S15_L007_R1_001.fastq.gz </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KO03 </td>
   <td style="text-align:left;"> KO03 </td>
   <td style="text-align:left;"> BK-AA-6_S16_L007_R355 </td>
   <td style="text-align:left;"> RF </td>
   <td style="text-align:left;"> Illumina </td>
   <td style="text-align:left;"> GCF </td>
   <td style="text-align:left;"> 2018-11-07 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> K00242 </td>
   <td style="text-align:left;"> 1x50 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 355 </td>
   <td style="text-align:left;"> hg38 </td>
   <td style="text-align:left;"> rnaseq </td>
   <td style="text-align:left;"> KO </td>
   <td style="text-align:left;"> CRI-BIO-646/data/180216_K00242_0355_AHT7NTBBXX-RKeenan-AA-RS10/FastQ </td>
   <td style="text-align:left;"> BK-AA-6_S16_L007_R1_001.fastq.gz </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WT01 </td>
   <td style="text-align:left;"> WT01 </td>
   <td style="text-align:left;"> BK-AA-1_S11_L007_R357 </td>
   <td style="text-align:left;"> RF </td>
   <td style="text-align:left;"> Illumina </td>
   <td style="text-align:left;"> GCF </td>
   <td style="text-align:left;"> 2018-11-07 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> K00242 </td>
   <td style="text-align:left;"> 1x50 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 357 </td>
   <td style="text-align:left;"> hg38 </td>
   <td style="text-align:left;"> rnaseq </td>
   <td style="text-align:left;"> WT </td>
   <td style="text-align:left;"> CRI-BIO-646/data/180223_K00242_0357_AHT37VBBXX-RKeenan-AA-RS10/FastQ </td>
   <td style="text-align:left;"> BK-AA-1_S11_L007_R1_001.fastq.gz </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WT02 </td>
   <td style="text-align:left;"> WT02 </td>
   <td style="text-align:left;"> BK-AA-2_S12_L007_R357 </td>
   <td style="text-align:left;"> RF </td>
   <td style="text-align:left;"> Illumina </td>
   <td style="text-align:left;"> GCF </td>
   <td style="text-align:left;"> 2018-11-07 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> K00242 </td>
   <td style="text-align:left;"> 1x50 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 357 </td>
   <td style="text-align:left;"> hg38 </td>
   <td style="text-align:left;"> rnaseq </td>
   <td style="text-align:left;"> WT </td>
   <td style="text-align:left;"> CRI-BIO-646/data/180223_K00242_0357_AHT37VBBXX-RKeenan-AA-RS10/FastQ </td>
   <td style="text-align:left;"> BK-AA-2_S12_L007_R1_001.fastq.gz </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WT03 </td>
   <td style="text-align:left;"> WT03 </td>
   <td style="text-align:left;"> BK-AA-3_S13_L007_R357 </td>
   <td style="text-align:left;"> RF </td>
   <td style="text-align:left;"> Illumina </td>
   <td style="text-align:left;"> GCF </td>
   <td style="text-align:left;"> 2018-11-07 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> K00242 </td>
   <td style="text-align:left;"> 1x50 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 357 </td>
   <td style="text-align:left;"> hg38 </td>
   <td style="text-align:left;"> rnaseq </td>
   <td style="text-align:left;"> WT </td>
   <td style="text-align:left;"> CRI-BIO-646/data/180223_K00242_0357_AHT37VBBXX-RKeenan-AA-RS10/FastQ </td>
   <td style="text-align:left;"> BK-AA-3_S13_L007_R1_001.fastq.gz </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KO01 </td>
   <td style="text-align:left;"> KO01 </td>
   <td style="text-align:left;"> BK-AA-4_S14_L007_R357 </td>
   <td style="text-align:left;"> RF </td>
   <td style="text-align:left;"> Illumina </td>
   <td style="text-align:left;"> GCF </td>
   <td style="text-align:left;"> 2018-11-07 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> K00242 </td>
   <td style="text-align:left;"> 1x50 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 357 </td>
   <td style="text-align:left;"> hg38 </td>
   <td style="text-align:left;"> rnaseq </td>
   <td style="text-align:left;"> KO </td>
   <td style="text-align:left;"> CRI-BIO-646/data/180223_K00242_0357_AHT37VBBXX-RKeenan-AA-RS10/FastQ </td>
   <td style="text-align:left;"> BK-AA-4_S14_L007_R1_001.fastq.gz </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KO02 </td>
   <td style="text-align:left;"> KO02 </td>
   <td style="text-align:left;"> BK-AA-5_S15_L007_R357 </td>
   <td style="text-align:left;"> RF </td>
   <td style="text-align:left;"> Illumina </td>
   <td style="text-align:left;"> GCF </td>
   <td style="text-align:left;"> 2018-11-07 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> K00242 </td>
   <td style="text-align:left;"> 1x50 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 357 </td>
   <td style="text-align:left;"> hg38 </td>
   <td style="text-align:left;"> rnaseq </td>
   <td style="text-align:left;"> KO </td>
   <td style="text-align:left;"> CRI-BIO-646/data/180223_K00242_0357_AHT37VBBXX-RKeenan-AA-RS10/FastQ </td>
   <td style="text-align:left;"> BK-AA-5_S15_L007_R1_001.fastq.gz </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KO03 </td>
   <td style="text-align:left;"> KO03 </td>
   <td style="text-align:left;"> BK-AA-6_S16_L007_R357 </td>
   <td style="text-align:left;"> RF </td>
   <td style="text-align:left;"> Illumina </td>
   <td style="text-align:left;"> GCF </td>
   <td style="text-align:left;"> 2018-11-07 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> K00242 </td>
   <td style="text-align:left;"> 1x50 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 357 </td>
   <td style="text-align:left;"> hg38 </td>
   <td style="text-align:left;"> rnaseq </td>
   <td style="text-align:left;"> KO </td>
   <td style="text-align:left;"> CRI-BIO-646/data/180223_K00242_0357_AHT37VBBXX-RKeenan-AA-RS10/FastQ </td>
   <td style="text-align:left;"> BK-AA-6_S16_L007_R1_001.fastq.gz </td>
  </tr>
</tbody>
</table></div>



## <a name="Pre"/>Prerequisites | [Top](#Top)


We will use SSH (Secure Shell) to connect to CRI's HPC. SSH now is included or can be installed in all standard operating systems (Windows, Linux, and OS X).



### <a name="Login"/>Login and Setup Tutorial Working Directory | [Top](#Top)


The login procedure varies slightly depending on whether you use a Mac/Unix/Linux computer or a Windows computer.

* Log into one of entry nodes in CRI HPC  
    1. Open a terminal session.  
    2. Connect to the login node of the CRI HPC cluster:  
        
        ```bash
        $ ssh -l username@gardner.cri.uchicago.edu
        ```
    3. If it's your first time to log in, you will be asked to accept the ssh key. Type "**yes**"  
    4. Type in the password when prompted  

        > Make sure that you replace **_<span style="color:red">`username`</span>_** with your login name.

>  
> - **<span style="color:red">CAUTION</span>**
>     - THIS PACKAGE IS LARGE, PLEASE DO NOT DOWNLOAD IT TO YOUR HOME DIRECTORY
>     - USE OTHER LOCATION LIKE /gpfs/data/bioinformatics/username
>  

* Set up a tutorial directory  
    1. You should be in your home directory after logging in  
        
        ```bash
        $ pwd
        /home/username
        ```
    2. Instead of downloading the pipeline package to your local home directory, use other location like /gpfs/data/bioinformatics/username
        
        ```bash
        $ cd /gpfs/data/bioinformatics/username; pwd
        /gpfs/data/bioinformatics/username
        ```

* Download the pipeline package  
    1. One way to download the pipeline package via `wget`  
        
        ```bash
        $ wget ftp://logia.cri.uchicago.edu/bioinformatics/tutorials/Nov2018/CRI-BIO-646-BMB-RKeenan.tgz
        ```
    2. Or, copy the package directly from HPC  
        
        ```bash
        $ cp /gpfs/data/bioinformatics/shared/tutorials/Nov2018/CRI-BIO-646-BMB-RKeenan.tgz .
        ```
    3. Uncompress the tarball file  
        
        ```bash
        $ tar -zxvf CRI-BIO-646-BMB-RKeenan.tgz
        ```
    4. Change working directory to pipeline dirctory  
        
        ```bash
        $ cd CRI-BIO-646-BMB-RKeenan
        $ tree -d -L 4
        ```
        
        ```
        ## ../
        ## ├── CRI-BIO-646
        ## │   ├── data
        ## │   │   ├── 180216_K00242_0355_AHT7NTBBXX-RKeenan-AA-RS10
        ## │   │   │   └── FastQ
        ## │   │   └── 180223_K00242_0357_AHT37VBBXX-RKeenan-AA-RS10
        ## │   │       └── FastQ
        ## │   └── references
        ## │       └── v28_92_GRCh38.p12
        ## │           └── STAR
        ## ├── DATA
        ## ├── SRC
        ## │   ├── Python
        ## │   │   ├── lib
        ## │   │   ├── module
        ## │   │   └── util
        ## │   └── R
        ## │       ├── lib
        ## │       │   └── 3.5
        ## │       ├── module
        ## │       └── util
        ## └── docs
        ##     ├── IMG
        ##     └── result
        ## 
        ## 23 directories
        ```

* File structure  
    - Raw sequencing data files (*.fastq.gz) are located at **`CRI-BIO-646/data/`**  
        
        ```bash
        $ tree CRI-BIO-646/data/
        |-- 180216_K00242_0355_AHT7NTBBXX-RKeenan-AA-RS10
        |   `-- FastQ
        |       |-- BK-AA-1_S11_L007_R1_001.fastq.gz
        |       |-- BK-AA-2_S12_L007_R1_001.fastq.gz
        |       |-- BK-AA-3_S13_L007_R1_001.fastq.gz
        |       |-- BK-AA-4_S14_L007_R1_001.fastq.gz
        |       |-- BK-AA-5_S15_L007_R1_001.fastq.gz
        |       `-- BK-AA-6_S16_L007_R1_001.fastq.gz
        `-- 180223_K00242_0357_AHT37VBBXX-RKeenan-AA-RS10
            `-- FastQ
                |-- BK-AA-1_S11_L007_R1_001.fastq.gz
                |-- BK-AA-2_S12_L007_R1_001.fastq.gz
                |-- BK-AA-3_S13_L007_R1_001.fastq.gz
                |-- BK-AA-4_S14_L007_R1_001.fastq.gz
                |-- BK-AA-5_S15_L007_R1_001.fastq.gz
                `-- BK-AA-6_S16_L007_R1_001.fastq.gz
        ```
    - Genome data are located at **`/gpfs/data/bioinformatics/ReferenceData/cri_rnaseq_2018/vM18_93_GRCm38.p6`**  
        
        ```bash
        $ tree /gpfs/data/bioinformatics/cri_rnaseq_2018_ex/example/references/v28_92_GRCh38.p12
        |-- GRCh38_rRNA.bed
        |-- GRCh38_rRNA.bed.interval_list
        |-- STAR
        |   |-- Genome
        |   |-- Log.out
        |   |-- SA
        |   |-- SAindex
        |   |-- chrLength.txt
        |   |-- chrName.txt
        |   |-- chrNameLength.txt
        |   |-- chrStart.txt
        |   |-- exonGeTrInfo.tab
        |   |-- exonInfo.tab
        |   |-- geneInfo.tab
        |   |-- genomeParameters.txt
        |   |-- run_genome_generate.logs
        |   |-- sjdbInfo.txt
        |   |-- sjdbList.fromGTF.out.tab
        |   |-- sjdbList.out.tab
        |   `-- transcriptInfo.tab
        |-- genes.gtf
        |-- genes.gtf.bed12
        |-- genes.refFlat.txt
        |-- genome.chrom.sizes
        |-- genome.dict
        `-- genome.fa
        ```

* Pipeline/project related files  
    - project related files (i.e., metadata & configuration file) as used in this tutorial are located under **`CRI-BIO-646/`**  
        
        ```bash
        $ ls -l CRI-BIO-646/CRI-BIO-646.*
        ```
        
        ```
        ## CRI-BIO-646/CRI-BIO-646.metadata.txt
        ## CRI-BIO-646/CRI-BIO-646.pipeline.yaml
        ```
        + Here are the first few lines in the configuration example file CRI-BIO-646/**`CRI-BIO-646.pipeline.yaml`**
            
            ```
            ---
            pipeline:
              flags:
                aligners:
                  run_star: True
                quantifiers:
                  run_featurecounts: True
                  run_rsem: False
                  run_kallisto: False
                callers:
                  run_edger: True
                  run_deseq2: True
                  run_limma: True
              software:
                main:
                  use_module: 0
                  adapter_pe: AGATCGGAAGAGCGGTTCAG,AGATCGGAAGAGCGTCGTGT
                  adapter_se: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
                  fastq_format: 33
                  genome_assembly: hg38
            ```
        
        >   
        > When running on another dataset, you will need to modify these <span style="color:red">`two files`</span> and the master pipeline script (i.e., <span style="color:red">`Build_RNAseq.CRI-BIO-646.sh`</span>) (as described below) accordingly.  
        >   
        > For instance, if you would like to turn off the DE analysis tool limma, you can set the respecitve paramter to 'False' in configuration file  
        >     - run_limma: False  
        >   
        
        >  
        > For metadata file, you might pay attendtion on the following settings
        > 1. **Single End (SE)** Library
        >     1. Set Flavor column as 1xReadLength (e.g., 1x50)
        >     2. Set Seqfile1 column as the file name of the repective sequencing file
        > 2. **Paired End (PE)** Library
        >     1. Set Flavor column as 2xReadLength (e.g., 2x50)
        >     2. Set Seqfile1 column as the file name of the repective read 1 (R1) sequencing file
        >     3. Set an additional column named '**Seqfile2**' as the file name of the repective read 2 (R2) sequencing file
        > 3. **Non strand-specific** Library
        >     1. Set LibType column to NS
        > 4. **Strand-specific** Library
        >     1. Inquire the library type from your seuqencing center and set LibType column to FR (the left-most end of the fragment (in transcript coordinates, or the first-strand synthesis) is the first sequenced) or RF (the right-most end of the fragment (in transcript coordinates) is the first sequenced, or the second-strand synthesis). You can read this [blog](http://onetipperday.sterding.com/2012/07/how-to-tell-which-library-type-to-use.html) for more details of strand-specific RNA-seq.
        >  
    
    - Master pipeline script  
        
        ```bash
        $ cat Build_RNAseq.CRI-BIO-646.sh
        ```
        
        ```
        ## 
        ## 
        ## ## build pipeline scripts
        ## 
        ## now=$(date +"%m-%d-%Y_%H:%M:%S")
        ## 
        ## ## project info
        ## project="CRI-BIO-646"
        ## SubmitRNAseqExe="Submit_${PWD##*/}.sh"
        ## padding="CRI-BIO-646/"
        ## 
        ## ## command
        ## echo "START" `date` " Running build_rnaseq.py"
        ## python3 SRC/Python/build_rnaseq.py \
        ## 	--projdir $PWD \
        ## 	--metadata $PWD/${padding}$project.metadata.txt \
        ## 	--config $PWD/${padding}$project.pipeline.yaml \
        ## 	--systype cluster \
        ## 	--threads 8 \
        ## 	--log_file $PWD/Build_RNAseq.$project.$now.log
        ## 
        ## ## submit pipeline master script
        ## echo "START" `date` " Running $SubmitRNAseqExe"
        ## echo "bash $SubmitRNAseqExe"
        ## 
        ## echo "END" `date`
        ```
        
        >   
        > Basically, when running on your own dataset, you will need to modify this master pipeline script (i.e., <span style="color:red">`Build_RNAseq.CRI-BIO-646.sh`</span>) accordingly.    
        >   
        > For instance, you can change respective parameters as follows.  
        >     - project="PROJECT_AS_PREFIX" (e.g., **CRI-BIO-646** which is used as a prefix of metadata file **CRI-BIO-646**.metadata.txt and configuration file **CRI-BIO-646**.pipeline.yaml)  
        >     - padding="DIRECTORY_NAME_CONTAINING_PROJECT_DATA" (e.g., **CRI-BIO-646** which is the folder name to accommodate metadata file, configuration file, sequencing data folder, and references folder)  
        >   


## <a name="Steps"/>Pipeline Steps | [Top](#Top)

* Before running, you need to know

    >
    > The master BigDataScript script can be <span style="color:red">`ONLY run on a head/entry node`</span> other than any other compuatation node.
    > But, you can step by step run individual sub-task bash scripts in any computation node interactively.
    >

* Generate sub-task scirpts of the pipepline
    
    ```bash
    # load modules
    $ module purge;module load gcc udunits R/3.5.0 python/3.6.0; module update
    
    # This step is optional but it will install all necessary R packages ahead.
    # In case the pipeline was terminated due to the failure of R package installation later when running the pipeline.
    # A successful execution result will show the package versions of all necessary packages from devtools to pheatmap
    $ Rscript --vanilla SRC/R/util/prerequisite.packages.R
    
    # create directories and generate all necessary scripts
    $ bash Build_RNAseq.CRI-BIO-646.sh
    ```

    + this step will execute **`SRC/Python/build_rnaseq.py`** using python3 to generate all sub-task bash scripts and directories according to the provided metadata and configuration files (i.e., CRI-BIO-646/**`CRI-BIO-646.metadata.txt`** and CRI-BIO-646/**`CRI-BIO-646.pipeline.yaml`** )
        
        ```bash
        $ tree -d RNAseq
        |-- Aln
        |   `-- star
        |       |-- KO01
        |       |   |-- BK-AA-4_S14_L007_R355
        |       |   `-- BK-AA-4_S14_L007_R357
        |       |-- KO02
        |       |   |-- BK-AA-5_S15_L007_R355
        |       |   `-- BK-AA-5_S15_L007_R357
        |       |-- KO03
        |       |   |-- BK-AA-6_S16_L007_R355
        |       |   `-- BK-AA-6_S16_L007_R357
        |       |-- WT01
        |       |   |-- BK-AA-1_S11_L007_R355
        |       |   `-- BK-AA-1_S11_L007_R357
        |       |-- WT02
        |       |   |-- BK-AA-2_S12_L007_R355
        |       |   `-- BK-AA-2_S12_L007_R357
        |       `-- WT03
        |           |-- BK-AA-3_S13_L007_R355
        |           `-- BK-AA-3_S13_L007_R357
        |-- AlnQC
        |   |-- picard
        |   |   `-- star
        |   |       |-- KO01
        |   |       |   `-- tmp
        |   |       |-- KO02
        |   |       |   `-- tmp
        |   |       |-- KO03
        |   |       |   `-- tmp
        |   |       |-- WT01
        |   |       |   `-- tmp
        |   |       |-- WT02
        |   |       |   `-- tmp
        |   |       `-- WT03
        |   |           `-- tmp
        |   `-- rseqc
        |       `-- star
        |           |-- KO01
        |           |   `-- tmp
        |           |-- KO02
        |           |   `-- tmp
        |           |-- KO03
        |           |   `-- tmp
        |           |-- WT01
        |           |   `-- tmp
        |           |-- WT02
        |           |   `-- tmp
        |           `-- WT03
        |               `-- tmp
        |-- DEG
        |   |-- deseq2
        |   |   `-- featurecounts
        |   |       `-- star
        |   |-- edger
        |   |   `-- featurecounts
        |   |       `-- star
        |   `-- limma
        |       `-- featurecounts
        |           `-- star
        |-- LociStat
        |   `-- featurecounts
        |       `-- star
        |           `-- CRI-BIO-646-BMB-RKeenan
        |-- PostAna
        |   |-- clusterprofiler
        |   |   `-- featurecounts
        |   |       `-- star
        |   |           `-- CRI-BIO-646-BMB-RKeenan
        |   `-- pheatmap
        |       `-- featurecounts
        |           `-- star
        |               `-- CRI-BIO-646-BMB-RKeenan
        |-- QuantQC
        |   `-- featurecounts
        |       `-- star
        |-- Quantification
        |   `-- featurecounts
        |       `-- star
        |           |-- KO01
        |           |-- KO02
        |           |-- KO03
        |           |-- WT01
        |           |-- WT02
        |           `-- WT03
        |-- RawReadQC
        |   |-- KO01
        |   |   |-- BK-AA-4_S14_L007_R355
        |   |   |   `-- BK-AA-4_S14_L007_R1_001_fastqc
        |   |   |       |-- Icons
        |   |   |       `-- Images
        |   |   `-- BK-AA-4_S14_L007_R357
        |   |       `-- BK-AA-4_S14_L007_R1_001_fastqc
        |   |           |-- Icons
        |   |           `-- Images
        |   |-- KO02
        |   |   |-- BK-AA-5_S15_L007_R355
        |   |   |   `-- BK-AA-5_S15_L007_R1_001_fastqc
        |   |   |       |-- Icons
        |   |   |       `-- Images
        |   |   `-- BK-AA-5_S15_L007_R357
        |   |       `-- BK-AA-5_S15_L007_R1_001_fastqc
        |   |           |-- Icons
        |   |           `-- Images
        |   |-- KO03
        |   |   |-- BK-AA-6_S16_L007_R355
        |   |   |   `-- BK-AA-6_S16_L007_R1_001_fastqc
        |   |   |       |-- Icons
        |   |   |       `-- Images
        |   |   `-- BK-AA-6_S16_L007_R357
        |   |       `-- BK-AA-6_S16_L007_R1_001_fastqc
        |   |           |-- Icons
        |   |           `-- Images
        |   |-- WT01
        |   |   |-- BK-AA-1_S11_L007_R355
        |   |   |   `-- BK-AA-1_S11_L007_R1_001_fastqc
        |   |   |       |-- Icons
        |   |   |       `-- Images
        |   |   `-- BK-AA-1_S11_L007_R357
        |   |       `-- BK-AA-1_S11_L007_R1_001_fastqc
        |   |           |-- Icons
        |   |           `-- Images
        |   |-- WT02
        |   |   |-- BK-AA-2_S12_L007_R355
        |   |   |   `-- BK-AA-2_S12_L007_R1_001_fastqc
        |   |   |       |-- Icons
        |   |   |       `-- Images
        |   |   `-- BK-AA-2_S12_L007_R357
        |   |       `-- BK-AA-2_S12_L007_R1_001_fastqc
        |   |           |-- Icons
        |   |           `-- Images
        |   `-- WT03
        |       |-- BK-AA-3_S13_L007_R355
        |       |   `-- BK-AA-3_S13_L007_R1_001_fastqc
        |       |       |-- Icons
        |       |       `-- Images
        |       `-- BK-AA-3_S13_L007_R357
        |           `-- BK-AA-3_S13_L007_R1_001_fastqc
        |               |-- Icons
        |               `-- Images
        `-- shell_scripts
        ```
* Execute the entire analysis with just one command
    
    ```bash
    # This will start to run the entire pipeline.
    # You can chekc teh BDS report to know the running status.
    $ bash Submit_CRI-BIO-646-BMB-RKeenan.sh
    ```
    + Again, this step will execute the master BigDataScript script **`Submit_CRI-BIO-646-BMB-RKeenan.bds`**, so you will need to run this command <span style="color:red">**on a head/entry node**</span>.
    
    > 
    > Check running status
    > 1. $ qstat -a
    > 2. tail *.bds.log
    > 
  



## <a name="Navi"/>Navigation Map | [Top](#Top)


![Navigation Map](IMG/CRI-RnaSeq-Workflow.png)



### <a name="ReadQC"/>Step 1: Quality Control | [Top](#Top)


For the first step, the pipeline will perform quality assessment on the raw fastq files.

The BDS code snippet for the sample KO01 will look like:


```bash
$ grep -A1 run.RawReadQC.FastQC.BK-AA-4_S14_L007_R355.sh Submit_*.bds
```

```
## dep( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/RawReadQC/KO01/BK-AA-4_S14_L007_R355/BK-AA-4_S14_L007_R1_001_fastqc.zip' ] <- [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/CRI-BIO-646/data/180216_K00242_0355_AHT7NTBBXX-RKeenan-AA-RS10/FastQ/BK-AA-4_S14_L007_R1_001.fastq.gz' ], cpus := 1, mem := 16*G, timeout := 72*hour, taskName := "FastQC.BK-AA-4_S14_L007_R355") sys bash /gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/shell_scripts/run.RawReadQC.FastQC.BK-AA-4_S14_L007_R355.sh; sleep 2
## goal( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/RawReadQC/KO01/BK-AA-4_S14_L007_R355/BK-AA-4_S14_L007_R1_001_fastqc.zip' ] )
```

This code chunk will invoke the bash script RNAseq/shell_scripts/**`run.RawReadQC.FastQC.BK-AA-4_S14_L007_R355.sh`** to execute FastQC on the KO01(BK-AA-1_S11_L007_R355) sequencing library.

After the completion of the entire pipeline, you can check FastQC report per individual libraries; for instance, the partial report of KO01 will be as follows or a full [report](result/BK-AA-4_S14_L007_R355_fastqc.html).

<img src="IMG/BK-AA-4_S14_L007_R1_001_fastqc.png" alt="KO01_FastQC"  width="800" height="1200">

You can check [FastQC Help](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/) for more details about how to interpret a FastQC report.

Or, compare your reports to the example reports provided by FastQC for a [Good Illumina Data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) or [Bad Illumina Data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html).



### <a name="Aln"/>Step 2.1: Read Alignment | [Top](#Top)


In this step, the pipeline will conduct read alignment on the raw fastq files.

The BDS code snippet for the sample KO01 will look like:


```bash
$ grep -A1 run.alignRead.star.BK-AA-4_S14_L007_R355.sh Submit_*.bds
```

```
## dep( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Aln/star/KO01/BK-AA-4_S14_L007_R355/BK-AA-4_S14_L007_R355.star.bam' ] <- [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/CRI-BIO-646/data/180216_K00242_0355_AHT7NTBBXX-RKeenan-AA-RS10/FastQ/BK-AA-4_S14_L007_R1_001.fastq.gz' ], cpus := 4, mem := 64*G, timeout := 72*hour, taskName := "star.BK-AA-4_S14_L007_R355") sys bash /gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/shell_scripts/run.alignRead.star.BK-AA-4_S14_L007_R355.sh; sleep 2
## goal( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Aln/star/KO01/BK-AA-4_S14_L007_R355/BK-AA-4_S14_L007_R355.star.bam' ] )
```

This code chunk will invoke the bash script RNAseq/shell_scripts/**`run.alignRead.star.BK-AA-4_S14_L007_R355.sh`** to execute STAR on the KO01(BK-AA-4_S14_L007_R355) sequencing library.

After the completion of the entire pipeline, you can check the alignment result of each individual libraries; for instance, the result of KO01(BK-AA-4_S14_L007_R355) will be as follows.


```bash
$ tree RNAseq/Aln/star/KO01/BK-AA-4_S14_L007_R355
RNAseq/Aln/star/KO01/BK-AA-4_S14_L007_R355
|-- BK-AA-4_S14_L007_R355.star.Aligned.sortedByCoord.out.bam
|-- BK-AA-4_S14_L007_R355.star.Log.final.out
|-- BK-AA-4_S14_L007_R355.star.Log.out
|-- BK-AA-4_S14_L007_R355.star.Log.progress.out
|-- BK-AA-4_S14_L007_R355.star.SJ.out.tab
|-- BK-AA-4_S14_L007_R355.star.Unmapped.out.mate1
|-- BK-AA-4_S14_L007_R355.star.bai
|-- BK-AA-4_S14_L007_R355.star.bam -> BK-AA-4_S14_L007_R355.star.Aligned.sortedByCoord.out.bam
`-- run.alignRead.star.BK-AA-4_S14_L007_R355.log
```

You can check a log file (e.g., **RNAseq/Aln/star/KO01/BK-AA-4_S14_L007_R355/`BK-AA-4_S14_L007_R355.star.Log.final.out`**) for more alignment information provided by STAR.


```bash
$ cat RNAseq/Aln/star/KO01/BK-AA-4_S14_L007_R355/BK-AA-4_S14_L007_R355.star.Log.final.out
```

```
##                                  Started job on |	Nov 07 12:19:19
##                              Started mapping on |	Nov 07 12:19:34
##                                     Finished on |	Nov 07 12:21:59
##        Mapping speed, Million of reads per hour |	430.40
## 
##                           Number of input reads |	17335560
##                       Average input read length |	50
##                                     UNIQUE READS:
##                    Uniquely mapped reads number |	14809963
##                         Uniquely mapped reads % |	85.43%
##                           Average mapped length |	49.84
##                        Number of splices: Total |	2272172
##             Number of splices: Annotated (sjdb) |	2259453
##                        Number of splices: GT/AG |	2249028
##                        Number of splices: GC/AG |	19534
##                        Number of splices: AT/AC |	2520
##                Number of splices: Non-canonical |	1090
##                       Mismatch rate per base, % |	0.20%
##                          Deletion rate per base |	0.00%
##                         Deletion average length |	1.44
##                         Insertion rate per base |	0.00%
##                        Insertion average length |	1.28
##                              MULTI-MAPPING READS:
##         Number of reads mapped to multiple loci |	0
##              % of reads mapped to multiple loci |	0.00%
##         Number of reads mapped to too many loci |	2403663
##              % of reads mapped to too many loci |	13.87%
##                                   UNMAPPED READS:
##        % of reads unmapped: too many mismatches |	0.00%
##                  % of reads unmapped: too short |	0.56%
##                      % of reads unmapped: other |	0.15%
##                                   CHIMERIC READS:
##                        Number of chimeric reads |	0
##                             % of chimeric reads |	0.00%
```



### <a name="AlnQC"/>Step 2.2: Alignment QC | [Top](#Top)


In this step, the pipeline will conduct a QC on alignment result.

The BDS code snippets for the sample KO01 will look like:


```bash
$ grep -A1 run.alnQC.*.star.KO01.*.sh Submit_*.bds
```

```
## dep( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/picard/star/KO01/KO01.star.picard.RNA_Metrics', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/picard/star/KO01/KO01.star.picard.RNA_Metrics.pdf' ] <- [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Aln/star/KO01/KO01.star.bai' ], cpus := 4, mem := 32*G, timeout := 72*hour, taskName := "picard.star.KO01") sys bash /gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/shell_scripts/run.alnQC.picard.star.KO01.CollectRnaSeqMetrics.sh; sleep 2
## goal( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/picard/star/KO01/KO01.star.picard.RNA_Metrics', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/picard/star/KO01/KO01.star.picard.RNA_Metrics.pdf' ] )
## --
## dep( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/rseqc/star/KO01/KO01.star.rseqc.clipping_profile.xls', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/rseqc/star/KO01/KO01.star.rseqc.clipping_profile.r', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/rseqc/star/KO01/KO01.star.rseqc.clipping_profile.pdf' ] <- [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Aln/star/KO01/KO01.star.bai' ], cpus := 8, mem := 64*G, timeout := 72*hour, taskName := "rseqc.star.KO01") sys bash /gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/shell_scripts/run.alnQC.rseqc.star.KO01.clipping_profile.py.sh; sleep 2
## goal( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/rseqc/star/KO01/KO01.star.rseqc.clipping_profile.xls', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/rseqc/star/KO01/KO01.star.rseqc.clipping_profile.r', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/rseqc/star/KO01/KO01.star.rseqc.clipping_profile.pdf' ] )
## --
## dep( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/rseqc/star/KO01/KO01.star.rseqc.infer_experiment.txt' ] <- [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Aln/star/KO01/KO01.star.bai' ], cpus := 8, mem := 64*G, timeout := 72*hour, taskName := "rseqc.star.KO01") sys bash /gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/shell_scripts/run.alnQC.rseqc.star.KO01.infer_experiment.py.sh; sleep 2
## goal( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/rseqc/star/KO01/KO01.star.rseqc.infer_experiment.txt' ] )
```

This code chunk will invoke few bash scripts (e.g., RNAseq/shell_scripts/**`run.alnQC.picard.star.KO01.CollectRnaSeqMetrics.sh`**, **`run.alnQC.rseqc.star.KO01.clipping_profile.py.sh`**, and **`run.alnQC.rseqc.star.KO01.infer_experiment.py.sh`**) to execute alignment QC tools (i.e., [Picard](https://broadinstitute.github.io/picard/) and [RSeQC](http://rseqc.sourceforge.net/)) on the sample KO01.

After the completion of the entire pipeline, you can check the alignment QC results of each individual samples; for instance, the results of KO01 will be as follows.


```bash
$ tree RNAseq/AlnQC/*/star/KO01
RNAseq/AlnQC/picard/star/KO01
|-- KO01.star.picard.RNA_Metrics
|-- KO01.star.picard.RNA_Metrics.pdf
|-- run.alnQC.picard.star.KO01.CollectRnaSeqMetrics.log
`-- tmp
RNAseq/AlnQC/rseqc/star/KO01
|-- KO01.star.rseqc.clipping_profile.pdf
|-- KO01.star.rseqc.clipping_profile.r
|-- KO01.star.rseqc.clipping_profile.xls
|-- KO01.star.rseqc.infer_experiment.txt
|-- run.alnQC.rseqc.star.KO01.clipping_profile.py.log
`-- tmp
```

You can check alignment statistics (e.g., **RNAseq/AlnQC/picard/star/KO01/`KO01.star.picard.RNA_Metrics`**) for more information provided by [Picard](https://broadinstitute.github.io/picard/).


```bash
$ head RNAseq/AlnQC/picard/star/KO01/KO01.star.picard.RNA_Metrics
```

```
## ## htsjdk.samtools.metrics.StringHeader
## # picard.analysis.CollectRnaSeqMetrics REF_FLAT=/gpfs/data/bioinformatics/cri_rnaseq_2018/example/references/v28_92_GRCh38.p12/genes.refFlat.txt RIBOSOMAL_INTERVALS=/gpfs/data/bioinformatics/cri_rnaseq_2018/example/references/v28_92_GRCh38.p12/GRCh38_rRNA.bed.interval_list STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND CHART_OUTPUT=/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/picard/star/KO01/KO01.star.picard.RNA_Metrics.pdf METRIC_ACCUMULATION_LEVEL=[SAMPLE, ALL_READS] INPUT=/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Aln/star/KO01/KO01.star.bam OUTPUT=/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/picard/star/KO01/KO01.star.picard.RNA_Metrics TMP_DIR=[/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/AlnQC/picard/star/KO01/tmp]    MINIMUM_LENGTH=500 RRNA_FRAGMENT_PERCENTAGE=0.8 ASSUME_SORTED=true STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json
## ## htsjdk.samtools.metrics.StringHeader
## # Started on: Wed Nov 07 12:32:15 CST 2018
## 
## ## METRICS CLASS	picard.analysis.RnaSeqMetrics
## PF_BASES	PF_ALIGNED_BASES	RIBOSOMAL_BASES	CODING_BASES	UTR_BASES	INTRONIC_BASES	INTERGENIC_BASES	IGNORED_READS	CORRECT_STRAND_READS	INCORRECT_STRAND_READS	PCT_RIBOSOMAL_BASES	PCT_CODING_BASES	PCT_UTR_BASES	PCT_INTRONIC_BASES	PCT_INTERGENIC_BASES	PCT_MRNA_BASES	PCT_USABLE_BASES	PCT_CORRECT_STRAND_READS	MEDIAN_CV_COVERAGE	MEDIAN_5PRIME_BIAS	MEDIAN_3PRIME_BIAS	MEDIAN_5PRIME_TO_3PRIME_BIAS	SAMPLE	LIBRARY	READ_GROUP
## 1582833150	1577842710	334200	876931136	608664079	71989244	19924828	0	4418837	102189	0.000212	0.555779	0.385757	0.045625	0.012628	0.941536	0.938567	0.977397	0.904811	0.026924	0.99656	0.041207			
## 1582833150	1577842710	334200	876931136	608664079	71989244	19924828	0	4418837	102189	0.000212	0.555779	0.385757	0.045625	0.012628	0.941536	0.938567	0.977397	0.904811	0.026924	0.99656	0.041207	unknown		
```

Or, the respective coverage plot of the sample KO01 produced by [Picard](https://broadinstitute.github.io/picard/) will be as follows.

<img src="IMG/KO01.star.picard.RNA_Metrics.png" alt="KO01_coverage"  width="600" height="600">

There are two alignment measurements performed using [RSeQC](http://rseqc.sourceforge.net).

1. [clipping_profile.py](http://rseqc.sourceforge.net/#clipping-profile-py)
    * Calculate the distributions of clipped nucleotides across reads
2. [infer_experiment.py](http://rseqc.sourceforge.net/#infer-experiment-py)
    * Use to “guess” how RNA-seq sequencing was configured, particularly how reads were stranded for strand-specific RNA-seq data, through comparing the “strandedness of reads” with the “strandedness of transcripts”.
<!-- 3. [RPKM_saturation.py](http://rseqc.sourceforge.net/#rpkm-saturation-py) -->
<!--     * Calculate RPKM value using a series of resampled subsets from total RNA reads to check if the current sequencing depth was saturated or not (or if the RPKM values were stable or not) in terms of genes’ expression estimation -->

The results will be as follows. Please check the [RSeQC](http://rseqc.sourceforge.net) website for more measurements and details.

<img src="IMG/KO01.star.rseqc.clipping_profile.png" alt="KO01_clipping_profile"  width="600" height="600">


```bash
$cat RNAseq/AlnQC/rseqc/star/KO01/KO01.star.rseqc.infer_experiment.txt
```

```
## 
## 
## This is SingleEnd Data
## Fraction of reads failed to determine: 0.0757
## Fraction of reads explained by "++,--": 0.0179
## Fraction of reads explained by "+-,-+": 0.9064
```

Here, you can confirm again the sample KO01 is a single-end, strand-specific library using the second-strand synthesis.



### <a name="Quant"/>Step 3: Expression Quantification | [Top](#Top)

In this step, the pipeline will conduct expression quantification over alignments.

The BDS code snippet for the sample KO01 will look like:


```bash
$ grep -A1 run.quant.featurecounts.star.KO01.sh Submit_*.bds
```

```
## dep( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/KO01/KO01.star.featurecounts.count' ] <- [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Aln/star/KO01/KO01.star.bai' ], cpus := 4, mem := 32*G, timeout := 72*hour, taskName := "featurecounts.star.KO01") sys bash /gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/shell_scripts/run.quant.featurecounts.star.KO01.sh; sleep 2
## goal( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/KO01/KO01.star.featurecounts.count' ] )
```

This code chunk will invoke the bash script (e.g., RNAseq/shell_scripts/**`run.quant.featurecounts.star.KO01.sh`**) to execute expression quantification tool (i.e., [Subread](http://subread.sourceforge.net/)::[featureCounts](http://bioinf.wehi.edu.au/featureCounts/) on the sample KO01.

After the completion of the entire pipeline, you can check the quantification results of each individual samples; for instance, the results of KO01 will be as follows.


```bash
$ tree RNAseq/Quantification/featurecounts/star/KO01
RNAseq/Quantification/featurecounts/star/KO01
|-- KO01.star.featurecounts.count
|-- KO01.star.featurecounts.count.jcounts
|-- KO01.star.featurecounts.count.summary
`-- run.quant.featurecounts.star.KO01.log
```

You can check quantification statistics (e.g., **RNAseq/Quantification/featurecounts/star/KO01`KO01.star.featurecounts.count.summary`**) for more information provided by [Subread](http://subread.sourceforge.net/)::[featureCounts](http://bioinf.wehi.edu.au/featureCounts/)


```bash
$ cat RNAseq/Quantification/featurecounts/star/KO01/KO01.star.featurecounts.count.summary
```

```
## Status	/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Aln/star/KO01/KO01.star.bam
## Assigned	28638683
## Unassigned_Unmapped	0
## Unassigned_MappingQuality	0
## Unassigned_Chimera	0
## Unassigned_FragmentLength	0
## Unassigned_Duplicate	0
## Unassigned_MultiMapping	0
## Unassigned_Secondary	0
## Unassigned_Nonjunction	0
## Unassigned_NoFeatures	2036590
## Unassigned_Overlapping_Length	0
## Unassigned_Ambiguity	981390
```

Or, the top 10 most abundant genes in the sample KO01 will be as follows.


```bash
$ cat <(head -n2 RNAseq/Quantification/featurecounts/star/KO01/KO01.star.featurecounts.count | tail -n+2 | cut -f1,7) <(cut -f1,7 RNAseq/Quantification/featurecounts/star/KO01/KO01.star.featurecounts.count | sort -k2,2nr | head)
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>Top 10 most abundant genes in KO01</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Geneid </th>
   <th style="text-align:left;"> Chr </th>
   <th style="text-align:left;"> Start </th>
   <th style="text-align:left;"> End </th>
   <th style="text-align:left;"> Strand </th>
   <th style="text-align:left;"> Length </th>
   <th style="text-align:left;"> KO01 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSG00000210082.2 </td>
   <td style="text-align:left;"> chrM </td>
   <td style="text-align:left;"> 1671 </td>
   <td style="text-align:left;"> 3229 </td>
   <td style="text-align:left;"> + </td>
   <td style="text-align:left;"> 1559 </td>
   <td style="text-align:left;"> 129707 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000198804.2 </td>
   <td style="text-align:left;"> chrM </td>
   <td style="text-align:left;"> 5904 </td>
   <td style="text-align:left;"> 7445 </td>
   <td style="text-align:left;"> + </td>
   <td style="text-align:left;"> 1542 </td>
   <td style="text-align:left;"> 129552 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000167658.15 </td>
   <td style="text-align:left;"> chr19 </td>
   <td style="text-align:left;"> 3976056 </td>
   <td style="text-align:left;"> 3985469 </td>
   <td style="text-align:left;"> - </td>
   <td style="text-align:left;"> 4027 </td>
   <td style="text-align:left;"> 105680 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000229807.10 </td>
   <td style="text-align:left;"> chrX </td>
   <td style="text-align:left;"> 73820651 </td>
   <td style="text-align:left;"> 73851091 </td>
   <td style="text-align:left;"> - </td>
   <td style="text-align:left;"> 19961 </td>
   <td style="text-align:left;"> 99939 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000074800.14 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:left;"> 8861000 </td>
   <td style="text-align:left;"> 8879190 </td>
   <td style="text-align:left;"> - </td>
   <td style="text-align:left;"> 5276 </td>
   <td style="text-align:left;"> 88419 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000198886.2 </td>
   <td style="text-align:left;"> chrM </td>
   <td style="text-align:left;"> 10760 </td>
   <td style="text-align:left;"> 12137 </td>
   <td style="text-align:left;"> + </td>
   <td style="text-align:left;"> 1378 </td>
   <td style="text-align:left;"> 86567 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000111640.14 </td>
   <td style="text-align:left;"> chr12 </td>
   <td style="text-align:left;"> 6533927 </td>
   <td style="text-align:left;"> 6538371 </td>
   <td style="text-align:left;"> + </td>
   <td style="text-align:left;"> 2981 </td>
   <td style="text-align:left;"> 79933 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000184009.10 </td>
   <td style="text-align:left;"> chr17 </td>
   <td style="text-align:left;"> 81509971 </td>
   <td style="text-align:left;"> 81523847 </td>
   <td style="text-align:left;"> - </td>
   <td style="text-align:left;"> 2993 </td>
   <td style="text-align:left;"> 75963 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000198712.1 </td>
   <td style="text-align:left;"> chrM </td>
   <td style="text-align:left;"> 7586 </td>
   <td style="text-align:left;"> 8269 </td>
   <td style="text-align:left;"> + </td>
   <td style="text-align:left;"> 684 </td>
   <td style="text-align:left;"> 70428 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSG00000067225.17 </td>
   <td style="text-align:left;"> chr15 </td>
   <td style="text-align:left;"> 72199029 </td>
   <td style="text-align:left;"> 72231822 </td>
   <td style="text-align:left;"> - </td>
   <td style="text-align:left;"> 8571 </td>
   <td style="text-align:left;"> 69910 </td>
  </tr>
</tbody>
</table></div>



### <a name="DEG"/>Step 4-1: Identify Differentially Expressed Genes (DEGs) | [Top](#Top)


In this step, the pipeline will identify differentially expressed genes (DEG) according to the alignment result files (i.e., BAM files) after the alignment step.

The BDS code snippets for the example dataset will look like:


```bash
$ grep -A1 run.call.*.featurecounts.star.*.sh Submit_*.bds
```

```
## dep( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/DEG/edger/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.edger.count.txt', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/DEG/edger/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.edger.test.DEG.txt' ] <- [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/WT01/WT01.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/WT02/WT02.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/WT03/WT03.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/KO01/KO01.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/KO02/KO02.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/KO03/KO03.star.featurecounts.count' ], cpus := 4, mem := 32*G, timeout := 72*hour, taskName := "edger.featurecounts.star.CRI-BIO-646-BMB-RKeenan") sys bash /gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/shell_scripts/run.call.edger.featurecounts.star.CRI-BIO-646-BMB-RKeenan.sh; sleep 2
## goal( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/DEG/edger/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.edger.count.txt' ] )
## --
## dep( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/DEG/deseq2/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.count.txt', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/DEG/deseq2/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.test.DEG.txt' ] <- [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/WT01/WT01.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/WT02/WT02.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/WT03/WT03.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/KO01/KO01.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/KO02/KO02.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/KO03/KO03.star.featurecounts.count' ], cpus := 4, mem := 32*G, timeout := 72*hour, taskName := "deseq2.featurecounts.star.CRI-BIO-646-BMB-RKeenan") sys bash /gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/shell_scripts/run.call.deseq2.featurecounts.star.CRI-BIO-646-BMB-RKeenan.sh; sleep 2
## goal( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/DEG/deseq2/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.count.txt' ] )
## --
## dep( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/DEG/limma/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.limma.count.txt', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/DEG/limma/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.limma.test.DEG.txt' ] <- [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/WT01/WT01.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/WT02/WT02.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/WT03/WT03.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/KO01/KO01.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/KO02/KO02.star.featurecounts.count', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/Quantification/featurecounts/star/KO03/KO03.star.featurecounts.count' ], cpus := 4, mem := 32*G, timeout := 72*hour, taskName := "limma.featurecounts.star.CRI-BIO-646-BMB-RKeenan") sys bash /gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/shell_scripts/run.call.limma.featurecounts.star.CRI-BIO-646-BMB-RKeenan.sh; sleep 2
## goal( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/DEG/limma/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.limma.count.txt' ] )
```

This code chunk will invoke few bash scripts (e.g., RNAseq/shell_scripts/**`run.call.edger.featurecounts.star.CRI-BIO-646-BMB-RKeenan.sh`**, **`run.call.deseq2.featurecounts.star.CRI-BIO-646-BMB-RKeenan.sh`**, and **`run.call.limma.featurecounts.star.CRI-BIO-646-BMB-RKeenan.sh`**) to execute differential expression (DE) analysis using three the state-of-the-art tools  (i.e., [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), and [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)) on the example dataset of six samples from KO01 to WT03.


There are three DE analysis tools used in the current pipeline, including

1. [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html): Empirical Analysis of Digital Gene Expression Data in R
2. [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html): Differential gene expression analysis based on the negative binomial distribution
3. [limma](https://bioconductor.org/packages/release/bioc/html/limma.html): Linear Models for Microarray Data

After the completion of the entire pipeline, you can check the calling results of each individual methods; for instance, the analysis results of the example dataset will be as follows.


```bash
$ tree RNAseq/DEG/*/featurecounts/star/
RNAseq/DEG/deseq2/featurecounts/star/
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.RData
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.count.ntd.meanSdPlot.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.count.ntd.txt
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.count.rld.meanSdPlot.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.count.rld.txt
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.count.txt
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.count.vst.meanSdPlot.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.count.vst.txt
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.plotDispEsts.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.plotMA.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.test.DEG.txt
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.test.txt
`-- run.call.deseq2.featurecounts.star.CRI-BIO-646-BMB-RKeenan.log
RNAseq/DEG/edger/featurecounts/star/
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.edger.RData
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.edger.count.txt
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.edger.plotBCV.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.edger.plotMA.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.edger.plotSmear.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.edger.test.DEG.txt
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.edger.test.txt
`-- run.call.edger.featurecounts.star.CRI-BIO-646-BMB-RKeenan.log
RNAseq/DEG/limma/featurecounts/star/
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.limma.RData
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.limma.count.txt
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.limma.count.voom.meanSdPlot.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.limma.plotMA.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.limma.test.DEG.txt
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.limma.test.txt
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.limma.voom.mean-variance.pdf
`-- run.call.limma.featurecounts.star.CRI-BIO-646-BMB-RKeenan.log
```

You can check statistical test results per gene (e.g., **RNAseq/DEG/deseq2/featurecounts/star/`CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.test.txt`**) for more information generated by each method.


```bash
$ head RNAseq/DEG/deseq2/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.test.txt
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>Statistical test result per gene using DESeq2</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Geneid </th>
   <th style="text-align:right;"> baseMean </th>
   <th style="text-align:right;"> log2FoldChange </th>
   <th style="text-align:right;"> lfcSE </th>
   <th style="text-align:right;"> stat </th>
   <th style="text-align:right;"> pvalue </th>
   <th style="text-align:right;"> FDR </th>
   <th style="text-align:right;"> DEG </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 8943 </td>
   <td style="text-align:left;"> ENSG00000139219.18 </td>
   <td style="text-align:right;"> 7593.6731 </td>
   <td style="text-align:right;"> 3.364988 </td>
   <td style="text-align:right;"> 0.0909209 </td>
   <td style="text-align:right;"> 37.01005 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 9845 </td>
   <td style="text-align:left;"> ENSG00000092067.5 </td>
   <td style="text-align:right;"> 3022.5103 </td>
   <td style="text-align:right;"> 3.543611 </td>
   <td style="text-align:right;"> 0.0958003 </td>
   <td style="text-align:right;"> 36.98957 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 14075 </td>
   <td style="text-align:left;"> ENSG00000101160.13 </td>
   <td style="text-align:right;"> 1270.4139 </td>
   <td style="text-align:right;"> 3.254167 </td>
   <td style="text-align:right;"> 0.1034734 </td>
   <td style="text-align:right;"> 31.44932 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 13557 </td>
   <td style="text-align:left;"> ENSG00000142530.10 </td>
   <td style="text-align:right;"> 1241.2705 </td>
   <td style="text-align:right;"> -3.737343 </td>
   <td style="text-align:right;"> 0.1221921 </td>
   <td style="text-align:right;"> -30.58581 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 9650 </td>
   <td style="text-align:left;"> ENSG00000136167.13 </td>
   <td style="text-align:right;"> 2428.1942 </td>
   <td style="text-align:right;"> 2.931033 </td>
   <td style="text-align:right;"> 0.0974567 </td>
   <td style="text-align:right;"> 30.07524 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 8362 </td>
   <td style="text-align:left;"> ENSG00000084207.16 </td>
   <td style="text-align:right;"> 1082.8763 </td>
   <td style="text-align:right;"> 3.350925 </td>
   <td style="text-align:right;"> 0.1188343 </td>
   <td style="text-align:right;"> 28.19830 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 9368 </td>
   <td style="text-align:left;"> ENSG00000111331.12 </td>
   <td style="text-align:right;"> 575.0466 </td>
   <td style="text-align:right;"> 3.131198 </td>
   <td style="text-align:right;"> 0.1228025 </td>
   <td style="text-align:right;"> 25.49785 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5104 </td>
   <td style="text-align:left;"> ENSG00000135346.8 </td>
   <td style="text-align:right;"> 374.1435 </td>
   <td style="text-align:right;"> -5.235472 </td>
   <td style="text-align:right;"> 0.2229213 </td>
   <td style="text-align:right;"> -23.48574 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 14910 </td>
   <td style="text-align:left;"> ENSG00000165194.15 </td>
   <td style="text-align:right;"> 1085.3843 </td>
   <td style="text-align:right;"> -2.145027 </td>
   <td style="text-align:right;"> 0.0973406 </td>
   <td style="text-align:right;"> -22.03630 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2440 </td>
   <td style="text-align:left;"> ENSG00000163053.10 </td>
   <td style="text-align:right;"> 2297.1436 </td>
   <td style="text-align:right;"> -1.763135 </td>
   <td style="text-align:right;"> 0.0818927 </td>
   <td style="text-align:right;"> -21.52983 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -1 </td>
  </tr>
</tbody>
</table></div>

>

Or, the exploratory plots of the example dataset produced by [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) will be as follows.

1. An exploratory plot of the per-gene dispersion estimates together with the fitted mean-dispersion relationship
    * <img src="IMG/CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.plotDispEsts.png" alt="KO01_DispEsts"  width="600" height="600">
2. An exploratory plot of row standard deviations versus row means using the normalized counts transformation (f(count + pc))
    * <img src="IMG/CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.count.ntd.meanSdPlot.png" alt="KO01_ntd_meanSd"  width="600" height="600">
3. An exploratory plot of row standard deviations versus row means using the variance stabilizing transformation
    * <img src="IMG/CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.count.rld.meanSdPlot.png" alt="KO01_rld_meanSd"  width="600" height="600">
4. An exploratory plot of row standard deviations versus row means using the 'regularized log' transformation
    * <img src="IMG/CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.count.vst.meanSdPlot.png" alt="KO01_vst_meanSd"  width="600" height="600">
5. A MA plot of log2 fold changes (on the y-axis) versus the mean of normalized counts (on the x-axis)
    * <img src="IMG/CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.plotMA2.png" alt="KO01_plotMA"  width="600" height="600">
6. A scatter plot of log2 fold changes (on the y-axis) versus the FDR (on the x-axis)
    * <img src="IMG/CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.plotMA.png" alt="KO01_plotMA"  width="600" height="600">



<!-- ## **<span style="color:red">RUN AS PRACTICE</span>** -->

<!-- >  -->
<!-- >  -->
<!-- > To demonstrate the full power of the following analyses, the pipeline will use the pre-run count tables and DEG lists from the full sequencing datasets.   -->
<!-- > Please make sure that there is no any changes on the metadata file (i.e., <span style="color:red">```CRI-BIO-646/CRI-BIO-646.metadata.txt```</span>).   -->
<!-- > Otherwise, it will be considered as not running for practice and should expect the following results not be the same as shown here.   -->
<!-- >  -->
<!-- >  -->



### <a name="LociStat"/>Step 4-2: DEG Statistics | [Top](#Top)


In this step, the pipeline will collect DEG statistics and identify the overlapping set of identified DEGs from the previous methods.

The BDS code snippet for the sample KO01 will look like:


```bash
$ grep -A1 run.lociStat.featurecounts.star.*.sh Submit_*.bds
```

```
## dep( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/LociStat/featurecounts/star/CRI-BIO-646-BMB-RKeenan/CRI-BIO-646-BMB-RKeenan.star.featurecounts.overlap.txt' ] <- [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/DEG/edger/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.edger.test.DEG.txt', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/DEG/deseq2/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.test.DEG.txt', '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/DEG/limma/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.limma.test.DEG.txt' ], cpus := 4, mem := 32*G, timeout := 72*hour, taskName := "lociStat.featurecounts.star.CRI-BIO-646-BMB-RKeenan") sys bash /gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/shell_scripts/run.lociStat.featurecounts.star.CRI-BIO-646-BMB-RKeenan.sh; sleep 2
## goal( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/LociStat/featurecounts/star/CRI-BIO-646-BMB-RKeenan/CRI-BIO-646-BMB-RKeenan.star.featurecounts.overlap.txt' ] )
```

This code chunk will invoke the bash script (e.g., RNAseq/shell_scripts/**`run.lociStat.featurecounts.star.CRI-BIO-646-BMB-RKeenan.sh`**) to collect DEG statistics and to make a Venn diagram plot.

After the completion of the entire pipeline, you can check the statistics result of DEGs per method; for instance, the example dataset CRI-BIO-646 will be as follows.


```bash
$ grep -A5 'Up/Down regulated DEGs per methods' RNAseq/LociStat/featurecounts/star/*/run.lociStat.featurecounts.star.*.log | tail -n+2
```

```
## grep: result/run.lociStat.featurecounts.star.*.log: No such file or directory
```



```bash
$ cut -f1,2,4 RNAseq/LociStat/featurecounts/star/*/*.star.featurecounts.VennList.txt
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>DEG Statistics</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Methods </th>
   <th style="text-align:right;"> Method.Num </th>
   <th style="text-align:right;"> ID.Num </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> edger </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> deseq2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> limma </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 68 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> edger&amp;deseq2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 45 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> edger&amp;limma </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 107 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> deseq2&amp;limma </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> edger&amp;deseq2&amp;limma </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 2669 </td>
  </tr>
</tbody>
</table></div>

  
There is a Venn diagram plot will be generated after this step.

<img src="IMG/CRI-BIO-646-BMB-RKeenan.star.featurecounts.VennDiagram.png" alt="Venn"  width="600" height="600">



### <a name="QuantQC"/>Step 5: Sample Correlation | [Top](#Top)


In this step, the pipeline will make a PCA plot based on the transcriptional profiling of all samples.

The BDS code snippet for the sample KO01 will look like:


```bash
$ grep -A1 run.quantQC.pca.featurecounts.star.*.sh Submit_*.bds
```

```
## dep( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/QuantQC/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.pca.pdf' ] <- [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/DEG/deseq2/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.deseq2.count.txt' ], cpus := 8, mem := 64*G, timeout := 72*hour, taskName := "pca.featurecounts.star.CRI-BIO-646-BMB-RKeenan") sys bash /gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/shell_scripts/run.quantQC.pca.featurecounts.star.CRI-BIO-646-BMB-RKeenan.sh; sleep 2
## goal( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/QuantQC/featurecounts/star/CRI-BIO-646-BMB-RKeenan.star.featurecounts.pca.pdf' ] )
```

This code chunk will invoke the bash script (e.g., RNAseq/shell_scripts/**`run.quantQC.pca.featurecounts.star.CRI-BIO-646-BMB-RKeenan.sh`**) to make a PCA plot based on the alignment quantification result generated by [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) or one of DE analysis tools.

After the completion of the entire pipeline, you can check the PCA plot under the folder of `QuantQC/`.

<img src="IMG/CRI-BIO-646-BMB-RKeenan.star.featurecounts.pca.png" alt="PCA"  width="600" height="600">



### <a name="HeatMap"/>Step 6: Heat Map | [Top](#Top)


In this step, the pipeline will make a heat map based on the overlapping set of DEGs identified across different DE analysis tools.

The BDS code snippet for the sample KO01 will look like:


```bash
$ grep -A1 run.postAna.pheatmap.featurecounts.star.*.sh Submit_*.bds
```

```
## dep( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/PostAna/pheatmap/featurecounts/star/CRI-BIO-646-BMB-RKeenan/CRI-BIO-646-BMB-RKeenan.star.featurecounts.heatmap.pdf' ] <- [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/LociStat/featurecounts/star/CRI-BIO-646-BMB-RKeenan/CRI-BIO-646-BMB-RKeenan.star.featurecounts.overlap.txt' ], cpus := 4, mem := 32*G, timeout := 72*hour, taskName := "postAna.pheatmap.featurecounts.star.CRI-BIO-646-BMB-RKeenan") sys bash /gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/shell_scripts/run.postAna.pheatmap.featurecounts.star.CRI-BIO-646-BMB-RKeenan.sh; sleep 2
## goal( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/PostAna/pheatmap/featurecounts/star/CRI-BIO-646-BMB-RKeenan/CRI-BIO-646-BMB-RKeenan.star.featurecounts.heatmap.pdf' ] )
```

This code chunk will invoke the bash script (e.g., RNAseq/shell_scripts/**`run.postAna.pheatmap.featurecounts.star.CRI-BIO-646-BMB-RKeenan.sh`**) to make a heat map plot based on the overlapping set of DEGs identified across different DE analysis tools (i.e., [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), and [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)).

After the completion of the entire pipeline, you can check the heat map under the folder of `PostAna/pheatmap`.

<img src="IMG/CRI-BIO-646-BMB-RKeenan.star.featurecounts.heatmap.png" alt="CRI-BIO-646_full_Venn"  width="600" height="600">



### <a name="GSEA"/>Step 7: Functional Enrichment Analysis | [Top](#Top)


In this step, the pipeline will conduct enrichment analysis and make several exploratory plots based on the overlapping set of DEGs identified across different DE analysis tools.

The BDS code snippet for the projet CRI-BIO-646-BMB-RKeenan will look like:


```bash
$ grep -A1 run.postAna.clusterprofiler.featurecounts.star.*.sh Submit_*.bds
```

```
## dep( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/PostAna/clusterprofiler/featurecounts/star/CRI-BIO-646-BMB-RKeenan/CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichGO.ALL.txt' ] <- [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/LociStat/featurecounts/star/CRI-BIO-646-BMB-RKeenan/CRI-BIO-646-BMB-RKeenan.star.featurecounts.overlap.txt' ], cpus := 4, mem := 32*G, timeout := 72*hour, taskName := "postAna.clusterprofiler.featurecounts.star.CRI-BIO-646-BMB-RKeenan") sys bash /gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/shell_scripts/run.postAna.clusterprofiler.featurecounts.star.CRI-BIO-646-BMB-RKeenan.sh; sleep 2
## goal( [ '/gpfs/data/bioinformatics/wchan10/CRI-BIO-646-BMB-RKeenan/RNAseq/PostAna/clusterprofiler/featurecounts/star/CRI-BIO-646-BMB-RKeenan/CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichGO.ALL.txt' ] )
```

This code chunk will invoke the bash script (e.g., RNAseq/shell_scripts/**`run.postAna.clusterprofiler.featurecounts.star.CRI-BIO-646-BMB-RKeenan.sh`**) to conduct enrichment analyses including [GO](http://www.geneontology.org/) and [KEGG](https://www.genome.jp/kegg/) pathway erichment analyses as well as gene set enrichment analysis (GSEA).

After the completion of the entire pipeline, you can check the heat map under the folder of `PostAna/pheatmap`.


```bash
$ tree RNAseq/PostAna/clusterprofiler/featurecounts/star/*/
RNAseq/PostAna/clusterprofiler/featurecounts/star/CRI-BIO-646-BMB-RKeenan/
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichGO.ALL.cnetplot.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichGO.ALL.dotplot.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichGO.ALL.emapplot.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichGO.ALL.txt
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichGSEAGO.ALL.pos001.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichGSEAGO.ALL.txt
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichKEGG.cnetplot.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichKEGG.dotplot.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichKEGG.emapplot.pdf
|-- CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichKEGG.txt
`-- run.postAna.clusterprofiler.featurecounts.star.CRI-BIO-646-BMB-RKeenan.log
```

Below, several plots will be generated based on the overlapping set of DEGs using [GO](http://www.geneontology.org/) database as example.

1. An exploratory dot plot for enrichment result
    * <img src="IMG/CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichGO.ALL.dotplot.png" alt="CRI-BIO-646_full_enrichGO.ALL.dotplot"  width="600" height="600">
2. An exploratory enrichment map for enrichment result of over-representation test
    * <img src="IMG/CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichGO.ALL.emapplot.png" alt="CRI-BIO-646_full_enrichGO.ALL.emapplot"  width="600" height="600">
3. An exploratory Gene-Concept Network plot of over-representation test
    * <img src="IMG/CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichGO.ALL.cnetplot.png" alt="CRI-BIO-646_full_enrichGO.ALL.cnetplot"  width="600" height="600">
4. An exploratory plot of GSEA result with NES great than zero
    * <img src="IMG/CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichGSEAGO.ALL.pos001.png" alt="CRI-BIO-646_full_enrichGSEAGO.ALL.pos001"  width="600" height="600">
<!-- 5. An exploratory plot of GSEA result with NES less than zero -->
<!--     * <img src="IMG/CRI-BIO-646-BMB-RKeenan.star.featurecounts.enrichGSEAGO.ALL.neg001.png" alt="CRI-BIO-646_full_enrichGSEAGO.ALL.neg001"  width="600" height="600"> -->


## <a name="BDS"/>BigDataScript Report | [Top](#Top)



Considering the environment setting in the CRI HPC system, [BigDataScript](https://pcingola.github.io/BigDataScript/) was used as a job management system in the current development to achieve an automatic pipeline. It can handle the execution dependency of all sub-task bash scripts and resume from a failed point, if any.

After the completion of the entire pipeline, you will see a BigDataScript report in HTML under the pipeline folder. For instance, this is the report from one test run. The graphic timeline will tell you the execution time per sub-task script.

<img src="IMG/Submit.PROJECT.bds.report.png" alt="Submit.RNAseq.CRI-BIO-646.bds"  width="800" height="600">



[Last Updated on 2018/11/07] | [Top](#Top)
