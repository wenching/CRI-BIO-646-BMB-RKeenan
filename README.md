# CRI RNAseq 2018

RNA-Seq Analysis Pipeline based on [CRI](http://cri.uchicago.edu/) HPC system

## Dataset

The RNA-seq data used in this tutorial are from the project CRI-BIO-646-BMB-RKeenan.

## File description

This repository contains the following items:
- ```docs/README.html``` - the main tutorial documentation
- ```CRI-BIO-646/``` - folder accommodating all data (i.e., metadata file, configuration file, sequencing data folder, and references folder) for running this pipeline
- ```README.md``` - this description file
- ```SRC/``` - automatic pipelines for RNA-seq analysis

### Prerequisites

* [R](https://www.r-project.org/)
* [CRAN](https://cran.r-project.org/)
* [Bioconductor](https://www.bioconductor.org/)
* [Python 3](https://www.python.org/download/releases/3.0/)
* [BigDataScript](https://pcingola.github.io/BigDataScript/)

### Installing

```bash
# download the package
wget ftp://logia.cri.uchicago.edu/bioinformatics/tutorials/Nov2018/CRI-BIO-646-BMB-RKeenan.tgz
# Or, copy directly from HPC
cp /gpfs/data/bioinformatics/shared/tutorials/Nov2018/CRI-BIO-646-BMB-RKeenan.tgz .

# uncompress the tarball file
tar -zxvf CRI-BIO-646-BMB-RKeenan.tgz

# change working directory to the package directory
cd CRI-BIO-646-BMB-RKeenan

# load modules
$ module purge;module load gcc udunits python/3.6.0 R/3.5.0; module update

# This step is optional but it will install all necessary R packages ahead.
# In case the pipeline was terminated due to the failure of R package installation later when running the pipeline.
$ Rscript --vanilla SRC/R/util/prerequisite.packages.R

# create directories and generate all necessary scripts
$ bash Build_RNAseq.CRI-BIO-646.sh

# run the entire pipeline with just this command
$ bash Submit_CRI-BIO-646-BMB-RKeenan.sh
```

## Contributing

Please read [CONTRIBUTING.md](https://github.com/wenching/cri_rnaseq_2018/blob/master/CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

For the versions available, see the [tags on this repository](https://github.com/wenching/cri_rnaseq_2018/tags).

## Authors

* **Wen-Ching Chan** - *Initial work* - [Wen-Ching](https://github.com/wenching)

See also the list of [contributors](https://github.com/wenching/cri_rnaseq_2018/graphs/contributors) who participated in this project.

## License

This project is licensed under the [LGPLv3](https://www.gnu.org/licenses/lgpl-3.0.en.html) License - see the [LICENSE](LICENSE) file with a copy for details

## Acknowledgments

* Thank [Kyle Hernandez](https://github.com/kmhernan) for providing a private Git repository in the begining of this project
* Thank [Riyue Bao](https://github.com/riyuebao) for providing tutorial dataset, figures, and the pipeline in the Perl version
* Thank [PurpleBooth](https://gist.github.com/PurpleBooth) for providing the templates of [README.md](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2) & [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426)


