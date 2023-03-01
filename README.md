## Welcome to DEcode!

The goal of this project is to enable you to utilize genomic big data in identifying regulatory mechanisms for differential expression (DE).

DEcode predicts inter-tissue variations and inter-person variations in gene expression levels from TF-promoter interactions, RNABP-mRNA interactions, and miRNA-mRNA interactions.

You can read more about this method in [this paper](https://doi.org/10.1038/s42256-020-0201-6) (full text is available at [https://rdcu.be/b5r3p](https://rdcu.be/b5r3p)) where we conducted a series of evaluation and applications by predicting transcript usage, drivers of aging DE, gene coexpression relationships on a genome-wide scale, and frequent DE in diverse conditions.

<img src="https://raw.githubusercontent.com/stasaki/DEcode/master/img/fig1.jpg" width="80%">

### Run DEcode on Code Ocean

You can run DEcode on [Code Ocean platform](https://doi.org/10.24433/CO.0084803.v1) without setting up a computational environment. Our Code Ocean capsule provides reproducible workflows, all processed data, and pre-trained models for tissue- and person-specific transcriptomes and DEprior, at gene- or transcript level.   

### Add RNA features
#### Prerequisites
Before running the code in this repository, you must have the following python libraries installed on your system:

- pandas
- scipy

Also you need to install the following R packages:

- GenomicFeatures
- tidyverse
- data.table
- rtracklayer
- plyranges
- optparse

#### Usage
- First, download the example files by running the following commands:
```bash
#GTF file from GTEXv7 (hg19)
mkdir gtf
wget https://storage.googleapis.com/gtex_analysis_v7/reference/gencode.v19.genes.v7.patched_contigs.gtf -P ./gtf/
#eCLIP-seq peaks from Encode (hg19)
mkdir bed
wget https://www.encodeproject.org/files/ENCFF039BKT/@@download/ENCFF039BKT.bed.gz -P ./bed/
wget https://www.encodeproject.org/files/ENCFF379UQU/@@download/ENCFF379UQU.bed.gz -P ./bed/
```
This will create directories gtf and bed, and download the GTF file and two eCLIP-seq bed files into them.

Our input data for the gene-level model was constructed based on the gencode.v19.transcripts.patched_contigs.gtf file from the GTEXv7 dataset. This file contains only one representative transcript for each gene for the human genome (hg19).

If you use a different gene model, please make sure to select a representative transcript for each gene before running the pipeline.

The input data for the transcript-level model was created based on  `https://storage.googleapis.com/gtex_analysis_v7/reference/gencode.v19.transcripts.patched_contigs.gtf`. It is not necessary to filter the GTF file for the transcript-level model.

1. Convert genome coordinates to RNA coordinates using the following command:
```bash
Rscript functions/bed_to_RNA_coord.R -b ./bed/ -n 100 -g gtf/gencode.v19.genes.v7.patched_contigs.gtf -o custom
```
This will convert bed files in the genome coordinates in the ./bed/ directory to RNA coordinates using the gencode.v19.genes.v7.patched_contigs.gtf file and output as custom.txt.

- To convert RNA-coordinate peaks to Pandas format, use the following command:
```bash
python functions/to_sparse.py custom.txt
```
This will convert the RNA-coordinate peaks in the custom.txt file to a sparse Pandas DataFrame (custom.pkl).

- Place custom.pkl in the directory where RNA features are located, for example: ./data/toy/RNA_features/.

- Modify the code (Run_DEcode_toy.ipynb) as follows:
```python
mRNA_data_loc = "./data/toy/RNA_features/"
mRNA_annotation_data = ["POSTAR","TargetScan","custom"]
```
This will update the location of the mRNA data and specify that the custom.pkl file should be used as part of the RNA annotation data.

### Add promoter features



#### If you find DEcode useful in your work, please cite our manuscript.

Tasaki, S., Gaiteri, C., Mostafavi, S. & Wang, Y. Deep learning decodes the principles of differential gene expression.  Nature Machine Intelligence (2020) [[link to paper](https://doi.org/10.1038/s42256-020-0201-6)] (full text is available at [https://rdcu.be/b5r3p](https://rdcu.be/b5r3p))


#### Source databases for traning data. 
* GTEx transcriptome data - [GTEx portal](https://www.gtexportal.org/home/)
* Transcription factor binding peaks - [GTRD](https://doi.org/10.1093/nar/gky1128)
* RNA binding protein binding peaks - [POSTAR2](https://doi.org/10.1093/nar/gky830)
* miRNA binding locations - [TargetScan](https://doi.org/10.7554/eLife.05005)


