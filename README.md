## Welcome to DEcode!

The goal of this project is to enable you to utilize genomic big data in identifying regulatory mechanisms for differential expression (DE).

DEcode predicts inter-tissue variations and inter-person variations in gene expression levels from TF-promoter interactions, RNABP-mRNA interactions, and miRNA-mRNA interactions.

You can read more about this method in [this paper](https://doi.org/10.1038/s42256-020-0201-6) (full text is available at https://rdcu.be/b5r3p) where we conducted a series of evaluation and applications by predicting transcript usage, drivers of aging DE, gene coexpression relationships on a genome-wide scale, and frequent DE in diverse conditions.

<img src="https://raw.githubusercontent.com/stasaki/DEcode/master/img/fig1.jpg" width="80%">

### Run DEcode on Google Colab

This tutorial shows you a way to run DEcode on [Google Colab](https://colab.research.google.com) that provides you free access to a ready-to-use machine learning environment with a high-end GPU.

1. Go to [Google Colab](https://colab.research.google.com) and sign in to your Google account.
2. Open Jupyter notebook.
    - Menu -> File -> Open notebook -> GITHUB tab
    - Search `https://github.com/stasaki/DEcode`
    - Select `Run_DEcode_toy.ipynb`
3. Run each block of code.

Gene level features for all genes are available at data/GTEx53_gene/DEcode_data.tar.gz. However, due to the memory limitation in Google Colab, you cannot train a model with all genes. Please consider setting up Keras with GPU in your environment or use Code Ocean platform.   

### Run DEcode on Code Ocean

You can run DEcode on [Code Ocean platform](https://doi.org/10.24433/CO.0084803.v1) without setting up a computational environment. Our Code Ocean capsule provides reproducible workflows, all processed data, and pre-trained models for tissue- and person-specific transcriptomes and DEprior, at gene- or transcript level.   



#### If you find DEcode useful in your work, please cite our manuscript.

Tasaki, S., Gaiteri, C., Mostafavi, S. & Wang, Y. Deep learning decodes the principles of differential gene expression.  Nature Machine Intelligence (2020) [[link to paper](https://doi.org/10.1038/s42256-020-0201-6)] (full text is available at https://rdcu.be/b5r3p)


#### Source databases for traning data. 
* GTEx transcriptome data - [GTEx portal](https://www.gtexportal.org/home/)
* Transcription factor binding peaks - [GTRD](https://doi.org/10.1093/nar/gky1128)
* RNA binding protein binding peaks - [POSTAR2](https://doi.org/10.1093/nar/gky830)
* miRNA binding locations - [TargetScan](https://doi.org/10.7554/eLife.05005)


