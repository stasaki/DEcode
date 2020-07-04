## Welcome to DEcode!

The goal of this project is to enable you to utilize genomic big data in identifying regulatory mechanisms for differential expression (DE).

DEcode predicts inter-tissue variations and inter-person variations in gene expression levels from TF-promoter interactions, RNABP-mRNA interactions, and miRNA-mRNA interactions.

You can read more about this method in [this paper](https://doi.org/10.1101/2020.01.10.894238) where we conducted a series of evaluation and applications by predicting transcript usage, drivers of aging DE, gene coexpression relationships on a genome-wide scale, and frequent DE in diverse conditions.

<img src="https://raw.githubusercontent.com/stasaki/DEcode/master/img/fig1.jpg" width="80%">

### Run DEcode on Google Colab

This tutorial shows you a way to run DEcode on [Google Colab](https://colab.research.google.com) that provides you free access to a ready-to-use machine learning environment with a high-end GPU.

1. Go to [Google Colab](https://colab.research.google.com) and sign in to your Google account.
2. Open Jupyter notebook.
    - Menu -> File -> Open notebook -> GITHUB tab
    - Search `https://github.com/stasaki/DEcode`
    - Select `Run_DEcode.ipynb`
3. Run each block of code.  

### Run DEcode on Code Ocean

You can run DEcode on [Code Ocean platform](https://doi.org/10.24433/CO.0084803.v1) without setting up computational environment. Our Code Ocean capsule providesr reproducible workflows, pre-processed data, and pre-trained models for tissue- and person-specific transcriptomes.   



### Citations

GTRD - [Yevshin, I., Sharipov, R., Kolmykov, S., Kondrakhin, Y. & Kolpakov, F. GTRD: a database on gene transcription regulation—2019 update. Nucleic Acids Res 47, D100-D105 (2019).](https://doi.org/10.1093/nar/gky1128)

POSTAR2 - [Zhu, Y. et al. POSTAR2: deciphering the post-transcriptional regulatory logics. Nucleic acids research 47, D203-D211 (2019).](https://doi.org/10.1093/nar/gky830)

TargetScan - [Agarwal, V., Bell, G. W., Nam, J. & Bartel, D. P. Predicting effective microRNA target sites in mammalian mRNAs. eLife 4 (2015).](https://doi.org/10.7554/eLife.05005)

DEcode - [Tasaki, S., Gaiteri, C., Mostafavi, S. & Wang, Y. Decoding differential gene expression. bioRxiv 2020.01.10.894238; doi: https://doi.org/10.1101/2020.01.10.894238](https://doi.org/10.1101/2020.01.10.894238)
