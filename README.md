# hiPSCchromatin-SVM
Repository for code associated with the paper "Signatures of digital Polycomb regulation in functional iPSC heterogeneity between individuals", Movilla Miangolarra et al., 2025 (https://www.biorxiv.org/content/10.1101/2025.07.25.666753v1)

## Contents

### Reproducing results of the paper
Differential expression analysis (for hiPSCs and pre-MEs) and annotation files within `Transcription_analysis` folder. To be run prior to the jupyter notebooks below, to obtain the differentially expressed genes to work with.
The Jupyter notebooks `PCA.ipynb`, `SVM.ipynb`, `SVM_prediction.ipynb` and `Mathematical_model.ipynb` provide the code to reproduce Figs. 2 to 6 in the paper. Note that data to run this code with will be uploaded upon publication.

### Main functions -- usage in other datasets
`main_funcs.py` contains the main functions used for the SVM-based analysis of hiPSC data (which enables downstream analysis).
The module is separate to ease the usage of similar workflows in other datasets. Consists of four functions, grouped in pairs:
- `binarize()`, which binarises transcriptomic data in order to feed it to `SVM_analysis()`, that will train an SVM for each gene and perform some post-processing of the resulting coefficients, to identify the regulatory logic of each gene.
- `binarize_test()`, which binarises transcriptomic data for training data and test data separately and can be used as a test of the classifier. `SVM_prediction()` trains the SVM on the training data and returns a prediction score (between 0 and 1, Platt scaling) and a prediction (0 or 1) for each gene, which can be used to check the accuracy of the classifier (e.g. AUROC).
