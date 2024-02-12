# Automated_IP

This repo provides the code related to the analysis pipelines reported the manuscript "Automated tumor immunophenotyping predicts clinical benefit to anti-PD-L1 immunotherapy" by Xiao Li, et.al. Specifically, 

* "forest_plot.R" generates figure.5
* "Fig2.R" generates figure.2
* "Supp.fig1.R" generates supplementary figure.1
* "plot_confusion_matrix.py" generates figure.4 panel (A)
* "POP_OAK_IMP130_pip4_publication" does the PFS and OS analysis of OAK and IMPassion130, for the propoased MICHA-BITE method (figure.3, figure.4 panel (B), supplementary figure.2, supplementary figure.3)
* "POP_OAK_IMP130_pip1-3_publication" does the PFS and OS analysis of OAK and IMPassion130, for Density cutoff, Binned CD8 density and MOCHA methods (supplementary figure.4-7). 
* "HKoeppen_IMP_cd8_10x_tiles_nonoverlap.m" creates the slide level and tile based measurements of CD8 staining localized to pan-CK positive and pan-CK negative regions on IMPassion130.
* "HKoeppen_OAK_cd8_10x_tiles_nonoverlap.m" creates the slide level and tile based measurements of CD8 staining localized to pan-CK positive and pan-CK negative regions on OAK.
* "HKoeppen_POPLAR_cd8_10x_tiles_nonoverlap.m" creates the slide level and tile based measurements of CD8 staining localized to pan-CK positive and pan-CK negative regions on Poplar.
* "HKoeppen_OAK_cd8_10x_tiles_noCK.m" creates the slide level and tile based measurements of CD8 staining irrespective of pan-CK positive or pan-CK negative regions on OAK.
* "HKoeppen_OAK_cd8_tsne.m" creates the various tSNE plot visualizations, trains a multiclass SVM on POPLAR datapoints, and generates predictions on OAK and IMPassion130 datapoints.
* "mil_Auto_IP.ipynb" trains the multi-instance learning model using POPLAR and test on OAK and IMPassion130.

