This directory contains the individual patient data files for the data set E-MTAB-6149, which has data on 5 patients with lung cancer.
The pateints are called 
Patient1
Patient2
Patient3
Patient4
Patient5

Every patient has both cancer and non-cancer cells.
While files with more genes are available, we put only the filtered version with 38/1279 genes in the repository.
These files have been preprocessed and some normalization has been applied.

The files with 38 genes are called:
<patient>_cancer_cells_thresh1_38genes.txt
<patient>_noncancer_cells_thresh1_38genes.txt
<patient>_allnontumor_cells_thresh1_38genes.txt

The files with 1269 genes are called:
<patient>_cancer_cells_thresh1_1269genes.txt
<patient>_noncancer_cells_thresh1_1269genes.txt
<patient>_allnontumor_cells_thresh1_1269genes.txt

In this data set 739/1269 genes encoding cell surface receptors were measured.

The cancer and non-cancer cells are separated and indicated in the file names.

There is one header row before the genes start on the second row.

There are two opening columns on the left of each file, and the first of the two columns contains the gene name.

Expression levels are on a log2 scale
