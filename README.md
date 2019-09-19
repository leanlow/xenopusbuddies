~~~~~~~~~~~~~~~~
Saha Lab Scripts
~~~~~~~~~~~~~~~~

----------------
shared_finder.py
----------------
Takes two gene list (csv file) and outputs the shared genes between the lists in a new csv file
Include the rank of the gene in its parent list, log2FoldChange, and padj 

USAGE: shared_finder.py <gene_list_dir_1> <gene_list_dir_2> <output_dir>
OUTPUT: output_dir_genelist1.csv_shared_results.csv

------------------
diff_homeo_bias.py
------------------
Takes two gene list of homeologs (csv file) and outputs the shared genes between the lists in a new csv file
Include the rank of the gene in its parent list, log2FoldChange, and padj. This does a THREE-WAY comparison.

USAGE: shared_finder.py <gene_list_dir_1> <gene_list_dir_2> <gene_list_dir_3> <output_dir>
OUTPUT: output_dir_genelist1.csv_shared_results.csv
