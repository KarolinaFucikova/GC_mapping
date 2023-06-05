# GC_mapping

Evolution of GC content in Green algae

The R project contains the data on a selection of Chlorophyte green algae. 
The table (GC_content.csv) contains information on 69 species of green algae: their habitat of origin (aquatic, soil, snow), the precipitation and min and max temperatures at the sites of the algae's origins (average in the warmest and coldest month, respectively), their GC content in chloroplast genomes (overall and coding-only), and GC in their nuclear ribosomal genes (18S and 28S), as well as GC of their chloroplast rrs gene (missing for several taxa due to the partial nature of the genome sequences in GenBank). 

The tree file contains a Bayesian consensus tree based on chloroplast protein-coding genes, onto which the GC content and temperature data (or any other traits from the table) can be mapped. The R script allows for drawing the tree, mapping the temperature and GC data on it (after dropping tips with no data for that variable), and also graphs correlations among selected variables using ggplot. Further, the code also calculates phylogenetically independent contrasts, allowing for a phylogeny-corrected regression between temperature and GC content. Several plotting options and color schemes are explored and can be recreated using the script. Example plots are saved in the plots folder.

