# MapUCEtoGenome

Scripts for mapping UCEs to genomes, plotting and analysing.

Mapping part generally following 
[integrating-functional-genomics-into-phylogenomics](https://github.com/matthewhvandam/integrating-functional-genomics-into-phylogenomics) 
and relying on  [calacademy's uce_types](https://github.com/calacademy-research/ccgutils/tree/master/uce_types).

Author of current repo: Hanaka Mera   mailto:hanaka.mera[@]my.jcu.edu.au

---

## Getting input data

1. For mapping UCEs: [1_MapUCEtoGenome_general.sh](https://github.com/mhanaka/MapUCEtoGenome/blob/main/1_MapUCEtoGenome_general.sh)
2. For ML tree: [2_phylogeny.sh](https://github.com/mhanaka/MapUCEtoGenome/blob/main/2_phylogeny.sh)

## Plotting

1. Plot phylogeny: [R_1_tree.R](https://github.com/mhanaka/MapUCEtoGenome/blob/main/R_1_tree.R)
2. Data wrangling mapped UCEs: [R_2_datawrangle.R](https://github.com/mhanaka/MapUCEtoGenome/blob/main/R_2_datawrangle.R)
3. Plot mapped UCEs (general overview): [R_3_plot.R](https://github.com/mhanaka/MapUCEtoGenome/blob/main/R_3_plot.R)
4. Plot mapped UCEs (using [chromoMap](https://lakshay-anand.github.io/chromoMap/docs.html)): [R_4_chromomap.R](https://github.com/mhanaka/MapUCEtoGenome/blob/main/R_4_chromomap.R)
5. Extract annotations/gene descriptions: [R_5_annotations.R](https://github.com/mhanaka/MapUCEtoGenome/blob/main/R_5_annotations.R)
