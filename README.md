# Phylogeography
R code to generate Figure 5, trees and accuracies.

## Code Layout
### analysis
* meannode.tree : the MCC tree generated from BEAST
* tabel.txt : text file that match tips to locations (states)
* implementation.R : R code for randomly sampling Guinea or Sierra Leone. Generating accuracies and rates.
* travelers.R : R code for sampling travelers. Tree generated by ace function contains parent nodes' states, use for loop to find each tips' parents state to define travellers. 


## package and software versions
* strap
* stringr
* ips
* ggtree
* treeio
* ggplot2
* reshape2
* ggstance
* diversitree
* phytools
* ggimage
* ggpubr
* dplyr
