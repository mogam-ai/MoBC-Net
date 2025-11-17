
# MoBCnet

<!-- badges: start -->
<!-- badges: end -->

The goal of MoBCnet is to ...

## Installation

You can install the development version of MoBCnet from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hanbi-99/MoBC-net")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MoBCnet)
library(dplyr)

## basic example code

getwd()
setwd("/Users/userid/test")
dir()

#----- load ppi and module info
edge_list.file <- system.file("data","toy_net1_A.txt", package='MoBCnet')
membership_info.file <- system.file("data","toy_net1_M.txt", package='MoBCnet')

ppi <- read.csv(edge_list.file,header = FALSE,sep='')
membership_info <- read.csv(membership_info.file,header = FALSE,sep='')
colnames(membership_info) = c('node','module')
modules=split(membership_info$node, membership_info$module)

#----- run MoBC
# run community distance
dist.res <- CommuinityDistance(ppi,modules,randomMethod='RandSDM',random=1000, ratio=0.5,nCore=3)
re1 = dist.res$distance
re2 = dist.res$filtered.modules
re3 = dist.res$graph

mst.net = dist.res$mst.net
plot(mst.net)

# run MoBC
mobc.res <- MoBC.genes(ppi,module1.gene='1',module2.gene = '2',randomMethod='RandSDM',random=1000, ratio=0.5,nCore=3)
mobc.res %>% head

```

