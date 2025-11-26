### installation
```r
install.packages("devtools")
devtools::install_github("hanbi-99/MoBC-net")
```



### Tutorial

```r
library(dplyr)
library(MoBCnet)
# tutorial


adjm.file =  ystem.file("data","A_N1000.txt", package='MoBCnet')
module.file =  ystem.file("data","M_N1000.txt", package='MoBCnet')

adjm <- read.table(adjm.file)
module <- read.table(module.file)
modules = split(rownames(module),module[,1])[-1]

#--- setting id
ppi = reshape2::melt(as.matrix(adjm))
ppi$Var2 = gsub('V','',ppi$Var2) %>% as.numeric
ppi = subset(ppi, value==1)[,1:2]

module.gene.list = modules


#--- load file
# run community distance
dist.res <- CommuinityDistance(ppi,modules,randomMethod='RandSDM',random=1000, ratio=0.1,nCore=3)
re1 = dist.res$distance
re2 = dist.res$filtered.modules
re3 = dist.res$graph

# run MoBC
mobc.res <- MoBC.genes(ppi, module.gene.list=modules, module1='1', module2 = '2', randomMethod='RandSDM',random=1000, ratio=0.07,nCore=3)
mobc.res %>% head

```
