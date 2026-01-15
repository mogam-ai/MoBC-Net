
## MoBC-Net

**MoBC-Net** is an R package for identifying statistically significant interacting module pairs and their key link genes in biological networks.

Our method selects interacting module pairs based on a **closest-distance measure** between modules, and identifies **link genes** using **module betweenness centrality (MoBC)**. To assess the robustness and significance of the detected interactions, MoBC-Net employs an **adaptive binning and randomization scheme** that preserves key network properties, including **module size**, **node degree distribution**, and **network modularity**.

By integrating network topology with rigorous statistical evaluation, MoBC-Net enables systematic discovery of biologically meaningful module–module interactions and their critical connector genes.

<br>

<p align="center">
  <img src="https://github.com/user-attachments/assets/b1fdc11d-6927-45af-9668-777d29b85bc3"
       alt="image"
       width="800" />
</p>

<br>


---
### installation
```r
install.packages("devtools")
devtools::install_github("hanbi-99/MoBC-net")
```


### Data load
```r
library(dplyr)
library(MoBCnet)
library(reshape2)

#--- load data & setting
adjm.file =  system.file("data","A_N1000.txt", package='MoBCnet')
module.file =  system.file("data","M_N1000.txt", package='MoBCnet')

adjm <- read.table(adjm.file)
module <- read.table(module.file)
modules = split(rownames(module),module[,1])[-1]

ppi = reshape2::melt(as.matrix(adjm))
ppi$Var2 = gsub('V','',ppi$Var2) %>% as.numeric
ppi = subset(ppi, value==1)[,1:2]

module.gene.list = modules
```


### Caculate community distance
```r
dist.res <- CommuinityDistance(ppi,modules,randomMethod='RandSDM',random=1000, ratio=0.1,nCore=3)
re1 = dist.res$distance
re2 = dist.res$filtered.modules
re3 = dist.res$graph
```


### Get plot dendrogram
```r
xre = re1; 
zmin = min(xre$z_score)
xre$weight = xre$z_score -zmin +0.01

xre1 = xre[,c(2,1,3:7)] %>% 'colnames<-'(colnames(xre))
xre = rbind(xre,xre1)

xrem = reshape2::dcast(xre, Module1~Module2, value.var='weight')
rownames(xrem) = xrem$Module1; xrem = xrem[,-1]
xrem = as.dist(xrem)

hc = hclust(xrem, method='single')
plot(hc, cex=2, hang=-1, yaxt = "n",main=' ',ylab=' ')

# weight 기준 y축 위치
yticks = pretty(hc$height)

# z_score로 변환한 label
zlabels = yticks + zmin - 0.01
axis(2, at = yticks, labels = round(zlabels, 2))
mtext("z_score", side = 2, line = 3)
```


### Get link genes
```r
mobc.res <- MoBC.genes(ppi, module.gene.list=modules, module1='1', module2 = '2', randomMethod='RandSDM',random=1000, ratio=0.07,nCore=3)
mobc.res %>% head
```


### Get plot for link path
```r
m1 = modules[['1']]
m2 = modules[['2']]
linkgene='200'
link.gene.path(ppi, m1, m2, linkgene, col1='pink',col2='lightblue',link.col='forestgreen')

```
