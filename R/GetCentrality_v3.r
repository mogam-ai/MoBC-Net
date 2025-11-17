


# #--------------------------- new


# get.freq <-function(g, snode, enode){
# 	edges = igraph::all_shortest_paths(g, snode, enode)
# 	edges = edges$res %>% lapply(function(xx) setdiff(names(xx), names(xx)[c(1,length(xx))]))
# 	# etab = edges %>% unlist
# 	return(edges) #!!
# }

# get.freq.v2 <-function(g, snode, enode){
# 	edges = igraph::all_shortest_paths(g, snode, enode)
# 	edges = edges$res %>% lapply(function(xx)  setdiff(names(xx), names(xx)[1]))
# 	# etab = edges %>% unlist
# 	return(edges) #!!
# }




# #' Calculate centrality between two modules from MoBC result 
# #' 
# #' 
# #' @title cal.MoBCgenes.values
# #' @param g graph
# #' @param community1 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
# #' @param community2 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
# #' @returns vector
# #' @export
# #' @examples
# #' cal.MoBCgenes.values(graph, 'module_1','module_2', allg)


# #--- new
# cal.MoBCgenes.values <- function(g, community1, community2, allg){
    
# 	scorevec = rep(0, length(igraph::V(g))) %>% 'names<-'(igraph::V(g)$name)
# 	shortestm = igraph::distances(g, community1, community2)
# 	rmin  = apply(shortestm,1,function(xx) colnames(shortestm)[which(xx %in% min(xx))])
# 	cmin  = apply(shortestm,2,function(xx) rownames(shortestm)[which(xx %in% min(xx))])


#     r.endNode = lengths(rmin) %>% sum
#     c.endNode = lengths(cmin) %>% sum

#     # comm1
#     r.sp.genel = sapply(names(rmin), function(start.node){
# 		end.node = rmin[[start.node]]
#         vv = sapply(end.node, function(endn){
#     		etab = get.freq(g, start.node, endn)
#             etab1 = etab %>% unlist %>% table
#             val = etab1/length(etab)
#     		return(val[allg])
#         }) %>% 'rownames<-'(allg)
#         vv[is.na(vv)] = 0
#         return(vv)
# 	})
#     rval = do.call(cbind, r.sp.genel)
#     if(ncol(rval)!=r.endNode) stop('Not matched column number (r)', call. = FALSE)
#     r.val = apply(rval,1,sum)

#     # comm2
#     c.sp.genel = sapply(names(cmin), function(start.node){
# 		end.node = cmin[[start.node]]
#         vv = sapply(end.node, function(endn){
#     		etab = get.freq(g, start.node, endn)
#             etab1 = etab %>% unlist %>% table
#             val = etab1/length(etab)
#     		return(val[allg])
#         }) %>% 'rownames<-'(allg)
#         vv[is.na(vv)] = 0
#         return(vv)
# 	})
#     cval = do.call(cbind, c.sp.genel)
#     if(ncol(cval)!=c.endNode) stop('Not matched column number (r)', call. = FALSE)
#     c.val = apply(cval,1,sum)



#     scorev = (r.val+c.val)/sum(r.endNode+c.endNode) #!!

# 	# score.df = data.frame(gene=allg,community1.score =as.numeric(r.val), community2.score=as.numeric(c.val), score=all.score)
#     # scorev = score.df$community1.score + score.df$community2.score
#     # names(scorev) = allg

# 	return(scorev)
# }



# # #---- old version
# # cal.MoBCgenes.values <- function(g, community1, community2, allg){
    
# # 	scorevec = rep(0, length(igraph::V(g))) %>% 'names<-'(igraph::V(g)$name)
# # 	shortestm = igraph::distances(g, community1, community2)
# # 	rmin  = apply(shortestm,1,function(xx) colnames(shortestm)[which(xx %in% min(xx))])

# #     # comm1
# #     r.sp.genel = sapply(names(rmin), function(start.node){
# # 		end.node = rmin[[start.node]]
# # 		etab = get.freq(g, start.node, end.node)
# # 		return(etab)
# # 	})
    
# #     r.pathn = sum(lengths(r.sp.genel))
# #     r.tab = unlist(r.sp.genel) %>% table %>% sort

# #     # comm2
# # 	cmin  = apply(shortestm,2,function(xx) rownames(shortestm)[which(xx %in% min(xx))])
# #     c.sp.genel = sapply(names(cmin), function(start.node){
# # 		end.node = cmin[[start.node]]
# # 		etab = get.freq(g, start.node, end.node)
# # 		return(etab)
# # 	})
    
# #     c.pathn = sum(lengths(c.sp.genel))
# #     c.tab = unlist(c.sp.genel) %>% table %>% sort


# #     r.num = length(community1)
# #     c.num = length(community2)

# #     r.score = r.tab/r.pathn*c.num/sum(r.num+c.num) #!!
# #     c.score = c.tab/c.pathn*r.num/sum(r.num+c.num) #!!
    
# #     r.score = r.score[allg] %>% 'names<-'(allg)
# #     c.score = c.score[allg] %>% 'names<-'(allg)
# #     r.score[is.na(r.score)] = 0
# #     c.score[is.na(c.score)] = 0

# # 	score.df = data.frame(gene=allg,community1.score =as.numeric(r.score), community2.score=as.numeric(c.score))
# #     scorev = score.df$community1.score + score.df$community2.score
# #     names(scorev) = allg

# # 	return(scorev)
# # }






# # g = res2@graph
# # community1 = res2@filtered.communities[[2]]
# # community2 = res2@filtered.communities[[3]]
# # random = 1000
# # ratio = 0.1

# # random part
# cal.MoBC.random <- function(g, community1, community2,random,ratio,cal.p,show.binning, nCore){

#     allg = igraph::V(g)$name %>% as.character()

#     cl1g = community1
#     cl2g = community2

#     deg = igraph::degree(g)
#     membership <- rep(0, length(deg))
#     membership[names(deg) %in% cl1g] <- 1
#     membership[names(deg) %in% cl2g] <- 2 
#     use.random = cal.p

#     cl1f = cl1g %>% sort %>% paste0(collapse='_')
#     cl1f = substr(cl1f,1,30)
#     cl2f = cl2g %>% sort %>% paste0(collapse='_')
#     cl2f = substr(cl2f,1,30)

#     kk = c(cl1f,cl2f) %>% sort
#     options=paste0(c(kk,random,ratio),collapse='_')

#     dirn1 = paste0('./MoBCtmp/MoBC/',options,'/',cl1f,'/',cal.p,'_',random)
#     dirn2 = paste0('./MoBCtmp/MoBC/',options,'/',cl2f,'/',cal.p,'_',random)
#     dir.flag = dir.exists(dirn1) & dir.exists(dirn2)
#     l.flag = ifelse(length(list.files(dirn1))==random & length(list.files(dirn2))==random,TRUE,FALSE)

#     if(dir.flag & l.flag){
#         tag1 = all(sapply(1:random, function(j) read.csv(paste0(dirn1,'/rand',j,'.csv'))[,1] %>% length)==length(cl1g))
#         tag2 = all(sapply(1:random, function(j) read.csv(paste0(dirn2,'/rand',j,'.csv'))[,1] %>% length)==length(cl2g))
#         length.flag = tag1 & tag2
#     } else{
#         length.flag = FALSE
#     }

#     if(dir.flag & l.flag & length.flag){
#         cat(paste0('You have tmp files for random sampling - ',cal.p,". We will use these files.\n"))
#         comm.distance.list = parallel::mclapply(1:random,mc.cores=nCore, function(j){
#             rs1 = read.csv(paste0(dirn1,'/rand',j,'.csv'))[,1]
#             rs2 = read.csv(paste0(dirn2,'/rand',j,'.csv'))[,1]
#             comm.distance = cal.MoBCgenes.values(g, rs1,rs2, allg) 
#             return(comm.distance)
#         })
#         comm.distance.list = do.call(cbind, comm.distance.list) %>% 'rownames<-'(allg)
#         return(comm.distance.list)
#     }


#     cat(paste0("You don't have tmp files for random sampling - ",cal.p,". Generating random samples will take time...  \n"))
#     dir.create(dirn1, recursive=TRUE,showWarnings = FALSE)
#     dir.create(dirn2, recursive=TRUE,showWarnings = FALSE)


#     # cal.MoBCgenes.values(g, cl1g.random, cl2g.random, allg) 
#     if(cal.p=='random1'){
#         # pb <- progress::progress_bar$new(total = random)
#         comm.distance.list = parallel::mclapply(1:random, mc.cores=nCore, function(j){
#             rsamplel = simple_sampling(deg, membership,random)                
#             write.csv(rsamplel[[1]],paste0(dirn1,'/rand',j,'.csv'), row.names=F)
#             write.csv(rsamplel[[2]],paste0(dirn2,'/rand',j,'.csv'), row.names=F)
#             comm.distance = cal.MoBCgenes.values(g, rsamplel[[1]], rsamplel[[2]], allg) 
#             return(comm.distance)
#         })
#         comm.distance.list = do.call(cbind, comm.distance.list) %>% 'rownames<-'(allg)

#     } else if(cal.p=='random2'){
#         hist.bin0 = estimate_deg_bag(deg, membership, ratiov=ratio, ncv=random) 
#         hist.bin = hist.bin0$node_bag
#         names(hist.bin) = 1:length(hist.bin)
#         comm.distance.list = parallel::mclapply(1:random, mc.cores=nCore, function(j){
#             # pb <- progress::progress_bar$new(total = random)
            
#             samplingN = sapply(hist.bin, function(xx) sum(names(xx) %in% cl1g)) %>% 'names<-'(names(hist.bin))
#             cl1g.random = lapply(names(samplingN), function(xn){
#                 use.bg = (hist.bin[[xn]])
#                 # set.seed(m+j)
#                 sample(use.bg, samplingN[[xn]], replace=FALSE)
#             }) %>% unlist %>% unique

#             samplingN = sapply(hist.bin, function(xx) sum(names(xx) %in% cl2g)) %>% 'names<-'(names(hist.bin))
#             cl2g.random = lapply(names(samplingN), function(xn){
#                 use.bg = setdiff(hist.bin[[xn]], c(cl1g.random)) 
#                 # set.seed(n+j)
#                 sample(use.bg, samplingN[[xn]], replace=FALSE)
#             }) %>% unlist %>% unique
#             write.csv(cl1g.random,paste0(dirn1,'/rand',j,'.csv'), row.names=F)
#             write.csv(cl2g.random,paste0(dirn2,'/rand',j,'.csv'), row.names=F)
#             comm.distance = cal.MoBCgenes.values(g, cl1g.random, cl2g.random, allg) 
#             # pb$tick()
#             return(comm.distance)
#         })
#         comm.distance.list = do.call(cbind, comm.distance.list) %>% 'rownames<-'(allg)

#     } else if(cal.p=='random4'){
#         S <- igraph::distances(g, algorithm = "unweighted")
#         re = estimate_deg_bag(deg, membership, ratiov=ratio, ncv = random)
#         # pb <- progress::progress_bar$new(total = random)
#         # comm.distance.list = c()
#         # for(j in 1:random){
#         comm.distance.list= parallel::mclapply(1:random, mc.cores=nCore, function(j){
#             rsamplel = modularity_sampling(re, deg, membership, S)                
#             write.csv(rsamplel[[1]],paste0(dirn1,'/rand',j,'.csv'), row.names=F)
#             write.csv(rsamplel[[2]],paste0(dirn2,'/rand',j,'.csv'), row.names=F)      
#             comm.distance = cal.MoBCgenes.values(g, rsamplel[[1]], rsamplel[[2]], allg) 
#             # comm.distance.list = cbind(comm.distance.list,as.matrix(comm.distance))
#             # pb$tick()
#             return(comm.distance)
#         })
#         comm.distance.list = do.call(cbind, comm.distance.list) %>% 'rownames<-'(allg)

#     } else if(cal.p=='random3'){
#         cat('Select randomization3 method')
#         S <- igraph::distances(g, algorithm = "unweighted")                                
#         m1 <- length(cl1g)
#         n1 <- length(cl2g)

#         tnodes <- seq_len(length(igraph::V(g)$name)) 

#         # clustermn 함수를 사용해 2개 그룹으로 나눈다고 가정 (사용자 정의 필요)
#         # 예) clusterAssignment <- clustermn(S[sind, sind], m, n)
#         #
#         # 여기서는 간단히, 유사 함수가 있다고 가정하고 아래처럼 호출합니다.
#         #
#         # 사용자는 clustermn 함수를 별도로 정의해서
#         # "m개와 n개로 구성되면서, 같은 모듈 내 노드끼리는 평균 거리가 짧도록"  
#         # 하는 방식을 구현해야 합니다.

        
#         # pb <- progress::progress_bar$new(total = random)
#         # comm.distance.list = c()
#         comm.distance.list= parallel::mclapply(1:random, mc.cores=nCore, function(j){
#             sind <- sample(tnodes, m1 + n1, replace = FALSE)
#             clusterAssignment <- clustermn(S[sind, sind], m1, n1)
#             Urand <- sind[clusterAssignment == 1]
#             Wrand <- sind[clusterAssignment == 2]      
#             write.csv(Urand,paste0(dirn1,'/rand',j,'.csv'), row.names=F)
#             write.csv(Wrand,paste0(dirn2,'/rand',j,'.csv'), row.names=F)   
#             comm.distance = cal.MoBCgenes.values(g, Urand, Wrand, allg) 
#             # comm.distance.list = cbind(comm.distance.list,as.matrix(comm.distance))

#             # pb$tick()
#         })
#         comm.distance.list = do.call(cbind, comm.distance.list) %>% 'rownames<-'(allg)

#     } else {
#         stop('Please enter the right method for randomization.', call.=FALSE)
#     }
#     cat('fin comm.dist\n')
#     # start.time-end.time
#     return(comm.distance.list)
# }






# #' Calculate centrality between two modules from MoBC result 
# #' 
# #' 
# #' @title cal.MoBCgenes
# #' @param g graph
# #' @param community1 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
# #' @param community2 The name of the module for which centrality is being calculated. This should be one of the communities provided as input
# #' @returns vector
# #' @export
# #' @examples
# #' cal.MoBCgenes(graph, 'community1','community2',random,ratio,cal.p, nCore)



# cal.MoBCgenes <- function(g, community1, community2,random,ratio,cal.p, nCore){

#     # cat('cal.MoBCgenes')

#     allg = igraph::V(g)$name %>% as.character()
#     scorev = cal.MoBCgenes.values(g, community1, community2, allg)

# 	score.df = data.frame(gene=allg,score=scorev)
# 	score.df$node_type = 'link'
# 	score.df$node_type[score.df$gene %in% c(community1, community2)] = 'community genes'
# 	score.df = score.df %>% dplyr::arrange(-score)
#     score.df = subset(score.df, score>0)

#     colix = c('gene','score','node_type')
#     # cat('test.\n')
#     if(cal.p!='None'){
#         cat('Normalized MoBC score will be provided.\n')
#         colix = c('gene','normalized_score','node_type')
#         random.mat = cal.MoBC.random(g, community1, community2,random,ratio,cal.p,show.binning=FALSE, nCore=nCore)
#         pval = sapply(score.df$gene, function(gn){
#             xval = score.df[match(gn, score.df$gene),'score']
#             pval = sum(random.mat[gn,]>xval)/random
#         })
#         nscore = sapply(score.df$gene, function(gn){
#             xval = score.df[match(gn, score.df$gene),'score']
#             zvalv = c(random.mat[gn,],xval)
#             (xval-mean(zvalv))/sd(zvalv)
#         })
#         score.df$normalized_score = nscore
#         score.df$pval = pval
#         pv = sapply(score.df$pval, function(vv) ifelse(vv>0.5, 1-vv,vv))
#         score.df$p.adj = p.adjust(pv,'BH')
#         colix = c('gene','score','normalized_score','node_type','pval','p.adj')

#     }
# 	return(score.df[,colix])
# }




# # cal.MoBCgenes_test <- function(g, community1, community2,random,ratio,cal.p, nCore){


# #     allg = igraph::V(g)$name %>% as.character()
# #     scorev = cal.MoBCgenes.values(g, community1, community2, allg)

# # 	score.df = data.frame(gene=allg,score=scorev)
# # 	score.df$node_type = 'link'
# # 	score.df$node_type[score.df$gene %in% c(community1, community2)] = 'community genes'
# # 	score.df = score.df %>% dplyr::arrange(-score)
# #     score.df = subset(score.df, score>0)

# #     colix = c('gene','score','node_type')

# #     if(cal.p!='None'){
# #         random.mat = cal.MoBC.random(g, community1, community2,random,ratio,cal.p,show.binning=FALSE, nCore=nCore)
# #         pval = sapply(score.df$gene, function(gn){
# #             xval = score.df[match(gn, score.df$gene),'score']
# #             pval = sum(random.mat[gn,]>xval)/random
# #         })
# #         score.df$pval = pval
# #         pv = sapply(score.df$pval, function(vv) ifelse(vv>0.5, 1-vv,vv))
# #         score.df$p.adj = p.adjust(pv,'BH')
# #         colix = c('gene','score','node_type','pval')

# #     }
# # 	return(list(res=score.df[,colix],random=random.mat))
# # }





# #' Calculate centrality between two modules from MoBC result 
# #' 
# #' 
# #' @title cal.FCgene
# #' @param g graph
# #' @param module1.name The name of the module for which centrality is being calculated. This should be one of the communities provided as input
# #' @param module2.name The name of the module for which centrality is being calculated. This should be one of the communities provided as input
# #' @returns data.frame
# #' @export
# #' @examples
# #' cal.FCgene(graph, 'module_1','module_2')




# cal.FCgene <- function(g, community1, community2){

#     gene.ix = igraph::V(g)$name

#     # re = igraph::all_shortest_paths(g, community1[1:3],community2[1:4])
#     # lapply(re$res, function(xx) names(xx)[1] ) %>% unlist %>% unique
#     # lapply(re$res, function(xx) names(xx)[length(xx)] ) %>% unlist %>% unique
    
#     re = sapply(community1, function(g1){
            
#         edges = igraph::all_shortest_paths(g, g1, community2)
#         end.ix = edges$res %>% sapply(function(xx) names(xx)[length(xx)])
#         end.ix = split(1:length(end.ix), end.ix)
#         resl = lapply(end.ix, function(ixix){
#             intg = edges$res[ixix] %>% lapply(function(xx) setdiff(names(xx), names(xx)[c(1,length(xx))]))
#             intg.tab = unlist(intg) %>% table
#             intg.tab = intg.tab/length(intg)
#             intg.tab[gene.ix] %>% 'names<-'(gene.ix)
#         })
#         res = do.call(rbind, resl) %>% as.matrix
#         res[is.na(res)] = 0
#         resv = apply(res,2,sum)
#     })
#     re1 = apply(re,1,sum)
#     re1 = re1/length(community1)/length(community2)

# 	score.df = data.frame(gene=gene.ix,score=re1[gene.ix])
#     score.df$node_type = 'link'
# 	score.df$node_type[score.df$gene %in% c(community1, community2)] = 'module genes'
# 	score.df = score.df %>% dplyr::arrange(-score)

# 	return(score.df)
# }


# #' Calculate centrality between two modules from MoBC result 
# #' 
# #' 
# #' @title Get.Centrality
# #' @param network results from CommDistFunction function
# #' @param module1.gene The name of the module for which centrality is being calculated. This should be one of the communities provided as input
# #' @param module2.gene The name of the module for which centrality is being calculated. This should be one of the communities provided as input
# #' @returns data.frame
# #' @export
# #' @examples
# #' Get.Centrality(MoBC.result, 'module_1','module_2')




# # MoBC.genes <- function(MoBC.result, module1.name, module2.name,random,ratio,cal.p=FALSE){
# # 	if(!is(MoBC.result, 'MoBCresult')){
# # 		stop("input should be MoBC class", call. = FALSE)
# # 	}
# # 	communities = MoBC.result@filtered.modules

# # 	if(!all(c(module1.name, module2.name) %in% names(communities))){
# # 		stop('module name should be included in name of pre-defined module', call. = FALSE)
# # 	}

# # 	cal.MoBCgenes(MoBC.result@graph, 
# # 					community1=MoBC.result@filtered.modules[[module1.name]], 
# # 					community2=MoBC.result@filtered.modules[[module2.name]],
# #                     random=random, ratio=ratio,cal.p=cal.p)
# # }


# MoBC.genes <- function(network,
#                              module1.gene, module2.gene,
#                              randomMethod=c('None','RandC','RandCD','RandCM','RandCDM'),
# 							 random = 1000,
#                              nCore=1,
#                              ratio = 0.1) {
#     overlap_filtering=TRUE
#     # cat(method,'\n')
#     if (is.character(randomMethod)){
#         randomMethod <- match.arg(randomMethod)
#         # cat(dist.function,'\n')
#         randomMethod <- switch(randomMethod,
#             None = 'None',
#             RandC = 'random1',
#             RandCD = 'random2',
#             RandCM = 'random3',
#             RandCDM = 'random4')
#     } else {
#         stop('Method function is wrong. Check the method function', call.=FALSE)
#     }
# 	if (is.null(module1.gene) | is.null(module2.gene)) {
# 		stop('Module gene list is empty', call.=FALSE)
# 	}
# 	if (any(is.na(module1.gene)) | any(is.na(module1.gene))) {
# 		stop('Please assign right node names', call.=FALSE)
# 	}
#     module.genelist = list(module1=module1.gene, module2=module2.gene)
# 	g.res  <- preprocessedNetwork(network)
#     comm.genelist <- CommunityGenelist(module.genelist, g.res, overlap_filtering = overlap_filtering)

# 	communities = comm.genelist

# 	x=cal.MoBCgenes(g.res, 
# 					community1=comm.genelist[['module1']], 
# 					community2=comm.genelist[['module2']],
#                     random=random, ratio=ratio,cal.p=randomMethod, nCore=nCore)
#     # x = subset(x, score>0)
# 	return(x)
# 	}



# # MoBC.genes_test <- function(network,
# #                              module1.gene, module2.gene,
# #                              randomMethod=c('None','RandC','RandCD','RandCM','RandCDM'),
# # 							 random = 1000,
# #                              nCore=1,
# #                              ratio = 0.1) {
# #     overlap_filtering=TRUE
# #     # cat(method,'\n')
# #     if (is.character(randomMethod)){
# #         randomMethod <- match.arg(randomMethod)
# #         # cat(dist.function,'\n')
# #         randomMethod <- switch(randomMethod,
# #             None = 'None',
# #             RandC = 'random1',
# #             RandCD = 'random2',
# #             RandCM = 'random3',
# #             RandCDM = 'random4')
# #     } else {
# #         stop('Method function is wrong. Check the method function', call.=FALSE)
# #     }
# # 	if (is.null(module1.gene) | is.null(module2.gene)) {
# # 		stop('Module gene list is empty', call.=FALSE)
# # 	}
# # 	if (any(is.na(module1.gene)) | any(is.na(module1.gene))) {
# # 		stop('Please assign right node names', call.=FALSE)
# # 	}
# #     module.genelist = list(module1=module1.gene, module2=module2.gene)
# # 	g.res  <- preprocessedNetwork(network)
# #     comm.genelist <- CommunityGenelist(module.genelist, g.res, overlap_filtering = overlap_filtering)

# # 	communities = comm.genelist

# # 	x=cal.MoBCgenes_test(g.res, 
# # 					community1=comm.genelist[['module1']], 
# # 					community2=comm.genelist[['module2']],
# #                     random=random, ratio=ratio,cal.p=randomMethod, nCore=nCore)
# #     # x = subset(x, score>0)
# # 	return(x)
# # 	}







# #' Calculate centrality between two modules from MoBC result 
# #' 
# #' 
# #' @title plotDist
# #' @param MoBC.result results from CommDistFunction function
# #' @param module1.name The name of the community for which centrality is being calculated. This should be one of the communities provided as input
# #' @param module2.name The name of the community for which centrality is being calculated. This should be one of the communities provided as input
# #' @param top 
# #' @param module1.color 
# #' @param module2.color 
# #' @returns plot
# #' @export
# #' @examples
# #' plot.MoBC.genes(MoBC.result, module1.name, module2.name, 
# #'                    top=10, module1.color='lightblue1',module2.color='lightpink')



# # plotMoBC_genes <- function(MoBC.result, module1.name, module2.name, 
# #                     top=10, module1.color='lightblue1',module2.color='lightpink'){
# # 	if(!is(MoBC.result, 'MoBCresult')){
# # 		stop("input should be MoBC class", call. = FALSE)
# # 	}
# # 	communities = MoBC.result@filtered.modules

# # 	if(!all(c(module1.name, module2.name) %in% names(communities))){
# # 		stop('community name should be included in name of pre-defined community', call. = FALSE)
# # 	}

# # 	re = cal.MoBCgenes(MoBC.result@graph, 
# # 					community1=MoBC.result@filtered.modules[[module1.name]], 
# # 					community2=MoBC.result@filtered.modules[[module2.name]],
# #                     random=1, ratio=1,cal.p=FALSE)
# #     re = re[1:top,]

# #     useg = unlist(MoBC.result@filtered.modules[c(module1.name, module2.name)])
# #     useg = c(useg, re$gene)

# #     g2 <- igraph::induced_subgraph(MoBC.result@graph, useg)

# # 	layout <- igraph::layout_with_fr(g2)
    
# #     vcolor = rep('grey', length(igraph::V(g2)$name))
# #     vcolor[igraph::V(g2)$name %in% MoBC.result@filtered.modules[[module1.name]]] = module1.color
# #     vcolor[igraph::V(g2)$name %in% MoBC.result@filtered.modules[[module2.name]]] = module2.color

# #     tcolor = rep('white', length(igraph::V(g2)$name))
# #     tcolor[igraph::V(g2)$name %in% re$gene] = 'red'

# # 	plre = plot(g2, 
# # 		layout = layout, 
# # 		# mark.groups = split(V(g)$name,clv),
# # 		# vertex.label = fgid1[match(V(g)$name, fgid1$EntrezID),'gene_name'],
# # 		# vertex.label = '', #vns
# # 		vertex.color=vcolor,
# # 		vertex.frame.width=5,
# # 		vertex.frame.color=tcolor,#'white',
# # 		edge.color ='grey', #adjustcolor('black', alpha=0.6),
# # 		# vertex.size= (cln[V(cl.g2)$name]^0.5)*4,
# # 		vertex.size=10,
# # 		# vertex.label.dist=1,
# # 		# vertex.frame.color = 'grey90',
# # 		vertex.label.color='black',
# # 		# vertex.label.font=ifelse(V(g)$name %in% np.gl[[pn]], 2,1),
# # 		vertex.label.size = 0.001
# # 		# edge.width=(igraph::E(g2)$weight)*2
# # 	)
# #     return(plre)
# # }



# #' Calculate centrality between two modules from MoBC result 
# #' 
# #' 
# #' @title plotDist
# #' @param MoBC.result results from CommDistFunction function
# #' @param pval cut-off for filtering edges between communities
# #' @returns plot
# #' @export
# #' @examples
# #' plotDist(MoBC.result, pval=0.05)




# plotDist <- function(MoBC.result, pval=0.05){
# 	if(!is(MoBC.result, 'MoBCresult')){
# 		stop("input should be MoBC class", call. = FALSE)
# 	}

# 	distm = MoBC.result@MoBCresults
# 	sig.dist = subset(distm, pvalue < pval)[,1:3]
# 	sig.dist$weight = -sig.dist$z_score
# 	ntkg = igraph::graph_from_data_frame(sig.dist[,c('Module1','Module2','weight')], directed=FALSE)
# 	ntkg = igraph::simplify(ntkg, remove.multiple = TRUE, remove.loops = TRUE)

#     maxn = max(lengths(MoBC.result@filtered.modules))
#     comm.col = colorspace::sequential_hcl(length(MoBC.result@filtered.modules), "Terrain") %>% 'names<-'(names(MoBC.result@filtered.modules))
    
#     commn = lengths(MoBC.result@filtered.modules)[igraph::V(ntkg)$name]
#     sizev = (commn-min(commn))/(max(commn)-min(commn))
#     sizev = sizev*20+20

# 	layout <- igraph::layout_with_fr(ntkg)

# 	plre = plot(ntkg, 
# 		layout = layout, 
# 		# mark.groups = split(V(g)$name,clv),
# 		# vertex.label = fgid1[match(V(g)$name, fgid1$EntrezID),'gene_name'],
# 		# vertex.label = '', #vns
# 		vertex.color=comm.col[igraph::V(ntkg)$name],
# 		vertex.frame.width=0.3,
# 		vertex.frame.color='white',
# 		edge.color ="grey",#adjustcolor('black', alpha=0.6),
# 		# vertex.size= (cln[V(cl.ntkg)$name]^0.5)*4,
# 		vertex.size=sizev,
# 		# vertex.label.dist=1,
# 		# vertex.frame.color = 'grey90',
# 		vertex.label.color='black',
# 		# vertex.label.font=ifelse(V(g)$name %in% np.gl[[pn]], 2,1),
# 		vertex.label.size = 0.1,
# 		edge.width=(rank(igraph::E(ntkg)$weight))
# 	)

#     legend("bottomright", col=comm.col, pch=19, legend=names(comm.col), title='Module')
# }



