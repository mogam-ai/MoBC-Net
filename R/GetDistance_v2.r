# #' module distance function
# #' 
# #' Calculate the distance and z_score of the communities
# #' use.f can be choosed in the DistFunction.r
# #' ex) get.shortest.dist, get.kernel.dist, get.centre.dist, get.separation.dist, get.closest.dist
# #' As a default, get.closest.dist function is used to measure the distance between communities
# #' User also can make the dist function and use it for calculating module distance
# #' This function is made to know the z-score of a measured distance from distances of degree-preserved random networks
# #' 
# #' @param network network dataframe
# #' @param module.genelist module genes list
# #' @param random the number of random network formation and distance calculation
# #' @param overlap_filtering overlap genes filtering (TRUE/FALSE)
# #' @param method distance measuring method
# #' @return distance results, module.genelist, network
# #' @export 
# #' 


# # 1.RandC, 2. RandCD, 3. RandCM, 4. RandCDM
# CommuinityDistance <- function(network,
# 							 module.genelist,
#                              randomMethod=c('None','RandC','RandCD','RandCM','RandCDM'),
# 							 random = 1000,
#                              ratio = 0.1,
#                              nCore=1,
#                              method = c('closest', 'shortest', 'kernel', 'centre', 'separation')) {
#     # cat(method,'\n')
#     overlap_filtering=TRUE
#     if (is.character(method)){
#         dist.function <- match.arg(method)
#         # cat(dist.function,'\n')
#         dist.function <- switch(dist.function,
#             closest = get.closest.dist,
#             shortest = get.shortest.dist,
#             kernel = get.kernel.dist,
#             centre = get.centre.dist,
#             separation = get.separation.dist)
#     } else if(is.function(method)){
#         dist.function = method
#     } else {
#         stop('Method function is wrong. Check the method function', call.=FALSE)
#     }
# 	if (is.null(names(module.genelist)) | any(is.na(names(module.genelist)))) {
# 		stop('Please assign names to the module list.', call.=FALSE)
# 	}
#     if (is.character(randomMethod)){
#         randomMethod <- match.arg(randomMethod)
#         # cat(dist.function,'\n')
#         randomMethod <- switch(randomMethod,
#             none = FALSE,
#             RandC = 'random1',
#             RandCD = 'random2',
#             RandCM = 'random3',
#             RandCDM = 'random4')
#     } else {
#         stop('Method function is wrong. Check the method function', call.=FALSE)
#     }
# 	g.res  <- preprocessedNetwork(network)
#     comm.genelist <- CommunityGenelist(module.genelist, g.res, overlap_filtering = overlap_filtering)
# 	distm <- igraph::distances(g.res, igraph::V(g.res), igraph::V(g.res))

#     cat('Dist matrix :', dim(distm)[1],'X',dim(distm)[2], 'is made','\n')
#     cat('Random distance measuring is going to be processed by', random, 'times','\n')

#     # 클러스터 초기화


#     binl = list()
#     # m=1, n=2
# 	results = lapply(1:(length(comm.genelist)-1), function(m){
# 		 dist.rel =lapply((m+1):length(comm.genelist), function(n){ #dist.rel =
# 			cat('Distance measuring :','module',names(comm.genelist)[m],' - ','module', names(comm.genelist)[n],'\n')
#             idv = paste0(m,'_',n)
#             cl1g = comm.genelist[[m]]
# 			cl2g = comm.genelist[[n]]

#             deg = igraph::degree(g.res)
#             membership <- rep(0, length(deg))
#             membership[names(deg) %in% cl1g] <- 1
#             membership[names(deg) %in% cl2g] <- 2 

#             kk = c(names(comm.genelist)[m],names(comm.genelist)[n]) %>% sort
#             options=paste0(c(kk,random,ratio),collapse='_')

#             dirn1 = paste0('./MoBCtmp/Dist/',options,'/',names(comm.genelist)[m],'_',randomMethod,'_',random)
#             dirn2 = paste0('./MoBCtmp/Dist/',options,'/',names(comm.genelist)[n],'_',randomMethod,'_',random)
#             dir.flag = dir.exists(dirn1) & dir.exists(dirn2)
#             l.flag = ifelse(length(list.files(dirn1))==random & length(list.files(dirn2))==random,TRUE,FALSE)

#             if(dir.flag & l.flag){
#                 tag1 = all(sapply(1:random, function(j) read.csv(paste0(dirn1,'/rand',j,'.csv'))[,1] %>% length)==length(cl1g))
#                 tag2 = all(sapply(1:random, function(j) read.csv(paste0(dirn2,'/rand',j,'.csv'))[,1] %>% length)==length(cl2g))
#                 length.flag = tag1 & tag2
#             } else{
#                 length.flag = FALSE
#             }

#             if(dir.flag & l.flag & length.flag){

#                 cat(paste0('You have tmp files for random sampling - ',randomMethod,". We will use these files.\n"))
#                 comm.distance.list = parallel::mclapply(1:random,mc.cores=nCore, function(j){
#                     rs1 = read.csv(paste0(dirn1,'/rand',j,'.csv'))[,1]
#                     rs2 = read.csv(paste0(dirn2,'/rand',j,'.csv'))[,1]
#                     comm.distance = dist.function(distm, rs1,rs2) 
#                     return(comm.distance)
#                 }) %>% unlist

#                 xval = dist.function(distm, cl1g, cl2g) 
#                 zval = (xval-mean(comm.distance.list))/sd(comm.distance.list)
#                 pval = sum(comm.distance.list<xval)/random

#                 df1 = data.frame(Module1=names(comm.genelist)[m], Module2=names(comm.genelist)[n], z_score=zval, distance_score=xval, pvalue=pval)
#                 pv = sapply(df1$pvalue, function(vv) ifelse(vv>0.5, 1-vv,vv))
#                 df1$p.adj = p.adjust(pv,'BH')
#                 return(df1)
#             }

#             cat(paste0("You don't have tmp files for random sampling - ",randomMethod,". Generating random samples will take time...  \n"))
#             dir.create(dirn1, recursive=TRUE,showWarnings = FALSE)
#             dir.create(dirn2, recursive=TRUE,showWarnings = FALSE)
            
#             if(randomMethod=='random1'){
                
#                 comm.distance.list = parallel::mclapply(1:random, mc.cores=nCore, function(j){
#                     rsamplel = simple_sampling(deg, membership,random)
#                     write.csv(rsamplel[[1]],paste0(dirn1,'/rand',j,'.csv'), row.names=F)
#                     write.csv(rsamplel[[2]],paste0(dirn2,'/rand',j,'.csv'), row.names=F)
#                     comm.distance = dist.function(distm, rsamplel[[1]], rsamplel[[2]]) 
#                     return(comm.distance)
#                 }) %>% unlist

#             } else if(randomMethod=='random2'){
#                 hist.bin0 = estimate_deg_bag(deg, membership, ratiov=ratio, ncv=random)
#                 re = hist.bin0
#                 hist.bin = hist.bin0$node_bag
#                 names(hist.bin) = 1:length(hist.bin)
#                 comm.distance.list = parallel::mclapply(1:random, mc.cores=nCore, function(j){
                    
#                     samplingN = sapply(hist.bin, function(xx) sum(names(xx) %in% cl1g)) %>% 'names<-'(names(hist.bin))
#                     cl1g.random = lapply(names(samplingN), function(xn){
#                         use.bg = (hist.bin[[xn]])
#                         # set.seed(m+j)
#                         sample(use.bg, samplingN[[xn]], replace=FALSE)
#                     }) %>% unlist %>% unique

#                     samplingN = sapply(hist.bin, function(xx) sum(names(xx) %in% cl2g)) %>% 'names<-'(names(hist.bin))
#                     cl2g.random = lapply(names(samplingN), function(xn){
#                         use.bg = setdiff(hist.bin[[xn]], c(cl1g.random)) 
#                         # set.seed(n+j)
#                         sample(use.bg, samplingN[[xn]], replace=FALSE)
#                     }) %>% unlist %>% unique
#                     write.csv(cl1g.random,paste0(dirn1,'/rand',j,'.csv'), row.names=F)
#                     write.csv(cl2g.random,paste0(dirn2,'/rand',j,'.csv'), row.names=F)
#                     comm.distance = dist.function(distm, cl1g.random, cl2g.random) 
#                     return(comm.distance)
#                 }) %>% unlist

#             } else if(randomMethod=='random4'){
#                 S <- igraph::distances(g.res, algorithm = "unweighted")
#                 # pb <- progress::progress_bar$new(total = random)
#                 re = estimate_deg_bag(deg, membership, ratiov=ratio, ncv = random)

#                 comm.distance.list = parallel::mclapply(1:random, mc.cores=nCore, function(j){
#                     rsamplel = modularity_sampling(re, deg, membership, S)     
#                     write.csv(rsamplel[[1]],paste0(dirn1,'/rand',j,'.csv'), row.names=F)
#                     write.csv(rsamplel[[2]],paste0(dirn2,'/rand',j,'.csv'), row.names=F)           
#                     comm.distance = dist.function(distm, rsamplel[[1]], rsamplel[[2]]) 
#                     return(comm.distance)
#                 }) %>% unlist

#                 # comm.distance.list = c()
#                 # for(j in 1:random){
#                 #     rsamplel = modularity_sampling(re, deg, membership, S)     
#                 #     write.csv(rsamplel[[1]],paste0(dirn1,'/rand',j,'.csv'), row.names=F)
#                 #     write.csv(rsamplel[[2]],paste0(dirn2,'/rand',j,'.csv'), row.names=F)           
#                 #     comm.distance = dist.function(distm, rsamplel[[1]], rsamplel[[2]]) 
#                 #     comm.distance.list = c(comm.distance.list,comm.distance)
#                 #     pb$tick()
#                 # }
#                 # comm.distance.list = sapply(1:random, function(j){
#                 #     rsamplel = modularity_sampling(deg, membership, S, random)                
#                 #     comm.distance = dist.function(distm, rsamplel[[1]], rsamplel[[2]]) 
#                 #     return(comm.distance)
#                 # }) %>% unlist
#             } else if(randomMethod=='random3'){
#                 cat('Select randomization3 method')
#                 S <- igraph::distances(g.res, algorithm = "unweighted")                                
#                 m1 <- length(cl1g)
#                 n1 <- length(cl2g)

#                 tnodes <- seq_len(length(igraph::V(g.res)$name)) 

#                 # clustermn 함수를 사용해 2개 그룹으로 나눈다고 가정 (사용자 정의 필요)
#                 # 예) clusterAssignment <- clustermn(S[sind, sind], m, n)
#                 #
#                 # 여기서는 간단히, 유사 함수가 있다고 가정하고 아래처럼 호출합니다.
#                 #
#                 # 사용자는 clustermn 함수를 별도로 정의해서
#                 # "m개와 n개로 구성되면서, 같은 모듈 내 노드끼리는 평균 거리가 짧도록"  
#                 # 하는 방식을 구현해야 합니다.

                
#                 # pb <- progress::progress_bar$new(total = random)

#                 comm.distance.list = parallel::mclapply(1:random,mc.cores=nCore, function(j){
#                     sind <- sample(tnodes, m1 + n1, replace = FALSE)
#                     clusterAssignment <- clustermn(S[sind, sind], m1, n1)
#                     Urand <- sind[clusterAssignment == 1]
#                     Wrand <- sind[clusterAssignment == 2]         
#                     write.csv(Urand,paste0(dirn1,'/rand',j,'.csv'), row.names=F)
#                     write.csv(Wrand,paste0(dirn2,'/rand',j,'.csv'), row.names=F)
#                     comm.distance = dist.function(distm, Urand, Wrand) 
#                     return(comm.distance)
#                 }) %>% unlist

#                 # comm.distance.list = c()
#                 # for(j in 1:random){
#                 #     sind <- sample(tnodes, m1 + n1, replace = FALSE)
#                 #     clusterAssignment <- clustermn(S[sind, sind], m1, n1)
#                 #     Urand <- sind[clusterAssignment == 1]
#                 #     Wrand <- sind[clusterAssignment == 2]         
#                 #     write.csv(Urand,paste0(dirn1,'/rand',j,'.csv'), row.names=F)
#                 #     write.csv(Wrand,paste0(dirn2,'/rand',j,'.csv'), row.names=F)
#                 #     comm.distance = dist.function(distm, Urand, Wrand) 
#                 #     comm.distance.list = c(comm.distance.list,comm.distance)
#                 #     pb$tick()
#                 # }

#             } else {
#             	stop('Please enter the right method for randomization.', call.=FALSE)
#             }

# 			xval = dist.function(distm, cl1g, cl2g) 
# 			zval = (xval-mean(comm.distance.list))/sd(comm.distance.list)
# 			pval = sum(comm.distance.list<xval)/random

# 			df1 = data.frame(Module1=names(comm.genelist)[m], Module2=names(comm.genelist)[n], z_score=zval, distance_score=xval, pvalue=pval)
#             pv = sapply(df1$pvalue, function(vv) ifelse(vv>0.5, 1-vv,vv))
#             df1$p.adj = p.adjust(pv,'BH')
#             df1$pvalue = pv
# 			# cat('-end\n')
#             fn1 = paste0('./MoBCtmp/Dist/',options,'/',idv,'_',randomMethod,'.RDS')
#             saveRDS(re, file=fn1)
# 			return(df1)
# 			})
# 		dist.rel = do.call(rbind, dist.rel)
# 		}) #%>% dplyr::bind_rows %>% as.data.frame
# 	results = do.call(rbind, results)
#     x = list(
#         distance = results,
#         filtered.modules = comm.genelist,
#         graph = g.res)

#     for(m in 1:(length(comm.genelist)-1)){
# 		for(n in (m+1):length(comm.genelist)){
#             idv = paste0(m,'_',n)
#             kk = c(names(comm.genelist)[m],names(comm.genelist)[n]) %>% sort
#             options=paste0(c(kk,random,ratio),collapse='_')
#             fn1 = paste0('./MoBCtmp/Dist/',options,'/',idv,'_',randomMethod,'.RDS')
#             re = readRDS(fn1)
#             binl[[idv]] = re
#         }
#     }

#     x$bag = binl

#     xre = results
#     xre$weight = xre$z_score -min(xre$z_score)
 
#     xre.ntkg = igraph::graph_from_data_frame(xre[,c('Module1','Module2','weight')], directed=FALSE)
#     # xre.ntkg1 = mst(xre.ntkg)
#     # x$mst.net = xre.ntkg1

#    # xre$weight = -max(x$weight)


# 	# x <- new("MoBCresult",
#     #     MoBCresults = results,
#     #     filtered.modules = comm.genelist,
#     #     graph = g.res)
# 	return(x)
# 	}

