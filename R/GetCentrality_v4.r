


#--------------------------- new


get.freq <-function(g, snode, enode){
	edges = igraph::all_shortest_paths(g, snode, enode)
	edges = edges$res %>% lapply(function(xx) setdiff(names(xx), names(xx)[c(1,length(xx))]))
	# etab = edges %>% unlist
	return(edges) #!!
}

get.freq.v2 <-function(g, snode, enode){
	edges = igraph::all_shortest_paths(g, snode, enode)
	edges = edges$res %>% lapply(function(xx)  setdiff(names(xx), names(xx)[1]))
	# etab = edges %>% unlist
	return(edges) #!!
}


#' Calculate link gene centrality between two modules
#'
#' @title MoBC.genes
#' @param network A data frame representing the input biological network, where rows define interactions between genes.
#' @param module.gene.list A list of gene sets corresponding to predefined network modules.
#' @param module1 A character string specifying the first module for which MoBC-based centrality is calculated.
#'   This must match one of the module names provided in \code{module.gene.list}.
#' @param module2 A character string specifying the second module for which MoBC-based centrality is calculated.
#'   This must match one of the module names provided in \code{module.gene.list}.
#' @param randomMethod A character string specifying the method used to generate randomized networks for significance assessment.
#' @param random An integer specifying the number of randomized network generations used to estimate centrality significance.
#' @param ratio A numeric value indicating the maximum allowed overlap ratio when sampling or generating randomized modules.
#' @param nCore An integer specifying the number of CPU cores to use for parallel computation.
#' @returns A data frame containing MoBC-based centrality scores and associated statistics for genes linking the two modules.
#' @export


MoBC.genes <- function(network,
                             module.gene.list,
                             module1, module2,
                             randomMethod=c('None','RandSD','RandSDM'),
							 random = 1000,
                             nCore=1,
                             ratio = 0.1) {
    overlap_filtering=TRUE
    # cat(method,'\n')
    if (is.character(randomMethod)){
        randomMethod <- match.arg(randomMethod)
        # cat(dist.function,'\n')
        randomMethod <- switch(randomMethod,
            None = 'None',
            # RandC = 'random1',
            RandSD = 'randSD',
            # RandCM = 'random3',
            RandSDM = 'RandSDM')
    } else {
        stop('Method function is wrong. Check the method function', call.=FALSE)
    }
	if (!any(names(module.gene.list) %in% module1) | !any(names(module.gene.list) %in% module2)) {
		stop('Module name is not included in module list. Please assign right module names.', call.=FALSE)
	}
	if (any(lengths(module.gene.list)==0)) {
		stop('Please assign right node names', call.=FALSE)
	}
	if (!is.character(module1) | !is.character(module2)) {
		stop('Please assign character name as module name', call.=FALSE)
	}

	g.res  <- preprocessedNetwork(network)
    comm.genelist <- CommunityGenelist(module.gene.list, g.res, overlap_filtering = overlap_filtering)


	x=cal.MoBCgenes(g.res, comm.genelist,
					community1n=module1, 
					community2n=module2,
                    random=random, ratio=ratio,randomMethod=randomMethod, nCore=nCore)
    # x = subset(x, score>0)
	return(x)
	}





cal.MoBCgenes <- function(g, comm.genelist, community1n, community2n,random,ratio,randomMethod, nCore){

    # cat('cal.MoBCgenes')

    community1 = comm.genelist[[community1n]]
    community2 = comm.genelist[[community2n]]

    allg = igraph::V(g)$name %>% as.character()
    scorev = cal.MoBCgenes.values(g, community1, community2, allg)

	score.df = data.frame(gene=allg,score=scorev)
	score.df$node_type = 'link'
	score.df$node_type[score.df$gene %in% c(community1, community2)] = 'community genes'
	score.df = score.df %>% dplyr::arrange(-score)
    score.df = subset(score.df, score>0)

    colix = c('gene','score','node_type')
    # cat('test.\n')
    if(randomMethod!='None'){
        cat('Normalized MoBC score will be provided.\n')
        colix = c('gene','normalized_score','node_type')
        random.mat = cal.MoBC.random(g, comm.genelist, community1n, community2n,random,ratio,randomMethod,show.binning=FALSE, nCore=nCore)
        pval = sapply(score.df$gene, function(gn){
            xval = score.df[match(gn, score.df$gene),'score']
            pval = sum(random.mat[gn,]>xval)/random
        })
        nscore = sapply(score.df$gene, function(gn){
            xval = score.df[match(gn, score.df$gene),'score']
            zvalv = c(random.mat[gn,],xval)
            (xval-mean(zvalv))/sd(zvalv)
        })
        score.df$normalized_score = nscore
        score.df$pval = pval
        pv = sapply(score.df$pval, function(vv) ifelse(vv>0.5, 1-vv,vv))
        score.df$p.adj = p.adjust(pv,'BH')
        colix = c('gene','score','normalized_score','node_type','pval','p.adj')

    }
	return(score.df[,colix])
}


#--- new
cal.MoBCgenes.values <- function(g, community1, community2, allg){
    
	scorevec = rep(0, length(igraph::V(g))) %>% 'names<-'(igraph::V(g)$name)
	shortestm = igraph::distances(g, community1, community2)
	rmin  = apply(shortestm,1,function(xx) colnames(shortestm)[which(xx %in% min(xx))])
	cmin  = apply(shortestm,2,function(xx) rownames(shortestm)[which(xx %in% min(xx))])


    r.endNode = lengths(rmin) %>% sum
    c.endNode = lengths(cmin) %>% sum

    # comm1
    r.sp.genel = sapply(names(rmin), function(start.node){
		end.node = rmin[[start.node]]
        vv = sapply(end.node, function(endn){
    		etab = get.freq(g, start.node, endn)
            etab1 = etab %>% unlist %>% table
            val = etab1/length(etab)
    		return(val[allg])
        }) %>% 'rownames<-'(allg)
        vv[is.na(vv)] = 0
        return(vv)
	})
    rval = do.call(cbind, r.sp.genel)
    if(ncol(rval)!=r.endNode) stop('Not matched column number (r)', call. = FALSE)
    r.val = apply(rval,1,sum)

    # comm2
    c.sp.genel = sapply(names(cmin), function(start.node){
		end.node = cmin[[start.node]]
        vv = sapply(end.node, function(endn){
    		etab = get.freq(g, start.node, endn)
            etab1 = etab %>% unlist %>% table
            val = etab1/length(etab)
    		return(val[allg])
        }) %>% 'rownames<-'(allg)
        vv[is.na(vv)] = 0
        return(vv)
	})
    cval = do.call(cbind, c.sp.genel)
    if(ncol(cval)!=c.endNode) stop('Not matched column number (r)', call. = FALSE)
    c.val = apply(cval,1,sum)



    scorev = (r.val+c.val)/sum(r.endNode+c.endNode) #!!

	# score.df = data.frame(gene=allg,community1.score =as.numeric(r.val), community2.score=as.numeric(c.val), score=all.score)
    # scorev = score.df$community1.score + score.df$community2.score
    # names(scorev) = allg

	return(scorev)
}


# g = res2@graph
# community1 = res2@filtered.communities[[2]]
# community2 = res2@filtered.communities[[3]]
# random = 1000
# ratio = 0.1

# cal.MoBC.random(g, comm.genelist, community1n, community2n,random,ratio,cal.p,show.binning=FALSE, nCore=nCore)
# community1 --> genes in community1
# community2 --> genes in community2
# comm.genelist


# random part
cal.MoBC.random <- function(g, comm.genelist, community1n, community2n,random,ratio,randomMethod,show.binning, nCore){
    # cat("cal.MoBC.radnom - T_T")
    allg = igraph::V(g)$name %>% as.character()

    cl1g = comm.genelist[[community1n]]
    cl2g = comm.genelist[[community2n]]


    fn = make_cache_key(g, modules=comm.genelist, params=paste0(c(random,ratio),collapse='_'))
    dirn = paste0('./MoBCtmp/',fn,'/',randomMethod)

    # ㅡ make membership
    deg = igraph::degree(g)
    membership <- rep(0, length(deg))
    for(ii in 1:length(comm.genelist)){
        membership[names(deg) %in% comm.genelist[[ii]]] <- ii #membership --> numeric
    }

    # - make file name
    module_names <- names(comm.genelist)
    files_valid <- check_module_files(dirn, module_names, random)
    # cat(files_valid,'\n')
    # cat(randomMethod,'\n')


    # make directory
    if(!all(files_valid)){
        cat(paste0("You don't have optimal tmp files for random sampling - ",randomMethod,". Generating random samples will take time...  \n"))

        if(file.exists(dirn)) unlink(dirn, recursive=TRUE)
        dir.create(dirn, recursive=TRUE,showWarnings = FALSE)
        hist.bin0 = abinning_estimate_deg_bag_prob(deg, membership, kappa=ratio, ncv = random)

        fn1 = paste0(dirn,'/hist_bin.RDS')
        saveRDS(hist.bin0, file=fn1)
    }

    flag = seq(1,10)/10*random
    # if not random --> make random file first
    if(!all(files_valid) & randomMethod=='randSD'){
        hist.bin = hist.bin0$node_bag
        names(hist.bin) = 1:length(hist.bin)
        make_random = parallel::mclapply(1:random, mc.cores=nCore, function(j){
            if(any(flag %in% j)) cat('Random ',j,' is generating...\n')
            sample.save = list()
            for(ii in 1:length(comm.genelist)){
                # cat('Random samples for ',names(comm.genelist)[ii],' module\n')
                samplingN = sapply(hist.bin, function(xx) sum(names(xx) %in% comm.genelist[[ii]])) %>% 'names<-'(names(hist.bin))
                clg.random = lapply(names(samplingN), function(xn){
                    use.bg = setdiff(hist.bin[[xn]], unlist(sample.save)) 
                    sample(use.bg, samplingN[[xn]], replace=FALSE)
                }) %>% unlist %>% unique

                sample.save[[ii]] = clg.random
                dir.create(paste0(dirn,'/',names(comm.genelist)[ii]), recursive=TRUE,showWarnings = FALSE)
                write.csv(clg.random,paste0(dirn,'/',names(comm.genelist)[ii],'/rand',j,'.csv'), row.names=F)
            }
            return(NULL)
        })

    } else if(!all(files_valid) & randomMethod=='RandSDM'){
        S <- igraph::distances(g, algorithm = "unweighted")

        comm.distance.list = parallel::mclapply(1:random, mc.cores=nCore, function(j){
            if(any(flag %in% j)) cat('Random ',j,' is generating...\n')
            rsamplel = modularity_sampling_multi(hist.bin0, deg, membership, S)     
            for(ii in 1:length(rsamplel)){
                # cat('Random samples for ',names(comm.genelist)[ii],' module\n')

                clg.random = rsamplel[[ii]]
                dir.create(paste0(dirn,'/',names(comm.genelist)[ii]), recursive=TRUE,showWarnings = FALSE)
                write.csv(clg.random,paste0(dirn,'/',names(comm.genelist)[ii],'/rand',j,'.csv'), row.names=F)
            }
            return(NULL)
        })

    } else if(randomMethod!='RandSD' & randomMethod!='RandSDM'){
        stop('Please enter the right method for randomization.', call.=FALSE)
    }


    # check again
    files_valid1 <- check_module_files(dirn, module_names, random)

    if(all(files_valid1)){### edit needed

        cat(paste0('You have tmp files for random sampling based on ',randomMethod,": ",dirn,". We will use these files.\n"))

        comm.distance.list = parallel::mclapply(1:random,mc.cores=nCore, function(j){
            if(any(flag %in% j)) cat('We load ',j,' random.\n')
            # cat('read ',j,'\n')
            rs1 = read.csv(paste0(dirn,'/',community1n,'/rand',j,'.csv'))[,1]
            rs2 = read.csv(paste0(dirn,'/',community2n,'/rand',j,'.csv'))[,1]
            # cat('cal ',j,'\n')
            
            comm.distance = cal.MoBCgenes.values(g, rs1,rs2, allg) 
            # cat('done ',j,'\n')
            return(comm.distance)
        })
        comm.distance.list = do.call(cbind, comm.distance.list) %>% 'rownames<-'(allg)
        return(comm.distance.list)
    }

    cat('fin comm.dist\n')
    # start.time-end.time
    return(comm.distance.list)
}








cal.FCgene <- function(g, community1, community2){

    gene.ix = igraph::V(g)$name

    # re = igraph::all_shortest_paths(g, community1[1:3],community2[1:4])
    # lapply(re$res, function(xx) names(xx)[1] ) %>% unlist %>% unique
    # lapply(re$res, function(xx) names(xx)[length(xx)] ) %>% unlist %>% unique
    
    re = sapply(community1, function(g1){
            
        edges = igraph::all_shortest_paths(g, g1, community2)
        end.ix = edges$res %>% sapply(function(xx) names(xx)[length(xx)])
        end.ix = split(1:length(end.ix), end.ix)
        resl = lapply(end.ix, function(ixix){
            intg = edges$res[ixix] %>% lapply(function(xx) setdiff(names(xx), names(xx)[c(1,length(xx))]))
            intg.tab = unlist(intg) %>% table
            intg.tab = intg.tab/length(intg)
            intg.tab[gene.ix] %>% 'names<-'(gene.ix)
        })
        res = do.call(rbind, resl) %>% as.matrix
        res[is.na(res)] = 0
        resv = apply(res,2,sum)
    })
    re1 = apply(re,1,sum)
    re1 = re1/length(community1)/length(community2)

	score.df = data.frame(gene=gene.ix,score=re1[gene.ix])
    score.df$node_type = 'link'
	score.df$node_type[score.df$gene %in% c(community1, community2)] = 'module genes'
	score.df = score.df %>% dplyr::arrange(-score)

	return(score.df)
}



#' Plot shortest paths through a link gene between two modules
#' 
#' @param g Network data frame
#' @param module1 First module genes
#' @param module2 Second module genes  
#' @param linkgene Target link gene to visualize paths through
#' @param col1 Color for first module genes
#' @param col2 Color for second module genes
#' @param link.col Color for link genes
#' @export
#'

link.gene.path<-function(g, x, y, linkgene,col1,col2,link.col) {

    g.graph = preprocessedNetwork(g)
    shortestm = igraph::distances(g.graph, x, y)
    rmin  = apply(shortestm,1,function(xx) colnames(shortestm)[which(xx %in% min(xx))])
    cmin  = apply(shortestm,2,function(xx) rownames(shortestm)[which(xx %in% min(xx))])
    shorteste = lapply(names(rmin), function(snode) {
        enode = rmin[[snode]]
        edges = igraph::all_shortest_paths(g.graph, snode, enode)
        edges$res
    })
    shorteste.2 = lapply(names(cmin), function(snode) {
        enode = cmin[[snode]]
        edges = igraph::all_shortest_paths(g.graph, snode, enode)
        edges$res
    })

    path.res = lapply(linkgene, function(gene) {
        # cat(gene,'\n')
        vv = lapply(shorteste, function(x1) {
            x1 = lapply(x1, names)
            tfv = sapply(x1, function(path) any(path %in% gene))
            if(!any(tfv)) return(NULL)
            x1[tfv]
            })
        vv1 = vv[!sapply(vv, is.null)]
        vv2 = list()
        for(ii in 1:length(vv1)){
            for(jj in 1:length(vv1[[ii]])){
                vv2 = c(vv2, vv1[[ii]][jj])
            }
        }
        
        return(vv2)
    }) %>% 'names<-'(linkgene)

    path.res.2 = lapply(linkgene, function(gene) {
        # cat(gene,'\n')
        vv = lapply(shorteste.2, function(x1) {
            x1 = lapply(x1, names)
            tfv = sapply(x1, function(path) any(path %in% gene))
            if(!any(tfv)) return(NULL)
            x1[tfv]
            })
        vv1 = vv[!sapply(vv, is.null)]
        if(length(vv1)==0) return(NULL)
        vv2 = list()
        for(ii in 1:length(vv1)){
            for(jj in 1:length(vv1[[ii]])){
                vv2 = c(vv2, vv1[[ii]][jj])
            }
        }
        
        return(vv2)
    }) %>% 'names<-'(linkgene)

    path.res = list(path1 = path.res, path2 =path.res.2)

    # plot function
    ppi.df = g

    tdf.a =  lapply(path.res[['path1']][[linkgene]], function(vv){
        tdf2 = lapply(2:length(vv), function(ii){
            data.frame(from=vv[ii-1],to=vv[ii])
        }) %>% bind_rows %>% as.data.frame
    }) %>% bind_rows %>% as.data.frame
    tdf.b =  lapply(path.res[['path2']][[linkgene]], function(vv){
        tdf2 = lapply(2:length(vv), function(ii){
            data.frame(from=vv[ii-1],to=vv[ii])
        }) %>% bind_rows %>% as.data.frame
    }) %>% bind_rows %>% as.data.frame


    tt = rbind(tdf.a, tdf.b)

    #-- key leftg

    ta =  lapply(path.res[['path1']][[linkgene]], function(vv){
        vv[1:(which(vv==linkgene)-1)]
        }) %>% unlist %>% unique

    tb =  lapply(path.res[['path2']][[linkgene]], function(vv){
        vv[(which(vv==linkgene)+1):length(vv)]
        }) %>% unlist %>% unique

    intersect(ta, tb)
    t.all = c(ta, tb)


    ntkg = igraph::graph_from_data_frame(tt, directed=TRUE)
    ntkg = igraph::simplify(ntkg, remove.multiple = TRUE, remove.loops = TRUE)


    shortestm = igraph::distances(ntkg, igraph::V(ntkg)$name, linkgene) %>% as.data.frame %>% 'colnames<-'(c('x'))
    shortestm = shortestm*2
    shortestm[which(rownames(shortestm) %in% t.all),1] = -shortestm[which(rownames(shortestm) %in% t.all),1]
    shortestm = shortestm %>% arrange(x)
    shl = split(shortestm, shortestm[,1]) %>% lapply(function(vv){

        if(nrow(vv)==1){
            val=0
        } else{
            val = seq(-2,2,length.out=nrow(vv))
        }
        vv$y = val
        return(vv)

    }) %>% bind_rows %>% as.matrix
    shl = shl[igraph::V(ntkg)$name,]

    # shortest path


    # plot(ntkg, layout=layout_on_grid)

    vv = igraph::V(ntkg)$name
    colv = rep('grey',length(vv))
    colv[vv %in% x] = col1
    colv[vv %in% y] = col2
    colv[vv %in% linkgene] = link.col


    if(length(unique(tt[,1]))>10){

        ids = shl %>% rownames()
        ids[shl[,1]<0] = ''
        shl[which(shl[,1]==4),1]=3

        # xx = norm_coords(shl, xmin = -1, xmax) --> test


        plot(ntkg, 
            layout = shl, 
            rescale=TRUE,
            vertex.size=ifelse(ids=='',4,12),
            edge.arrow.size=0.5,
            vertex.label = ids, #felse(shl[,1]<0,'',),#labels,
            vertex.color=colv,#'orange',
            vertex.frame.width=2,#ifelse(tfv,2,1),
            vertex.frame.color='black',#border.col, #rep('white', length(tcolor)),#tcolor,#'white',
            edge.color ='black', #adjustcolor('black', alpha=0.6),
            vertex.label.family='Arial',
            vertex.label.color='black',
            vertex.label.cex = 0.8,
            edge.width=1
        )

        xx = norm_coords(shl)
        textl = xx[ids=='',]
        textl[,1]=textl[,1]-0.15

        text(textl[,1],textl[,2], rownames(textl), col='black',cex=0.9)


    } else{

        plot(ntkg, 
            layout = shl, 
            rescale=TRUE,
            vertex.color=colv,#'orange',
            vertex.frame.width=2,#ifelse(tfv,2,1),
            vertex.frame.color='black',#border.col, #rep('white', length(tcolor)),#tcolor,#'white',
            edge.color ='black', #adjustcolor('black', alpha=0.6),
            vertex.label.family='Arial',
            vertex.label.color='black',
            vertex.label.cex = 0.8,
            edge.width=1
        )

    }

}
