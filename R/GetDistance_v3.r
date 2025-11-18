#' module distance function
#' 
#' Calculate the distance and z_score of the communities
#' use.f can be choosed in the DistFunction.r
#' ex) get.shortest.dist, get.kernel.dist, get.centre.dist, get.separation.dist, get.closest.dist
#' As a default, get.closest.dist function is used to measure the distance between communities
#' User also can make the dist function and use it for calculating module distance
#' This function is made to know the z-score of a measured distance from distances of degree-preserved random networks
#' 
#' @param network network dataframe
#' @param module.genelist module genes list
#' @param random the number of random network formation and distance calculation
#' @param overlap_filtering overlap genes filtering (TRUE/FALSE)
#' @param method distance measuring method
#' @return distance results, module.genelist, network
#' @export 
#' 


# 1.RandC, 2. RandCD, 3. RandCM, 4. RandCDM
CommuinityDistance <- function(network,
							 module.genelist,
                             randomMethod=c('None','RandSD','RandSDM'),
							 random = 1000,
                             ratio = 0.1,
                             nCore=1,
                             method = c('closest', 'shortest', 'kernel', 'centre', 'separation')) {
    # cat(method,'\n')
    overlap_filtering=TRUE
    if (is.character(method)){
        dist.function <- match.arg(method)
        # cat(dist.function,'\n')
        dist.function <- switch(dist.function,
            closest = get.closest.dist,
            shortest = get.shortest.dist,
            kernel = get.kernel.dist,
            centre = get.centre.dist,
            separation = get.separation.dist)
    } else if(is.function(method)){
        dist.function = method
    } else {
        stop('Method function is wrong. Check the method function', call.=FALSE)
    }
	if (is.null(names(module.genelist)) | any(is.na(names(module.genelist)))) {
		stop('Please assign names to the module list.', call.=FALSE)
	}
    if (is.character(randomMethod)){
        randomMethod <- match.arg(randomMethod)
        # cat(dist.function,'\n')
        randomMethod <- switch(randomMethod,
            none = FALSE,
            RandSD = 'RandSD',
            RandSDM = 'RandSDM')
    } else {
        stop('Method function is wrong. Check the method function', call.=FALSE)
    }
	g.res  <- preprocessedNetwork(network)
    comm.genelist <- CommunityGenelist(module.genelist, g.res, overlap_filtering = overlap_filtering)
	distm <- igraph::distances(g.res, igraph::V(g.res), igraph::V(g.res))

    cat('Dist matrix :', dim(distm)[1],'X',dim(distm)[2], 'is made','\n')
    cat('Random distance measuring is going to be processed by', random, 'times','\n')

    # 클러스터 초기화
    # 기존 데이터 있는지 확인



    fn = make_cache_key(g.res, modules=comm.genelist, params=paste0(c(random,ratio),collapse='_'))
    dirn = paste0('./MoBCtmp/',fn,'/',randomMethod)

    # ㅡ make membership
    deg = igraph::degree(g.res)
    membership <- rep(0, length(deg))
    for(ii in 1:length(comm.genelist)){
        membership[names(deg) %in% comm.genelist[[ii]]] <- ii #membership --> numeric
    }

    # - make file name
    module_names <- names(comm.genelist)
    files_valid <- check_module_files(dirn, module_names, random)
    cat(files_valid,'\n')
    cat(randomMethod,'\n')

    # make directory
    if(!all(files_valid)){
        cat(paste0("You don't have optimal tmp files for random sampling - ",randomMethod,". Generating random samples will take time...  \n"))

        if(file.exists(dirn)) unlink(dirn, recursive=TRUE)
        dir.create(dirn, recursive=TRUE,showWarnings = FALSE)
        hist.bin0 = abinning_estimate_deg_bag_prob(deg, membership, kappa=ratio, ncv = random)

        fn1 = paste0(dirn,'/hist_bin.RDS')
        saveRDS(hist.bin0, file=fn1)
        cat('Generated A bin.\n')
    }

    flag = seq(1,10)/10*random
    # if not random --> make random file first
    if(!all(files_valid) & randomMethod=='RandSD'){
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
                    base::sample(use.bg, samplingN[[xn]], replace=FALSE)
                }) %>% unlist %>% unique

                sample.save[[ii]] = clg.random
                dir.create(paste0(dirn,'/',names(comm.genelist)[ii]), recursive=TRUE,showWarnings = FALSE)
                write.csv(clg.random,paste0(dirn,'/',names(comm.genelist)[ii],'/rand',j,'.csv'), row.names=F)
            }
            return(NULL)
        })

    } else if(!all(files_valid) & randomMethod=='RandSDM'){
        cat('hi?ㅇㅅㅇ\n')
        S <- igraph::distances(g.res, algorithm = "unweighted")
        cat('hi?ㅎㅅㅎ\n')
        comm.distance.list = parallel::mclapply(1:random, mc.cores=nCore, function(j){
        # comm.distance.list = lapply(1:random,  function(j){
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

        cat(paste0('You have tmp files for random sampling - ',randomMethod,". We will use these files.\n"))
        
        binl = list()
        results = lapply(1:(length(comm.genelist)-1), function(m){
            dist.rel =lapply((m+1):length(comm.genelist), function(n){ #dist.rel =
                cat('Distance measuring :','module',names(comm.genelist)[m],' - ','module', names(comm.genelist)[n],'\n')

                cl1g = comm.genelist[[m]]
                cl2g = comm.genelist[[n]]

                comm.distance.list = parallel::mclapply(1:random,mc.cores=nCore, function(j){
                    if(any(flag %in% j)) cat('We load ',j,' random.\n')
                    rs1 = read.csv(paste0(dirn,'/',names(comm.genelist)[m],'/rand',j,'.csv'))[,1]
                    rs2 = read.csv(paste0(dirn,'/',names(comm.genelist)[n],'/rand',j,'.csv'))[,1]
                    
                    comm.distance = dist.function(distm, rs1,rs2) 
                    return(comm.distance)
                }) %>% unlist

                xval = dist.function(distm, cl1g, cl2g) 
                zval = (xval-mean(comm.distance.list))/sd(comm.distance.list)
                pval = sum(comm.distance.list<xval)/random

                df1 = data.frame(Module1=names(comm.genelist)[m], Module2=names(comm.genelist)[n], z_score=zval, distance_score=xval, pvalue=pval)
                pv = sapply(df1$pvalue, function(vv) ifelse(vv>0.5, 1-vv,vv))
                df1$p.adj = p.adjust(pv,'BH')
                return(df1)
            })
            dist.rel = do.call(rbind, dist.rel)
        })
        results = do.call(rbind, results)

        x = list(
            distance = results,
            filtered.modules = comm.genelist,
            graph = g.res)

    }

    fn1 = paste0(dirn,'/hist_bin.RDS')
    hist.bin0 = readRDS(fn1)
        

    x$bag = hist.bin0


	return(x)
    }





make_cache_key <- function(g, modules, params) {
    mod_sorted = lapply(modules, function(g) sort(unique(g)))
    mod_sorted = mod_sorted[order(names(mod_sorted))]
    digest::digest(list(
        nodes = sort(igraph::V(g)$name),
        edges = nrow(igraph::as_data_frame(g, what = "edges")),
        modules = mod_sorted,
        params = params  # 예: list(n_random = 1000, method = "degree")
    ))
}


check_module_files <- function(base_dir, module_names, random_count) {
    all_valid = rep(TRUE, length(module_names)) %>% 'names<-'(module_names)
    
    for(module_name in module_names) {
        module_dir <- file.path(base_dir, module_name)
        
        # Directory 존재 확인
        if(!dir.exists(module_dir)) {
            cat("Directory not found:", module_dir, "\n")
            all_valid[module_name] <- FALSE
            next
        }
        
        # 파일 개수 확인
        files <- list.files(module_dir, pattern = "^rand[0-9]+\\.csv$")
        if(length(files) != random_count) {
            cat("File count mismatch for", module_name, ":", length(files), "!=", random_count, "\n")
            all_valid[module_name] <- FALSE
        }
    }
    
    return(all_valid)
}
