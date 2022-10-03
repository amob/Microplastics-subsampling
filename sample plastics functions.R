library(coda)
library(iNEXT)#coverage and asymptotic richness estimators

range01=function(x){
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
}

#we have data with a known number of particles in different categories; where all particles where tested for polymer type

#we want to test what a perfect sampler would find following sampling protocol for proportion of polymer types in dataset
# we want to generate a curve for # of particles sampled from each category

sampleplastics <- function(plastics,numpercat,scorabletypecolnum,truetypecolnum, random=F) {
	 #dataset of scored plastics, the number per category sampling (or random must =T and this is the total subsample size), 
	 #and the column that defines the category by which subsampling will be structured
	colnames(plastics)[truetypecolnum] <- "truetype"
	colnames(plastics)[scorabletypecolnum] <- "scorabletype"  #rename column to be used for subsampling
	unfilteredB<- plastics[sample(1:nrow(plastics),nrow(plastics),replace=F), ]  #randomize order
	if(random==FALSE){
		partskeptB<- list()
		for(i in 1:length(unique(unfilteredB$scorabletype))){  #for each category 
			partsincat <- which(unfilteredB$scorabletype == unique(unfilteredB$scorabletype)[i]) #which particles are in that category in the pre-randomized data
			partskeptB[[i]] <- partsincat[1:numpercat] #subsample row numbers
		}
		partskeptB.1 <- unlist(partskeptB)[!is.na(unlist(partskeptB))] #unlist
		filteredB <- unfilteredB[partskeptB.1,] #extract rows from the dataset
	} else { #if random -- i.e. no category to assist randomizing
			filteredB <- unfilteredB[1:numpercat,] #just sample the first number of rows equal to the sample size
		}
	
	if(random==FALSE){ 
		summarizeplastics <- table(paste(plastics$scorabletype)) / nrow(plastics)
		sortUscoretype <- sort(unique(paste(plastics$scorabletype)))
		proptable <- matrix(NA, nrow=length(sortUscoretype),ncol=length(sort(unique(plastics$truetype))))
		for(category in (1:length(sortUscoretype)) ){
			 catsample <-  filteredB[ paste(filteredB$scorabletype) == as.character(  sortUscoretype[category]  ) ,] 
			 proptable[category,] <- sapply(  sort(unique(plastics$truetype)),  function(z)  sum(as.character(catsample$truetype) == as.character(z))  )  / nrow(catsample)			
			}
		sampledrates <- sapply(1:ncol(proptable), function(type) sum(proptable[,type]*summarizeplastics) ) 
	} else{
		summarizeB <- sapply(1:length(sort(unique(plastics$truetype))), function(type)  
			length(which(as.character(filteredB$truetype)==as.character(sort(unique(plastics$truetype))[type]) ) ))
		sampledrates <- summarizeB/sum(summarizeB)
	}
	
	summarizeAll <- table(plastics$truetype) 		
	names(sampledrates) <- paste( sort(unique(plastics$truetype)),"Sample",sep="_") 
	 truerates <- summarizeAll/sum(summarizeAll)
	names(truerates) <- paste( sort(unique(plastics$truetype)),"Input",sep="_") 
	sumerrrates <- sum (abs(  sampledrates -  truerates  )  )/2
	avgerrrates <- mean (abs( sampledrates -  truerates   )  )

	S.obs		<- DataInfo(as.vector(table(filteredB$truetype)),datatype="abundance")$S.obs 
	AsyChaoRaw <- ChaoRichness(as.vector(table(filteredB$truetype)),datatype="abundance")$Estimator 
	Coverage <- DataInfo(as.vector(table(filteredB$truetype)),datatype="abundance")$SC 

	return(c(  sumerrrates = sumerrrates, avgerrrates=avgerrrates,
				nobs = nrow(filteredB[!is.na(filteredB$truetype),]),  
				S.obs=S.obs,Coverage=Coverage, AsyChaoRaw = AsyChaoRaw)) 
}




####this next function is nearly identical to the previous, in the first part
####### but now we have multiple samples from it (e.g. individual birds, individual raingardens) -- 
############ -- and we want to retain separate information across our samples

StructuredData_depthper <- function(plastics,numpercat,sampleidcolnum, scorabletypecolnum,truetypecolnum,random=F) { 
	 #dataset of scored plastics, the number per category sampling (or random must =T and this is the total subsample size), 
	 #and the column that defines the category by which subsampling will be structured
	colnames(plastics)[truetypecolnum] <- "truetype"
	colnames(plastics)[scorabletypecolnum] <- "scorabletype"  #rename column to be used for subsampling
	colnames(plastics)[sampleidcolnum] <- "sampleid"  #rename column to be used for subsampling
	listsamps <- lapply(unique(plastics$sampleid), function(id) plastics[plastics$sampleid == id,]   ) #make individual samples into separate dataframes in a list
	unfilteredB<- lapply(listsamps, function(id) id[sample(1:nrow(id),nrow(id),replace=F), ] ) #randomize order, now for each in list
	
	
	filteredB.1<- list() #an empty list 
	if(random==FALSE){
		partskeptB<- vector(mode="list",length = length(listsamps)) #an empty list of length samples
			partskeptB<-lapply(partskeptB, function(x) list()) #each element is now itself an initialized list of null length
		partskeptB.1<- vector(mode="list",length = length(listsamps)) #an empty list of length samples
			partskeptB.1 <- lapply(partskeptB.1, function(x) list()) #each element is now itself an initialized list of null length
		for(id in 1:length(unfilteredB)){
			for(i in 1:length(unique(unfilteredB[[id]]$scorabletype))){  #for each category 
				partsincat <- which(unfilteredB[[id]]$scorabletype == unique(unfilteredB[[id]]$scorabletype)[i]) #which particles are in that category in the pre-randomized data
				partskeptB[[id]][[i]] <- partsincat[1:numpercat] #subsample row numbers
			}
			partskeptB.1[[id]] <- unlist(partskeptB[[id]])[!is.na(unlist(partskeptB[[id]]))] #unlist
			filteredB.1[[id]] <- unfilteredB[[id]][partskeptB.1[[id]],] #extract rows from the dataset
		}
	} else { #if random -- i.e. no category to assist randomizing
			for(id in 1:length(unfilteredB)){
				rows2sampl <- ifelse(numpercat<nrow(unfilteredB[[id]]), numpercat, nrow(unfilteredB[[id]]) )
				filteredB.1[[id]] <- unfilteredB[[id]][1:rows2sampl,] #just sample the first number of rows equal to the numpercat, or length of the sample
			}
		}
	filteredB <- data.frame(matrix(ncol=ncol(plastics),nrow=0))
	colnames(filteredB) <- colnames(plastics)
	for(i in 1:length(filteredB.1)){
		filteredB <- rbind(filteredB,filteredB.1[[i]])
	}
	
 
	if(random==FALSE){
		summarizeplastics <- table(paste(plastics$scorabletype)) / nrow(plastics)
		sortUscoretype <- sort(unique(paste(plastics$scorabletype)))
		proptable <- matrix(NA, nrow=length(sortUscoretype),ncol=length(sort(unique(plastics$truetype))))
		for(category in (1:length(sortUscoretype)) ){
			 catsample <-  filteredB[ paste(filteredB$scorabletype) == as.character(  sortUscoretype[category]  ) ,] 
			 proptable[category,] <- sapply(  sort(unique(plastics$truetype)),  function(z)  sum(as.character(catsample$truetype) == as.character(z))  )  / nrow(catsample)			
			}
		sampledrates <- sapply(1:ncol(proptable), function(type) sum(proptable[,type]*summarizeplastics) ) 

	} else{
		summarizeB <- sapply(1:length(sort(unique(plastics$truetype))), function(type)  
			length(which(as.character(filteredB$truetype)==as.character(sort(unique(plastics$truetype))[type]) ) ))
		sampledrates <-  summarizeB/sum(summarizeB) #
	}
	
	summarizeAll <- table(plastics$truetype)  #we summarize at the full dataset level for the true summary
	 truerates <- summarizeAll/sum(summarizeAll)
		#it is the goal to summarize the full matrix	
	names(sampledrates) <- paste( sort(unique(plastics$truetype)),"Sample",sep="_") 
	names(truerates) <- paste( sort(unique(plastics$truetype)),"Input",sep="_") 
	sumerrrates <- sum (abs(  sampledrates -  truerates  )  )/2
	avgerrrates <- mean (abs( sampledrates -  truerates   )  )
	
	#error, when averages consider stratified structure
	incTmat <- sapply(1:length(unique(plastics$truetype)), function(z)  tapply(plastics$truetype==sort(unique(plastics$truetype))[z] ,plastics$sampleid,sum) )
	prpTmat <- incTmat/rowSums(incTmat)
	trueStrat <- colMeans(prpTmat)
	incSmat <- matrix(0, ncol = length(unique(plastics$truetype)), nrow=length(unique(filteredB$sampleid)))#matrix with columns for each plastic truetype, rows for each sample in filteredB
	smplIDs <- sort(unique(filteredB$sampleid))
	for(i in 1:nrow(incSmat)){
		incSmat[i,]<- sapply(1:length(unique(plastics$truetype)), function(z) sum( filteredB$truetype==sort(unique(plastics$truetype))[z] & filteredB$sampleid==smplIDs[i]  ) )
	}
	prpSmat <- incSmat/rowSums(incSmat)
	sampleStrat <- colMeans(prpSmat)
	sumErrStrat <- sum(abs(sampleStrat - trueStrat))/2

	S.obs		<- DataInfo(as.vector(table(filteredB$truetype)),datatype="abundance")$S.obs #on RAW breakdown, not incidence freq. to be comparable across methods
	AsyChaoRaw <- ChaoRichness(as.vector(table(filteredB$truetype)),datatype="abundance")$Estimator 
	Coverage <- DataInfo(as.vector(table(filteredB$truetype)),datatype="abundance")$SC 


	return(   c(sumerrrates = sumerrrates, avgerrrates=avgerrrates, sumErrStrat = sumErrStrat,
				nobs = nrow(filteredB[!is.na(filteredB$truetype),]),  
				S.obs=S.obs,Coverage=Coverage, AsyChaoRaw = AsyChaoRaw) )
}


####this next function is nearly identical to StructuredData_depthper()
############ -- but now we don't care about retaining information across individual samples

StructuredData_pool <- function(plastics,numpercat,sampleidcolnum, scorabletypecolnum,truetypecolnum,random=F) { 
	 #dataset of scored plastics, the number per category sampling (or random must =T and this is the total subsample size), 
	 #and the column that defines the category by which subsampling will be structured
	colnames(plastics)[truetypecolnum] <- "truetype"
	colnames(plastics)[scorabletypecolnum] <- "scorabletype"  #rename column to be used for subsampling
	colnames(plastics)[sampleidcolnum] <- "sampleid"  #rename column to be used for subsampling
	unfilteredB<- plastics[sample(1:nrow(plastics),nrow(plastics),replace=F), ]  #randomize order
	if(random==FALSE){
		partskeptB<- list()
		for(i in 1:length(unique(unfilteredB$scorabletype))){  #for each category 
			partsincat <- which(unfilteredB$scorabletype == unique(unfilteredB$scorabletype)[i]) #which particles are in that category in the pre-randomized data
			partskeptB[[i]] <- partsincat[1:numpercat] #subsample row numbers
		}
		partskeptB.1 <- unlist(partskeptB)[!is.na(unlist(partskeptB))] #unlist
		filteredB <- unfilteredB[partskeptB.1,] #extract rows from the dataset
	} else { #if random -- i.e. no category to assist randomizing
			filteredB <- unfilteredB[1:numpercat,] #just sample the first number of rows equal to the sample size
		}
	
	if(random==FALSE){
		colnames(plastics)[scorabletypecolnum] <- "scorabletype"  
		summarizeplastics <- table(paste(plastics$scorabletype)) / nrow(plastics)
		sortUscoretype <- sort(unique(paste(plastics$scorabletype)))
		proptable <- matrix(NA, nrow=length(sortUscoretype),ncol=length(sort(unique(plastics$truetype))))
		for(category in (1:length(sortUscoretype)) ){
			 catsample <-  filteredB[ paste(filteredB$scorabletype) == as.character(  sortUscoretype[category]  ) ,] 
			 proptable[category,] <- sapply(  sort(unique(plastics$truetype)),  function(z)  sum(as.character(catsample$truetype) == as.character(z))  )  / nrow(catsample)			
			}
		sampledrates <- sapply(1:ncol(proptable), function(type) sum(proptable[,type]*summarizeplastics) ) 
	} else{
		summarizeB <- sapply(1:length(sort(unique(plastics$truetype))), function(type)  
			length(which(as.character(filteredB$truetype)==as.character(sort(unique(plastics$truetype))[type]) ) ))
		sampledrates <- summarizeB/sum(summarizeB)
	}
	
	summarizeAll <- table(plastics$truetype) 		
	names(sampledrates) <- paste( sort(unique(plastics$truetype)),"Sample",sep="_") 
	 truerates <- summarizeAll/sum(summarizeAll)
	names(truerates) <- paste( sort(unique(plastics$truetype)),"Input",sep="_") 
	sumerrrates <- sum (abs(  sampledrates -  truerates  )  )/2
	avgerrrates <- mean (abs( sampledrates -  truerates   )  )
	
	#error, when averages consider stratified structure
	incTmat <- sapply(1:length(unique(plastics$truetype)), function(z)  tapply(plastics$truetype==sort(unique(plastics$truetype))[z] ,plastics$sampleid,sum) )
	prpTmat <- incTmat/rowSums(incTmat)
	trueStrat <- colMeans(prpTmat)
	incSmat <- matrix(0, ncol = length(unique(plastics$truetype)), nrow=length(unique(filteredB$sampleid)))#matrix with columns for each plastic truetype, rows for each sample in filteredB
	smplIDs <- sort(unique(filteredB$sampleid))
	for(i in 1:nrow(incSmat)){
		incSmat[i,]<- sapply(1:length(unique(plastics$truetype)), function(z) sum( filteredB$truetype==sort(unique(plastics$truetype))[z] & filteredB$sampleid==smplIDs[i]  ) )
	}
	prpSmat <- incSmat/rowSums(incSmat)
	sampleStrat <- colMeans(prpSmat)
	sumErrStrat <- sum(abs(sampleStrat - trueStrat))/2


	S.obs		<- DataInfo(as.vector(table(filteredB$truetype)),datatype="abundance")$S.obs #on RAW breakdown
	AsyChaoRaw <- ChaoRichness(as.vector(table(filteredB$truetype)),datatype="abundance")$Estimator ###ON *RAW* breakdown!
	Coverage <- DataInfo(as.vector(table(filteredB$truetype)),datatype="abundance")$SC #on RAW breakdown

	return(c(  sumerrrates = sumerrrates, avgerrrates=avgerrrates,  sumErrStrat = sumErrStrat,
				nobs = nrow(filteredB[!is.na(filteredB$truetype),]),  
				S.obs=S.obs,Coverage=Coverage, AsyChaoRaw = AsyChaoRaw)) 
}


#short function for processing output of sampleplastics()
div.process <- function(divout,nrows=NULL){   #must supply output of sampleplastics(), and number of rows (the values over which the datasets were simulated)
	N <- 		lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) divout[[x]][[smplsz]]$nobs[1] ) )
	averr.m <- 	lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) mean(divout[[x]][[smplsz]]$avgerrrates,na.rm=T) ) )
	smerr.m <- 	lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) mean(divout[[x]][[smplsz]]$sumerrrates,na.rm=T) ) )
	Sobs.m <- 	lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) mean(divout[[x]][[smplsz]]$S.obs,na.rm=T) ) )
	cov.m <- 	lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) mean(divout[[x]][[smplsz]]$Coverage,na.rm=T) ) )
	asychR.m <- 	lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) mean(divout[[x]][[smplsz]]$AsyChaoRaw,na.rm=T) ) )
	averr.hpdi<-lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) HPDinterval(as.mcmc(divout[[x]][[smplsz]]$avgerrrates),0.95) ) )
	smerr.hpdi<-lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) HPDinterval(as.mcmc(divout[[x]][[smplsz]]$sumerrrates),0.95) ) )
	Sobs.hpdi <-lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) HPDinterval(as.mcmc(divout[[x]][[smplsz]]$S.obs),0.95) ) )
	cov.hpdi<-lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) HPDinterval(as.mcmc(divout[[x]][[smplsz]]$Coverage),0.95) ) )
	asychR.hpdi<-lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) HPDinterval(as.mcmc(divout[[x]][[smplsz]]$AsyChaoRaw),0.95) ) )	
	return(list(N=N,averr.m = averr.m, averr.hpdi=averr.hpdi,smerr.m = smerr.m, smerr.hpdi=smerr.hpdi, Sobs.m=Sobs.m,Sobs.hpdi=Sobs.hpdi,
				 asychR.m = asychR.m, asychR.hpdi = asychR.hpdi, cov.m = cov.m, cov.hpdi = cov.hpdi))
}

#short function for processing output of StructuredData_depthper() or StructuredData_pool(), same as above, but for datasets with multiple samples
div.processM <- function(divout,nrows=NULL) { #must supply output of StructuredData_depthper() or StructuredData_pool(), and number of rows (the values over which the datasets were simulated)
	N <- 		lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) mean(divout[[x]][[smplsz]]$nobs) ) )
	smSerr.m <- lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) mean(divout[[x]][[smplsz]]$sumErrStrat,na.rm=T) ) )
	smerr.m <- 	lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) mean(divout[[x]][[smplsz]]$sumerrrates,na.rm=T) ) )
	Sobs.m <- 	lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) mean(divout[[x]][[smplsz]]$S.obs,na.rm=T) ) )
	cov.m <- 	lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) mean(divout[[x]][[smplsz]]$Coverage,na.rm=T) ) )
	asychR.m <- 	lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) mean(divout[[x]][[smplsz]]$AsyChaoRaw,na.rm=T) ) )
	smSerr.hpdi<-lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) HPDinterval(as.mcmc(divout[[x]][[smplsz]]$sumErrStrat),0.95) ) )
	smerr.hpdi<-lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) HPDinterval(as.mcmc(divout[[x]][[smplsz]]$sumerrrates),0.95) ) )
	Sobs.hpdi <-lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) HPDinterval(as.mcmc(divout[[x]][[smplsz]]$S.obs),0.95) ) )
	cov.hpdi<-lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) HPDinterval(as.mcmc(divout[[x]][[smplsz]]$Coverage),0.95) ) )
	asychR.hpdi<-lapply(1:length(divout), function(x) sapply(1:(nrows[x]), function(smplsz) HPDinterval(as.mcmc(divout[[x]][[smplsz]]$AsyChaoRaw),0.95) ) )	
		divstats <- list(asychR.m = asychR.m, asychR.hpdi = asychR.hpdi, cov.m = cov.m, cov.hpdi = cov.hpdi)
	return(c(list(N=N,smerr.m = smerr.m, smerr.hpdi=smerr.hpdi,smSerr.m = smSerr.m, smSerr.hpdi=smSerr.hpdi, Sobs.m=Sobs.m,Sobs.hpdi=Sobs.hpdi),
				divstats ))
}