#####THIS SCRIPT IS AN EXAMPLE.
#####IT IS RUN WITH A SUBSET OF THE DATASETS USED IN THE ASSOCIATED PUBLISHED ANALYSIS
#####IT IS REPETITIVE IN SOME PLACES TO RETAIN THE LINE-BY-LINE STRUCTURE USED BY THE SCRIPTS FOR THE FULL ANALYSIS

#goals:
#estimating %particle types in full set of samples from subsample / set of subsamples
#this is for datasets with multiple samples ("many jars")
#success is getting same % in each category as when estimated by full sample 
#another success metric is the extent to which the subsample (and also full sample) is expected to characterize the environment
#This second metric is evaluated with sample coverage and Chao's estimator
#There are two sets of categories here, "broad" types (plastic, anthropogenic, natural) and all material types ("narrow"", all chemically identified types)
#Because these are multi-sample datasets, there are two primary types of subsampling, and two ways the subsample is evaluated for success: 
#	-drawing evenly from samples, in which case the summed error comparison is weighted by sample of origin  ("dp"/"equal")
#	-pooling all particles and drawing from the pool, in which case the summed error comparison is not weighted ("pool")
#check success of each subsample method across datasets

library(coda)
library(iNEXT) #diversity metrics


range01=function(x){
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
}

dropNA <- function(vect){return(vect[!is.na(vect)])}

#read in functions for sampling plastics
source('~/Microplastics-subsampling/sample\ plastics\ functions.R') 
#source location depends on local file structure! might need to change path, or set working directory

# read in the data; checked categories with e.g. table(plasticsLF$Particle.type); if problems noticed, code replaces some names; same for color and morphology (no typos in morphology)
plasticsRG <- read.csv("~/Microplastics-subsampling/Rain_garden_combined.csv",header=T,stringsAsFactors=F)
	#white_grey and grey_white are presumed non-equivalent. they are restricted to one subsample each. likewise, blue_white and white_blue
	#particle type typo corrections:****NOTE****these are either dealt with (glass) or consistent within subdatasets (paint, polyethylene), so in the other script, there is no need to change (because they are separate datasets)
	plasticsRG$Particle.type[which(plasticsRG$Particle.type=="glass ")] <- "glass"	
	plasticsRG$Particle.type[which(plasticsRG$Particle.type=="paint ")] <- "paint"	
	plasticsRG$Particle.type[which(plasticsRG$Particle.type=="polethylene")] <- "polyethylene"	

datasets <- list(plasticsRG)
datasetsbtype <- list(plasticsRG)

datnames <- c("rain garden")
datnamesbtype <- c("lorenz plankton", "lorenz sediment","bourdages birds", "teboul birds","CJ sediment", "CJ inverts", "CJ inverts time", "rain garden","vidal")

for(i in 1:length(datasets)){datasets[[i]]$ColorMorph <- paste(datasets[[i]]$Color, datasets[[i]]$Morphology,sep="") } #create a pasted color and morphology category
for(i in 1:length(datasetsbtype)){datasetsbtype[[i]]$ColorMorph <- paste(datasetsbtype[[i]]$Color, datasetsbtype[[i]]$Morphology,sep="") } #create a pasted color and morphology category


data_titles <- c("Storm water")
databtype_titles <- c("Storm water")

#get number of particles in each dataset (for both groups of datasets)
nrows<- unlist(lapply(datasets,nrow))
nrowsbtype<- unlist(lapply(datasetsbtype,nrow))
nsamps <- unlist(lapply(datasets,function(z) length(unique(z$Subsample))))

#extract measures of diversity using metrics to estimate what the sampling site might truly contain (beyond even what was sampled, chao);
	#and a simple tally of what was sampled: "actual" = observed polymer richness
	# -- only for datasets with narrow type characterization
	# pooling all particles from all samples
richness <- unlist(lapply(datasets, function(x) length(table(x$Particle.type))))
coverage <- unlist(lapply(datasets, function(x) DataInfo( as.vector(table(x$Particle.type)),datatype="abundance")$SC ) )
chao <- unlist(lapply(datasets, function(x) ChaoRichness(as.vector(table(x$Particle.type)),datatype="abundance")$Estimator ) )#asymptotic chao estimate of richness
#descriptive table of coverage and richness -- only for datasets with narrow type characterization
DescStats <- data.frame(nsamps, nrows,coverage, 1-coverage, richness, chao, chao/richness, chao - richness)
DescStats <- DescStats[order(nrows),]
colnames(DescStats) <- c("Samples","SampleSize", "Coverage", "1MinusCoverage","ObservedTypes","ChaoEstimator","Chao/ObsTypes","ChaoMinusObsTypes")
write.csv(round(DescStats,digits=4),"~/Microplastics-subsampling/descriptive stats_narrow_manyjars_example.csv",row.names=F)

#columns for necessary information are not identical in order/number across datasets, this exctracts the needed column numbers
#for datasets where narrow characterization is possible
colorcol <- lapply(datasets,function(z) which(colnames(z)=="Color"))
morphcol <- lapply(datasets,function(z) which(colnames(z)=="Morphology"))
colmorphcol <- lapply(datasets,function(z) which(colnames(z)=="ColorMorph"))
subsampcol <- lapply(datasets,function(z) which(colnames(z)=="Subsample"))
PTcol <- lapply(datasets,function(z) which(colnames(z)=="Particle.type"))
#for datasets where only broad characterization is possible
colorcolb <- lapply(datasetsbtype,function(z) which(colnames(z)=="Color"))
morphcolb <- lapply(datasetsbtype,function(z) which(colnames(z)=="Morphology"))
colmorphcolb <- lapply(datasetsbtype,function(z) which(colnames(z)=="ColorMorph"))
subsampcolb <- lapply(datasetsbtype,function(z) which(colnames(z)=="Subsample"))
PTBcolb <- lapply(datasetsbtype,function(z) which(colnames(z)=="Particle.type.broad"))

######NUMBERS needed for subsampling when goal is to characterize particles at the NARROW level
#get max number of particles in any category for stratified sampling (e.g. in shpmax, most particles are fibers, #fibers)
#in this case the max is also limited by the largest sample 
#when subsampling evenly across samples, all datest particles will be sampled when particles/sample in subsample = size of the largest sample
nrowsmax <- unlist(lapply(datasets, function(z) max(table(z$Subsample))  ) )
shpmax <- unlist(lapply(datasets, function(z) max(table(paste(z$Subsample,z$Morphology))) ) )
colmax <-  unlist(lapply(datasets, function(z) max(table(paste(z$Subsample, z$Color)))) )
shpcolmax <-  unlist(lapply(datasets, function(z) max(table(paste(z$Subsample, z$Morphology,z$Color)))) )
#get the maximum subsample particles when subsampling this way (for narrow characterization)
subsamps <-  unlist(lapply(datasets, function(z) length(table(z$Subsample)) ) )
#get max number of particles in any category for stratified sampling (e.g. in shpmax, most particles are fibers, #fibers)
#in this case, these are for subsampling from the particles when they are pooled across dataset samples
shppmax <- unlist(lapply(datasets, function(z) max(table(z$Morphology))) )
colpmax <-  unlist(lapply(datasets, function(z) max(table(z$Color))) )
shpcolpmax <-  unlist(lapply(datasets, function(z) max(table(paste(z$Morphology,z$Color)))) )

######NUMBERS needed for subsampling when goal is to characterize particles at the BROAD level
######This is the same as the previous section, but with one additional dataset
#get max number of particles in any category for stratified sampling (e.g. in shpmax, most particles are fibers, #fibers)
#in this case the max is also limited by the largest sample 
#when subsampling evenly across samples, all datest particles will be sampled when particles/sample in subsample = size of the largest sample
nrowsmaxb <- unlist(lapply(datasetsbtype, function(z) max(table(z$Subsample))  ) )
shpmaxb <- unlist(lapply(datasetsbtype, function(z) max(table(paste(z$Subsample,z$Morphology))) ) )
colmaxb <-  unlist(lapply(datasetsbtype, function(z) max(table(paste(z$Subsample, z$Color)))) )
shpcolmaxb <-  unlist(lapply(datasetsbtype, function(z) max(table(paste(z$Subsample, z$Morphology,z$Color)))) )
#get the maximum subsample particles when subsampling this way (for broad characterization)
subsampsbtype <-  unlist(lapply(datasetsbtype, function(z) length(table(z$Subsample)) ) )
#get max number of particles in any category for stratified sampling (e.g. in shpmax, most particles are fibers, #fibers)
#in this case, these are for subsampling from the particles when they are pooled across dataset samples
shppmaxb <- unlist(lapply(datasetsbtype, function(z) max(table(z$Morphology))) )
colpmaxb <-  unlist(lapply(datasetsbtype, function(z) max(table(z$Color))) )
shpcolpmaxb <-  unlist(lapply(datasetsbtype, function(z) max(table(paste(z$Morphology,z$Color)))) )



#######################
##The next 55 lines have code that takes a long time to run

######Simulations for subsampling when goal is to characterize particles at the NARROW level
###Subsampling evenly across samples
Rand_str_dp_narrow <- 	lapply(1:length(datasets), function(x) lapply(1:(nrowsmax[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_depthper(datasets[[x]],smplsz,subsampcol[[x]],1,PTcol[[x]],random=TRUE) ) )) ))# sampling randomly, note the scoreable type column number is not used and can be any number other than the truetype column
save(Rand_str_dp_narrow,file="~/Microplastics-subsampling/Rand_str_dp_narrow_example.Rdata")
Shp_str_dp_narrow <- 		lapply(1:length(datasets), function(x) lapply(1:(shpmax[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_depthper(datasets[[x]],smplsz,subsampcol[[x]],morphcol[[x]],PTcol[[x]],random=FALSE) ) )) ))#
save(Shp_str_dp_narrow,file="~/Microplastics-subsampling/Shp_str_dp_narrow_example.Rdata")
Col_str_dp_narrow <- 		lapply(1:length(datasets), function(x) lapply(1:(colmax[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_depthper(datasets[[x]],smplsz,subsampcol[[x]],colorcol[[x]],PTcol[[x]],random=FALSE) ) )) ))#
save(Col_str_dp_narrow,file="~/Microplastics-subsampling/Col_str_dp_narrow_example.Rdata")
ShpCol_str_dp_narrow <- 		lapply(1:length(datasets), function(x) lapply(1:(shpcolmax[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_depthper(datasets[[x]],smplsz,subsampcol[[x]],colmorphcol[[x]],PTcol[[x]],random=FALSE) ) )) ))#
save(ShpCol_str_dp_narrow,file="~/Microplastics-subsampling/ShpCol_str_dp_narrow_example.Rdata")
###Subsampling from the pool of particles from all samples in the dataset
Rand_str_pool_narrow <- 	lapply(1:length(datasets), function(x) lapply(1:(nrows[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_pool(datasets[[x]],smplsz,subsampcol[[x]],1,PTcol[[x]],random=TRUE) ) )) ))# sampling randomly, note the scoreable type column number is not used and can be any number other than the truetype column
save(Rand_str_pool_narrow,file="~/Microplastics-subsampling/Rand_str_pool_narrow_example.Rdata")
Shp_str_pool_narrow <- 		lapply(1:length(datasets), function(x) lapply(1:(shppmax[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_pool(datasets[[x]],smplsz,subsampcol[[x]],morphcol[[x]],PTcol[[x]],random=FALSE) ) )) ))#
save(Shp_str_pool_narrow,file="~/Microplastics-subsampling/Shp_str_pool_narrow_example.Rdata")
Col_str_pool_narrow <- 		lapply(1:length(datasets), function(x) lapply(1:(colpmax[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_pool(datasets[[x]],smplsz,subsampcol[[x]],colorcol[[x]],PTcol[[x]],random=FALSE) ) )) ))#
save(Col_str_pool_narrow,file="~/Microplastics-subsampling/Col_str_pool_narrow_example.Rdata")
ShpCol_str_pool_narrow <- 		lapply(1:length(datasets), function(x) lapply(1:(shpcolpmax[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_pool(datasets[[x]],smplsz,subsampcol[[x]],colmorphcol[[x]],PTcol[[x]],random=FALSE) ) )) ))#
save(ShpCol_str_pool_narrow,file="~/Microplastics-subsampling/ShpCol_str_pool_narrow_example.Rdata")

######Simulations for subsampling when goal is to characterize particles at the BROAD level
###Subsampling evenly across samples
Rand_str_dp_b <- 	lapply(1:length(datasetsbtype), function(x) lapply(1:(nrowsmaxb[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_depthper(datasetsbtype[[x]],smplsz,subsampcolb[[x]],1,PTBcolb[[x]],random=TRUE) ) )) ))# sampling randomly, note the scoreable type column number is not used and can be any number other than the true type column
save(Rand_str_dp_b,file="~/Microplastics-subsampling/Rand_str_dp_b_example.Rdata")
Shp_str_dp_b <- 		lapply(1:length(datasetsbtype), function(x) lapply(1:(shpmaxb[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_depthper(datasetsbtype[[x]],smplsz,subsampcolb[[x]],morphcolb[[x]],PTBcolb[[x]],random=FALSE) ) )) ))#
save(Shp_str_dp_b,file="~/Microplastics-subsampling/Shp_str_dp_b_example.Rdata")
Col_str_dp_b <- 		lapply(1:length(datasetsbtype), function(x) lapply(1:(colmaxb[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_depthper(datasetsbtype[[x]],smplsz,subsampcolb[[x]],colorcolb[[x]],PTBcolb[[x]],random=FALSE) ) )) ))#
save(Col_str_dp_b,file="~/Microplastics-subsampling/Col_str_dp_b_example.Rdata")
ShpCol_str_dp_b <- 		lapply(1:length(datasetsbtype), function(x) lapply(1:(shpcolmaxb[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_depthper(datasetsbtype[[x]],smplsz,subsampcolb[[x]],colmorphcolb[[x]],PTBcolb[[x]],random=FALSE) ) )) ))#
save(ShpCol_str_dp_b,file="~/Microplastics-subsampling/ShpCol_str_dp_b_example.Rdata")
###Subsampling from the pool of particles from all samples in the dataset
Rand_str_pool_b <- 	lapply(1:length(datasetsbtype), function(x) lapply(1:(nrowsbtype[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_pool(datasetsbtype[[x]],smplsz,subsampcolb[[x]],1,PTBcolb[[x]],random=TRUE) ) )) ))# sampling randomly, note the scoreable type column number is not used and can be any number other than the true type column
save(Rand_str_pool_b,file="~/Microplastics-subsampling/Rand_str_pool_b_example.Rdata")
Shp_str_pool_b <- 		lapply(1:length(datasetsbtype), function(x) lapply(1:(shppmaxb[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_pool(datasetsbtype[[x]],smplsz,subsampcolb[[x]],morphcolb[[x]],PTBcolb[[x]],random=FALSE) ) )) ))#
save(Shp_str_pool_b,file="~/Microplastics-subsampling/Shp_str_pool_b_example.Rdata")
Col_str_pool_b <- 		lapply(1:length(datasetsbtype), function(x) lapply(1:(colpmaxb[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_pool(datasetsbtype[[x]],smplsz,subsampcolb[[x]],colorcolb[[x]],PTBcolb[[x]],random=FALSE) ) )) ))#
save(Col_str_pool_b,file="~/Microplastics-subsampling/Col_str_pool_b_example.Rdata")
ShpCol_str_pool_b <- 		lapply(1:length(datasetsbtype), function(x) lapply(1:(shpcolpmaxb[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) 
						StructuredData_pool(datasetsbtype[[x]],smplsz,subsampcolb[[x]],colmorphcolb[[x]],PTBcolb[[x]],random=FALSE) ) )) ))#
save(ShpCol_str_pool_b,file="~/Microplastics-subsampling/ShpCol_str_pool_b_example.Rdata")

#########################end section of code that takes a long time to run

#load pre-run (or re-run in previous lines) results

load("~/Microplastics-subsampling/Rand_str_dp_narrow_example.Rdata")
load("~/Microplastics-subsampling/Shp_str_dp_narrow_example.Rdata")
load("~/Microplastics-subsampling/Col_str_dp_narrow_example.Rdata")
load("~/Microplastics-subsampling/ShpCol_str_dp_narrow_example.Rdata")
RandSdpN <- div.processM(Rand_str_dp_narrow, nrows=nrowsmax)
ShpSdpN <- div.processM(Shp_str_dp_narrow, nrows=shpmax)
ColSdpN <- div.processM(Col_str_dp_narrow, nrows=colmax)
ShpColSdpN <- div.processM(ShpCol_str_dp_narrow, nrows=shpcolmax)

load("~/Microplastics-subsampling/Rand_str_pool_narrow_example.Rdata")
load("~/Microplastics-subsampling/Shp_str_pool_narrow_example.Rdata")
load("~/Microplastics-subsampling/Col_str_pool_narrow_example.Rdata")
load("~/Microplastics-subsampling/ShpCol_str_pool_narrow_example.Rdata")
RandSpoolN <- div.processM(Rand_str_pool_narrow, nrows=nrows)
ShpSpoolN <- div.processM(Shp_str_pool_narrow, nrows=shppmax)
ColSpoolN <- div.processM(Col_str_pool_narrow, nrows=colpmax)
ShpColSpoolN <- div.processM(ShpCol_str_pool_narrow, nrows=shpcolpmax)

load("~/Microplastics-subsampling/Rand_str_dp_b_example.Rdata")
load("~/Microplastics-subsampling/Shp_str_dp_b_example.Rdata")
load("~/Microplastics-subsampling/Col_str_dp_b_example.Rdata")
load("~/Microplastics-subsampling/ShpCol_str_dp_b_example.Rdata")
RandSdpB <- div.processM(Rand_str_dp_b, nrows=nrowsmaxb)
ShpSdpB <- div.processM(Shp_str_dp_b, nrows=shpmaxb)
ColSdpB <- div.processM(Col_str_dp_b, nrows=colmaxb)
ShpColSdpB <- div.processM(ShpCol_str_dp_b, nrows=shpcolmaxb)

load("~/Microplastics-subsampling/Rand_str_pool_b_example.Rdata")
load("~/Microplastics-subsampling/Shp_str_pool_b_example.Rdata")
load("~/Microplastics-subsampling/Col_str_pool_b_example.Rdata")
load("~/Microplastics-subsampling/ShpCol_str_pool_b_example.Rdata")
RandSpoolB <- div.processM(Rand_str_pool_b, nrows=nrowsbtype)
ShpSpoolB <- div.processM(Shp_str_pool_b, nrows=shppmaxb)
ColSpoolB <- div.processM(Col_str_pool_b, nrows=colpmaxb)
ShpColSpoolB <- div.processM(ShpCol_str_pool_b, nrows=shpcolpmaxb)


#########MININUM SUBSAMPLING CALCULATIONS FROM WHICH TO DRAW RECOMMENDATIONS

#####Error in the proportion of particles estimated in each material type
###subsampling equally across samples, comparing to structured error
##NARROW characterization of material types
minimumsNSerr <- matrix(ncol=8,nrow=length(datasets))
for(i in 1:length(datasets)) {
	minimumsNSerr[i,] <- c( RandSdpN$N[[i]] [min(which(  RandSdpN$smSerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  RandSdpN$smSerr.hpdi[[i]][2,] <0.2)),   
					   ColSdpN$N[[i]] [min(which( ColSdpN$smSerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ColSdpN$smSerr.hpdi[[i]][2,] <0.2)),   
					   ShpSdpN$N[[i]] [min(which(  ShpSdpN$smSerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ShpSdpN$smSerr.hpdi[[i]][2,] <0.2)),
					   ShpColSdpN$N[[i]] [min(which( ShpColSdpN$smSerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ShpColSdpN$smSerr.hpdi[[i]][2,]<0.2))    )
					}   
minimumsNSerr <- data.frame(minimumsNSerr)
colnames(minimumsNSerr) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimumsNSerr <- cbind( data.frame(dataset = datnames, totparticles = nrows, totsubsamps = subsamps), minimumsNSerr)
write.csv(minimumsNSerr,"~/Microplastics-subsampling/minimums_SummedStrErrorLessThan20_narrow_manyjarsEqual_example.csv",row.names=F)
##BROAD characterization of material types
minimumsSerr <- matrix(ncol=8,nrow=length(datasetsbtype))
for(i in 1:length(datasetsbtype)) {
	minimumsSerr[i,] <- c( RandSdpB$N[[i]] [min(which(  RandSdpB$smSerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  RandSdpB$smSerr.hpdi[[i]][2,] <0.2)),   
					   ColSdpB$N[[i]] [min(which( ColSdpB$smSerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ColSdpB$smSerr.hpdi[[i]][2,] <0.2)),   
					   ShpSdpB$N[[i]] [min(which(  ShpSdpB$smSerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ShpSdpB$smSerr.hpdi[[i]][2,] <0.2)),
					   ShpColSdpB$N[[i]] [min(which( ShpColSdpB$smSerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ShpColSdpB$smSerr.hpdi[[i]][2,]<0.2))    )
					}   
minimumsSerr <- data.frame(minimumsSerr)
colnames(minimumsSerr) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimumsSerr <- cbind( data.frame(dataset = datnamesbtype, totparticles = nrowsbtype, totsubsamps = subsampsbtype), minimumsSerr)
write.csv(minimumsSerr,"~/Microplastics-subsampling/minimums_SummedStrErrorLessThan20_broad_manyjarsEqual_example.csv",row.names=F)
###subsampling from pooled samples compared to pooled sample error
##NARROW characterization of material types
minimumsNpool <- matrix(ncol=8,nrow=length(datasets))
for(i in 1:length(datasets)) {
	minimumsNpool[i,] <- c( RandSpoolN$N[[i]] [min(which(  RandSpoolN$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  RandSpoolN$smerr.hpdi[[i]][2,] <0.2)),   
					   ColSpoolN$N[[i]] [min(which( ColSpoolN$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ColSpoolN$smerr.hpdi[[i]][2,] <0.2)),   
					   ShpSpoolN$N[[i]] [min(which(  ShpSpoolN$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ShpSpoolN$smerr.hpdi[[i]][2,] <0.2)),
					   ShpColSpoolN$N[[i]] [min(which( ShpColSpoolN$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ShpColSpoolN$smerr.hpdi[[i]][2,]<0.2))    )
					}   
minimumsNpool <- data.frame(minimumsNpool)
colnames(minimumsNpool) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimumsNpool <- cbind( data.frame(dataset = datnames, totparticles = nrows, totsubsamps = subsamps), minimumsNpool)
write.csv(minimumsNpool,"~/Microplastics-subsampling/minimums_SummedErrorLessThan20_narrow_manyjarsPool_example.csv",row.names=F)
##BROAD characterization of material types
minimumsBpool <- matrix(ncol=8,nrow=length(datasetsbtype))
for(i in 1:length(datasetsbtype)) {
	minimumsBpool[i,] <- c( RandSpoolB$N[[i]] [min(which(  RandSpoolB$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  RandSpoolB$smerr.hpdi[[i]][2,] <0.2)),   
					   ColSpoolB$N[[i]] [min(which( ColSpoolB$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ColSpoolB$smerr.hpdi[[i]][2,] <0.2)),   
					   ShpSpoolB$N[[i]] [min(which(  ShpSpoolB$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ShpSpoolB$smerr.hpdi[[i]][2,] <0.2)),
					   ShpColSpoolB$N[[i]] [min(which( ShpColSpoolB$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ShpColSpoolB$smerr.hpdi[[i]][2,]<0.2))    )
					}   
minimumsBpool <- data.frame(minimumsBpool)
colnames(minimumsBpool) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimumsBpool <- cbind( data.frame(dataset = datnamesbtype, totparticles = nrowsbtype, totsubsamps = subsampsbtype), minimumsBpool)
write.csv(minimumsBpool,"~/Microplastics-subsampling/minimums_SummedErrorLessThan20_broad_manyjarsPool_example.csv",row.names=F)

#####Recapitulating number of particle types in the full set of samples (richness)
#####ALL NARROW characterization of material types
###subsampling equally across samples
minimumsR <- matrix(ncol=8,nrow=length(datasets))
for(i in 1:length(datasets)) {
	minimumsR[i,] <- c(RandSdpN$N[[i]] [min(which(  RandSdpN$Sobs.hpdi[[i]][1,] >= 0.8*richness[i]))],   
					   min(which(  RandSdpN$Sobs.hpdi[[i]][1,] >= 0.8*richness[i])),   
					   ColSdpN$N[[i]] [min(which(  ColSdpN$Sobs.hpdi[[i]][1,]>= 0.8*richness[i]))],   
					   min(which(  ColSdpN$Sobs.hpdi[[i]][1,]>= 0.8*richness[i])),   
					   ShpSdpN$N[[i]] [min(which(  ShpSdpN$Sobs.hpdi[[i]][1,]>= 0.8*richness[i]))],   
					   min(which(  ShpSdpN$Sobs.hpdi[[i]][1,]>= 0.8*richness[i])),   
					   ShpColSdpN$N[[i]] [min(which(  ShpColSdpN$Sobs.hpdi[[i]][1,]>= 0.8*richness[i]))],
					   min(which(  ShpColSdpN$Sobs.hpdi[[i]][1,]>= 0.8*richness[i]))   )   
					}   
minimumsR <- data.frame(minimumsR)
colnames(minimumsR) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimumsR <- cbind( data.frame(dataset = datnames, totparticles = nrows, totsubsamps = subsamps), minimumsR)
write.csv(minimumsR,"~/Microplastics-subsampling/minimums_EstSammpleRichnessGreaterThan80_narrow_manyjarsEqual_example.csv",row.names=F)
###subsampling from pooled samples
minimumsRP <- matrix(ncol=8,nrow=length(datasets))
for(i in 1:length(datasets)) {
	minimumsRP[i,] <- c(RandSpoolN$N[[i]] [min(which(  RandSpoolN$Sobs.hpdi[[i]][1,] >= 0.8*richness[i]))],   
					   min(which(  RandSpoolN$Sobs.hpdi[[i]][1,] >= 0.8*richness[i])),   
					   ColSpoolN$N[[i]] [min(which(  ColSpoolN$Sobs.hpdi[[i]][1,]>= 0.8*richness[i]))],   
					   min(which(  ColSpoolN$Sobs.hpdi[[i]][1,]>= 0.8*richness[i])),   
					   ShpSpoolN$N[[i]] [min(which(  ShpSpoolN$Sobs.hpdi[[i]][1,]>= 0.8*richness[i]))],   
					   min(which(  ShpSpoolN$Sobs.hpdi[[i]][1,]>= 0.8*richness[i])),   
					   ShpColSpoolN$N[[i]] [min(which(  ShpColSpoolN$Sobs.hpdi[[i]][1,]>= 0.8*richness[i]))],
					   min(which(  ShpColSpoolN$Sobs.hpdi[[i]][1,]>= 0.8*richness[i]))   )   
					}   
minimumsRP <- data.frame(minimumsRP)
colnames(minimumsRP) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimumsRP <- cbind( data.frame(dataset = datnames, totparticles = nrows, totsubsamps = subsamps), minimumsRP)
write.csv(minimumsRP,"~/Microplastics-subsampling/minimums_EstSammpleRichnessGreaterThan80_narrow_manyjarsPool_example.csv",row.names=F)

#####Recapitulating number of particle types estimated to be in the environment by the full set of samples
#####Chao's estimator
#####ALL NARROW characterization of material types
###subsampling equally across samples
minimumsAE <- matrix(ncol=8,nrow=length(datasets))
for(i in 1:length(datasets)) {
	minimumsAE[i,] <- c( RandSdpN$N[[i]] [min(which(  RandSdpN$asychR.hpdi[[i]][1,] >= 0.8*chao[i]))],   
					   min(which(  RandSdpN$asychR.hpdi[[i]][1,] >= 0.8*chao[i])),   
					   ColSdpN$N[[i]] [min(which(  ColSdpN$asychR.hpdi[[i]][1,]>= 0.8*chao[i]))],
					   min(which(  ColSdpN$asychR.hpdi[[i]][1,]>= 0.8*chao[i])),
					   ShpSdpN$N[[i]] [min(which(  ShpSdpN$asychR.hpdi[[i]][1,]>= 0.8*chao[i]))],
					   min(which(  ShpSdpN$asychR.hpdi[[i]][1,]>= 0.8*chao[i])),
					   ShpColSdpN$N[[i]] [min(which(  ShpColSdpN$asychR.hpdi[[i]][1,]>= 0.8*chao[i]))],
					   min(which(  ShpColSdpN$asychR.hpdi[[i]][1,]>= 0.8*chao[i])))   
					}   
minimumsAE <- data.frame(minimumsAE)
colnames(minimumsAE) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimumsAE <- cbind( data.frame(dataset = datnames, totparticles = nrows, totsubsamps = subsamps), minimumsAE)
write.csv(minimumsAE,"~/Microplastics-subsampling/minimums_EstEnvRichnessgreaterthan80_narrow_manyjarsEqual_example.csv",row.names=F)
###subsampling from pooled samples
minimumsAEP <- matrix(ncol=8,nrow=length(datasets))
for(i in 1:length(datasets)) {
	minimumsAEP[i,] <- c( RandSpoolN$N[[i]] [min(which(  RandSpoolN$asychR.hpdi[[i]][1,] >= 0.8*chao[i]))],   
					   min(which(  RandSpoolN$asychR.hpdi[[i]][1,] >= 0.8*chao[i])),   
					   ColSpoolN$N[[i]] [min(which(  ColSpoolN$asychR.hpdi[[i]][1,]>= 0.8*chao[i]))],
					   min(which(  ColSpoolN$asychR.hpdi[[i]][1,]>= 0.8*chao[i])),
					   ShpSpoolN$N[[i]] [min(which(  ShpSpoolN$asychR.hpdi[[i]][1,]>= 0.8*chao[i]))],
					   min(which(  ShpSpoolN$asychR.hpdi[[i]][1,]>= 0.8*chao[i])),
					   ShpColSpoolN$N[[i]] [min(which(  ShpColSpoolN$asychR.hpdi[[i]][1,]>= 0.8*chao[i]))],
					   min(which(  ShpColSpoolN$asychR.hpdi[[i]][1,]>= 0.8*chao[i])))   
					}   
minimumsAEP <- data.frame(minimumsAEP)
colnames(minimumsAEP) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimumsAEP <- cbind( data.frame(dataset = datnames, totparticles = nrows, totsubsamps = subsamps), minimumsAEP)
write.csv(minimumsAEP,"~/Microplastics-subsampling/minimums_EstEnvRichnessgreaterthan80_narrow_manyjarsPool_example.csv",row.names=F)

#####Recapitulating number of particle types estimated to be in the environment by the full set of samples
#####Chao's estimator
#####ALL NARROW characterization of material types
###subsampling equally across samples
minimumsC <- matrix(ncol=8,nrow=length(datasets))
for(i in 1:length(datasets)) {
	minimumsC[i,] <- c( RandSdpN$N[[i]] [min(which(  RandSdpN$cov.hpdi[[i]][1,] >= 0.8)) ],   
					   min(which(  RandSdpN$cov.hpdi[[i]][1,] >= 0.8)),   
					   ColSdpN$N[[i]] [min(which(  ColSdpN$cov.hpdi[[i]][1,]>= 0.8))],
					   min(which(  ColSdpN$cov.hpdi[[i]][1,]>= 0.8)),
					   ShpSdpN$N[[i]] [min(which(  ShpSdpN$cov.hpdi[[i]][1,]>= 0.8))],
					   min(which(  ShpSdpN$cov.hpdi[[i]][1,]>= 0.8)),
					   ShpColSdpN$N[[i]] [min(which(  ShpColSdpN$cov.hpdi[[i]][1,]>= 0.8))],
					   min(which(  ShpColSdpN$cov.hpdi[[i]][1,]>= 0.8)))   
					}   
minimumsC <- data.frame(minimumsC)
colnames(minimumsC) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimumsC <- cbind( data.frame(dataset = datnames, totparticles = nrows, totsubsamps = subsamps), minimumsC)
write.csv(minimumsC,"~/Microplastics-subsampling/minimums_EstEnvCoverageGreaterThan80_narrow_manyjarsEqual_example.csv",row.names=F)
##
minimumsCP <- matrix(ncol=8,nrow=length(datasets))
for(i in 1:length(datasets)) {
	minimumsCP[i,] <- c( RandSpoolN$N[[i]] [min(which(  RandSpoolN$cov.hpdi[[i]][1,-1] >= 0.8)) +1],   #coverage cannot be estimated from a sample size 1, must remove first obs
					   min(which(  RandSpoolN$cov.hpdi[[i]][1,-1] >= 0.8))+1,    #coverage cannot be estimated from a sample size 1, must remove first obs
					   ColSpoolN$N[[i]] [min(which(  ColSpoolN$cov.hpdi[[i]][1,]>= 0.8))],
					   min(which(  ColSpoolN$cov.hpdi[[i]][1,]>= 0.8)),
					   ShpSpoolN$N[[i]] [min(which(  ShpSpoolN$cov.hpdi[[i]][1,]>= 0.8))],
					   min(which(  ShpSpoolN$cov.hpdi[[i]][1,]>= 0.8)),
					   ShpColSpoolN$N[[i]] [min(which(  ShpColSpoolN$cov.hpdi[[i]][1,]>= 0.8))],
					   min(which(  ShpColSpoolN$cov.hpdi[[i]][1,]>= 0.8)))   
					}   
minimumsCP <- data.frame(minimumsCP)
colnames(minimumsCP) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimumsCP <- cbind( data.frame(dataset = datnames, totparticles = nrows, totsubsamps = subsamps), minimumsCP)
write.csv(minimumsCP,"~/Microplastics-subsampling/minimums_EstEnvCoverageGreaterThan80_narrow_manyjarsPool_example.csv",row.names=F)



#####Make Figures###

bysize <- order(nrows)
bysizeb <- order(nrowsbtype)

#summed structured error, subsampling equally
pdf("~/Microplastics-subsampling/SummedStrError_narrow_manyjarsEqual_example.pdf",height=3,width=3.5)
par(oma=c(1,4,1,2))
par(mar=c(3,0,2,0))
i=1
plot(1~1, pch=NA, xlim=c(0,nrows[i]),ylim=c(0,1),ylab="",xlab="",yaxt="n")
	polygon( c(RandSdpN$N[[i]], rev(RandSdpN$N[[i]])), y= c(RandSdpN$smSerr.hpdi[[i]][1,],rev(RandSdpN$smSerr.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(RandSdpN$smSerr.m[[i]] ~ I(RandSdpN$N[[i]]) )
	polygon( c(ColSdpN$N[[i]], rev(ColSdpN$N[[i]])), y= c(ColSdpN$smSerr.hpdi[[i]][1,],rev(ColSdpN$smSerr.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(ColSdpN$smSerr.m[[i]] ~ I(ColSdpN$N[[i]]),  col=rgb(1,0,0))
	polygon( c(ShpSdpN$N[[i]], rev(ShpSdpN$N[[i]])), y= c(ShpSdpN$smSerr.hpdi[[i]][1,],rev(ShpSdpN$smSerr.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(ShpSdpN$smSerr.m[[i]] ~ I(ShpSdpN$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColSdpN$N[[i]], rev(ShpColSdpN$N[[i]])), y= c(ShpColSdpN$smSerr.hpdi[[i]][1,],rev(ShpColSdpN$smSerr.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColSdpN$smSerr.m[[i]] ~ I(ShpColSdpN$N[[i]]),col=rgb(1,0,1))
	abline(h=0.2,lty=2,lwd=1.5)
	points( cbind(minimumsNSerr$randnum,minimumsNSerr$colnum,minimumsNSerr$shpnum,minimumsNSerr$shpcolnum)[i,], jitter(rep(0.2, times = 4),factor=1),col=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1) ),cex=1.5,pch=16 )
	mtext(data_titles[i], side=3, line =0.25, cex=1)
	mtext("Summed error",line=2.5,side=2)
 	axis(side=2)
	legend(150,y=0.55,c("20%Error"),lty=2,bty="n")
	legend(150,y=1.1,c("Rand","Col","Morph","MorphCol"),fill=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1)),bty="n")
	mtext("Particles chemically identified",line=2.5,side=1)
dev.off()
pdf("~/Microplastics-subsampling/SummedStrError_broad_manyjarsEqual_example.pdf",height=3,width=3.5)
par(oma=c(1,4,1,2))
par(mar=c(3,0,2,0))
i=1
plot(1~1, pch=NA, xlim=c(0,nrows[i]),ylim=c(0,1),ylab="",xlab="",yaxt="n")
	polygon( c(RandSdpB$N[[i]], rev(RandSdpB$N[[i]])), y= c(RandSdpB$smSerr.hpdi[[i]][1,],rev(RandSdpB$smSerr.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(RandSdpB$smSerr.m[[i]] ~ I(RandSdpB$N[[i]]) )
	polygon( c(ColSdpB$N[[i]], rev(ColSdpB$N[[i]])), y= c(ColSdpB$smSerr.hpdi[[i]][1,],rev(ColSdpB$smSerr.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(ColSdpB$smSerr.m[[i]] ~ I(ColSdpB$N[[i]]),  col=rgb(1,0,0))
	polygon( c(ShpSdpB$N[[i]], rev(ShpSdpB$N[[i]])), y= c(ShpSdpB$smSerr.hpdi[[i]][1,],rev(ShpSdpB$smSerr.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(ShpSdpB$smSerr.m[[i]] ~ I(ShpSdpB$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColSdpB$N[[i]], rev(ShpColSdpB$N[[i]])), y= c(ShpColSdpB$smSerr.hpdi[[i]][1,],rev(ShpColSdpB$smSerr.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColSdpB$smSerr.m[[i]] ~ I(ShpColSdpB$N[[i]]),col=rgb(1,0,1))
	abline(h=0.2,lty=2,lwd=1.5)
	points( cbind(minimumsSerr$randnum,minimumsSerr$colnum,minimumsSerr$shpnum,minimumsSerr$shpcolnum)[i,], jitter(rep(0.2, times = 4),factor=1),col=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1) ),cex=1.5,pch=16 )
	mtext(databtype_titles[i], side=3, line =0.25, cex=1)
 	axis(side=2)
	mtext("Summed error",line=2.5,side=2)
	legend(150,y=0.55,c("20%Error"),lty=2,bty="n")
	legend(150,y=1.1,c("Rand","Col","Morph","MorphCol"),fill=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1)),bty="n")
	mtext("Particles chemically identified",line=2.5,side=1)
dev.off()

#Figures for summed error (when subsampling from pooled samples) are only in supplemental figures.
#see "sample plastics grouped figures_example.R"


#see "sample plastics grouped figures_example.R" for figures depicting how well subsamples capture the number of particle types in the sample
#and how well subsamples/datasets represent particles in the environment
