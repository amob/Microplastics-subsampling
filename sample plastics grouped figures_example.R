#####THIS SCRIPT IS AN EXAMPLE.
#####IT IS RUN WITH A SUBSET OF THE DATASETS USED IN THE ASSOCIATED PUBLISHED ANALYSIS
#####IT IS REPETITIVE IN SOME PLACES TO RETAIN THE LINE-BY-LINE STRUCTURE USED BY THE SCRIPTS FOR THE FULL ANALYSIS

#goal:
#use previous analyses, and figure code, to make multipanel figures for paper supplement

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

####################CODE SECTION BEGINNING HERE SAME AS IN OTHER SCRIPTS, notes in this script abbreviated/absent
### proceed to line 245

# read in the data from group samples, correct as in other analysis
plasticsRG <- read.csv("~/Microplastics-subsampling/Rain_garden_combined.csv",header=T,stringsAsFactors=F)
	#white_grey and grey_white are presumed non-equivalent. they are restricted to one subsample each. likewise, blue_white and white_blue
	#particle type typo corrections:****NOTE****these are either dealt with (glass) or consistent within subdatasets (paint, polyethylene), so in the other script, there is no need to change (because they are separate datasets)
	plasticsRG$Particle.type[which(plasticsRG$Particle.type=="glass ")] <- "glass"	
	plasticsRG$Particle.type[which(plasticsRG$Particle.type=="paint ")] <- "paint"	
	plasticsRG$Particle.type[which(plasticsRG$Particle.type=="polethylene")] <- "polyethylene"	

datasets <- list(plasticsRG)
datasetsbtype <- list(plasticsRG)
datnames <- c("rain garden")
datnamesbtype <- c("rain garden")
for(i in 1:length(datasets)){datasets[[i]]$ColorMorph <- paste(datasets[[i]]$Color, datasets[[i]]$Morphology,sep="") } #create a pasted color and morphology category
for(i in 1:length(datasetsbtype)){datasetsbtype[[i]]$ColorMorph <- paste(datasetsbtype[[i]]$Color, datasetsbtype[[i]]$Morphology,sep="") } #create a pasted color and morphology category
data_titles <- c("Storm water")
databtype_titles <- c("Storm water")

nrows<- unlist(lapply(datasets,nrow))
nrowsbtype<- unlist(lapply(datasetsbtype,nrow))

nrowsmax <- unlist(lapply(datasets, function(z) max(table(z$Subsample))  ) )
shpmax <- unlist(lapply(datasets, function(z) max(table(paste(z$Subsample,z$Morphology))) ) )
colmax <-  unlist(lapply(datasets, function(z) max(table(paste(z$Subsample, z$Color)))) )
shpcolmax <-  unlist(lapply(datasets, function(z) max(table(paste(z$Subsample, z$Morphology,z$Color)))) )

subsamps <-  unlist(lapply(datasets, function(z) length(table(z$Subsample)) ) )
subsampsbtype <-  unlist(lapply(datasetsbtype, function(z) length(table(z$Subsample)) ) )

shppmax <- unlist(lapply(datasets, function(z) max(table(z$Morphology))) )
colpmax <-  unlist(lapply(datasets, function(z) max(table(z$Color))) )
shpcolpmax <-  unlist(lapply(datasets, function(z) max(table(paste(z$Morphology,z$Color)))) )

nrowsmaxb <- unlist(lapply(datasetsbtype, function(z) max(table(z$Subsample))  ) )
shpmaxb <- unlist(lapply(datasetsbtype, function(z) max(table(paste(z$Subsample,z$Morphology))) ) )
colmaxb <-  unlist(lapply(datasetsbtype, function(z) max(table(paste(z$Subsample, z$Color)))) )
shpcolmaxb <-  unlist(lapply(datasetsbtype, function(z) max(table(paste(z$Subsample, z$Morphology,z$Color)))) )

shppmaxb <- unlist(lapply(datasetsbtype, function(z) max(table(z$Morphology))) )
colpmaxb <-  unlist(lapply(datasetsbtype, function(z) max(table(z$Color))) )
shpcolpmaxb <-  unlist(lapply(datasetsbtype, function(z) max(table(paste(z$Morphology,z$Color)))) )


#read in data from individual samples, correct as in analysis
plasticsRG1 <- read.csv("~/Microplastics-subsampling/raingarden.csv",header=T,stringsAsFactors=F) 
plasticsRG2 <- read.csv("~/Microplastics-subsampling/raingarden2.csv",header=T,stringsAsFactors=F) 
	plasticsRG2$Particle.type[plasticsRG2$Particle.type=="glass "] <- "glass"
plasticsRG3 <- read.csv("~/Microplastics-subsampling/raingarden3.csv",header=T,stringsAsFactors=F) 
datasetsI <- list(plasticsRG1,plasticsRG2,plasticsRG3)
for(i in 1:length(datasetsI)){datasetsI[[i]]$ColorMorph <- paste(datasetsI[[i]]$Color, datasetsI[[i]]$Morphology,sep="") } #create a pasted color and morphology category
data_titlesI <- c("Storm water", "Storm water","Storm water")
datnamesI <- c("Raingarden1","Raingarden2","Raingarden3")
nrowsI <- unlist(lapply(datasetsI,nrow))
shpmaxI <- unlist(lapply(datasetsI, function(z) max(table(z$Morphology))) )
colmaxI <-  unlist(lapply(datasetsI, function(z) max(table(z$Color))) )
shpcolmaxI <-  unlist(lapply(datasetsI, function(z) max(table(paste(z$Morphology,z$Color)))) )

datasetsBI <- list(plasticsRG1,plasticsRG2,plasticsRG3)
for(i in 1:length(datasetsBI)){datasetsBI[[i]]$ColorMorph <- paste(datasetsBI[[i]]$Color, datasetsBI[[i]]$Morphology,sep="") } #create a pasted color and morphology category
data_titlesBI <- c("Storm water", "Storm water","Storm water")
datnamesBI <- c("Raingarden1","Raingarden2","Raingarden3")
nrowsBI <- unlist(lapply(datasetsBI,nrow))
shpmaxBI <- unlist(lapply(datasetsBI, function(z) max(table(z$Morphology))) )
colmaxBI <-  unlist(lapply(datasetsBI, function(z) max(table(z$Color))) )
shpcolmaxBI <-  unlist(lapply(datasetsBI, function(z) max(table(paste(z$Morphology,z$Color)))) )


#load  results for grouped samples, process as in analysis
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


#load results for individual samples, process as in analysis.
load("~/Microplastics-subsampling/Rand_i_example.Rdata")
load("~/Microplastics-subsampling/Col_i_example.Rdata")
load("~/Microplastics-subsampling/Shp_i_example.Rdata")
load("~/Microplastics-subsampling/ShpCol_i_example.Rdata")
Randp <- div.process(Rand, nrows=nrowsI)
Colp <- div.process(Col, nrows=colmaxI) 
Shpp <- div.process(Shp, nrows=shpmaxI)
ShpColp <- div.process(ShpCol, nrows=shpcolmaxI)
load("~/Microplastics-subsampling/Rand_i_b_example.Rdata")
load("~/Microplastics-subsampling/Col_i_b_example.Rdata")
load("~/Microplastics-subsampling/Shp_i_b_example.Rdata")
load("~/Microplastics-subsampling/ShpCol_i_b_example.Rdata")
RandpB <- div.process(Randb, nrows=nrowsBI)
ColpB <- div.process(Colb, nrows=colmaxBI) 
ShppB <- div.process(Shpb, nrows=shpmaxBI)
ShpColpB <- div.process(ShpColb, nrows=shpcolmaxBI)


##minimum sampling calculations, first for grouped samples
#equal sampling comparing to structured error
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
#pooled sampling compared to pooled sample error
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

###minimum sampling calculations for individual samples
minimums <- matrix(ncol=8,nrow=length(datasetsI))
for(i in 1:length(datasetsI)) {
	minimums[i,] <- c( Randp$N[[i]] [min(which(  Randp$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  Randp$smerr.hpdi[[i]][2,] <0.2)),   
					   Colp$N[[i]] [min(which( Colp$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  Colp$smerr.hpdi[[i]][2,] <0.2)),   
					   Shpp$N[[i]] [min(which(  Shpp$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  Shpp$smerr.hpdi[[i]][2,] <0.2)),
					   ShpColp$N[[i]] [min(which( ShpColp$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ShpColp$smerr.hpdi[[i]][2,]<0.2))    )
					}   
minimums <- data.frame(minimums)
colnames(minimums) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimums <- cbind( data.frame(dataset = datnamesI, totparticles = nrowsI), minimums)
#
minimumsBI <- matrix(ncol=8,nrow=length(datasetsBI))
for(i in 1:length(datasetsBI)) {
	minimumsBI[i,] <- c( RandpB$N[[i]] [min(which(  RandpB$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  RandpB$smerr.hpdi[[i]][2,] <0.2)),   
					   ColpB$N[[i]] [min(which( ColpB$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ColpB$smerr.hpdi[[i]][2,] <0.2)),   
					   ShppB$N[[i]] [min(which(  ShppB$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ShppB$smerr.hpdi[[i]][2,] <0.2)),
					   ShpColpB$N[[i]] [min(which( ShpColpB$smerr.hpdi[[i]][2,] <0.2))],   
					   min(which(  ShpColpB$smerr.hpdi[[i]][2,]<0.2))  
	  )
					}   
minimumsBI <- data.frame(minimumsBI)
colnames(minimumsBI) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimums <- cbind( data.frame(dataset = datnamesBI, totparticles = nrowsBI), minimumsBI)


####################CODE SECTION ENDING HERE SAME AS IN OTHER SCRIPTS, notes in this script abbreviated


#####Make Figures###

##dataset order when sorted by total number of particles. datasets with single/multiple samples, and broad/narrow characterization
bysize <- order(nrows)
bysizeb <- order(nrowsbtype)
bysizeI <- order(nrowsI)
bysizeBI <- order(nrowsBI)

#######Summed Error across all dataset types, sampling types, panels a, b, and c for supplemental MS figures
#####For BROAD material types
pdf("~/Microplastics-subsampling/Broad_all_panel_a_indSumE_example.pdf",height=3.2,width=8)
layout(matrix(c(1:3),ncol=3,byrow=T))
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
for(i in bysizeBI){
plot(1~1, pch=NA, xlim=c(0,nrowsBI[i]),ylim=c(0,1),ylab="",xlab="")
	polygon( c(RandpB$N[[i]], rev(RandpB$N[[i]])), y= c(RandpB$smerr.hpdi[[i]][1,],rev(RandpB$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(RandpB$smerr.m[[i]] ~ I(RandpB$N[[i]]) )
	polygon( c(ColpB$N[[i]], rev(ColpB$N[[i]])), y= c(ColpB$smerr.hpdi[[i]][1,],rev(ColpB$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(ColpB$smerr.m[[i]] ~ I(ColpB$N[[i]]),  col=rgb(1,0,0))
	polygon( c(ShppB$N[[i]], rev(ShppB$N[[i]])), y= c(ShppB$smerr.hpdi[[i]][1,],rev(ShppB$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(ShppB$smerr.m[[i]] ~ I(ShppB$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColpB$N[[i]], rev(ShpColpB$N[[i]])), y= c(ShpColpB$smerr.hpdi[[i]][1,],rev(ShpColpB$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColpB$smerr.m[[i]] ~ I(ShpColpB$N[[i]]),col=rgb(1,0,1))
	abline(h=0.2,lty=2,lwd=1.5)
	points( cbind(minimumsBI$randnum,minimumsBI$colnum,minimumsBI$shpnum,minimumsBI$shpcolnum)[i,], jitter(rep(0.2, times = 4),factor=3),col=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1) ),cex=1.5,pch=16 )
	mtext(data_titlesBI[i], side=3, line =0.25, cex=1)
	if(i==2){mtext("Summed error",line=2.5,side=2)}
	if(i==1){mtext("Particles chemically identified",line=2.5,side=1)}
	if(i==3){legend(150,y=0.65,c("20%Error"),lty=2,bty="n")}
	if(i==3){legend(150,y=1,c("Rand","Col","Morph","MorphCol"),fill=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1)),bty="n")}
}
dev.off()
#
pdf("~/Microplastics-subsampling/Broad_all_panel_b_equalSumStrE_example.pdf",height=4.9,width=4.5)
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
i=1
plot(1~1, pch=NA, xlim=c(0,nrowsbtype[i]),ylim=c(0,1),ylab="",xlab="")
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
	mtext("Summed error",line=2.5,side=2)
	mtext("Particles chemically identified",line=2.5,side=1)
	legend(300,y=0.65,c("20%Error"),lty=2,bty="n")
dev.off()
#
pdf("~/Microplastics-subsampling/Broad_all_panel_c_pooledSumE_example.pdf",height=4.9,width=4.5)
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
i=1
plot(1~1, pch=NA, xlim=c(0,nrowsbtype[i]),ylim=c(0,1),ylab="",xlab="")
	polygon( c(RandSpoolB$N[[i]], rev(RandSpoolB$N[[i]])), y= c(RandSpoolB$smerr.hpdi[[i]][1,],rev(RandSpoolB$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(RandSpoolB$smerr.m[[i]] ~ I(RandSpoolB$N[[i]]) )
	polygon( c(ColSpoolB$N[[i]], rev(ColSpoolB$N[[i]])), y= c(ColSpoolB$smerr.hpdi[[i]][1,],rev(ColSpoolB$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(ColSpoolB$smerr.m[[i]] ~ I(ColSpoolB$N[[i]]),  col=rgb(1,0,0))
	polygon( c(ShpSpoolB$N[[i]], rev(ShpSpoolB$N[[i]])), y= c(ShpSpoolB$smerr.hpdi[[i]][1,],rev(ShpSpoolB$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(ShpSpoolB$smerr.m[[i]] ~ I(ShpSpoolB$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColSpoolB$N[[i]], rev(ShpColSpoolB$N[[i]])), y= c(ShpColSpoolB$smerr.hpdi[[i]][1,],rev(ShpColSpoolB$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColSpoolB$smerr.m[[i]] ~ I(ShpColSpoolB$N[[i]]),col=rgb(1,0,1))
	abline(h=0.2,lty=2,lwd=1.5)
	points( cbind(minimumsBpool$randnum,minimumsBpool$colnum,minimumsBpool$shpnum,minimumsBpool$shpcolnum)[i,], jitter(rep(0.2, times = 4),factor=1),col=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1) ),cex=1.5,pch=16 )
	mtext(databtype_titles[i], side=3, line =0.25, cex=1)
	mtext("Summed error",line=2.5,side=2)
	mtext("Particles chemically identified",line=2.5,side=1)
	legend(300,y=1,c("Rand","Col","Morph","MorphCol"),fill=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1)),bty="n")
dev.off()

#####For NARROW material types
pdf("~/Microplastics-subsampling/Narrow_all_panel_a_indSumE_example.pdf",height=3.2,width=8)
layout(matrix(c(1:4),ncol=4,byrow=T))
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
for(i in bysizeI){
plot(1~1, pch=NA, xlim=c(0,nrowsI[i]),ylim=c(0,1),ylab="",xlab="")
	polygon( c(Randp$N[[i]], rev(Randp$N[[i]])), y= c(Randp$smerr.hpdi[[i]][1,],rev(Randp$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(Randp$smerr.m[[i]] ~ I(Randp$N[[i]]) )
	polygon( c(Colp$N[[i]], rev(Colp$N[[i]])), y= c(Colp$smerr.hpdi[[i]][1,],rev(Colp$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(Colp$smerr.m[[i]] ~ I(Colp$N[[i]]),  col=rgb(1,0,0))
	polygon( c(Shpp$N[[i]], rev(Shpp$N[[i]])), y= c(Shpp$smerr.hpdi[[i]][1,],rev(Shpp$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(Shpp$smerr.m[[i]] ~ I(Shpp$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColp$N[[i]], rev(ShpColp$N[[i]])), y= c(ShpColp$smerr.hpdi[[i]][1,],rev(ShpColp$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColp$smerr.m[[i]] ~ I(ShpColp$N[[i]]),col=rgb(1,0,1))
	abline(h=0.2,lty=2,lwd=1.5)
	points( cbind(minimums$randnum,minimums$colnum,minimums$shpnum,minimums$shpcolnum)[i,], jitter(rep(0.2, times = 4),factor=3),col=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1) ),cex=1.5,pch=16 )
	mtext(data_titlesI[i], side=3, line =0.25, cex=1)
	if(i==2){mtext("Summed error",line=2.5,side=2)}
	if(i==1){mtext("Particles chemically identified",line=2.5,side=1)}
}
plot(1~1, pch=NA, xlim=c(0,nrowsI[i]),ylim=c(0,1),ylab="",xlab="",bty="n",xaxt="n",yaxt="n")
	legend(0,y=0.65,c("20%Error"),lty=2,bty="n")
	legend(0,y=1,c("Rand","Col","Morph","MorphCol"),fill=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1)),bty="n")
dev.off()
#
pdf("~/Microplastics-subsampling/Narrow_all_panel_b_equalSumStrE_example.pdf",height=4.9,width=4.5)
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
i=1
plot(1~1, pch=NA, xlim=c(0,nrows[i]),ylim=c(0,1),ylab="",xlab="")
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
	mtext("Particles chemically identified",line=2.5,side=1)
dev.off()
pdf("~/Microplastics-subsampling/Narrow_all_panel_c_poolSumE_example.pdf",height=4.9,width=4.5)
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
i=1
plot(1~1, pch=NA, xlim=c(0,nrows[i]),ylim=c(0,1),ylab="",xlab="")
	polygon( c(RandSpoolN$N[[i]], rev(RandSpoolN$N[[i]])), y= c(RandSpoolN$smerr.hpdi[[i]][1,],rev(RandSpoolN$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(RandSpoolN$smerr.m[[i]] ~ I(RandSpoolN$N[[i]]) )
	polygon( c(ColSpoolN$N[[i]], rev(ColSpoolN$N[[i]])), y= c(ColSpoolN$smerr.hpdi[[i]][1,],rev(ColSpoolN$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(ColSpoolN$smerr.m[[i]] ~ I(ColSpoolN$N[[i]]),  col=rgb(1,0,0))
	polygon( c(ShpSpoolN$N[[i]], rev(ShpSpoolN$N[[i]])), y= c(ShpSpoolN$smerr.hpdi[[i]][1,],rev(ShpSpoolN$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(ShpSpoolN$smerr.m[[i]] ~ I(ShpSpoolN$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColSpoolN$N[[i]], rev(ShpColSpoolN$N[[i]])), y= c(ShpColSpoolN$smerr.hpdi[[i]][1,],rev(ShpColSpoolN$smerr.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColSpoolN$smerr.m[[i]] ~ I(ShpColSpoolN$N[[i]]),col=rgb(1,0,1))
	abline(h=0.2,lty=2,lwd=1.5)
	points( cbind(minimumsNpool$randnum,minimumsNpool$colnum,minimumsNpool$shpnum,minimumsNpool$shpcolnum)[i,], jitter(rep(0.2, times = 4),factor=1),col=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1) ),cex=1.5,pch=16 )
	mtext(data_titles[i], side=3, line =0.25, cex=1)
	mtext("Summed error",line=2.5,side=2)
	mtext("Particles chemically identified",line=2.5,side=1)
dev.off()


#####How well the subsample captures the number of material types observed in the full sample
#####across all dataset types, sampling types, panels a, b, and c for supplemental MS figure
#####For NARROW material types only
richnessI <- unlist(lapply(datasetsI, function(x) length(table(x$Particle.type))))
richness <- unlist(lapply(datasets, function(x) length(table(x$Particle.type))))
pdf("~/Microplastics-subsampling/Narrow_panel_a_indEstSampleRichness_example.pdf",height=3.2,width=8)
layout(matrix(c(1:3,4),ncol=4,byrow=T))
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
for(i in bysizeI){
plot(1~1, pch=NA, xlim=c(0,nrowsI[i]),ylim=c(0,20),ylab="",xlab="")
	polygon( c(Randp$N[[i]], rev(Randp$N[[i]])), y= c(Randp$Sobs.hpdi[[i]][1,],rev(Randp$Sobs.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(Randp$Sobs.m[[i]] ~ I(Randp$N[[i]]) )
	polygon( c(Colp$N[[i]], rev(Colp$N[[i]])), y= c(Colp$Sobs.hpdi[[i]][1,],rev(Colp$Sobs.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(Colp$Sobs.m[[i]] ~ I(Colp$N[[i]]),  col=rgb(1,0,0))
	polygon( c(Shpp$N[[i]], rev(Shpp$N[[i]])), y= c(Shpp$Sobs.hpdi[[i]][1,],rev(Shpp$Sobs.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(Shpp$Sobs.m[[i]] ~ I(Shpp$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColp$N[[i]], rev(ShpColp$N[[i]])), y= c(ShpColp$Sobs.hpdi[[i]][1,],rev(ShpColp$Sobs.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColp$Sobs.m[[i]] ~ I(ShpColp$N[[i]]),col=rgb(1,0,1))
	mtext(data_titlesI[i], side=3, line =0.25, cex=1)
	abline(h=richnessI[i],lty=3,lwd=1.5)
	if(i==2){mtext("# Material types observed",line=2.5,side=2)}
	if(i==1){mtext("Particles chemically identified",line=2.5,side=1)}
}
plot(1~1, pch=NA, xlim=c(0,nrowsI[i]),ylim=c(0,1),ylab="",xlab="",bty="n",xaxt="n",yaxt="n")
	legend(0,y=0.3,c("Full Sample/Dataset # obs"),lty=3,lwd=1.5,bty="n")
	legend(0,y=1,c("Rand","Col","Morph","MorphCol"),fill=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1)),bty="n")
dev.off()
#
pdf("~/Microplastics-subsampling/Narrow_panel_b_equalEstDatasetRichnes_example.pdf",height=4.9,width=4.5)
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
i=1
plot(1~1, pch=NA, xlim=c(0,nrows[i]),ylim=c(0,40),ylab="",xlab="")
	polygon( c(RandSdpN$N[[i]], rev(RandSdpN$N[[i]])), y= c(RandSdpN$Sobs.hpdi[[i]][1,],rev(RandSdpN$Sobs.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(RandSdpN$Sobs.m[[i]] ~ I(RandSdpN$N[[i]]) )
	polygon( c(ColSdpN$N[[i]], rev(ColSdpN$N[[i]])), y= c(ColSdpN$Sobs.hpdi[[i]][1,],rev(ColSdpN$Sobs.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(ColSdpN$Sobs.m[[i]] ~ I(ColSdpN$N[[i]]),  col=rgb(1,0,0))
	polygon( c(ShpSdpN$N[[i]], rev(ShpSdpN$N[[i]])), y= c(ShpSdpN$Sobs.hpdi[[i]][1,],rev(ShpSdpN$Sobs.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(ShpSdpN$Sobs.m[[i]] ~ I(ShpSdpN$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColSdpN$N[[i]], rev(ShpColSdpN$N[[i]])), y= c(ShpColSdpN$Sobs.hpdi[[i]][1,],rev(ShpColSdpN$Sobs.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColSdpN$Sobs.m[[i]] ~ I(ShpColSdpN$N[[i]]),col=rgb(1,0,1))
	abline(h=richness[i],lty=3,lwd=1.5,col=rgb(0.5,0.5,0.5))
	mtext(data_titles[i], side=3, line =0.25, cex=1)
	mtext("# Material types observed",line=2.5,side=2)
	mtext("Particles chemically identified",line=2.5,side=1)
dev.off()
#
pdf("~/Microplastics-subsampling/Narrow_panel_c_poolEstDatasetRichness_example.pdf",height=4.9,width=4.5)
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
i=1
plot(1~1, pch=NA, xlim=c(0,nrows[i]),ylim=c(0,40),ylab="",xlab="")
	polygon( c(RandSpoolN$N[[i]], rev(RandSpoolN$N[[i]])), y= c(RandSpoolN$Sobs.hpdi[[i]][1,],rev(RandSpoolN$Sobs.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(RandSpoolN$Sobs.m[[i]] ~ I(RandSpoolN$N[[i]]) )
	polygon( c(ColSpoolN$N[[i]], rev(ColSpoolN$N[[i]])), y= c(ColSpoolN$Sobs.hpdi[[i]][1,],rev(ColSpoolN$Sobs.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(ColSpoolN$Sobs.m[[i]] ~ I(ColSpoolN$N[[i]]),  col=rgb(1,0,0))
	polygon( c(ShpSpoolN$N[[i]], rev(ShpSpoolN$N[[i]])), y= c(ShpSpoolN$Sobs.hpdi[[i]][1,],rev(ShpSpoolN$Sobs.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(ShpSpoolN$Sobs.m[[i]] ~ I(ShpSpoolN$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColSpoolN$N[[i]], rev(ShpColSpoolN$N[[i]])), y= c(ShpColSpoolN$Sobs.hpdi[[i]][1,],rev(ShpColSpoolN$Sobs.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColSpoolN$Sobs.m[[i]] ~ I(ShpColSpoolN$N[[i]]),col=rgb(1,0,1))
	abline(h=richness[i],lty=3,lwd=1.5,col=rgb(0.5,0.5,0.5))
	mtext(data_titles[i], side=3, line =0.25, cex=1)
	mtext("# Material types observed",line=2.5,side=2)
	mtext("Particles chemically identified",line=2.5,side=1)
dev.off()


#####How well how well the (sub)sample captures the number of material types estimated to exist in the environment from which the full sample was drawn
#####across all dataset types, sampling types, panels a, b, and c for supplemental MS figure
#####For NARROW material types only
chaoI <- unlist(lapply(datasetsI, function(x) ChaoRichness(as.vector(table(x$Particle.type)),datatype="abundance")$Estimator ) )#asymptotic chao estimate of richness
chao <- unlist(lapply(datasets, function(x) ChaoRichness(as.vector(table(x$Particle.type)),datatype="abundance")$Estimator ) )#asymptotic chao estimate of richness
pdf("~/Microplastics-subsampling/Narrow_panel_a_indEstEnvRichness_example.pdf",height=3.2,width=8)
layout(matrix(c(1:3,4),ncol=4,byrow=T))
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
for(i in bysizeI){
plot(1~1, pch=NA, xlim=c(0,nrowsI[i]),ylim=c(0,65),ylab="",xlab="")
	polygon( c(Randp$N[[i]], rev(Randp$N[[i]])), y= c(Randp$asychR.hpdi[[i]][1,],rev(Randp$asychR.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(Randp$asychR.m[[i]] ~ I(Randp$N[[i]]) )
	polygon( c(Colp$N[[i]], rev(Colp$N[[i]])), y= c(Colp$asychR.hpdi[[i]][1,],rev(Colp$asychR.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(Colp$asychR.m[[i]] ~ I(Colp$N[[i]]),  col=rgb(1,0,0))
	polygon( c(Shpp$N[[i]], rev(Shpp$N[[i]])), y= c(Shpp$asychR.hpdi[[i]][1,],rev(Shpp$asychR.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(Shpp$asychR.m[[i]] ~ I(Shpp$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColp$N[[i]], rev(ShpColp$N[[i]])), y= c(ShpColp$asychR.hpdi[[i]][1,],rev(ShpColp$asychR.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColp$asychR.m[[i]] ~ I(ShpColp$N[[i]]),col=rgb(1,0,1))
	abline(h=chaoI[i],lty=3,lwd=1.5)
	abline(h=richnessI[i],lty=3,lwd=1.5,col=rgb(0.5,0.5,0.5))
	mtext(data_titlesI[i], side=3, line =0.25, cex=1)
	if(i==2){mtext("# Material types estimated (Chao)",line=2.5,side=2,adj=1)}#using chao from subsample
	if(i==1){mtext("Particles chemically identified",line=2.5,side=1,adj=-0.05)}
}
plot(1~1, pch=NA, xlim=c(0,nrowsI[i]),ylim=c(0,1),ylab="",xlab="",bty="n",xaxt="n",yaxt="n")
	legend(-25,y=0.4,c("Full Sample/Dataset Chao Estimated #","Full Sample/Dataset # obs"),lty=3,lwd=1.5,bty="n",col=c(rgb(0,0,0),rgb(0.5,0.5,0.5) ) )
	legend(-10,y=1,c("Rand","Col","Morph","MorphCol"),fill=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1)),bty="n")
dev.off()
#
pdf("~/Microplastics-subsampling/Narrow_panel_b_equalEstEnvRichness_example.pdf",height=4.9,width=4.5)
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
i=1
plot(1~1, pch=NA, xlim=c(0,nrows[i]),ylim=c(0,110),ylab="",xlab="")
	polygon( c(RandSdpN$N[[i]], rev(RandSdpN$N[[i]])), y= c(RandSdpN$asychR.hpdi[[i]][1,],rev(RandSdpN$asychR.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(RandSdpN$asychR.m[[i]] ~ I(RandSdpN$N[[i]]) )
	polygon( c(ColSdpN$N[[i]], rev(ColSdpN$N[[i]])), y= c(ColSdpN$asychR.hpdi[[i]][1,],rev(ColSdpN$asychR.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(ColSdpN$asychR.m[[i]] ~ I(ColSdpN$N[[i]]),  col=rgb(1,0,0))
	polygon( c(ShpSdpN$N[[i]], rev(ShpSdpN$N[[i]])), y= c(ShpSdpN$asychR.hpdi[[i]][1,],rev(ShpSdpN$asychR.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(ShpSdpN$asychR.m[[i]] ~ I(ShpSdpN$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColSdpN$N[[i]], rev(ShpColSdpN$N[[i]])), y= c(ShpColSdpN$asychR.hpdi[[i]][1,],rev(ShpColSdpN$asychR.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColSdpN$asychR.m[[i]] ~ I(ShpColSdpN$N[[i]]),col=rgb(1,0,1))
	abline(h=chao[i],lty=3,lwd=1.5)
	abline(h=richness[i],lty=3,lwd=1.5,col=rgb(0.5,0.5,0.5))
	mtext(data_titles[i], side=3, line =0.25, cex=1)
	mtext("# Material types estimated (Chao)",line=2.5,side=2)#using chao from subsample
	mtext("Particles chemically identified",line=2.5,side=1)
dev.off()
#
pdf("~/Microplastics-subsampling/Narrow_panel_c_poolEstEnvRichness_example.pdf",height=4.9,width=4.5)
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
i=1
plot(1~1, pch=NA, xlim=c(0,nrows[i]),ylim=c(0,110),ylab="",xlab="")
	polygon( c(RandSpoolN$N[[i]], rev(RandSpoolN$N[[i]])), y= c(RandSpoolN$asychR.hpdi[[i]][1,],rev(RandSpoolN$asychR.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(RandSpoolN$asychR.m[[i]] ~ I(RandSpoolN$N[[i]]) )
	polygon( c(ColSpoolN$N[[i]], rev(ColSpoolN$N[[i]])), y= c(ColSpoolN$asychR.hpdi[[i]][1,],rev(ColSpoolN$asychR.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(ColSpoolN$asychR.m[[i]] ~ I(ColSpoolN$N[[i]]),  col=rgb(1,0,0))
	polygon( c(ShpSpoolN$N[[i]], rev(ShpSpoolN$N[[i]])), y= c(ShpSpoolN$asychR.hpdi[[i]][1,],rev(ShpSpoolN$asychR.hpdi[[i]][2,]) ), border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(ShpSpoolN$asychR.m[[i]] ~ I(ShpSpoolN$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColSpoolN$N[[i]], rev(ShpColSpoolN$N[[i]])), y= c(ShpColSpoolN$asychR.hpdi[[i]][1,],rev(ShpColSpoolN$asychR.hpdi[[i]][2,]) ), border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColSpoolN$asychR.m[[i]] ~ I(ShpColSpoolN$N[[i]]),col=rgb(1,0,1))
	abline(h=chao[i],lty=3,lwd=1.5)
	abline(h=richness[i],lty=3,lwd=1.5,col=rgb(0.5,0.5,0.5))
	mtext(data_titles[i], side=3, line =0.25, cex=1)
	mtext("# Material types estimated (Chao)",line=2.5,side=2)#using chao from subsample
	mtext("Particles chemically identified",line=2.5,side=1)
dev.off()


#####How well the sample/subsample captures the estimated percent of particles in the environment
#####--proportion of environmental particles that are the same material type as a (sub)sampled particle
#####across all dataset types, sampling types, panels a, b, and c for supplemental MS figure
#####For NARROW material types only
coverageI <- unlist(lapply(datasetsI, function(x) DataInfo(as.vector(table(x$Particle.type)),datatype="abundance")$SC ) )
coverage <- unlist(lapply(datasets, function(x) DataInfo( as.vector(table(x$Particle.type)),datatype="abundance")$SC ) )
pdf("~/Microplastics-subsampling/Narrow_panel_a_indEstEnvCoverage_example.pdf",height=3.2,width=8)
layout(matrix(c(1:3,4),ncol=4,byrow=T))
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
for(i in bysizeI){
plot(1~1, pch=NA, xlim=c(0,nrowsI[i]),ylim=c(40,102),ylab="",xlab="")
	polygon( c(Randp$N[[i]][-1], rev(Randp$N[[i]][-1])), y= c(Randp$cov.hpdi[[i]][1,-1],rev(Randp$cov.hpdi[[i]][2,-1]) )*100, border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(Randp$cov.m[[i]][-1]*100 ~ I(Randp$N[[i]][-1]) )
	polygon( c(Colp$N[[i]], rev(Colp$N[[i]])), y= c(Colp$cov.hpdi[[i]][1,],rev(Colp$cov.hpdi[[i]][2,]) )*100, border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(Colp$cov.m[[i]]*100 ~ I(Colp$N[[i]]),  col=rgb(1,0,0))
	polygon( c(Shpp$N[[i]], rev(Shpp$N[[i]])), y= c(Shpp$cov.hpdi[[i]][1,],rev(Shpp$cov.hpdi[[i]][2,]) )*100, border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(Shpp$cov.m[[i]]*100 ~ I(Shpp$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColp$N[[i]], rev(ShpColp$N[[i]])), y= c(ShpColp$cov.hpdi[[i]][1,],rev(ShpColp$cov.hpdi[[i]][2,]) )*100, border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColp$cov.m[[i]]*100 ~ I(ShpColp$N[[i]]),col=rgb(1,0,1))
	abline(h=100,lty=2,lwd=1.5)
	mtext(data_titlesI[i], side=3, line =0.25, cex=1)
	if(i==2){mtext("Estimated %coverage, material types",line=2.5,side=2,adj=1,cex=0.8)}
	if(i==1){mtext("Particles chemically identified",line=2.5,side=1,adj=-0.05)}
}
plot(1~1, pch=NA, xlim=c(0,nrowsI[i]),ylim=c(0,1),ylab="",xlab="",bty="n",xaxt="n",yaxt="n")
	legend(0,y=0.3,c("100%"),lty=2,bty="n")
	legend(0,y=1,c("Rand","Col","Morph","MorphCol"),fill=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1)),bty="n")
dev.off()
#
pdf("~/Microplastics-subsampling/Narrow_panel_b_equalEstEnvCoverage_example.pdf",height=4.9,width=4.5)
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
i=1
plot(1~1, pch=NA, xlim=c(0,nrows[i]),ylim=c(40,102),ylab="",xlab="")
	polygon( c(RandSdpN$N[[i]][-1], rev(RandSdpN$N[[i]][-1])), y= c(RandSdpN$cov.hpdi[[i]][1,-1],rev(RandSdpN$cov.hpdi[[i]][2,-1]) )*100, border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(RandSdpN$cov.m[[i]][-1]*100 ~ I(RandSdpN$N[[i]][-1]) )
	polygon( c(ColSdpN$N[[i]], rev(ColSdpN$N[[i]])), y= c(ColSdpN$cov.hpdi[[i]][1,],rev(ColSdpN$cov.hpdi[[i]][2,]) )*100, border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(ColSdpN$cov.m[[i]]*100 ~ I(ColSdpN$N[[i]]),  col=rgb(1,0,0))
	polygon( c(ShpSdpN$N[[i]], rev(ShpSdpN$N[[i]])), y= c(ShpSdpN$cov.hpdi[[i]][1,],rev(ShpSdpN$cov.hpdi[[i]][2,]) )*100, border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(ShpSdpN$cov.m[[i]]*100 ~ I(ShpSdpN$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColSdpN$N[[i]], rev(ShpColSdpN$N[[i]])), y= c(ShpColSdpN$cov.hpdi[[i]][1,],rev(ShpColSdpN$cov.hpdi[[i]][2,]) )*100, border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColSdpN$cov.m[[i]]*100 ~ I(ShpColSdpN$N[[i]]),col=rgb(1,0,1))
	abline(h=100,lty=2,lwd=1.5)
	mtext(data_titles[i], side=3, line =0.25, cex=1)
	mtext("Estimated %coverage, material types",line=2.5,side=2)
	mtext("Particles chemically identified",line=2.5,side=1)
dev.off()
#
pdf("~/Microplastics-subsampling/Narrow_panel_c_poolEstEnvCoverage_example.pdf",height=4.9,width=4.5)
par(mar=c(2,1,2,1))
par(oma=c(2,3,1,0))
i=1
plot(1~1, pch=NA, xlim=c(0,nrows[i]),ylim=c(40,102),ylab="",xlab="")
	polygon( c(RandSpoolN$N[[i]][-1], rev(RandSpoolN$N[[i]][-1])), y= c(RandSpoolN$cov.hpdi[[i]][1,-1],rev(RandSpoolN$cov.hpdi[[i]][2,-1]) )*100, border=NA, col=rgb(0,0,0,alpha=0.25) )
	lines(RandSpoolN$cov.m[[i]][-1]*100 ~ I(RandSpoolN$N[[i]][-1]) )
	polygon( c(ColSpoolN$N[[i]], rev(ColSpoolN$N[[i]])), y= c(ColSpoolN$cov.hpdi[[i]][1,],rev(ColSpoolN$cov.hpdi[[i]][2,]) )*100, border=NA, col=rgb(1,0,0,alpha=0.25) )
	lines(ColSpoolN$cov.m[[i]]*100 ~ I(ColSpoolN$N[[i]]),  col=rgb(1,0,0))
	polygon( c(ShpSpoolN$N[[i]], rev(ShpSpoolN$N[[i]])), y= c(ShpSpoolN$cov.hpdi[[i]][1,],rev(ShpSpoolN$cov.hpdi[[i]][2,]) )*100, border=NA, col=rgb(0,0,1,alpha=0.25) )
	lines(ShpSpoolN$cov.m[[i]]*100 ~ I(ShpSpoolN$N[[i]]),col=rgb(0,0,1))
	polygon( c(ShpColSpoolN$N[[i]], rev(ShpColSpoolN$N[[i]])), y= c(ShpColSpoolN$cov.hpdi[[i]][1,],rev(ShpColSpoolN$cov.hpdi[[i]][2,]) )*100, border=NA, col=rgb(1,0,1,alpha=0.25) )
	lines(ShpColSpoolN$cov.m[[i]]*100 ~ I(ShpColSpoolN$N[[i]]),col=rgb(1,0,1))
	abline(h=100,lty=2,lwd=1.5)
	mtext(data_titles[i], side=3, line =0.25, cex=1)
	mtext("Estimated %coverage, material types",line=2.5,side=2)
	mtext("Particles chemically identified",line=2.5,side=1)
dev.off()
