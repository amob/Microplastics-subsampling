#goals:
#estimating %particle types in full sample from subsample
#this is for datasets with a single sample ("one jar")
#success is getting same % in each category as when estimated by full sample 
#another success metric is the extent to which the subsample (and also full sample) is expected to characterize the environment
#This second metric is evaluated with sample coverage and Chao's estimator
#categories here are all material types ("narrow"), meaning split into all chemically identified types
#check success of each subsample method across datasets

library(coda)
library(iNEXT) #diversity metrics

range01=function(x){
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
}

dropNA <- function(vect){return(vect[!is.na(vect)])}

#read in functions for sampling plastics
source('~/Sample-microplastics/sample\ plastics\ functions.R') 
#source location depends on local file structure! might need to change path, or set working directorys


# read in the data; checked categories with e.g. table(plasticsLF$Particle.type); if problems noticed, code replaces some names; same for color and morphology (no typos in morphology)
plasticsLF <- read.csv("~/Sample-microplastics/Lusher_Fish.csv",header=T,stringsAsFactors=F)
	plasticsLF$Color[plasticsLF$Color=="white "] <- "white"
plasticsLS <- read.csv("~/Sample-microplastics/Lusher_Sediment.csv",header=T,stringsAsFactors=F)
plasticsLM <- read.csv("~/Sample-microplastics/Lusher_Mussels.csv",header=T,stringsAsFactors=F)
plasticsBM <- read.csv("~/Sample-microplastics/Brate_Mussels.csv",header=T,stringsAsFactors=F) 
plasticsP1 <- read.csv("~/Sample-microplastics/Primpke_1.csv",header=T,stringsAsFactors=F)
plasticsP2 <- read.csv("~/Sample-microplastics/Primpke_2.csv",header=T,stringsAsFactors=F)
plasticsP3 <- read.csv("~/Sample-microplastics/Primpke_3.csv",header=T,stringsAsFactors=F)
plasticsRG1 <- read.csv("~/Sample-microplastics/raingarden.csv",header=T,stringsAsFactors=F) 
plasticsRG2 <- read.csv("~/Sample-microplastics/raingarden2.csv",header=T,stringsAsFactors=F) 
	plasticsRG2$Particle.type[plasticsRG2$Particle.type=="glass "] <- "glass"
plasticsRG3 <- read.csv("~/Sample-microplastics/raingarden3.csv",header=T,stringsAsFactors=F) 
datasets <- list(plasticsLF,plasticsLS,plasticsLM,plasticsBM,plasticsP1,plasticsP2,plasticsP3,plasticsRG1,plasticsRG2,plasticsRG3)

for(i in 1:length(datasets)){datasets[[i]]$ColorMorph <- paste(datasets[[i]]$Color, datasets[[i]]$Morphology,sep="") } #create a pasted color and morphology category

data_titles <- c("Animal tissue","Sediment", "Animal tissue","Animal tissue", "Surface water", "Surface water", "Surface water", "Storm water", "Storm water","Storm water")

#extract measures of diversity using metrics to estimate what the sampling site might truly contain (beyond even what was sampled, chao); 
	#and a simple tally of what was sampled: "richness" = observed polymer richness in full sample
richness <- unlist(lapply(datasets, function(x) length(table(x$Particle.type))))
coverage <- unlist(lapply(datasets, function(x) DataInfo(as.vector(table(x$Particle.type)),datatype="abundance")$SC ) )
chao <- unlist(lapply(datasets, function(x) ChaoRichness(as.vector(table(x$Particle.type)),datatype="abundance")$Estimator ) )#asymptotic chao estimate of richness

datnames <- c("LusherFish","LusherSediment","LusherMussels","BrateMussels","Primpke1","Primpke2","Primpke3","Raingarden1","Raingarden2","Raingarden3")

#get number of particles in each dataset
nrows<- unlist(lapply(datasets,nrow))
#get max number of particles in any category for stratified sampling (e.g. in shpmax, most particles are fibers, #fibers)
shpmax <- unlist(lapply(datasets, function(z) max(table(z$Morphology))) )
colmax <-  unlist(lapply(datasets, function(z) max(table(z$Color))) )
shpcolmax <-  unlist(lapply(datasets, function(z) max(table(paste(z$Morphology,z$Color)))) )

#descriptive table of coverage and richness
DescStats <- data.frame(nrows,coverage, 1-coverage, richness, chao, chao/richness, chao - richness)
DescStats <- DescStats[order(nrows),]
colnames(DescStats) <- c("SampleSize", "Coverage", "1MinusCoverage","ObservedTypes","ChaoEstimator","Chao/ObsTypes","ChaoMinusObsTypes")
write.csv(DescStats,"~/Sample-microplastics/descriptive stats_narrow_onejar.csv")



#next 9 lines take a long time to run, can skip ahead to load pre-run results
#take 200 replicate subsamples at every single subsample size, compute success metrics
Rand <- 	lapply(1:length(datasets), function(x) lapply(1:(nrows[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) sampleplastics(datasets[[x]],smplsz,5,5,random=TRUE) ) )) ))# sampling randomly, note the scoreable type column number is not used and can be any number
save(Rand,file="~/Sample-microplastics/Rand_i.Rdata")
Shp <- 		lapply(1:length(datasets), function(x) lapply(1:(shpmax[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) sampleplastics(datasets[[x]],smplsz,4,5,random=FALSE) ) )) ))# sampling by 4th column (Morphology)
save(Shp,file="~/Sample-microplastics/Shp_i.Rdata")
Col <- 		lapply(1:length(datasets), function(x) lapply(1:(colmax[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) sampleplastics(datasets[[x]],smplsz,3,5,random=FALSE) ) )) ))# sampling by 3rd column (Color)
save(Col,file="~/Sample-microplastics/Col_i.Rdata")
ShpCol<- 	lapply(1:length(datasets), function(x) lapply(1:(shpcolmax[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) sampleplastics(datasets[[x]],smplsz,7,5,random=FALSE) ) )) ))# sampling by 6th column (ColorMorph)
save(ShpCol,file="~/Sample-microplastics/ShpCol_i.Rdata")

#load pre-run (or re-run in previous lines) results
load("~/Sample-microplastics/Rand_i.Rdata")
load("~/Sample-microplastics/Col_i.Rdata")
load("~/Sample-microplastics/Shp_i.Rdata")
load("~/Sample-microplastics/ShpCol_i.Rdata")
Randp <- div.process(Rand, nrows=nrows)
Colp <- div.process(Col, nrows=colmax) 
Shpp <- div.process(Shp, nrows=shpmax)
ShpColp <- div.process(ShpCol, nrows=shpcolmax)

###minimum subsampling calculations from which to draw recommendations - summed error
###uses a 20% threshold, 0.2 as a proportion
minimums <- matrix(ncol=8,nrow=length(datasets))
for(i in 1:length(datasets)) {
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
minimums <- cbind( data.frame(dataset = datnames, totparticles = nrows), minimums)
write.csv(minimums,"~/Sample-microplastics/minimums_SummedErrorLessThan20_narrow_onejar.csv",row.names=F)
#

###minimum subsampling calculations from which to draw recommendations - compare observed material types to those observed in the full sample
#uses a 20% threshold (i.e. so coverage should be 80%, or 0.8 proportion)
minimumsR <- matrix(ncol=8,nrow=length(datasets))
for(i in 1:length(datasets)) {
#[2,] >= richness[i]
#[1,] >= richness[i]*.8
	minimumsR[i,] <- c(Randp$N[[i]] [min(which(  Randp$Sobs.hpdi[[i]][1,] >= 0.8*richness[i]))],   
					   min(which(  Randp$Sobs.hpdi[[i]][1,] >= 0.8*richness[i])),   
					   Colp$N[[i]] [min(which(  Colp$Sobs.hpdi[[i]][1,]>= 0.8*richness[i]))],   
					   min(which(  Colp$Sobs.hpdi[[i]][1,]>= 0.8*richness[i])),   
					   Shpp$N[[i]] [min(which(  Shpp$Sobs.hpdi[[i]][1,]>= 0.8*richness[i]))],   
					   min(which(  Shpp$Sobs.hpdi[[i]][1,]>= 0.8*richness[i])),   
					   ShpColp$N[[i]] [min(which(  ShpColp$Sobs.hpdi[[i]][1,]>= 0.8*richness[i]))],
					   min(which(  ShpColp$Sobs.hpdi[[i]][1,]>= 0.8*richness[i]))   )   
					}   
minimumsR <- data.frame(minimumsR)
colnames(minimumsR) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimumsR<- cbind( data.frame(dataset = datnames, totparticles = nrows), minimumsR)
write.csv(minimumsR,"~/Sample-microplastics/minimums_EstSampleRichnessGreaterThan80_narrow_onejar.csv",row.names=F)

###minimum subsampling calculations from which to draw recommendations - compare observed material types to those expected in the environment
minimumsAE <- matrix(ncol=8,nrow=length(datasets))
for(i in 1:length(datasets)) {
	minimumsAE[i,] <- c( Randp$N[[i]] [min(which(  Randp$asychR.hpdi[[i]][1,] >= 0.8*chao[i]))],   
					   min(which(  Randp$asychR.hpdi[[i]][1,] >= 0.8*chao[i])),   
					   Colp$N[[i]] [min(which(  Colp$asychR.hpdi[[i]][1,]>= 0.8*chao[i]))],
					   min(which(  Colp$asychR.hpdi[[i]][1,]>= 0.8*chao[i])),
					   Shpp$N[[i]] [min(which(  Shpp$asychR.hpdi[[i]][1,]>= 0.8*chao[i]))],
					   min(which(  Shpp$asychR.hpdi[[i]][1,]>= 0.8*chao[i])),
					   ShpColp$N[[i]] [min(which(  ShpColp$asychR.hpdi[[i]][1,]>= 0.8*chao[i]))],
					   min(which(  ShpColp$asychR.hpdi[[i]][1,]>= 0.8*chao[i])))   
					}   
minimumsAE <- data.frame(minimumsAE)
colnames(minimumsAE) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimumsAE <- cbind( data.frame(dataset = datnames, totparticles = nrows), minimumsAE)
write.csv(minimumsAE,"~/Sample-microplastics/minimums_EstEnvRichnessGreaterThan80_narrow_onejar.csv",row.names=F)


###minimum subsampling calculations from which to draw recommendations - estimate the percent of particles in the environment that are the same material type as a (sub)sampled particle 
minimumsC <- matrix(ncol=8,nrow=length(datasets))
for(i in 1:length(datasets)) {
	minimumsC[i,] <- c( Randp$N[[i]] [min(which(  Randp$cov.hpdi[[i]][1,-1] >= 0.8))+1 ],    #coverage cannot be estimated from a sample size 1, must remove first obs 
					   min(which(  Randp$cov.hpdi[[i]][1,-1] >= 0.8))+1,    #coverage cannot be estimated from a sample size 1, must remove first obs 
					   Colp$N[[i]] [min(which(  Colp$cov.hpdi[[i]][1,]>= 0.8))],
					   min(which(  Colp$cov.hpdi[[i]][1,]>= 0.8)),
					   Shpp$N[[i]] [min(which(  Shpp$cov.hpdi[[i]][1,]>= 0.8))],
					   min(which(  Shpp$cov.hpdi[[i]][1,]>= 0.8)),
					   ShpColp$N[[i]] [min(which(  ShpColp$cov.hpdi[[i]][1,]>= 0.8))],
					   min(which(  ShpColp$cov.hpdi[[i]][1,]>= 0.8)))   
					}   
minimumsC <- data.frame(minimumsC)
colnames(minimumsC) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimumsC <- cbind( data.frame(dataset = datnames, totparticles = nrows), minimumsC)
write.csv(minimumsC,"~/Sample-microplastics/minimums_EstEnvCoverageGreaterThan80_narrow_onejar.csv",row.names=F)


##

####Make Figures

bysize <- order(nrows)

###this figure graphically show results for the summed error metric - for how well the proportion of types in the subsample represents the proportion in the sample
pdf("~/Sample-microplastics/SummedError3Panel_narrow_onejar.pdf",height=3,width=6)
layout(matrix(1:3,ncol=3,byrow=T))
par(oma=c(1,4,1,2))
par(mar=c(3,0,2,0))
for(i in c(7,8,2)){ #datasets are: Primpke 3 raingarden 1, Lusher sediment
plot(1~1, pch=NA, xlim=c(0,nrows[i]),ylim=c(0,1),ylab="",xlab="",yaxt="n")
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
	mtext(data_titles[i], side=3, line =0.25, cex=1)
	if(i==7){mtext("Summed error",line=2.5,side=2)}
	if(i==7){axis(side=2)}
	if(i==2){legend(175,y=1,c("Rand","Col","Morph","MorphCol"),fill=c(rgb(0,0,0),rgb(1,0,0),rgb(0,0,1),rgb(1,0,1)),bty="n")}
	if(i==2){legend(175,y=0.65,c("20%Error"),lty=2,bty="n")}
 	if(i==8){mtext("Particles chemically identified",line=2.5,side=1)}
}
dev.off()
#this figure uses only 3 datasets
#see "sample plastics grouped figures.R" for figure of all datasets used in supplemental figure

#see "sample plastics grouped figures.R" for figures depicting how well subsamples capture the number of particle types in the sample
#and how well sub/samples represent particles in the environment



