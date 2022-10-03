#goals:
#estimating %particle types in full sample from subsample
#this is for datasets with a single sample ("one jar")
#success is getting same % in each category as when estimated by full sample 
#categories here are broad material types ("btype"), meaning split into plastic, anthropogenic, and natural
#check success of each subsample method across datasets

library(coda) #hpdi
library(iNEXT) #diversity metrics

range01=function(x){
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
}

dropNA <- function(vect){return(vect[!is.na(vect)])}

#read in functions for sampling plastics
source('~/Sample-microplastics/sample\ plastics\ functions.R') 
#source location depends on local file structure! might need to change path, or set working directory


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
plasticsVT1 <- read.csv("~/Sample-microplastics/Vidal_Tsukada_beach1.csv",header=T,stringsAsFactors=F) 
plasticsVT2 <- read.csv("~/Sample-microplastics/Vidal_Tsukada_beach2.csv",header=T,stringsAsFactors=F) 

datasets <- list(plasticsLF,plasticsLS,plasticsLM,plasticsBM,plasticsP1,plasticsP2,plasticsP3,plasticsRG1,plasticsRG2,plasticsRG3,plasticsVT1,plasticsVT2)

for(i in 1:length(datasets)){datasets[[i]]$ColorMorph <- paste(datasets[[i]]$Color, datasets[[i]]$Morphology,sep="") } #create a pasted color and morphology category

data_titles <- c("Animal tissue","Sediment", "Animal tissue","Animal tissue", "Surface water", "Surface water", "Surface water", "Storm water", "Storm water","Storm water","Sediment","Sediment")
datnames <- c("LusherFish","LusherSediment","LusherMussels","BrateMussels","Primpke1","Primpke2","Primpke3","Raingarden1","Raingarden2","Raingarden3","VidalBeach1","VidalBeach2")

#get number of particles in each dataset
nrows<- unlist(lapply(datasets,nrow))
#get max number of particles in any category for stratified sampling (e.g. in shpmax, most particles are fibers, #fibers)
shpmax <- unlist(lapply(datasets, function(z) max(table(z$Morphology))) )
colmax <-  unlist(lapply(datasets, function(z) max(table(z$Color))) )
shpcolmax <-  unlist(lapply(datasets, function(z) max(table(paste(z$Morphology,z$Color)))) )

#next 9 lines take a long time to run, can skip ahead to load pre-run results
#take 200 replicate subsamples at every single subsample size, compute success metrics
Randb <- 	lapply(1:length(datasets), function(x) lapply(1:(nrows[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) sampleplastics(datasets[[x]],smplsz,5,6,random=TRUE) ) )) ))# sampling randomly, note the scoreable type column number is not used when random and can be any number other than the true type column number
save(Randb,file="~/Sample-microplastics/Rand_i_b.Rdata")
Shpb <- 		lapply(1:length(datasets), function(x) lapply(1:(shpmax[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) sampleplastics(datasets[[x]],smplsz,4,6,random=FALSE) ) )) ))# sampling by 4th column (Morphology)
save(Shpb,file="~/Sample-microplastics/Shp_i_b.Rdata")
Colb <- 		lapply(1:length(datasets), function(x) lapply(1:(colmax[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) sampleplastics(datasets[[x]],smplsz,3,6,random=FALSE) ) )) ))# sampling by 3rd column (Color)
save(Colb,file="~/Sample-microplastics/Col_i_b.Rdata")
ShpColb<- 	lapply(1:length(datasets), function(x) lapply(1:(shpcolmax[x]), function(smplsz) as.data.frame(t(sapply(1:200, function(z) sampleplastics(datasets[[x]],smplsz,7,6,random=FALSE) ) )) ))# sampling by 6th column (ColorMorph)
save(ShpColb,file="~/Sample-microplastics/ShpCol_i_b.Rdata")

#load pre-run results (or re-run in previous lines)
load("~/Sample-microplastics/Rand_i_b.Rdata")
load("~/Sample-microplastics/Col_i_b.Rdata")
load("~/Sample-microplastics/Shp_i_b.Rdata")
load("~/Sample-microplastics/ShpCol_i_b.Rdata")
Randp <- div.process(Randb, nrows=nrows)
Colp <- div.process(Colb, nrows=colmax) 
Shpp <- div.process(Shpb, nrows=shpmax)
ShpColp <- div.process(ShpColb, nrows=shpcolmax)

###minimum subsampling calculations from which to draw recommendations
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
					   min(which(  ShpColp$smerr.hpdi[[i]][2,]<0.2))  
	  )
					}   
minimums <- data.frame(minimums)
colnames(minimums) <- c("randnum","randper","colnum","colper","shpnum","shpper","shpcolnum","shpcolper")
minimums <- cbind( data.frame(dataset = datnames, totparticles = nrows), minimums)
write.csv(minimums,"~/Sample-microplastics/minimums_SummedErrorLessThan20_broad_onejar.csv",row.names=F)
####Make Figures

bysize <- order(nrows)

#figures
pdf("~/Sample-microplastics/SummedError3panel_broad_onejar.pdf",height=3,width=6)
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

#see "sample plastics grouped figures.R" for figure of all datasets used in supplemental figure

