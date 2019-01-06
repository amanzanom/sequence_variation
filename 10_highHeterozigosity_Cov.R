#===============================================================================
#   Author: Alejandro Manzano-Marin,
#
#   File: plotShannonSeq
#   Date: 25-11-2012
#   Version: 0.1
#
#   Usage:
#      R --args shannonSAM_out_refID.relNucFreq windowSize < plotShannonSeq.R
#
#    Description: plotShannonSeq is a script designed to calculate the Shannon entropy for positions in a SAM
#                 alignment considering A,T,G,C and - extracted using the shannonSAM.pl perl script. The
#                 calculations are done as described in Zagordi et al Read length versus Depth Coverage for 
#                 Viral Quasispecies Recnstruction. PLoS ONE 7(10) e47046.
#
#    Contact: Contact the author at alejandro.manzano@uv.es using 'shannonSAM: ' as
#             as begining for subject for bug reporting, feature request or whatever
#             (of course realted to the software).
#
#===============================================================================

# Add libraries
library(seqinr);

# Declare variables

## Argument variables
inFile<- NULL;
prefix<- NULL;
wSize<- NULL;

## Other Variables
nucAbsFreq<- NULL;
mHeterozigosity<- NULL;
posHeterozigosity<- NULL;
scoreHeterozigosity<- NULL;
covPerSite<- NULL;


# Main Program

## Read input file and window size
inFile<-as.character(commandArgs(trailingOnly=T)[1]);
prefix<-as.character(commandArgs(trailingOnly=T)[2]);
wSize<-as.numeric(commandArgs(trailingOnly=T)[3]);

nucAbsFreq<- read.table(inFile, header=TRUE, row.names=1, sep="\t");

## Calculate shannon entropy per site and store in a vector
for (i in 1:nrow(nucAbsFreq)){
	posHeterozigosity[i]<- 0;
	for (j in 1:length(nucAbsFreq[i, ])){
#		if (sum(nucAbsFreq[i, ])>0){
			posHeterozigosity[i]<- posHeterozigosity[i] + (nucAbsFreq[i, j]/sum(nucAbsFreq[i, ]))^2
#		}
#		else {
#			posHeterozigosity[i]<- 0;
#		}
	}
	mHeterozigosity[i]<- length(nucAbsFreq[i, nucAbsFreq[i, ]>0]);
	posHeterozigosity[i]<- 1 - posHeterozigosity[i];
	covPerSite[i]<- sum(nucAbsFreq[i, ]);
}

#posHeterozigosity<- as.numeric(posHeterozigosity);

## Make matrix with x and y values for means of shannon entropy per overlapping wSize windows
#scoreHeterozigosity<- matrix(nrow=(length(posHeterozigosity)-wSize+1), ncol=2);
#for (i in 1:(length(posHeterozigosity)-wSize+1)){
#	scoreHeterozigosity[i,1]<- i+(wSize/2);
#	scoreHeterozigosity[i,2]<- mean(posHeterozigosity[i:(i+wSize-1)]);
#}

## Make matrix with x and y values for means of coverage per overlapping wSize windows
meanCov<- matrix(nrow=(length(covPerSite)-wSize+1), ncol=4);
for (i in 1:(length(covPerSite)-wSize+1)){
	meanCov[i,1]<- i+(wSize/2);
	meanCov[i,2]<- mean(covPerSite[i:(i+wSize-1)]);
}

for (i in 1:nrow(meanCov)){
	Zscore<- (meanCov[i,2]-mean(meanCov[,2]))/sd(meanCov[,2]);
	meanCov[i,3]<- Zscore;
	meanCov[i,4]<- 2*pnorm(-abs(Zscore));
}

write.table(meanCov[meanCov[,4]<.05 & meanCov[,3]>0,], file=paste(prefix, '_signDiffCov_up.txt', sep=''), sep="\t", row.names=FALSE, col.names=FALSE);
write.table(meanCov[meanCov[,4]<.05 & meanCov[,3]<0,], file=paste(prefix, '_signDiffCov_down.txt', sep=''), sep="\t", row.names=FALSE, col.names=FALSE);

## Shannon entropy statistically significant higher entropy
scoreHeterozigosity<- matrix(nrow=length(posHeterozigosity), ncol=4);
for (i in 1:length(posHeterozigosity)){
	scoreHeterozigosity[i,1]<- i;
	scoreHeterozigosity[i,2]<- posHeterozigosity[i];
	Zscore<- (posHeterozigosity[i]-mean(posHeterozigosity, na.rm=TRUE))/sd(posHeterozigosity, na.rm=TRUE);
	scoreHeterozigosity[i,3]<- Zscore;
	scoreHeterozigosity[i,4]<- 2*pnorm(-abs(Zscore));
}

#for (i in 1:nrow(meanCov)){
#	Zscore<- (meanCov[i,2]-mean(meanCov[,2]))/sd(meanCov[,2]);
#	meanCov[i,3]<- Zscore;
#	meanCov[i,4]<- 2*pnorm(-abs(Zscore));
#}

write.table(subset(scoreHeterozigosity, scoreHeterozigosity[, 4]< .05, na.rm=TRUE), file=paste(prefix, '_signHeterozigosity_up.txt', sep=''), sep="\t", row.names=FALSE, col.names=FALSE);
#write.table(meanCov[meanCov[,4]<.05 & meanCov[,3]<0,], file=paste(prefix, '_signShannon_down.txt', sep=''), sep="\t", row.names=FALSE, col.names=FALSE);
#print(scoreHeterozigosity);
#print(subset(scoreHeterozigosity, scoreHeterozigosity[, 4]< .05)[,1]);

## Make barplot
pdf(file=paste(prefix, ".pdf", sep=""));
	plot(posHeterozigosity, main=paste("Heterozigosity per site in", prefix, sep=" "), xlim=c(0,9677), ylim=c(0,0.5), xlab="Sequence Position", ylab="Heterozigosity", type="h", bty="n");
#	lines(which(posHeterozigosity == 0), rep.int(max(posHeterozigosity), length(which(posHeterozigosity == 0))), type="h", col=col2alpha("yellow", 0.05)); # This line gives color to zero-coverage zones
#	lines(subset(scoreHeterozigosity, scoreHeterozigosity[, 4]< .05)[,1], subset(scoreHeterozigosity, scoreHeterozigosity[, 4]< .05)[,2], type="p", col="red", lwd=10, cex=10, pch="*");
	
	text(subset(scoreHeterozigosity, scoreHeterozigosity[, 4]< .05, na.rm=TRUE)[,1], y=subset(scoreHeterozigosity, scoreHeterozigosity[, 4]< .05, na.rm=TRUE)[,2], labels="*", col="red");
	axis(side=1, at=axisTicks(c(0, length(posHeterozigosity)), log=FALSE, nint = 10), labels=FALSE, tick=TRUE);
	legend("topright", paste("mean states=", mean(mHeterozigosity), sep=" "));
dev.off();
pdf(file=paste(prefix, "_mHeterozigosity.pdf", sep=""), paper="a4r");
	hist(mHeterozigosity, main=paste("m param Heterozigosity per site in", prefix, sep=" "), xlab="Sequence Position", ylab="m");
dev.off();
#print (maxposHeterozigosity);

