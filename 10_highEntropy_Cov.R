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
aComp<- NULL;
tComp<- NULL;
gComp<- NULL;
cComp<- NULL;
gapComp<- NULL;
posShannon<- NULL;
meanShannon<- NULL;
covPerSite<- NULL;


# Main Program

## Read input file and window size
inFile<-as.character(commandArgs(trailingOnly=T)[1]);
prefix<-as.character(commandArgs(trailingOnly=T)[2]);
wSize<-as.numeric(commandArgs(trailingOnly=T)[3]);

nucAbsFreq<- read.table(inFile, header=TRUE, row.names=1, sep="\t");

## Calculate shannon entropy per site and store in a vector
for (i in 1:nrow(nucAbsFreq)){
	sumNuc<- nucAbsFreq$A[i] + nucAbsFreq$T[i] + nucAbsFreq$G[i] + nucAbsFreq$C[i] + nucAbsFreq$GAP[i];
	if (nucAbsFreq$A[i] == 0){
		aComp<- 0;
	}
	else {
		aComp<- (nucAbsFreq$A[i]/sumNuc)*log((nucAbsFreq$A[i]/sumNuc), 10);
	}
	if (nucAbsFreq$T[i] == 0){
		tComp<- 0;
	}
	else {
		tComp<- (nucAbsFreq$T[i]/sumNuc)*log((nucAbsFreq$T[i]/sumNuc), 10);
	}
	if (nucAbsFreq$G[i] == 0){
		gComp<- 0;
	}
	else {
		gComp<- (nucAbsFreq$G[i]/sumNuc)*log((nucAbsFreq$G[i]/sumNuc), 10);
	}
	if (nucAbsFreq$C[i] == 0){
		cComp<- 0;
	}
	else {
		cComp<- (nucAbsFreq$C[i]/sumNuc)*log((nucAbsFreq$C[i]/sumNuc), 10);
	}
	if (nucAbsFreq$GAP[i] == 0){
		gapComp<- 0;
	}
	else {
		gapComp<- (nucAbsFreq$GAP[i]/sumNuc)*log((nucAbsFreq$GAP[i]/sumNuc), 10);
	}
	posShannon[i]<- (aComp + tComp + gComp + cComp + gapComp)*-1;
	covPerSite[i]<- sumNuc;
}

#posShannon<- as.numeric(posShannon);

## Make matrix with x and y values for means of shannon entropy per overlapping wSize windows
#meanShannon<- matrix(nrow=(length(posShannon)-wSize+1), ncol=2);
#for (i in 1:(length(posShannon)-wSize+1)){
#	meanShannon[i,1]<- i+(wSize/2);
#	meanShannon[i,2]<- mean(posShannon[i:(i+wSize-1)]);
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
meanShannon<- matrix(nrow=length(posShannon), ncol=4);
for (i in 1:length(posShannon)){
	meanShannon[i,1]<- i;
	meanShannon[i,2]<- posShannon[i];
	Zscore<- (posShannon[i]-mean(posShannon))/sd(posShannon);
	meanShannon[i,3]<- Zscore;
	meanShannon[i,4]<- 2*pnorm(-abs(Zscore));;
}

#for (i in 1:nrow(meanCov)){
#	Zscore<- (meanCov[i,2]-mean(meanCov[,2]))/sd(meanCov[,2]);
#	meanCov[i,3]<- Zscore;
#	meanCov[i,4]<- 2*pnorm(-abs(Zscore));
#}

write.table(meanShannon[meanShannon[,4]<.05 & meanShannon[,3]>0,], file=paste(prefix, '_signShannon_up.txt', sep=''), sep="\t", row.names=FALSE, col.names=FALSE);
#write.table(meanCov[meanCov[,4]<.05 & meanCov[,3]<0,], file=paste(prefix, '_signShannon_down.txt', sep=''), sep="\t", row.names=FALSE, col.names=FALSE);


### Make barplot
#pdf(file=paste(prefix, ".pdf", sep=""), paper="a4r");
#	plot(rownames(nucAbsFreq), posShannon, main=paste("Shannon Entropy in", prefix, sep=" "), xlab="Sequence Position", ylab="entropy", type="h", bty="n");
#	lines(which(posShannon == 0), rep.int(max(posShannon), length(which(posShannon == 0))), type="h", col=col2alpha("yellow", 0.05)); # This line gives color to zero-coverage zones
#	lines(meanShannon, type="l", col="red");
#	axis(side=1, at=axisTicks(c(0, length(posShannon)), log=FALSE, nint = 10), labels=FALSE, tick=TRUE);
#dev.off();

