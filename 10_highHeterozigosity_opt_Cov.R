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

## Calculate Heterozigosity per site and store in a vector, skip all non-covered positions and all states with a ratio less than .01 
for (i in 1:nrow(nucAbsFreq)){
	## Avoid all non-covered sites
	if (sum(nucAbsFreq[i, ])!=0){
		## Avoid all states with a ratio less than .01
		currPosCharAbsFreq<- subset(nucAbsFreq[i, ], subset=TRUE, prop.table(nucAbsFreq[i, ])>=.01);
		posHeterozigosity[i]<- 0;
		for (j in 1:length(currPosCharAbsFreq)){
			posHeterozigosity[i]<- posHeterozigosity[i] + (currPosCharAbsFreq[j]/sum(currPosCharAbsFreq))^2
		}
		mHeterozigosity[i]<- length(currPosCharAbsFreq);
		posHeterozigosity[i]<- 1 - posHeterozigosity[i];
		covPerSite[i]<- sum(currPosCharAbsFreq);
	}
	else {
		posHeterozigosity[i]<- NA;
		mHeterozigosity[i]<- NA;
		posHeterozigosity[i]<- NA;
		covPerSite[i]<- 0;
	}
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

write.table(round(meanCov[meanCov[,4]<.05 & meanCov[,3]>0,], digits=4), file=paste(prefix, '_signDiffCov_up.txt', sep=''), sep="\t", row.names=FALSE, col.names=FALSE);
write.table(round(meanCov[meanCov[,4]<.05 & meanCov[,3]<0,], digits=4), file=paste(prefix, '_signDiffCov_down.txt', sep=''), sep="\t", row.names=FALSE, col.names=FALSE);

## Shannon entropy statistically significant higher entropy
scoreHeterozigosity<- matrix(nrow=length(posHeterozigosity), ncol=4);
for (i in 1:length(posHeterozigosity)){
	if (!is.na(posHeterozigosity[i])){
		scoreHeterozigosity[i,1]<- i;
		scoreHeterozigosity[i,2]<- posHeterozigosity[i];
		Zscore<- (posHeterozigosity[i]-mean(posHeterozigosity, na.rm=TRUE))/sd(posHeterozigosity, na.rm=TRUE);
		scoreHeterozigosity[i,3]<- Zscore;
		scoreHeterozigosity[i,4]<- 2*pnorm(-abs(Zscore));
	}
	else {
		scoreHeterozigosity[i]<- c(i, NA, NA, NA);
	}
}

#for (i in 1:nrow(meanCov)){
#	Zscore<- (meanCov[i,2]-mean(meanCov[,2]))/sd(meanCov[,2]);
#	meanCov[i,3]<- Zscore;
#	meanCov[i,4]<- 2*pnorm(-abs(Zscore));
#}

write.table(round(scoreHeterozigosity, digits=4), file=paste(prefix, '_signHeterozigosity.txt', sep=''), sep="\t", row.names=FALSE, col.names=FALSE);
write.table(round(subset(scoreHeterozigosity, scoreHeterozigosity[, 4] < .05 & scoreHeterozigosity[, 3] > 0, na.rm=TRUE), digits=4), file=paste(prefix, '_signHeterozigosity_up.txt', sep=''), sep="\t", row.names=FALSE, col.names=FALSE);
write.table(round(subset(scoreHeterozigosity, scoreHeterozigosity[, 4] < .05 & scoreHeterozigosity[, 3] < 0, na.rm=TRUE), digits=4), file=paste(prefix, '_signHeterozigosity_down.txt', sep=''), sep="\t", row.names=FALSE, col.names=FALSE);

plotLegend<- c(paste("mean states=", round(mean(mHeterozigosity, na.rm=TRUE), digits=2), sep=" "), paste("sd states=", round(sd(mHeterozigosity, na.rm=TRUE), digits=2), sep=" "), paste("mean H=", round(mean(posHeterozigosity, na.rm=TRUE), digits=2), sep=" "), paste("sd H=", round(sd(posHeterozigosity, na.rm=TRUE), digits=2), sep=" "));
## Make barplot
pdf(file=paste(prefix, ".pdf", sep=""));
	plot(posHeterozigosity, main=paste("Heterozigosity per site in", inFile, sep=" "), xlim=c(0,9677), ylim=c(0,0.5), xlab="Sequence Position", ylab="Heterozigosity", type="h", bty="n");
	text(subset(scoreHeterozigosity, scoreHeterozigosity[, 4]< .05  & scoreHeterozigosity[, 3] > 0, na.rm=TRUE)[,1], y=subset(scoreHeterozigosity, scoreHeterozigosity[, 4]< .05, na.rm=TRUE)[,2], labels="*", col="red");
	text(subset(scoreHeterozigosity, scoreHeterozigosity[, 4]< .05  & scoreHeterozigosity[, 3] < 0, na.rm=TRUE)[,1], y=subset(scoreHeterozigosity, scoreHeterozigosity[, 4]< .05, na.rm=TRUE)[,2], labels="*", col="blue");
	axis(side=1, at=axisTicks(c(0, length(posHeterozigosity)), log=FALSE, nint = 10), labels=FALSE, tick=TRUE);
	legend("topright", plotLegend);
dev.off();

quit(status=0);

