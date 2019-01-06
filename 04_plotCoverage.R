#===============================================================================
#   Author: Alejandro Manzano-Marin,
#
#   File: plotCoverage
#   Date: 25-11-2012
#   Version: 0.1
#
#   Usage:
#      R --args samSummary.absFreq windowSize < plotShannonSeq.R
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

absFreq<- read.table(inFile, header=TRUE, row.names=1, sep="\t", comment.char='#');
abc<- colnames(absFreq);

## Calculate coverage per-site
for (i in 1:nrow(absFreq)){
	covPerSite[i]<- 0;
	for (j in 1:length(abc)){
		covPerSite[i]<- covPerSite[i] + absFreq[i,j];
	}
}

meanCovPerSite<- round(mean(covPerSite[covPerSite!=0]), digits=2);
medianCovPerSite<- round(median(covPerSite[covPerSite!=0]), digits=2);
sdCovPerSite<- round(sd(covPerSite[covPerSite!=0]), digits=2);
minCovPerSite<- min(covPerSite[covPerSite!=0]);
maxCovPerSite<- max(covPerSite[covPerSite!=0]);


## Plot Coverage
pdf(file=paste(prefix, ".pdf", sep=''));
	hist(covPerSite, main="Per-site mapping coverage", xlab="Coverage");
	abline(v=meanCovPerSite, col="red", lty=4, lwd=2);
	abline(v=meanCovPerSite+sdCovPerSite*-1, col="blue", lty=4, lwd=2);
	abline(v=meanCovPerSite+sdCovPerSite, col="blue", lty=4, lwd=2);
	abline(v=meanCovPerSite+sdCovPerSite*-2, col="grey", lty=4, lwd=2);
	abline(v=meanCovPerSite+sdCovPerSite*2, col="grey", lty=4, lwd=2);
	legendText <- c(paste("Mean: ", meanCovPerSite, sep=''), paste("Median: ", medianCovPerSite, sep=''), paste("StDev: ", sdCovPerSite, sep= ''), paste("Min: ", minCovPerSite, sep= ''), paste("Max: ", maxCovPerSite, sep= ''));
	legend(x="topright", legendText, cex=1, bty="n");
dev.off();

