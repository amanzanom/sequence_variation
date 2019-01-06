## input file from sam file: 
## cut -f 9 readsGood_PE_vs_backbone_pre.sam | grep -P "\d+" | grep -v "^-" > insSize.txt

inFile<-as.character(commandArgs(trailingOnly=T)[1]);
outPrefix<-as.character(commandArgs(trailingOnly=T)[2]);

insSize<- read.table(inFile, header=FALSE, sep="\t", colClasses=c("numeric"));

## This step is required to avoid duplicate zero-sized insert sizes ##
insSizeZero<- insSize[insSize==0];
insSizeZero<- insSizeZero[1:(length(insSizeZero)/2)];
## This step is required to avoid negative insert sizes ##
insSizeNonZero<- insSize[insSize>0];

insSize<- c(insSizeZero[!is.na(insSizeZero)], insSizeNonZero);

insSizeMean<- round(mean(insSize), digits=2);
insSizeMedian<- round(median(insSize), digits=2);
insSizeSd<- round(sd(insSize), digits=2);
insSizeMin<- min(insSize);
insSizeMax<- max(insSize);

pdf(file=paste(outPrefix, '.pdf', sep=''));
	hist(insSize, main="Paired-Ends Insert Size", xlab="Insert Size");
	abline(v=insSizeMean, col="red", lty=4, lwd=2);
	abline(v=insSizeMean+insSizeSd*-1, col="blue", lty=4, lwd=2);
	abline(v=insSizeMean+insSizeSd, col="blue", lty=4, lwd=2);
	abline(v=insSizeMean+insSizeSd*-2, col="grey", lty=4, lwd=2);
	abline(v=insSizeMean+insSizeSd*2, col="grey", lty=4, lwd=2);
	legendText <- c(paste("Mean: ", insSizeMean, sep=''), paste("Median: ", insSizeMedian, sep=''), paste("StDev: ", insSizeSd, sep= ''), paste("Min: ", insSizeMin, sep= ''), paste("Max: ", insSizeMax, sep= ''));
	legend(x="topright", legendText, cex=1, bty="n");
dev.off();

