#===============================================================================
#   Author: Alejandro Manzano Marin,
#
#   Usage:
#      R --vanilla --args post_Hzyg_signHeterozigosity.txt outPrefix < 13.averageHeterozygosity.R
#
#===============================================================================

## Argument variables
inFile <- NULL;
outPrefix <- NULL;

## Read input file
inFile <- as.character(commandArgs(trailingOnly=T)[1]);
outPrefix <- as.character(commandArgs(trailingOnly=T)[2]);

data <- read.table(inFile);

## remove rows with NA values
comp.data <-data[complete.cases(data), ];
## remove rows with H = 0
data.sub <- subset(comp.data, comp.data[,2] != 0);
## calculate average H, correcting for the number of polymorphic sites
average.h <- (1/length(data.sub[,2])) * mean(data.sub[,2]);

## write average H to file
write.table(average.h, file=outPrefix, row.names=FALSE, col.names=FALSE);

