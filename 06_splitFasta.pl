#!/usr/local/bin/perl

#########################################################
# Perl script that takes a list of contigs produced by  #
# either wgs-assembler of MIRA and translates it to     #
# contigs with the corresponding assemlber outfiles     #
# describing this relations (MIRA *_contigreadlist.txt) #
# (wgs-assembler *.posmap.frgctg)			#
#########################################################


use Getopt::Long;
use Data::Dumper;

GetOptions(\%opts, "infile|i=s", "parts|p=i", "help|h!");
        sub usage(){
                die "USAGE :: perl analyze_seq.pl -i file -image name [-s|-summary] [-gc] [-p] [-v|-verbose] [-help|-h]\n\n
	-infile or -i\n\tFile(s) containing contig(s) to be converted to reads [String]\n\n
	-parts or -p\n\tNumber of files to split the inFasta to [Integer]\n\n
	-help or -help\n\tGet help for using this script\n\n";
        }
        
#main program
{

	if ($opts{'help'} || !$opts{'infile'}){
		if (!$opts{'infile'}){
			print STDERR "Arguments -i is missing please check the usage manual\n\n";
		}
		&usage;
	}
	$parts=$opts{'parts'};
	@splittedFiles=();
	$c="1";
	$numRead=0;
	open (IN, "<$opts{'infile'}") || die "Unable to open file for reading: $opts{'infile'}\n$!\n";
	$cmd='grep -c -P "^>" ' . $opts{'infile'};
	print STDERR "Running command $cmd\n";
	$numSeq=`$cmd`;
	chomp $numSeq;
	print STDERR "Found $numSeq sequences\n";
	$numSeq/=$parts;
	if ($numSeq=~ m/\./){
		$numSeq=~ s/\..+$//g;
		$numSeq++;
	}
	$fileNum=1;
	$splitFile=$opts{'infile'};
	$splitFile=~ s/\./_$fileNum\./;
	print STDERR "Working on file $splitFile\n";
	push(@splittedFiles, $splitFile);
	open (SPLITFILE, ">$splitFile") || die "Unable to open file for writting: $splitFile\n$!\n";
	while ($line=<IN>){
		if ($line=~/^>/){
			$numRead++;
			if ($numRead / ($numSeq+1) == 1){
				close (SPLITFILE);
				$numRead=1;
				$fileNum++;
				$splitFile=~ s/_\d+\./_$fileNum\./;
				print STDERR "Working on file $splitFile\n";
				push(@splittedFiles, $splitFile);
				open (SPLITFILE, ">$splitFile") || die "Unable to open file for reading: $splitFile\n$!\n";
			}
		}
		print SPLITFILE $line;
	}
	close (SPLITFILE);
	print STDERR "Finished generating splitted input\n";
	close (IN);		
}
