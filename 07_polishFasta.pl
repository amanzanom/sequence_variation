#!/usr/bin/perl

use Getopt::Long;
use Data::Dumper;

GetOptions(\%opts, "ref=s", "polish|p=s", "verbose|v!", "help|h!");
        sub usage(){
                die "USAGE :: perl polishFASTA.pl -ref reference.fasta -polish|p polisher.corrections [-help|-h]\n\n
	-ref\n\tFASTA containing contig reference sequence(s) [String]\n\n
	-polish|p\n\tPolisher corrections file [String]\n\n
	-verbose|v\n\tSet verbose [Boolean]\n\n
	-help or -help\n\tGet help for using this script [Boolean]\n\n";
        }

if ($opts{'help'} || !$opts{'ref'} || !$opts{'polish'}){
	if (!$opts{'ref'}){
		print "FASTA reference file missing, please check usage manual\n\n";
	}
	if (!$opts{'polish'}){
		print "Polisher corrections file not specified, check usage manual\n\n";
	}
	&usage;
}

# Define subroutines

sub printFastaFile {
	my $outFile_FH= $_[0];
	my $header= $_[1];
	my $seq= $_[2];
	my $charNum=$_[3];
	my $i=0;	
	print $outFile_FH ">$header\n";
	for ($i=0; $i<length($$seq); $i+=$charNum){
		if (length(substr($$seq, $i, length($$seq)-$i))< $charNum){
			print $outFile_FH substr($$seq, $i, length($$seq)-$i)."\n";
		}
		else {
			print $outFile_FH substr($$seq, $i, $charNum)."\n";
		}
	}
	return (1);
}

# Capture options
my $refFasta=$opts{'ref'};
my $polishFile=$opts{'polish'};
my $verbose=$opts{'verbose'};

# Define variables to use in the script
my $charPerLine=70;
my %refCorrections=();
my @tempArray='';
my $line='';
my $seqId='';
my $tempSeq='';
my $newSeq='';
my $i=0;

# Read fasta file and write corrected fasta
open (POLISH, "$polishFile") || die "Unable to open file for reading $polishFile\n$!\n";
while ($line=<POLISH>){
	if ($line && $line!~ m/^#/){
		@tempArray= split (/\t/, $line);
		$tempArray[6]=~ m/^(\w)\w+\:(\w*)\:/;
		$refCorrections{$tempArray[0]}{$tempArray[3]}=$1.$2;
	}
}
close (POLISH);

# Read fasta file and write corrected fasta
#open (POLISHFASTA, ">$refFasta.corrected") || die "Unable to open file for writting $refFasta.corrected\n$!\n";
open (REFFASTA, "$refFasta") || die "Unable to open file for reading $refFasta\n$!\n";
while ($line=<REFFASTA>){
	chomp $line;
	if ($line=~ m/^>/ || eof){
		if (eof){
			$tempSeq.= $line;
		}
		if ($tempSeq){
			if ($verbose){
				print ">Correcting sequence $seqId\n";
			}
			$newSeq="";
			for ($i=0; $i< length($tempSeq); $i++){
				if ($refCorrections{$seqId}{$i+1}){
					if ($refCorrections{$seqId}{$i+1}=~ m/^S(\w+)/){ # Case substitution
						$newSeq.= $1;
						if ($verbose){
							print "Substituted position " . ($i+1) . " from " . substr($tempSeq, $i, 1) . " to $1\n";
						}
					}
					elsif ($refCorrections{$seqId}{$i+1}=~ m/^I(\w+)/){ # Case insertion
						$newSeq.= substr($tempSeq, $i, 1) . $1;
						if ($verbose){
							print "Inserted string $1 in position " . ($i+1) . " after " . substr($tempSeq, $i, 1) . "\n";
						}
					}
				}
				else {
					$newSeq.= substr($tempSeq, $i, 1);
				}
			}
			$newSeq=uc($newSeq);
			&printFastaFile(\*STDOUT, $seqId, \$newSeq, $charPerLine);
		}
		$line=~ m/^>(\S+)/;
		$seqId= $1;
		$tempSeq= "";
		next;
	}
	$tempSeq.= $line;
}
close (REFFASTA);
#close (OUT);

