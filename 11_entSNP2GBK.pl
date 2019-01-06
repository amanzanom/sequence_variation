#!/usr/bin/perl

use Getopt::Long;
use Data::Dumper;
use Cwd 'abs_path';


GetOptions(\%opts, "shEnt:s", "SNPSAMtools:s", "SNPVarScan:s", "help|h!");
        sub usage(){
                die "USAGE :: perl SNP2GBK.pl [--shEnt shEntout.tab] [--SNPSAMtools SAMtools.vcf] [--SNPVarScan SNPVarScan.vcf] [-help|-h]\n\n
	--shEnt\n\tHighly entropic sites tabular file [String]\n\n
	--SNPSAMtools\n\tSNP file from SAMtools [String]\n\n
	--SNPVarScan\n\t*.vcf file from VarScan [String]\n\n
	--help or -help\n\tGet help for using this script\n\n";
        }

if ($opts{'help'}){
#	if (!$opts{'infile'}){
#		print "infile missing, please check usage manual\n\n";
#	}
	&usage;
}

# Variable definition

## Capture options
my $shEntFile= $opts{'shEnt'};
my $SNP_SAMtools_File= $opts{'SNPSAMtools'};
my $SNP_VarScan_File= $opts{'SNPVarScan'};

## Other variables
my %shEntSites= ();
my %SNPSites= ();
my $line='';
my $cigarString='';
my $start=0;
my $strand=0;
my $sequence='';
my $id='';
my $skippedLeft=0;
my $skippedRight=0;


# Main program
if ($shEntFile){
	open (SHENTFILE, "<$shEntFile") or die ("Unable to open file for reading: $shEntFile\n$!\n");
	while ($line= <SHENTFILE>){
		chomp ($line);
		if ($line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)$/){
			$pos= $1;
			$shannonScore= $2;
			$zVal= $3;
			$pVal= $4;
			$shEntSites{$pos}[0]= $shannonScore;
			$shEntSites{$pos}[1]= $zVal;
			$shEntSites{$pos}[2]= $pVal;
		}
	}
	close (SHENTFILE);
	
	foreach $pos (sort {$a <=> $b} keys %shEntSites){
		print STDOUT '     misc_feature    ' . $pos . '..' . $pos . "\n";
		$shEntSites{$pos}[2]=~ m/^([\d\.]+)(\S+)*$/;;
		$num= $1;
		$extra= $2;
		print STDOUT '                     /inference="ab initio prediction: Shannon entropy normal p-val <0.05"' . "\n";
		print STDOUT '                     /note="Statistically significant high entropy site, p-val ';
		print STDOUT sprintf("%.4f", $num);
		print STDOUT $extra;
		print STDOUT '"'. "\n";
	}
}

if ($SNP_SAMtools_File){
	open (SNPFILE, "<$SNP_SAMtools_File") or die ("Unable to open file for reading: $SNP_SAMtools_File\n$!\n");
	while ($line= <SNPFILE>){
		chomp ($line);
		if ($line=~ m/^##samtoolsVersion=(\S+)/){
			$samVersion= $1;
			next;
		}
		if ($line=~ m/^\S+\t(\S+)\t\S+\t\S+\t(\S+)\t.+DP=(\d+);\S+;DP4=\d+,\d+,(\d+),(\d+)/){
			$pos= $1;
			$replace= lc($2);
			$depthHighQual= $3;
			$altFwd= $4;
			$altRev= $5;
			$SNPSites{$pos}{$replace}= ($altFwd+$altRev)/$depthHighQual;
		}
	}
	close (SNPFILE);
	
	foreach $pos (sort {$a <=> $b} keys %SNPSites){
		foreach $replace (sort {lc($a) cmp lc($b)} keys %{$SNPSites{$pos}}){
			print STDOUT '     variation       ' . $pos . '..' . $pos . "\n";		
			print STDOUT '                     /inference="ab initio prediction:mpileup:' . "$samVersion\"\n";
			print STDOUT '                     /note="SNP"' . "\n";
			print STDOUT '                     /replace="' . $replace . '"' . "\n";
			print STDOUT '                     /frequency="';
			print STDOUT sprintf("%.4f", $SNPSites{$pos}{$replace});
			print STDOUT '"' . "\n";
		}
	}
}

if ($SNP_VarScan_File){
	open (SNPFILE, "<$SNP_VarScan_File") or die ("Unable to open file for reading: $SNP_VarScan_File\n$!\n");
	while ($line= <SNPFILE>){
		chomp ($line);
		if ($line=~ m/^\S+\t(\S+)\t\S+\t\S+\t(\S+)\t.+\S+:\S+:\S+:(\S+):\S+:\S+:\S+:\S+:\S+:\S+:\S+:\S+:(\S+):(\S+)$/){ # Real frequency of variant
#		if ($line=~ m/^\S+\t(\S+)\t\S+\t\S+\t(\S+)\t.+\S+:\S+:\S+:\S+:\S+:\S+:\S+:\S+:\S+:\S+:(\S+):(\S+):(\S+):(\S+)$/){ # Relative freq between two variants
			# Columns GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR
			$pos= $1;
			$replace= lc($2);
			
			## Real frequency of variant
			$depthHighQual= $3;
			$altFwd= $4;
			$altRev= $5;
			$SNPSites{$pos}{$replace}= ($altFwd+$altRev)/$depthHighQual;
			#####
			
			## Relative frequency between two variants
#			$refFwd= $3;
#			$refRev= $4;
#			$altFwd= $5;
#			$altRev= $6;
#			$SNPSites{$pos}{$replace}= ($altFwd+$altRev)/($altFwd+$altRev+$refFwd+$refRev);
			#####
		}
	}
	close (SNPFILE);
	
	foreach $pos (sort {$a <=> $b} keys %SNPSites){
		foreach $replace (sort {lc($a) cmp lc($b)} keys %{$SNPSites{$pos}}){
			print STDOUT '     variation       ' . $pos . '..' . $pos . "\n";		
			print STDOUT '                     /inference="ab initio prediction:VarScan:2.3.9"' . "\n";
			print STDOUT '                     /note="SNP"' . "\n";
			print STDOUT '                     /replace="' . $replace . '"' . "\n";
			print STDOUT '                     /frequency="';
			print STDOUT sprintf("%.4f", $SNPSites{$pos}{$replace});
			print STDOUT '"' . "\n";
		}
	}
}


