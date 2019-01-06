#!/usr/bin/perl

use Getopt::Long;
use Data::Dumper;
use Cwd 'abs_path';


GetOptions(\%opts, "shEnt:s", "SNPsam:s", "SNPvarscan:s", "outformat:s", "help|h!");
        sub usage(){
                die "USAGE :: perl cleanSamBWT.pl -i inFile -b backboneFile.fasta -minOvlp percentage100Scale [-help|-h]\n\n
	-shEnt\n\tHighly entropic sites tabular file [String]\n\n
	-SNPsam\n\tSNP vcf file from SAMtools [String]\n\n
	-SNPvarscan\n\tSNP vcf file from VarScan [String]\n\n
	-help or -help\n\tGet help for using this script\n\n";
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
my $SNPFileSam= $opts{'SNPsam'};
my $SNPFileVarscan= $opts{'SNPvarscan'};
my $outFormat= $opts{'outformat'};


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
my %varHash= ();
my @varArray= ('DP', 'VDB', 'RPB', 'AF1', 'AC1', 'DP4', 'MQ', 'FQ', 'PV4');


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

if ($SNPFileSam){
	if ($opts{'outformat'} eq 'tab'){
		print "CHROM\tPOS\tREF\tALT\tQUAL\tDP\tVDB\tRPB\tAF1\tAC1\tDP4_1\tDP4_2\tDP4_3\tDP4_4\tMQ\tFQ\tPV4_1\tPV4_2\tPV4_3\tPV4_4\tFORMAT\n";
	}
	open (SNPFILE, "<$SNPFileSam") or die ("Unable to open file for reading: $SNPFileSam\n$!\n");
	while ($line= <SNPFILE>){
		chomp ($line);
		# CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tbam
		#DP=18904;VDB=4.757385e-01;RPB=-3.395079e+00;AF1=0.5;AC1=1;DP4=2239,8592,1896,6048;MQ=44;FQ=225;PV4=1.9e-07,4.7e-09,0.29,1
		if ($line=~ m/^([^\s#]+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/){
			%varHash= ('DP' => 'NA', 'VDB' => 'NA', 'RPB' => 'NA', 'AF1' => 'NA', 'AC1' => 'NA', 'DP4' => ['NA','NA','NA','NA'], 'MQ' => 'NA', 'FQ' => 'NA', 'PV4' => ['NA','NA','NA','NA']);
			$chr= $1;
			$pos= $2;
			$ID= $3;
			$refChar= $4;
			$altChar= $5;
			$snpQual= $6;
			$filter= $7;
			$snpInfo= $8;
			$format= $9;
			if ($opts{'outformat'} eq "tab"){
				print STDOUT "$chr\t$pos\t$refChar\t$altChar\t$snpQual\t";
				while ($snpInfo=~ /(\w+)=([^;]+)/g){
					$varName= $1;
					$varString= $2;
					if ($varName eq 'DP4' || $varName eq 'PV4'){
						$c= 0;
						while ($varString=~ m/([^,]+)/g){
							$varHash{$varName}[$c]= $1;
							$c++;
						}
					}
					else {
						$varHash{$varName}= $varString;
					}
				}
				print STDOUT $varHash{'DP'} . "\t" . $varHash{'VDB'} . "\t" . $varHash{'RPB'} . "\t" . $varHash{'AF1'} . "\t" . $varHash{'AC1'} . "\t" . join("\t", @{$varHash{'DP4'}}) . "\t" . $varHash{'MQ'} . "\t" . $varHash{'FQ'} . "\t" . join("\t", @{$varHash{'PV4'}});
			
				print STDOUT "\t$format\n";
			}
			elsif ($opts{'outformat'} eq 'gbk' && $snpInfo=~ m/DP4=(\d+),(\d+),(\d+),(\d+)/){
				$refFwdFreq= $1;
				$refRevFreq= $2;
				$altFwdFreq= $3;
				$altRevFreq= $4;
				$SNPSites{$pos}{$altChar}= ($altFwdFreq+$altRevFreq)/($altFwdFreq+$altRevFreq+$refFwdFreq+$refRevFreq);
			}
		}
	}
	close (SNPFILE);
	if ($opts{'outformat'} eq 'gbk'){
		foreach $pos (sort {$a <=> $b} keys %SNPSites){
			foreach $replace (sort {lc($a) cmp lc($b)} keys %{$SNPSites{$pos}}){
				print STDOUT '     variation       ' . $pos . '..' . $pos . "\n";		
				print STDOUT '                     /frequency="';
				print STDOUT sprintf("%.4f", $SNPSites{$pos}{$replace});
				print STDOUT '"' . "\n";
				print STDOUT '                     /inference="ab initio prediction: SAMtools mpileup"' . "\n";
				print STDOUT '                     /note="SNP"' . "\n";
				print STDOUT '                     /replace="' . $replace . '"' . "\n";
			}
		}
	}
}

if ($SNPFileVarscan){
	if ($opts{'outformat'} eq 'tab'){
		print "CHROM\tPOS\tREF\tALT\tQUAL\tSDP\tFREQ\tRDF\tRDR\tADF\tADR\n";
	}
	open (SNPFILE, "<$SNPFileVarscan") or die ("Unable to open file for reading: $SNPFileVarscan\n$!\n");
	while ($line= <SNPFILE>){
		chomp ($line);
		# CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tbam
		#DP=18904;VDB=4.757385e-01;RPB=-3.395079e+00;AF1=0.5;AC1=1;DP4=2239,8592,1896,6048;MQ=44;FQ=225;PV4=1.9e-07,4.7e-09,0.29,1
		if ($line=~ m/^([^\s#]+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/){
			%varHash= ('SDP' => 'NA', 'FREQ' => 'NA', 'RDF' => 'NA', 'RDR' => 'NA', 'ADF' => 'NA', 'ADR' => 'NA');
			$chr= $1;
			$pos= $2;
			$ID= $3;
			$refChar= $4;
			$altChar= $5;
			$snpQual= $6;
			$filter= $7;
			$snpInfo= $8;
			$format= $9;
			$format_values= $10;
			if ($opts{'outformat'} eq "tab"){
				print STDOUT "$chr\t$pos\t$refChar\t$altChar\t$snpQual\t";
				@temp_val= split(/:/, $format_values);
				$varHash{'SDP'}= $temp_val[2];
				$varHash{'FREQ'}= $temp_val[6];
				$varHash{'RDF'}= $temp_val[10];
				$varHash{'RDR'}= $temp_val[11];
				$varHash{'ADF'}= $temp_val[12];
				$varHash{'ADR'}= $temp_val[13];
				print STDOUT $varHash{'SDP'} . "\t" . $varHash{'FREQ'} . "\t" . $varHash{'RDF'} . "\t" . $varHash{'RDR'} . "\t" . $varHash{'ADF'} . "\t" .$varHash{'ADR'} . "\n";
			}
			elsif ($opts{'outformat'} eq 'gbk' && $snpInfo=~ m/DP4=(\d+),(\d+),(\d+),(\d+)/){
				$refFwdFreq= $1;
				$refRevFreq= $2;
				$altFwdFreq= $3;
				$altRevFreq= $4;
				$SNPSites{$pos}{$altChar}= ($altFwdFreq+$altRevFreq)/($altFwdFreq+$altRevFreq+$refFwdFreq+$refRevFreq);
			}
		}
	}
	close (SNPFILE);
	if ($opts{'outformat'} eq 'gbk'){
		foreach $pos (sort {$a <=> $b} keys %SNPSites){
			foreach $replace (sort {lc($a) cmp lc($b)} keys %{$SNPSites{$pos}}){
				print STDOUT '     variation       ' . $pos . '..' . $pos . "\n";		
				print STDOUT '                     /frequency="';
				print STDOUT sprintf("%.4f", $SNPSites{$pos}{$replace});
				print STDOUT '"' . "\n";
				print STDOUT '                     /inference="ab initio prediction: SAMtools mpileup"' . "\n";
				print STDOUT '                     /note="SNP"' . "\n";
				print STDOUT '                     /replace="' . $replace . '"' . "\n";
			}
		}
	}
}


