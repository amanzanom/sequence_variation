#!/usr/bin/perl

#===============================================================================
#   Author: Alejandro Manzano Marin
#
#   File: maskFasta.pl
#   Date: 01-04-2013
#   Version: 1.0
#
#   Usage:
#      perl maskFasta.pl -infile|i inFile.fasta -maskRange 'int-int' [-maskChar 'X'] [options]
#
#      Check out 'perl maskFasta.pl -h' for short usage manual and info on the software.
#
#    Description: This program is intended as a masking tool for fast files. It will take a FASTA formatted
#                 file and mask the desired region with a desired masking charater (default 'X').
#                 
#
#    Contact: Contact the author at alejandro.manzano@uv.es using 'maskFasta: ' as
#             as begining for subject for bug reporting, feature request or whatever
#             (of course realted to the software).
#
#    COPYRIGHT: Copyright (C) 2013  Alejandro Manzano-Marin.
#
#    LICENCE: This program is free software: you can redistribute it and/or modify it under the terms
#             of the GNU General Public License as published by the Free Software Foundation, either
#             version 3 of the License, or (at your option) any later version.
#             This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#             without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#             See the GNU General Public License for more details.
#             You should have received a copy of the GNU General Public License along with this program.
#             If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


# Load modules
use Getopt::Long;
use Pod::Usage;


# Define subroutines
sub printVersion {
	my $software= $_[0];
	my $version= $_[1];
	my $fileHandle= $_[2];
	print $fileHandle "$software v$version\n";
	exit (0);
}

sub printFasta {
	my $header= $_[0];
	my $sequence= $_[1];
	my $charPerLine= $_[2];
	my $fileHandle= $_[3];
	
	my $i= 0;
	
	print $fileHandle ">" . $header . "\n";
	if ($charPerLine < 0){
		print STDERR "printFasta ERROR: Illegal number of characters per line, number must be >=0\n";
		return (1);
	}
	if ($charPerLine == 0){
		print $fileHandle $sequence . "\n";
	}
	else {
		for ($i=0; $i<length($sequence); $i+=70){
			if (length(substr($sequence, $i, length($sequence)-$i)) < $charPerLine){
				print $fileHandle substr($sequence, $i, length($sequence)-$i)."\n";
			}
			else {
				print $fileHandle substr($sequence, $i, $charPerLine)."\n";
			}
		}
	}
	return (0);
}


# Variable definition

## Define other variables
my $i= 0;
my $line= '';
my $id= '';
my %sequences= ();
my $range= '';
my $start= 0;
my $end= 0;
my $tempString= '';

## General variables
my $PROGRAMNAME= 'maskFasta';
my $VERSION= '1.0';

## Define options default values
my @opt_inFile= '';
my $opt_maskRange= '';

my $opt_maskChar= 'X';

my $opt_verbose= 0;
my $opt_man= 0;
my $opt_help= 0;
my $opt_printVersion= 0;

## Define options hash
GetOptions(\%opts, 
	'infile|i=s' => \$opt_inFile, 
	'maskRange=s' => \$opt_maskRange, 
	'maskChar:s'=> \$opt_maskChar, 
	'verbose|v!' => \$opt_verbose, 
	'help|h!' => \$opt_help, 
	'man!'  => \$opt_man, 
	'version!' => \$opt_printVersion) || pod2usage(-exitval => 1,  -verbose => 2);

if ($opt_help){
	pod2usage(-exitval => 1,  -verbose => 1);
}
if ($opt_man){
	pod2usage(-exitval => 0, -verbose => 2);
}
if ($opt_printVersion){
	&printVersion($PROGRAMNAME, $VERSION, \*STDERR);
}

# Script documetation

=pod

=head1 NAME

maskFasta

=head1 VERSION

maskFasta v1.0

=head1 SYNOPSIS

perl maskFasta.pl -infile|i inFile.fasta -maskRange 'int-int' [-maskChar 'X'] [-verbose|v] [-help|h] [-man] [-version]

=head1 DESCRIPTION

This program is intended as a masking tool for fast files. It will take a FASTA formatted file and mask the desired region with a desired masking charater (default 'X').

=head1 OPTIONS

=head2 INPUT

=over 8

=item B<-i> | B<-infile> <string> (mandatory)

File containing sequences to be filtered in FAST(A/Q) format.

=item B<-maskRange> <string> (mandatory)

Desired renge(s) to mask from sequence in 1-based index (format 'start-end[,start-end,..]').

=back

=head2 FORMAT

=over 8

=item B<-maskChar> <string> (default: 'X')

Desired character to use for masking.

=back

=head2 INFO AND HELP

=over 8

=item B<-v> | B<-verbose> <boolean> (default: 0)

Prints status and info messages while processing.

=item B<-h> | B<-help> <boolean>

Print useful help on using this script.

=item B<-man> <boolean>

Print the full documentation.

=item B<-version> <boolean>

Print program version.

=back

=head1 AUTHOR

Alejandro Manzano-Marin, C<< <alejandro_dot_manzano_at_uv_dot_es> >>

=head1 BUGS

If you find a bug please email me at C<< <alejandro_dot_manzano_at_uv_dot_es> >> so that I can keep making maskFasta better.

=head1 COPYRIGHT

Copyright (C) 2013  Alejandro Manzano-Marin.

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut


# Check options
if (!$opt_inFile || ($opt_maskRange && $opt_maskRange!~ m/^\d+\-\d+(,\d+\-\d+)*$/) || !$opt_maskChar){
	print STDERR "ERROR:\n";
	if (!$opt_inFile){
		print STDERR "FASTA infile missing\n";
	}
	if ($opt_maskRange && $opt_maskRange!~ m/^\d+\-\d+(,\d+\-\d+)*$/){
		print STDERR "Masking range(s) not in correct format\n";
	}
	if (!$opt_maskChar){
		print STDERR "Mask character not specified\n";
	}
	print STDERR "Please check usage manual\n";
	pod2usage(-exitval => 1,  -verbose => 1);
}


# Assign values to variables dependant on options


# Main program

## Read list of accepted reads


if ($opt_verbose){
	print STDERR "Reading FASTA file and masking sequences " . $opt_inFile . "..";
}

open (FASTAIN, "<$opt_inFile") || die "Unable to open file $opt_inFile for reading\n$!\n";;
while ($line=<FASTAIN>){
	if ($line=~ m/^>(\S+)/ || eof){
		if (eof){
			$sequences{$id}.= $line;
		}
		if (length($sequences{$id})>0){
			foreach $range (split(/,/, $opt_maskRange)){
				$range=~ m/(\d+)-(\d+)/;
				$start= $1;
				$end= $2;
				$tempString= '';
				for ($i=0; $i<($end-$start+1); $i++){
					$tempString.= $opt_maskChar;
				}
				substr($sequences{$id}, $start-1, $end-$start+1, $tempString);
			}
			print STDOUT '>' . $id . "\n" . $sequences{$id} . "\n";
		}
		$id= $1;
		$sequences{$id}= '';
	}
	else {
		$sequences{$id}.= $line;
	}
}
close (FASTAIN);

if ($opt_verbose){
	print STDERR "DONE\n";
}


exit (0);
