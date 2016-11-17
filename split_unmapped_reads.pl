#!/usr/bin/env perl

use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;

# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Takes as input a FASTQ file and splits the sequence of each entry based on specified sizes. For each new entry an integer corresponding to the break point is appended on the header. For each input entry two new entries with the same header are created (a small and a large) and are output in different files."],
	[],
	['ifile=s',
		'Input FASTQ file. Use - for STDIN',
		{required => 1}],
	['out-small=s',
		'Output FASTA file for the smaller split reads',
		{required => 1}],
	['out-large=s',
		'Output FASTA file for the larger split reads',
		{required => 1}],
	['min-small=i',
		'Minimum length for small read',
		{default => 0}],
	['min-large=i',
		'Minimum length for large read',
		{default => 0}],
	['max-small=i',
		'Maximum length for small read',
		{default => 100}],
	['max-large=i',
		'Maximum length for large read',
		{default => 100}],
	['verbose|v', 'Print progress'],
	['help|h', 'Print usage and exit',
		{shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

my $min_small = $opt->min_small;
my $min_large = $opt->min_large;
my $max_small = $opt->max_small;
my $max_large = $opt->max_large;

warn "opening input file\n" if $opt->verbose;
my $IN = filehandle_for($opt->ifile);

my ($SOUT, $LOUT);
open ($SOUT, ">", $opt->out_small);
open ($LOUT, ">", $opt->out_large);

warn "splitting reads\n" if $opt->verbose;
while (my $line = $IN->getline()){
	if ($line =~ /^@/){ #first line of 4 line block in fastq
		chomp $line;
		my $name = $line;
		$name =~ s/^@//g; #removes first '@'
		my $seq = $IN->getline();
		chomp $seq;
		$IN->getline(); #skip empty line
		$IN->getline(); #skip quality line
		
		my $read_length = length($seq);
		for (my $i = 1; $i < $read_length; $i++){
			my $length1 = $i;
			my $length2 = $read_length-$i;
			my $read1 = substr($seq, 0, $length1); #from zero to length i
			my $read2 = substr($seq, $i, $length2); #from i to the remainder
			
			if (
				(length($read1) >= $min_small) and
				(length($read1) <= $max_small) and
				(length($read2) >= $min_large) and
				(length($read2) <= $max_large) and
				(length($read1) <= length($read2))
			){
				print $SOUT ">".$name."-".$i."\n".$read1."\n";
				print $LOUT ">".$name."-".$i."\n".$read2."\n";
			}
			elsif (
				(length($read2) >= $min_small) and
				(length($read2) <= $max_small) and
				(length($read1) >= $min_large) and
				(length($read1) <= $max_large) and
				(length($read2) <= length($read1))
			){
				print $SOUT ">".$name."-".$i."\n".$read2."\n";
				print $LOUT ">".$name."-".$i."\n".$read1."\n";
			}
		}
	}
}

exit;

sub filehandle_for {
	my ($file) = @_;

	if ($file eq '-'){
		open(my $IN, "<-");
		return $IN;
	}
	else {
		open(my $IN, "<", $file);
		return $IN
	}
}
