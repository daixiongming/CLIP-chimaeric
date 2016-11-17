#!/usr/bin/env perl

use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;

# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Cut selected columns from table"],
	[],
	['ifile1=s',
		'SAM file of pairs. Use - for STDIN',
		{required => 1}],
	['ifile2=s',
		'SAM file of pairs. Use - for STDIN',
		{required => 1}],
	['ofile1=s',
		'output file for ifile1.',
		{required => 1}],
	['ofile2=s',
		'output file for ifile2.',
		{required => 1}],
	['removen',
		'flag to remove pairs containing N nucleotides',
		{default => 0}],
	['removeintrons',
		'flag to remove pairs containing intron',
		{default => 0}],
	['verbose|v', 'Print progress'],
	['help|h', 'Print usage and exit',
		{shortcircuit => 1}],
);
print($usage->text), exit if $opt->help; 


warn "opening input file\n" if $opt->verbose;
my $IN1 = filehandle_for($opt->ifile1);
my $IN2 = filehandle_for($opt->ifile2);

my %read;
my %best_read_index;

while (my $line = $IN1->getline()){
	chomp $line;
	my $bitflag = (split("\t", $line))[1];
	if ($bitflag !~ /^[0-9]+$/ ){next;} #skip header lines
	if ($bitflag & 256){next;}; #skip secondary alignments
	
	my $cigar = (split("\t", $line))[5];
	if ($opt->removeintrons){
		if ( $cigar =~ /N/){next;}
	}
	
	if ($opt->removen){
		my $seq = (split("\t", $line))[9];
		if ($seq =~ s/N//gi){next;}
	}
	
	my $qname = (split("\t", $line))[0];
	my @namefields = split("-", $qname);
	my $read_index = $namefields[-1];
	my $read_name = $namefields[-2];
	
	$read{$read_name}{$read_index}{'score'} += matchscore($cigar);
}
close $IN1;

while (my $line = $IN2->getline()){
	chomp $line;
	my $bitflag = (split("\t", $line))[1];
	if ($bitflag !~ /^[0-9]+$/ ){next;} #skip header lines
	if ($bitflag & 256){next;}; #skip secondary alignments
	my $cigar = (split("\t", $line))[5];
	if ($opt->removeintrons){
		if ( $cigar =~ /N/){next;}
	}
	
	if ($opt->removen){
		my $seq = (split("\t", $line))[9];
		if ($seq =~ s/N//gi){next;}
	}
	
	my $qname = (split("\t", $line))[0];
	my @namefields = split("-", $qname);
	my $read_index = $namefields[-1];
	my $read_name = $namefields[-2];
	
	if (exists $read{$read_name}{$read_index}){
		$read{$read_name}{$read_index}{'score'} += matchscore($cigar);
		if (exists $best_read_index{$read_name}){
			my $best_score_to_now = $read{$read_name}{$best_read_index{$read_name}}{'score'};
			if ($read{$read_name}{$read_index}{'score'} > $best_score_to_now){
				$best_read_index{$read_name} = $read_index;
			}
		}
		else {
			$best_read_index{$read_name} = $read_index;
		}
	}
}
close $IN2;

### we have picked the best pairs ###

$IN1 = filehandle_for($opt->ifile1);
$IN2 = filehandle_for($opt->ifile2);
open (my $OUT1, ">", $opt->ofile1) or die;
open (my $OUT2, ">", $opt->ofile2) or die;

while (my $line = $IN1->getline()){
	chomp $line;
	my $bitflag = (split("\t", $line))[1];
	if ($bitflag !~ /^[0-9]+$/ ){next;} #skip header lines
	if ($bitflag & 256){next;}; #skip secondary alignments
	my $qname = (split("\t", $line))[0];
	my @namefields = split("-", $qname);
	my $read_index = $namefields[-1];
	my $read_name = $namefields[-2];
	
	if (exists $best_read_index{$read_name}){
		if ($best_read_index{$read_name} eq $read_index){
			print $OUT1 $line."\n";
		}
	}
}
close $IN1;



while (my $line = $IN2->getline()){
	chomp $line;
	my $bitflag = (split("\t", $line))[1];
	if ($bitflag !~ /^[0-9]+$/ ){next;} #skip header lines
	if ($bitflag & 256){next;}; #skip secondary alignments
	my $qname = (split("\t", $line))[0];
	my @namefields = split("-", $qname);
	my $read_index = $namefields[-1];
	my $read_name = $namefields[-2];
	
	if (exists $best_read_index{$read_name}){
		if ($best_read_index{$read_name} eq $read_index){
			print $OUT2 $line."\n";
		}
	}
}
close $IN2;

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

sub matchscore {
	my ($cigar) = @_;
	my $count = 0;
	while ($cigar =~ /(\d+)M/g) {
		$count += $1;
	}
	return $count;
}