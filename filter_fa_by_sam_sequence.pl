#!/usr/bin/env perl

use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;

# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Print FASTA entries whose sequence match any of the query sequences defined in a SAM file."],
	[],
	['ifile=s',
		'Input FASTA file with sequences to be filtered. Use - for STDIN',
		{required => 1}],
	['ffile=s',
		'Input SAM file to filter with. Use - for STDIN',
		{required => 1}],
	['min=i',
		'Minimum length of sequences to be considered',
	],
	['max=i',
		'Maximum length of sequences to be considered',
	],
	['verbose|v', 'Print progress'],
	['help|h', 'Print usage and exit',
		{shortcircuit => 1}],
);
print($usage->text), exit if $opt->help; 

if ($opt->ifile eq '-' and $opt->ffile eq '-') {
	die "cannot use STDIN for both input files\n";
}

warn "opening filter file\n" if $opt->verbose;
my $FN = filehandle_for($opt->ffile);

my %accepted_seq;

while (my $line = <$FN>){
	chomp $line;
	my $bitflag = (split("\t", $line))[1];
	if ($bitflag !~ /^[0-9]+$/ ){next;} #skip header lines
	if ($bitflag & 256){next;} #skip secondary alignments
	
	my $seq = (split("\t", $line))[9];
	if ($bitflag & 16){$seq = revcomp($seq);} #reads aligned on the minus strand
	
	if (defined $opt->min and length($seq) < $opt->min){next;}
	if (defined $opt->max and length($seq) > $opt->max){next;}
	$accepted_seq{$seq} = 1;
}
close $FN;

warn "opening input file\n" if $opt->verbose;
my $IN = filehandle_for($opt->ifile);

my %names_to_pass;
while (my $line = $IN->getline()){
	
	if ($line =~ /^>/){
		my $header = $line;		
		chomp $header;
		my $seqline = $IN->getline();
		chomp $seqline;
		if ((exists $accepted_seq{$seqline})){
			print $header."\n".$seqline."\n";
		}
	}
}
close $IN;

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

sub revcomp {
	my ($seq) = @_;
	$seq = reverse($seq);
	$seq =~ tr/ACTGactg/TGACtgac/;
	return $seq;
}
