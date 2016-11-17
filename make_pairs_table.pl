#!/usr/bin/env perl

use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;

# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Join entries from two FASTQ files by name. Output a table file where the first column corresponds to the entry name and the next two columns are the sequences from the input FASTQ files corresponding to that name."],
	[],
	['ifile1=s',
		'Input FASTA file 1. Use - for STDIN',
		{required => 1}],
	['ifile2=s',
		'Input FASTA file 2. Use - for STDIN',
		{required => 1}],
	['shuffle', 'If set, a random join that preserves sequence sizes is performed'],
	['verbose|v', 'Print progress'],
	['help|h', 'Print usage and exit',
		{shortcircuit => 1}],
);
print($usage->text), exit if $opt->help; 

if ($opt->ifile1 eq '-' and $opt->ifile2 eq '-') {
	die "cannot use STDIN for both input files\n";
}

my %sequence1;
my %random_seq_by_size;
warn "reading file 1\n" if $opt->verbose;
my $FN = filehandle_for($opt->ifile1);
while (my $line = <$FN>){
	if ($line =~ /^>/){
		my $header = $line;
		chomp $header;
		my $seqline = $FN->getline();
		chomp $seqline;
		$sequence1{$header} = $seqline;
		
		if ($opt->shuffle){
			push @{$random_seq_by_size{length($seqline)}}, $seqline;
		}
	}
}
close $FN;

warn "reading file 2 and joining\n" if $opt->verbose;
my $IN = filehandle_for($opt->ifile2);
while (my $line = $IN->getline()){
	if ($line =~ /^>/){
		my $header = $line;
		chomp $header;
		my $seqline = $IN->getline();
		chomp $seqline;
		if (exists $sequence1{$header}){
			my $seq1 = $sequence1{$header};
			if ($opt->shuffle){
				my $length1 = length($seq1);
				my $number_of_random_reads = @{$random_seq_by_size{$length1}};
				my $randomseq1 = $random_seq_by_size{$length1}->[int(rand($number_of_random_reads))];
				$seq1 = $randomseq1;
			}
			$header =~ s/^>//g;
			print $header."\t".$seq1."\t".$seqline."\n";
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
