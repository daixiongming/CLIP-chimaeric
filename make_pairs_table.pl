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
		'fasta file of sequences. Use - for STDIN',
		{required => 1}],
	['ifile2=s',
		'fasta file of sequences. Use - for STDIN',
		{required => 1}],
	['shuffle', 'shuffle reads keeping size'],
	['verbose|v', 'Print progress'],
	['help|h', 'Print usage and exit',
		{shortcircuit => 1}],
);
print($usage->text), exit if $opt->help; 

warn "opening first file\n" if $opt->verbose;
my $FN = filehandle_for($opt->ifile1);

my %sequence1;
my %random_seq_by_size;

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

warn "opening second file\n" if $opt->verbose;
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
