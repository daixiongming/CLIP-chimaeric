#!/usr/bin/env perl

use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;

# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Filter alignments"],
	[],
	['ifile=s',
		'Alignment file. Use - for STDIN',
		{required => 1}],
	['min-score=i',
		'Keep reads equal or above this score',
	],
	['max-score=i',
		'Keep reads equal or below this score',
	],
	['min-midpoint=i',
		'Keep reads equal or above this midpoint',
	],
	['max-midpoint=i',
		'Keep reads equal or below this midpoint',
	],
	['verbose|v', 'Print progress'],
	['help|h', 'Print usage and exit',
		{shortcircuit => 1}],
);
print($usage->text), exit if $opt->help; 


warn "opening input file\n" if $opt->verbose;
my $IN = filehandle_for($opt->ifile); 

while (my $line = <$IN>){
	chomp $line;
	my ($name, $seqA, $seqB, $lengthA, $lengthB, $alignment_score, $center_of_binding_loc_on_lgclip, $bindingB, $binding_revcomp_seqA, $revcomp_seqA) = split("\t", $line);
	
	next if (defined $opt->min_score) and ($alignment_score < $opt->min_score);
	next if (defined $opt->max_score) and ($alignment_score > $opt->max_score);
	next if (defined $opt->min_midpoint) and ($center_of_binding_loc_on_lgclip < $opt->min_midpoint);
	next if (defined $opt->max_midpoint) and ($center_of_binding_loc_on_lgclip > $opt->max_midpoint);
	
	print join("\t", ($name, $seqA, $seqB, $lengthA, $lengthB, $alignment_score, $center_of_binding_loc_on_lgclip, $bindingB, $binding_revcomp_seqA, $revcomp_seqA))."\n";
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
