#!/usr/bin/env perl

use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;

# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Cut selected columns from table"],
	[],
	['ifile=s',
		'file of pairs. Use - for STDIN',
		{required => 1}],
	['verbose|v', 'Print progress'],
	['test=i', 'only run for few reads'],
	['help|h', 'Print usage and exit',
		{shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

warn "opening input file\n" if $opt->verbose;
my $IN = filehandle_for($opt->ifile);

my $counter = 0;
while (my $line = $IN->getline()){
	chomp $line;
	my ($name, $seqA, $seqB) = split("\t", $line);
	my $revcomp_seqA = revcomp($seqA);
	my ($alignment_score, $bindA, $bindB, $align1, $align2) = align($revcomp_seqA, $seqB);
	my $number_of_matches = add_matches($bindB);
	my $center_of_binding_loc_on_lgclip = center_of_binding_loc_on_lgclip($bindB, $number_of_matches);
	
	if ($opt->verbose){ warn join("\n", $alignment_score, $revcomp_seqA, $seqB, $align1, $align2)."\n";}
	
	print join("\t", $line, length($seqA), length($seqB), $alignment_score, $center_of_binding_loc_on_lgclip, join("", @{$bindB}), join("", @{$bindA}), $revcomp_seqA)."\n";
	
	if (defined $opt->test){
		$counter++;
		if ($counter >= $opt->test){
			last;
		}
	}
}

exit;

### Subroutines

sub revcomp {
	my ($seq) = @_;
	$seq = reverse($seq);
	$seq =~ tr/ACTGactg/TGACtgac/;
	return $seq;
}

sub align {

	my ($seq1, $seq2) = @_;

	# scoring scheme
	my $MATCH    =  1;
	my $MISMATCH = -1;
	my $GAP      = -2;

	# initialization
	my @matrix;
	$matrix[0][0]{score}   = 0;
	$matrix[0][0]{pointer} = "none";
	for(my $j = 1; $j <= length($seq1); $j++) {
		$matrix[0][$j]{score}   = 0;
		$matrix[0][$j]{pointer} = "none";
	}
	for (my $i = 1; $i <= length($seq2); $i++) {
		$matrix[$i][0]{score}   = 0;
		$matrix[$i][0]{pointer} = "none";
	}

	# fill
	my @max_i;
	my @max_j;
	my $max_score = 0;

	for(my $i = 1; $i <= length($seq2); $i++) {
		for(my $j = 1; $j <= length($seq1); $j++) {
			my ($diagonal_score, $left_score, $up_score);
			
			# calculate match score
			my $letter1 = substr($seq1, $j-1, 1);
			my $letter2 = substr($seq2, $i-1, 1);       
			if ($letter1 eq $letter2) {
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
			}
			else {
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
			}
			
			# calculate gap scores
			$up_score   = $matrix[$i-1][$j]{score} + $GAP;
			$left_score = $matrix[$i][$j-1]{score} + $GAP;
			
			if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
				$matrix[$i][$j]{score}   = 0;
				$matrix[$i][$j]{pointer} = "none";
				next; # terminate this iteration of the loop
			}
			
			# choose best score
			if ($diagonal_score >= $up_score) {
				if ($diagonal_score >= $left_score) {
					$matrix[$i][$j]{score}   = $diagonal_score;
					$matrix[$i][$j]{pointer} = "diagonal";
				}
				else {
					$matrix[$i][$j]{score}   = $left_score;
					$matrix[$i][$j]{pointer} = "left";
				}
			} else {
				if ($up_score >= $left_score) {
					$matrix[$i][$j]{score}   = $up_score;
					$matrix[$i][$j]{pointer} = "up";
				}
				else {
					$matrix[$i][$j]{score}   = $left_score;
					$matrix[$i][$j]{pointer} = "left";
				}
			}
			
			# set maximum score
			if ($matrix[$i][$j]{score} > $max_score) {
				@max_i=();
				@max_j=();
			}
			if ($matrix[$i][$j]{score} >= $max_score) {
				$max_score = $matrix[$i][$j]{score};
				push @max_i, $i;
				push @max_j, $j;
			}
		}
	}

	# trace-back

	my $align1 = "";
	my $align2 = "";
	my @match_cnt_j = dots(length($seq1));
	my @match_cnt_i = dots(length($seq2));
	
	my $random_max_index = int(rand @max_i);
	my $ori_j = $max_j[$random_max_index];
	my $ori_i = $max_i[$random_max_index];
		
	### get end i,j
	my $j = $ori_j;
	my $i = $ori_i;
	while ($matrix[$i][$j]{pointer} ne "none"){
		if ($matrix[$i][$j]{pointer} eq "diagonal") {$i--; $j--;}
		elsif ($matrix[$i][$j]{pointer} eq "left") {$j--;}
		elsif ($matrix[$i][$j]{pointer} eq "up") {$i--;}
	}
	
	### add S (soft-clip) to positions upstream
	my $end_j = $j;
	my $end_i = $i;
	while (($i >= 0) and ($j >= 0)){
		if (substr($seq1, $j, 1) eq substr($seq2, $i, 1)){
			$match_cnt_j[$j] = 'm';
			$match_cnt_i[$i] = 'm';
		}
		else {
			$match_cnt_j[$j] = 's';
			$match_cnt_i[$i] = 's';		
		}
		$i--;
		$j--;
	}
	
	### add S (soft-clip) to positions downstream
	$j = $ori_j;
	$i = $ori_i;
	while (($i < length($seq2)) and ($j < length($seq1))){
		if (substr($seq1, $j, 1) eq substr($seq2, $i, 1)){
			$match_cnt_j[$j] = 'm';
			$match_cnt_i[$i] = 'm';
		}
		else {
			$match_cnt_j[$j] = 's';
			$match_cnt_i[$i] = 's';		
		}
		
		$i++;
		$j++;
	}
	
	
	
	warn join("\t", $ori_i, $ori_j, $end_i, $end_j)."\n";
# 	print $i.",".$j."\t".$max_score."\n";
	


	### TRACKBACK
	$i = $ori_i;
	$j = $ori_j;
	while (1) {
		last if $matrix[$i][$j]{pointer} eq "none";
		
		if ($matrix[$i][$j]{pointer} eq "diagonal") {
			$align1 .= substr($seq1, $j-1, 1);
			$align2 .= substr($seq2, $i-1, 1);
			
			if (substr($seq1, $j-1, 1) eq substr($seq2, $i-1, 1)){
				$match_cnt_j[$j-1] = 1;
				$match_cnt_i[$i-1] = 1;
			}
			else {
				$match_cnt_j[$j-1] = '0';
				$match_cnt_i[$i-1] = '0';
			}
			
			$i--; $j--;
		}
		elsif ($matrix[$i][$j]{pointer} eq "left") {
			$align1 .= substr($seq1, $j-1, 1);
			$match_cnt_j[$j-1] = 'g';
			$align2 .= "-";
			$j--;
		}
		elsif ($matrix[$i][$j]{pointer} eq "up") {
			$align1 .= "-";
			$align2 .= substr($seq2, $i-1, 1);
			$match_cnt_i[$i-1] = 'g';
			$i--;
		}   
	}

	$align1 = reverse $align1;
	$align2 = reverse $align2;
	
	warn join("", @match_cnt_i)."\n";
	warn join("", @match_cnt_j)."\n";
	warn $seq1."\n";
	
	
	return ($max_score, \@match_cnt_j, \@match_cnt_i, $align1, $align2);
}

sub center_of_binding_loc_on_lgclip{
	my ($binding, $totalmatches) = @_;
	my $goal = int($totalmatches/2);
	my @binding = @{$binding};
	my $count=0;
	for (my $i = 0; $i < @binding; $i++){
		if ((defined $binding[$i]) and ($binding[$i] eq '1')){
			$count++;
			if ($count == $goal){
				return $i;
			}
		}
	}
}

sub add_matches{
	my ($inref) = @_;
	my $total = 0;
	foreach my $val (@$inref){
		if ((defined $val) and ($val eq '1')){
			$total++;
		}
	}
	return $total;
}

sub dots {
	my ($N) = @_;
	my @array;
	for (my $i = 0; $i < $N; $i++){
		push @array, '.';
	}
	return @array;
}


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