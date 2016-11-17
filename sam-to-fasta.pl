#!/usr/bin/env perl

# Load modules
use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;

# Load GenOO library
use GenOO::Data::File::FASTA;
use GenOO::RegionCollection::Factory;

# Define and read command line options
my ($opt, $usage) = describe_options(
	'%c %o',
	["Print in FASTA format the genomic sequence of each alignment in the input SAM file."],
	[],
	['sam=s', 'Input SAM file', { required => 1}],
	['chr_dir=s', 'Directory with chromosome fasta files', { required => 1}],
	['out-length=i', 'If set, output sequences will be extended on both sides up to INT'],
	['max-length=i', 'If set, reads longer than INT are discarded'],
	['help|h', 'Print usage message and exit'],
);
print($usage->text), exit if $opt->help;

warn "Reading regions\n";
my $feats_col = GenOO::RegionCollection::Factory->create('SAM', {
	file => $opt->sam
})->read_collection;
my @feats = sort {$a->rname cmp $b->rname} $feats_col->all_records;

warn "Extracting sequences\n";
my ($open_rname_file, $rname_seq) = ('', '');
foreach my $r (@feats) {
	my $rname_file = $opt->chr_dir . $r->rname . '.fa';
	if ($rname_file ne $open_rname_file) {
		warn "Opening file $rname_file\n";
		if (!-e $rname_file and !-e $rname_file.'.gz') {
			warn "skipping $rname_file\n" if ! -e $rname_file;
			next;
		}
		my $fp;
		if (-e $rname_file.'.gz') {
			$fp = GenOO::Data::File::FASTA->new(file => $rname_file.'.gz');
		} else {
			$fp = GenOO::Data::File::FASTA->new(file => $rname_file);
		}
		my $rec = $fp->next_record;
		if ($rec->header eq $r->rname) {
			$rname_seq = $rec->sequence;
			$open_rname_file = $rname_file;
		}
		else {
			die;
		}
	}

	my $flank = 0;
	
	if (defined $opt->out_length){
		my $original_read_length = $r->stop - $r->start + 1;
		
		my $total_seq_length_toadd = $opt->out_length - $original_read_length;
		if ($total_seq_length_toadd < 0){warn "read longer than out length - trimming read\n";}
		$flank = int($total_seq_length_toadd/2)+1; #adding one because rounding down makes smaller sequences - will be trimmed to size at the end
	}
	
	
	my $r_seq = region_sequence_from_seq(
		\$rname_seq, $r->strand, $r->rname, $r->start, $r->stop, $flank) ||	die;
	
	if (defined $opt->out_length){
		$r_seq = trim_to_size($r_seq, $opt->out_length);
		if (!defined $r_seq){
			warn "seq: $r_seq length: ".length($r_seq)." shorter than needed - skipping\n";
			next;
		}
	}
	
	if ((defined $opt->max_length) and (length($r_seq) > $opt->max_length)){
		next;
	}
	
	say '>'.$r->qname."\n".uc($r_seq);
}

##############################
##############################
sub region_sequence_from_seq {
	my ($seq_ref, $strand, $rname, $start, $stop, $flank) = @_;

	#out of bounds
	return if ($start - $flank < 0);
	return if ($stop + $flank > length($$seq_ref) - 1);

	if ($strand == 1) {
		return substr($$seq_ref, $start-$flank, 2*$flank+$stop-$start+1);
	}
	else {
		my $seq = reverse(
			substr($$seq_ref, $start-$flank, 2*$flank+$stop-$start+1));
		if ($seq =~ /U/i) {
			$seq =~ tr/ATGCUatgcu/UACGAuacga/;
		}
		else {
			$seq =~ tr/ATGCUatgcu/TACGAtacga/;
		}
		return $seq;
	}
}

sub trim_to_size {
	my ($seq, $length) = @_;
	
	return undef if (length($seq) < $length);
	
	while (length($seq) >= $length+2){
		#trim on both sides
		my @seq = split('', $seq);
		pop @seq;
		shift @seq;
		$seq = join('', @seq);		
	}
	
	if (length($seq) > $length){
		# if there is 1 nt left to trim we will do it randomly either from the beginning or end
		my @seq = split('', $seq);
		if (rand() > 0.5){pop @seq;}
		else {shift @seq;}
		$seq = join('', @seq);	
	}
	
	return $seq;
}
