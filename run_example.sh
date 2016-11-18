species='SPECIES_ASSEMBLY'

array=( \
	'samplename1' \
	'samplename2' \
	'samplename3' \
)
for ((i=0;i<${#array[@]};++i));
do
	samplename="${array[i]}"
	mkdir results/$samplename/

	# SPLIT UNMAPPED LONG READS

	perl dev/split_unmapped_reads.pl \
		--ifile data/$samplename/aligned/reads.Unmapped.out.mate1 \
		--min-small 23 \
		--max-small 29 \
		--min-large 30 \
		--out-small results/$samplename/fragments.small.fa \
		--out-large results/$samplename/fragments.large.fa

	## FILTER SMALL FRAGMENTS WITH LIBRARY OF KNOWN SMALL RNAs

	cat data/$samplename/aligned/$samplename-smallreads.Aligned.sortedByCoord.out.bam \
		| samtools-1.2 view \
			-h \
			- \
		| perl dev/filter_fa_by_sam_sequence.pl \
			--ifile results/$samplename/fragments.small.fa \
			--ffile - \
		> results/$samplename/fragments.small.filtered.fa

	## FILTER LARGE FRAGMENTS WITH SMALL FRAGMENTS

	perl dev/filter_by_name.pl \
		--ifile results/$samplename/fragments.large.fa \
		--ref-file results/$samplename/fragments.small.filtered.fa \
		> results/$samplename/fragments.large.filtered.fa

	### ALIGN LARGE FRAGMENTS

	STAR-2.5.2b \
		--genomeDir data/genome/$species/ \
		--readFilesIn results/$samplename/fragments.large.filtered.fa \
		--readFilesCommand cat \
		--runThreadN 12 \
		--outSAMtype SAM \
		--outSAMattributes All \
		--outFilterMultimapScoreRange 0 \
		--alignIntronMax 50000 \
		--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
		--outFilterMatchNmin 8 \
		--outFilterMatchNminOverLread 0.9 \
		--outFileNamePrefix results/$samplename/fragments.large.filtered. \
		--sjdbOverhang 100 \
		--sjdbGTFfile /store/data/species/$species/annotation/UCSC_gene_parts_genename.gtf

	### ALIGN SMALL FRAGMENTS

	STAR-2.5.2b \
		--genomeDir data/genome/$species/ \
		--readFilesIn results/$samplename/fragments.small.filtered.fa \
		--readFilesCommand cat \
		--runThreadN 12 \
		--outSAMtype SAM \
		--outSAMattributes All \
		--outFilterMultimapScoreRange 0 \
		--alignIntronMax 50000 \
		--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
		--outFilterMatchNmin 8 \
		--outFilterMatchNminOverLread 0.9 \
		--outFileNamePrefix results/$samplename/fragments.small.filtered. \
		--sjdbOverhang 100 \
		--sjdbGTFfile /store/data/species/$species/annotation/UCSC_gene_parts_genename.gtf
		
	### 
	### FILTERS 
	### 

	### FILTER LARGE FRAGMENTS ON GENES

	cat results/$samplename/fragments.large.filtered.Aligned.out.sam \
		| samtools-1.2 view \
			-b \
			-T data/genome/$species/genome.fa \
			- \
		> results/$samplename/fragments.large.filtered.Aligned.bam

	bedtools intersect \
		-a results/$samplename/fragments.large.filtered.Aligned.bam \
		-b /store/data/species/$species/annotation/UCSC_gene_parts_genename.gtf \
		-split \
		-u \
	| samtools-1.2 view \
		- \
		> results/$samplename/fragments.large.filtered.Aligned.OnGenes.sam

	### KEEP UNIQUE PAIR

	perl dev/filter_best_pairs.pl \
		--ifile1 results/$samplename/fragments.large.filtered.Aligned.OnGenes.sam \
		--ifile2 results/$samplename/fragments.small.filtered.Aligned.out.sam \
		--ofile1 results/$samplename/fragments.large.filtered.Aligned.OnGenes.paired.sam \
		--ofile2 results/$samplename/fragments.small.filtered.Aligned.paired.sam

	### EXTEND LARGE FRAGMENT TO 200 nt > FASTA

	perl dev/sam-to-fasta.pl \
		--sam results/$samplename/fragments.large.filtered.Aligned.OnGenes.paired.sam \
		--chr_dir /store/data/UCSC/$species/chromosomes/ \
		--out-length 200 \
		> results/$samplename/fragments.large.filtered.Aligned.OnGenes.paired.200nt.fa
		
	### SMALL FRAGMENT TO FASTA NO EXTENSION

	perl dev/sam-to-fasta.pl \
		--sam results/$samplename/fragments.small.filtered.Aligned.paired.sam \
		--chr_dir /store/data/UCSC/$species/chromosomes/ \
		--max-length 29 \
		> results/$samplename/fragments.small.filtered.Aligned.paired.fa

	### MAKE ALIGNMENT FILE -> SMALL AND LARGE FRAGMENTS TOGETHER

	perl dev/make_pairs_table.pl \
		--ifile1 results/$samplename/fragments.small.filtered.Aligned.paired.fa \
		--ifile2 results/$samplename/fragments.large.filtered.Aligned.OnGenes.paired.200nt.fa \
		> results/$samplename/pairs.tab

	perl dev/make_pairs_table.pl \
		--ifile1 results/$samplename/fragments.small.filtered.Aligned.paired.fa \
		--ifile2 results/$samplename/fragments.large.filtered.Aligned.OnGenes.paired.200nt.fa \
		--shuffle \
		> results/$samplename/pairs.shuffled.tab

	### ALIGN
		
	perl dev/align.pl \
		--ifile results/$samplename/pairs.tab \
		> results/$samplename/pairs.aligned.tab
		
	perl dev/align.pl \
		--ifile results/$samplename/pairs.shuffled.tab \
		> results/$samplename/pairs.shuffled.aligned.tab
		
	### MAKE PLOTS

	Rscript dev/plots.R \
		--ifile results/$samplename/pairs.aligned.tab \
		--sfile results/$samplename/pairs.shuffled.aligned.tab \
		--ofile results/$samplename/pairs.aligned.pdf

done
