# nano_tools

nano_tools is a collection of several python scripts I have generated that I find useful for analysing nanopore data.

# demultiplex_nanopore_barcodes.py

This tool demultiplexes nanopore libraries based on barcodes in the 5' or 3' ends. In particular, it was designed to work with libraries generated via PCR, in which the PCR primers contain barcodes. The barcodes can be in the forward or reverse primer, or can be in both primers (in which case combinatorial demultiplexing is performed).

This enables one to, for example, generate 96 uniquely barcoded PCR products using 8 forward and 12 reverse primers.

The algorithm accounts for Nanopore's relatively high error rate and lack of positional precision. For each read it:
1. Takes the first and last N bases (default 100). The first part is used to search for forward primers/barcodes and the last part is used for reverse primers/barcodes
2. Searches for (partial) matches to each barcode using the rapidfuzz library - this allows for insertions/deletions and mismatches. 
3. Repeats steps 1 and 2 for the reverse complement of the read (this can be disabled)
4. Filters for matches that confidently match only a single unique barcode or barcode pair
5. Writes each out to a separate fastq.gz file

It performs a quick check using randomly generated sequences to ensure the parameters used do not lead to significant false discovery rates.

Note that no trimming of the barcodes is performed.

Performance is not especially rapid: roughly 500-1,000 reads are demultiplexed per second using standard parameters (around 3,000,000 per hour)

## Usage

```
--fastq/-f: This is the fastq file that you wish to demultiplex
--primers/-p: This is the csv of primers/barcodes used to generate the library (see below for formatting details).
--output/-o: This is the prefix (including directory) that files are written to

--min_score/-s (default 90): This is the minimum partial ratio score from rapidfuzz for a match to be detected
--max_ambiguity (default 82): If a second match is gains this score or more, the match is rejected because it is considered ambiguous
--length/-l (default 100): This is the length at the start and end of each read that is searched for barcodes
--ignore_rc (default False): This option ignores the reverse complement of the read (i.e. only looks at the forward read)
```

## Formatting the primers/barcodes csv

The csv must have 3 columns:
1. The first column is the primer name (eg "forward_1")
2. The second column is the barcode sequence from the primer (this should be the 5'-3' sequence in the primer itself)
3. The third column must either be "F" or "R" depending on whether it was a forward or reverse barcode

Here is an example:
```
R1,GACATCAATTCGAACAATCC,R
R2,TAACTTACCTGTATAGCACA,R
R3,TTATATACATCTTCCTGGCT,R
F1,AGAATATCTTAGACACTTGC,F
F2,TGTTCCGTATTGCTTAACAA,F
F3,ACTAATCCAGAGTTCTCAAG,F
```

If you only have barcodes on the forward primers, or only on the reverse primers, you still need to put an F or R (respectively) in the third column of each entry.

# extract_splice_junctions_from_bam.py

This tool outputs a csv table with the number of reads that have a given:
1. first and last alignment position
2. set of splice junctions
3. mapping quality
4. sam flag
5. strand

The purpose of this function is to summarise data from a (targeted) RNA sequencing experiment, for easier downstream analysis of splicing. 

It will produce reasonably small files for targeted RNA seq, but files will likely be huge for general RNA seq.

Note that in its current form, it does not take into account read pairings. As such, it is suitable for Nanopore or single-end Illumina, but not for paired-end Illumina sequencing.

## Usage
```
--bam/-b: This is the bam file from which splice junctions etc will be extracted
--output/-o: This is the prefix (including directory) that files are written to

--min_intron_length (default 50): This is the minimum length for a gap in the alignment to be registered as an intron
--early_stop (default -1): An option to terminate the function after a given number of records (when negative this is disabled, hence the default of -1)
```


# perform_enhanced_pileup.py

This tool performs an "enhanced" pileup on a bam file and exports a compressed CSV. 

It is "enhanced" in the sense that, unlike some alternative pileup tools, it includes information about insertions.

There is currently no option by which to filter for specific genomic regions, but this could be added if requested.

## Usage
```
--bam/-b: This is the bam file
--output/-o: This is the output filename

--min_aligned_length/-l (Default 20): Records with alignments shorter than this are ignored. For Nanopore, a much higher value (eg 500) may sometimes be applicable.
--primary (Default False): Add this to only consider primary alignments (often secondary alignments are ignored anyway as they do not have sequence information)
--min_intron_length (Default 50): 'Deletions' of this length or more are ignored in the pileup. These are common in RNA seq where they correspond to introns.
