# nano_tools

nano_tools is a collection of several python scripts I have generated that I find useful for analysing nanopore data.

# demultiplex_nanopore_barcodes.py

This tool demultiplexes nanopore libraries based on barcodes in the 5' or 3' ends. In particular, it was designed to work with libraries generated via PCR, in which the PCR primers contain barcodes. The barcodes can be in the forward or reverse primer, or can be in both primers (in which case combinatorial demultiplexing is performed).

This enables one to, for example, generate 96 uniquely barcoded PCR products using 8 forward and 12 reverse primers.

The algorithm accounts for Nanopore's relatively high error rate and lack of positional precision. For each read it:
1. Takes the first and last N bases (default 100)
2. Searches for (partial) matches to each barcode using the rapidfuzz library - this allows for insertions/deletions and mismatches. 
3. (Optionally) does the same for the reverse complement
4. Filters for matches that confidently match only a single unique barcode or barcode pair
5. Writes each out to a separate fastq.gz file

It performs a quick check using randomly generated sequences to ensure the parameters used do not lead to significant false discovery rates.

Note that no trimming of the barcodes is performed.

Performance is not especially rapid: roughly 500-1,000 reads are demultiplexed per second using standard parameters (around 3,000,000 per hour)

## Usage

```
--fastq/-f: This is the fastq file that you wish to demultiplex
--primers/-p: This is the csv of primers/barcodes used to generate the library (see below for formatting details)
--output/-o: This is the prefix (including directory) that files are written to

--min_score/-s (default 90): This is the minimum partial ratio score from rapidfuzz for a match to be detected
--max_ambiguity (default 82): If a second match is gains this score or more, the match is rejected because it is considered ambiguous
--length/-l (default 100): This is the length at the start and end of each read that is searched for barcodes
--ignore_rc (default False): This option ignores the reverse complement of the read (i.e. only looks at the forward read)
```

## Formatting the primers/barcodes csv

The csv must have 3 columns:
1. The first column is the primer name (eg "forward_1")
2. The second column is the primer sequence (this should be the 5'-3' sequence in the primer itself) (TODO check this)
3. The third column must either be "F" or "R" depending on whether it was a forward or reverse barcode

Here is an example:
```R1,GACATCAATTCGAACAATCC,R
R2,TAACTTACCTGTATAGCACA,R
R3,TTATATACATCTTCCTGGCT,R
R4,AATATGTGGATCGGACTCTA,R
R5,ACCTAACTCCTATTCATACG,R
R6,CACACAAGCATGTTAGTCTT,R
Rn,AGCTAGACATCGCGCAATAG,R
R7,ATTATCGGTAACATCTGCAC,R
F1,AGAATATCTTAGACACTTGC,F
F2,TGTTCCGTATTGCTTAACAA,F
F3,ACTAATCCAGAGTTCTCAAG,F
F4,AGGTAACATGTGATCCGATA,F
F5,ATCATTGGTTGTGTGGAATC,F
F6,TCAAGTATGAGGAACACCAA,F
Fn,ATTCAGTCGCGGACTAATCG,F
```
