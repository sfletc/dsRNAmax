# dsRNAmax: a multi-target double-stranded RNA design package 

dsRNAmax (dsRNAmaximizer) uses a kmer-based approach for multi-target dsRNA design and off-target avoidance. The package maximizes the number of kmers (default = 21nt, a common length of DICER-processed siRNAs) that perfectly match each input target sequence while avoiding any contiguous match of a specified length to any off-target sequence (default = 21nt, but this off-target kmer length can be less than the target kmer length).
    

![Alt text](./bioinf_github.jpg "dsRNA design")


# Installation

- Precompiled binaries for Windows, Linux and MacOS (Arm and Intel) can be downloaded from the [Releases page](https://github.com/sfletc/dsRNAmax/releases). For Linux and MacOS, ```chmod +x {executable}``` may be required.  

- To compile from source, [install Go](https://go.dev/doc/install), clone this repository, and build with ```go build``` - this will generate an executable for the operating system it was built on, with build dependencies downloaded automatically.      


# Usage

- The only required input is a FASTA file containing the target sequences the dsRNA should be optimised for. This will generate a 300nt dsRNA sense arm sequence that attempts to maximise the number of dsRNA-derived 21nt kmers that exactly match each input sequence.

### Help (-h) 

Command:
```
 .\dsRNAmax_Win.exe -h
```
Output:
```
dsRNAmax - dsRNA maximizer
Version:         1.1.1.

Usage of .\dsRNAmax_Win.exe:
  -biasHeader string
        Header of target sequence to bias toward
  -biasLvl int
        Level of bias to apply
  -constructLen int
        dsRNA sense arm length (default 300)
  -iterations int
        No. of iterations (default 100)
  -kmerLen int
        Kmer length (default 21)
  -offTargets string
        Path to off-target FASTA file
  -targets string
        Path to target FASTA file (required)
```
-----

### dsRNA design for multiple target sequences

An example command and command line output is shown for the western corn rootworm and southern corn rootworm vATPase-A transcripts using the ```-targets``` flag.  Statistics for the output dsRNA include the number of sense-arm derived kmers perfectly matching each target sequence (maximum = dsRNA length - kmer length +1), the Smith-Waterman-Gotoh similarity of the sense arm to each target sequence, mean kmer GC content, and the percentage of kmers (in both orientations) with a 5'U, 5'A and 5'C.  The geometric mean of kmers matching to each target sequences along with the sense arm GC content are also shown.  

Command:
```
.\dsRNAmax_Win.exe -targets .\cr.fa
```
Output:
```

dsRNAmax - dsRNA maximizer
Version:         1.0.1

.\cr.fa
Loading target sequences
     ---> 2 sequences loaded
Getting target sequence kmers
Counting kmers
Finding best construct

Results:
+--------------------------------+--------------+--------------------+------------------+---------+---------+---------+
|     TARGET SEQUENCE HEADER     | 21NT MATCHES | SWG SIMILARITY (%) | KMER MEAN GC (%) | 5'U (%) | 5'A (%) | 5'C (%) |
+--------------------------------+--------------+--------------------+------------------+---------+---------+---------+
| XM_028294206.1 PREDICTED:      | 232          | 97.0               | 41.9             | 29.3    | 30.2    | 25.0    |
| Diabrotica virgifera virgifera |              |                    |                  |         |         |         |
| V-type proton ATPase catalytic |              |                    |                  |         |         |         |
| subunit A (LOC114343389), mRNA |              |                    |                  |         |         |         |
| KX982002.1 Diabrotica          | 238          | 96.0               | 43.9             | 26.1    | 30.3    | 26.5    |
| undecimpunctata howardi        |              |                    |                  |         |         |         |
| vacuolar ATPase subunit A      |              |                    |                  |         |         |         |
| mRNA, complete cds             |              |                    |                  |         |         |         |
+--------------------------------+--------------+--------------------+------------------+---------+---------+---------+

Geometric mean of kmer hits to each target sequence: 234.98085028359225

dsRNA sense-arm sequence - 43.3% GC content
ATCGGAGATGAAGAGAAGGAAGGGCAGTATGGTTACGTCCATGCTGTCTCAGGTCCAGTCGTTACTGCTGAGAAAATGTCTGGTTCTGCTATGTACGAACTGGTGCGTGTAGGATACTATGAGCTGGTAGGAGAAATCATTAGATTGGAAGGTGACATGGCTACTATTCAGGTATACGAAGAAACATCAGGTGTAACTGTTGGTGATCCAGTATTAAGAACTGGTAAACCACTTTCAGTAGAACTTGGACCTGGTATTATGGGTTCCATTTTTGATGGTATCCAACGTCCATTGAAAGAC
```
-----

## Addition of off-target sequences to avoid

Additionally, off-target sequences can be added as a single FASTA file.  In this example, the transcriptome of the beneficial seven-spotted ladybeetle is added using the ```-offTargets``` flag.  No 21nt match from the generated dsRNA (in either orientation) will perfectly match any off-target sequence.

Command:
```
.\dsRNAmax_Win.exe -targets .\cr.fa -offTargets .\GCF_907165205.1_icCocSept1.1_rna.fna
```
Output:
```

dsRNAmax - dsRNA maximizer
Version:         1.0.1

.\cr.fa
Loading target sequences
     ---> 2 sequences loaded
Getting target sequence kmers
Loading and removing off-target sequences
     ---> 25480 sequences loaded
Counting kmers
Finding best construct

Results:
+--------------------------------+--------------+--------------------+------------------+---------+---------+---------+
|     TARGET SEQUENCE HEADER     | 21NT MATCHES | SWG SIMILARITY (%) | KMER MEAN GC (%) | 5'U (%) | 5'A (%) | 5'C (%) |
+--------------------------------+--------------+--------------------+------------------+---------+---------+---------+
| XM_028294206.1 PREDICTED:      | 228          | 97.0               | 34.6             | 35.5    | 30.3    | 19.3    |
| Diabrotica virgifera virgifera |              |                    |                  |         |         |         |
| V-type proton ATPase catalytic |              |                    |                  |         |         |         |
| subunit A (LOC114343389), mRNA |              |                    |                  |         |         |         |
| KX982002.1 Diabrotica          | 226          | 95.0               | 37.3             | 37.2    | 27.0    | 19.0    |
| undecimpunctata howardi        |              |                    |                  |         |         |         |
| vacuolar ATPase subunit A      |              |                    |                  |         |         |         |
| mRNA, complete cds             |              |                    |                  |         |         |         |
+--------------------------------+--------------+--------------------+------------------+---------+---------+---------+

Geometric mean of kmer hits to each target sequence: 226.99779734614165

dsRNA sense-arm sequence - 37.3% GC content
GCAGAAACGGACAAAATCACCTTGGAAATTGCCAGGCTTCTTAAAGAAGATTTCTTGCAACAAAACTCATACTCTTCTTATGACAGATTCTGTCCATTCTATAAAACTGTCGGTATGTTGAGAAACATGATCGGTTTGTACGATATGGCGAGACACGCCGTAGAATCAACCGCACAATCAGAAAATAAGATCACTTGGAACGTAATAAGAGATTCAATGAGTGGAATTTTATATCAACTTAGCAGTATGAAATTTAAGGATCCCGTAAAAGATGGTGAAGCTAAGATCAAGGCAGATTTT
```
----

### Bias toward a particular sequence

In some cases, it's desirable to maximise the number of kmers matching a particular sequence, while still maintaining effectiveness against other input targets.  This can be achieved by using ```-biasLvL``` and ```-biasHeader```.  For ```-biasHeader```, the full header (excluding ">") should be entered - use quotes if there are spaces (see below for example).  For ```-biasLvl```, input an integer for the degree of bias to apply.  The integer used will add additional copies of the selected sequence to the design process, with its effect depended on the total number of input target sequences, so it's worth trialling different degrees of bias (starting at 1).  

```
.\dsRNAmax_Win.exe -targets .\software\dsRNAConsDesigner\cr.fa -biasLvl 1 -biasHeader "XM_028294206.1 PREDICTED: Diabrotica virgifera virgifera V-type proton ATPase catalytic subunit A (LOC114343389), mRNA"
```

Output:

```
dsRNAmax - dsRNA maximizer
Version:         1.0.1

.\cr.fa
Loading target sequences
     ---> 2 sequences loaded
Getting target sequence kmers
Counting kmers
Finding best construct

Results:
+--------------------------------+--------------+--------------------+------------------+---------+---------+---------+
|     TARGET SEQUENCE HEADER     | 21NT MATCHES | SWG SIMILARITY (%) | KMER MEAN GC (%) | 5'U (%) | 5'A (%) | 5'C (%) |
+--------------------------------+--------------+--------------------+------------------+---------+---------+---------+
| XM_028294206.1 PREDICTED:      | 280          | 100.0              | 42.9             | 27.1    | 30.0    | 26.1    |
| Diabrotica virgifera virgifera |              |                    |                  |         |         |         |
| V-type proton ATPase catalytic |              |                    |                  |         |         |         |
| subunit A (LOC114343389), mRNA |              |                    |                  |         |         |         |
| KX982002.1 Diabrotica          | 190          | 93.0               | 42.4             | 28.4    | 31.1    | 24.7    |
| undecimpunctata howardi        |              |                    |                  |         |         |         |
| vacuolar ATPase subunit A      |              |                    |                  |         |         |         |
| mRNA, complete cds             |              |                    |                  |         |         |         |
+--------------------------------+--------------+--------------------+------------------+---------+---------+---------+

Geometric mean of kmer hits to each target sequence: 230.65125189341592

dsRNA sense-arm sequence - 43.0% GC content
ATCGGAGATGAAGAGAAGGAAGGGCAGTATGGTTATGTCCATGCTGTCTCAGGTCCAGTCGTTACTGCTGAGAAAATGTCTGGTTCTGCTATGTACGAACTGGTACGTGTCGGATACTATGAGCTGGTAGGAGAAATCATTAGATTGGAAGGTGACATGGCTACTATTCAGGTATACGAAGAAACATCAGGTGTAACTGTTGGTGATCCAGTATTAAGAACTGGTAAACCACTTTCAGTAGAACTTGGACCTGGTATTATGGGTTCCATTTTTGATGGTATCCAACGTCCATTGAAAGAC
```

# Troubleshooting

Why was no dsRNA generated?
- If the off-target sequences are too similar to the target sequences, no construct of the specified sense arm length can be generated.  Try reducing the kmer length and/or the construct length.  


# Cite

Forthcoming.  
