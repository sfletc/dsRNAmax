# dsRNAmax: a multi-target double-stranded RNA design package 

dsRNAmax (dsRNAmaximizer) uses a kmer-based approach for multi-target dsRNA design and off-target avoidance. The package maximizes the number of kmers (default = 21nt, a common length of DICER-processed siRNAs) that perfectly match each input target sequence while avoiding any contiguous match of a specified length to any off-target sequence (default = 21nt, but this off-target kmer length can be less than the target kmer length).
    

![Alt text](./bioinf_github.jpg "dsRNA design")


# Installation

- Precompiled binaries for Windows, Linux and MacOS (Arm and Intel) can be downloaded from the [Releases page](https://github.com/sfletc/dsRNAmax/releases). 

- For Linux and MacOS, ```chmod +x {executable}``` may be required.  

- To compile from source, [install Go](https://go.dev/doc/install), clone this repository, and build with ```go build``` - this will generate an executable for the operating system it was built on, with build dependencies downloaded automatically.      


# Usage

- The only required input is a FASTA file containing the target sequences the dsRNA should be optimised for. This will generate a 300nt dsRNA sense arm sequence that attempts to maximise the number of dsRNA-derived 21nt kmers that exactly match each input sequence.

### Help (-h) 

Command:
```
dsRNAmax -h
```
Output:
```
dsRNAmax - dsRNA maximizer
Version:         1.1.12

Usage of dsRNAmax:
  -biasHeader string
    	Header of target sequence to bias toward
  -biasLvl int
    	Level of bias to apply
  -constructLen int
    	dsRNA sense arm length (default 300)
  -csv string
    	CSV file name (optional)
  -iterations int
    	No. of iterations (default 100)
  -kmerLen int
    	Kmer length (default 21)
  -offTargetKmers string
    	Path to off-target kmer file (optional)
  -offTargets string
    	Comma-separated list of off-target FASTA file/s
  -otKmerLen int
    	Off-target Kmer length (must be <= kmer length) (default 21)
  -targets string
    	Path to target FASTA file (required)

```
-----

### dsRNA design for multiple target sequences

An example command and command line output is shown for the western corn rootworm and southern corn rootworm vATPase-A transcripts using the ```-targets``` flag.  Statistics for the output dsRNA include the number of sense-arm derived kmers perfectly matching each target sequence (maximum = dsRNA length - kmer length +1), the Smith-Waterman-Gotoh similarity of the sense arm to each target sequence, mean kmer GC content, and the percentage of kmers (in both orientations) with a 5'U, 5'A and 5'C.  The geometric mean of kmers matching to each target sequences along with the sense arm GC content are also shown.  

Command:
```
dsRNAmax -targets wstrn_sthrn_corn_rootowrm_vATPaseA.fa 

```
Output:
```

dsRNAmax - dsRNA maximizer (Version: 1.1.12)
Target FASTA File: wstrn_sthrn_corn_rootowrm_vATPaseA.fa
Loading target sequences...
   ---> 2 sequences loaded
Getting target sequence kmers...
3,358 target kmers loaded
Finding best construct...
Identified optimal construct of length 300 with a median kmer match of 235

Results:
+------------------------+--------------+--------------------+------------------+---------+---------+---------+
| TARGET SEQUENCE HEADER | 21NT MATCHES | SWG SIMILARITY (%) | KMER MEAN GC (%) | 5'U (%) | 5'A (%) | 5'C (%) |
+------------------------+--------------+--------------------+------------------+---------+---------+---------+
| WCR_vATPase_A          | 217          | 95.0               | 43.2             | 27.6    | 30.4    | 26.3    |
| SCR_vATPase_A          | 253          | 98.0               | 42.7             | 26.9    | 30.8    | 24.9    |
+------------------------+--------------+--------------------+------------------+---------+---------+---------+


Median of kmer hits to each target sequence: 235

dsRNA sense-arm sequence - 43.3% GC content
ATCGGAGATGAAGAGAAGGAAGGGCAGTATGGTTACGTCCATGCTGTCTCAGGTCCAGTCGTTACTGCTGAGAAAATGTCTGGTTCTGCTATGTACGAACTGGTACGTGTCGGATACTATGAGCTGGTAGGAGAAATCATTAGATTGGAAGGTGACATGGCTACTATTCAGGTATATGAAGAAACTTCTGGTGTAACCGTTGGTGATCCAGTATTAAGAACTGGTAAACCACTTTCAGTAGAACTTGGACCTGGTATTATGGGTTCCATTTTTGATGGTATCCAACGTCCATTGAAAGAC

```
-----

## Addition of off-target sequences to avoid

Additionally, off-target sequences can be added as a single FASTA file.  In this example, the transcriptome of the beneficial seven-spotted ladybeetle is added using the ```-offTargets``` flag.  No 21nt match from the generated dsRNA (in either orientation) will perfectly match any off-target sequence.

Command:
```
dsRNAmax -targets wstrn_sthrn_corn_rootowrm_vATPaseA.fa -offTargets 7_spotted_ladybird.fa 
```
Output:
```
dsRNAmax - dsRNA maximizer (Version: 1.1.12)
Target FASTA File: wstrn_sthrn_corn_rootowrm_vATPaseA.fa
Off-target FASTA File: 7_spotted_ladybird.fa
Loading target sequences...
   ---> 2 sequences loaded
Getting target sequence kmers...
3,358 target kmers loaded
Removing off-target kmers from FASTA files...
Total off-target-matching kmers removed: 35

2024/04/17 16:44:06 Finding best construct...
2024/04/17 16:44:07 Identified optimal construct of length 300 with a median kmer match of 227.00

Results:
+------------------------+--------------+--------------------+------------------+---------+---------+---------+
| TARGET SEQUENCE HEADER | 21NT MATCHES | SWG SIMILARITY (%) | KMER MEAN GC (%) | 5'U (%) | 5'A (%) | 5'C (%) |
+------------------------+--------------+--------------------+------------------+---------+---------+---------+
| WCR_vATPase_A          | 226          | 97.0               | 34.6             | 35.5    | 30.3    | 19.3    |
| SCR_vATPase_A          | 228          | 95.0               | 37.5             | 38.9    | 25.2    | 18.1    |
+------------------------+--------------+--------------------+------------------+---------+---------+---------+

Median of kmer hits to each target sequence: 227

dsRNA sense-arm sequence - 37.3% GC content
AGGTAAAGCATCTCTAGCAGAAACGGACAAAATCACCTTGGAAATTGCCAGGCTTCTTAAAGAAGATTTCTTGCAACAAAACTCATACTCTTCTTATGACAGATTCTGTCCATTCTATAAAACTGTCGGTATGTTGAGAAACATGATCGGTTTGTACGATATGGCGAGACACGCCGTAGAATCAACCGCACAATCAGAAAATAAGATCACTTGGAACGTAATAAGAGATTCAATGAGTGGAATTTTATATCAACTTAGCAGTATGAAATTTAAGGATCCCGTAAAAGATGGTGAAGCTAA

```
----

### Bias toward a particular sequence

In some cases, it's desirable to maximise the number of kmers matching a particular sequence, while still maintaining effectiveness against other input targets.  This can be achieved by using ```-biasLvL``` and ```-biasHeader```.  For ```-biasHeader```, the full header (excluding ">") should be entered - use quotes if there are spaces.  For ```-biasLvl```, input an integer for the degree of bias to apply.  The integer used will add additional copies of the selected sequence to the design process, with its effect depended on the total number of input target sequences, so it's worth trialling different degrees of bias (starting at 1).  

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
+------------------------+--------------+--------------------+------------------+---------+---------+---------+
| TARGET SEQUENCE HEADER | 21NT MATCHES | SWG SIMILARITY (%) | KMER MEAN GC (%) | 5'U (%) | 5'A (%) | 5'C (%) |
+------------------------+--------------+--------------------+------------------+---------+---------+---------+
| WCR_vATPase_A          | 217          | 95.0               | 43.2             | 27.6    | 30.4    | 26.3    |
| SCR_vATPase_A          | 253          | 98.0               | 42.7             | 26.9    | 30.8    | 24.9    |
+------------------------+--------------+--------------------+------------------+---------+---------+---------+


Geometric mean of kmer hits to each target sequence: 230.65125189341592

dsRNA sense-arm sequence - 43.0% GC content
ATCGGAGATGAAGAGAAGGAAGGGCAGTATGGTTATGTCCATGCTGTCTCAGGTCCAGTCGTTACTGCTGAGAAAATGTCTGGTTCTGCTATGTACGAACTGGTACGTGTCGGATACTATGAGCTGGTAGGAGAAATCATTAGATTGGAAGGTGACATGGCTACTATTCAGGTATACGAAGAAACATCAGGTGTAACTGTTGGTGATCCAGTATTAAGAACTGGTAAACCACTTTCAGTAGAACTTGGACCTGGTATTATGGGTTCCATTTTTGATGGTATCCAACGTCCATTGAAAGAC
```

# Troubleshooting

Why was no dsRNA generated?
- If the off-target sequences are too similar to the target sequences, no construct of the specified sense arm length can be generated.  Try reducing the kmer length and/or the construct length.  


# Cite

Forthcoming.  
