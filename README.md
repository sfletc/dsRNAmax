# dsRNAd: a multi-target double-stranded RNA design package for enhanced crop protection

Crop protection is undergoing significant evolution, with a move toward sustainable approaches that do not adversely impact the environment or human health.  RNA interference (RNAi) is being used via application or expression of double-stranded RNA (dsRNA) to target viral, fungal and insect pests and pathogens.  

dsRNAd (dsRNA designer) is a software package using a biologically relevant k-mer based approach for multi-target dsRNA design while mitigating unintended impacts to beneficial organisms.  

Where a specific target is considered a priority, the design algorithm can be biased toward that target.  This package aims to increase the simplicity, efficacy and cost-effectiveness RNAi-based crop protection approaches, aiding translation to market.    

![Alt text](./bioinf_github.jpg "Title")


# Installation

- Precompiled binaries for Windows, Linux and MacOS (Arm and Intel) can be downloaded from the [Releases page](https://github.com/sfletc/dsRNAd/releases)
- To compile from source, [install Go](https://go.dev/doc/install), clone this repository, and build with ```go build``` - this will generate an executable for the operating system it was built on, with build dependencices downloaded automatically.      


# Usage

- The only required input is a FASTA file containing the target sequences the dsRNA should be optimised for. This will generate a 300nt dsRNA sense arm sequence that attempts to maximise the number of dsRNA-derived 21nt kmers that exactly match each input sequence.

Command:
```
 .\dsRNAd_Win.exe -h
```
Output:
```
dsRNAd - dsRNA designer
Version:         1.0.1

Usage of .\dsRNAd_Win.exe:
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

An example command and command line output is shown for the western corn rootworm and southern corn rootworm vATPase-A transcripts using the ```-targets``` flag.

Command:
```
.\dsRNAd_Win.exe -targets .\cr.fa
```
Output:
```

dsRNAd - dsRNA designer
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

Additionally, off-target sequences can be added as a single FASTA file.  In this example, the transcriptome of the beneficial seven-spotted ladybeetle is added using the ```-offTargets``` flag.  No 21nt match from the generated dsRNA (in either orientation) will perfectly match any off-target sequence.

Command:
```
.\dsRNAd_Win.exe -targets .\cr.fa -offTargets .\GCF_907165205.1_icCocSept1.1_rna.fna
```
Output:
```

dsRNAd - dsRNA designer
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


# Troubleshooting

Why was no dsRNA generated?
- if the off-target sequences are too similar to the target sequences, no construct of the specified sense arm length can be generated.  Try reducing the kmer length and/or the construct length.  


# Cite

