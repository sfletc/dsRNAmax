package main

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"os"
	"strings"
)

// HeaderRef is a struct comprising a reference sequence header, seques and reverse complement
type HeaderRef struct {
	Header     string
	Seq        string
	ReverseSeq string
}

// RefLoad loads a reference sequence DNA file (FASTA format).
// It returns a slice of HeaderRef structs (individual reference header, sequence and reverse complement).
// Lower case nucleotides are converted to uppercase.
func RefLoad(refFile string) []*HeaderRef {
	var refSlice []*HeaderRef
	var singleHeaderRef *HeaderRef
	var header string
	var refSeq bytes.Buffer
	f, err := os.Open(refFile)
	if err != nil {
		fmt.Println("Problem opening fasta reference file " + refFile)
		errorShutdown()
	}
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		fastaLine := scanner.Text()
		switch {
		case strings.HasPrefix(fastaLine, ">"):
			seq := refSeq.String()
			singleHeaderRef = &HeaderRef{header, seq, reverseComplement(seq)}
			refSlice = append(refSlice, singleHeaderRef)
			header = fastaLine[1:]
			refSeq.Reset()
		case len(fastaLine) != 0:
			refSeq.WriteString(strings.ToUpper(fastaLine))
		}
	}
	seq := refSeq.String()
	singleHeaderRef = &HeaderRef{header, seq, reverseComplement(seq)}
	refSlice = append(refSlice, singleHeaderRef)
	refSlice = refSlice[1:]
	fmt.Println("     --->", len(refSlice), "sequences loaded")
	f.Close()
	return refSlice
}

// biasMod adds additional copies of the selected HeaderRef to a new HeaderREf slice.
// This biases kmer selection toward that HeaderRef at the cost of the overall geomean
func biasMod(ref []*HeaderRef, header string, bias int) ([]*HeaderRef, error) {
	var biasRef []*HeaderRef
	headerPresent := false
	for _, hr := range ref {
		biasRef = append(biasRef, hr)
		if hr.Header == header {
			headerPresent = true
			for i := 0; i < bias; i++ {
				biasRef = append(biasRef, hr)
			}
		}
	}
	if headerPresent {
		return biasRef, nil
	} else {
		return biasRef, errors.New("bias reference header not present in input file")
	}
}

// Reverse complementary DNA sequence
// Nucleotides must be upper case
// Only complements to ACGT are substituted.  Others remain the same.
func reverseComplement(s string) string {
	var dnaComplement = strings.NewReplacer(
		"A", "T", "T", "A", "G", "C", "C", "G",
	)
	complement := dnaComplement.Replace(s)
	reverseComplement := make([]byte, len(complement))
	for i, j := 0, len(reverseComplement)-1; i < len(reverseComplement); i, j = i+1, j-1 {
		reverseComplement[i] = complement[j]
	}
	return string(reverseComplement)
}
