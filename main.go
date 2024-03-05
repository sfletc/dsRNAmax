package main

import (
	"errors"
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
)

var Version = "1.0.3"

func errorShutdown() {
	fmt.Println("\nExiting program")
	os.Exit(1)
}

// reverseSlice reverses a slice of strings. Needed for the intWithCommas function.
func reverseSlice(s []string) []string {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
	return s
}

// intWithCommas converts an integer to a comma-separated string.
func intWithCommas(i int) string {
	in := strconv.Itoa(i) // Convert int to string
	out := make([]string, 0, len(in)/3+1)

	// Loop through the string, inserting commas every three digits.
	for len(in) > 3 {
		threeDigits := in[len(in)-3:]
		out = append(out, threeDigits)
		in = in[:len(in)-3]
	}
	out = append(out, in)

	// Since we built the out slice from right to left, reverse it.
	out = reverseSlice(out)

	// Join the slice back into a single string.
	return strings.Join(out, ",")
}

func clInput() (*string, *string, *int, *string, *int, *int, *int, *string, *int, *string, error) {
	refFile := flag.String("targets", "", "Path to target FASTA file (required)")
	otRefFiles := flag.String("offTargets", "", "Comma-separated list of off-target FASTA file/s")
	otKmerFile := flag.String("offTargetKmers", "", "Path to off-target kmer file (optional)")
	kmerLength := flag.Int("kmerLen", 21, "Kmer length")
	otKmerLength := flag.Int("otKmerLen", *kmerLength, "Off-target Kmer length (must be <= kmer length)")
	consLength := flag.Int("constructLen", 300, "dsRNA sense arm length")
	iterations := flag.Int("iterations", 100, "No. of iterations")
	biasHeader := flag.String("biasHeader", "", "Header of target sequence to bias toward")
	biasLvl := flag.Int("biasLvl", 0, "Level of bias to apply")
	csv := flag.String("csv", "", "CSV file name (optional)")
	flag.Parse()
	if *refFile == "" {
		return refFile, otRefFiles, kmerLength, otKmerFile, otKmerLength, consLength, iterations, biasHeader, biasLvl, csv, errors.New("error: no target FASTA file was specificed")
	}
	return refFile, otRefFiles, kmerLength, otKmerFile, otKmerLength, consLength, iterations, biasHeader, biasLvl, csv, nil
}

func main() {
	fmt.Println("\ndsRNAmax - dsRNA maximizer")
	fmt.Println("Version:\t", Version)
	fmt.Println("")
	refFile, otRefFiles, kmerLength, otKmerFile, otKmerLength, consLength, iterations, biasHeader, biasLvl, csv, err := clInput()
	if err != nil {
		log.Fatal(err)
	}
	//TODO: check if otKmer length <= kmerLength
	fmt.Println("Target FASTA File:\t", *refFile)
	if *otRefFiles != "" {
		fmt.Println("Off-target FASTA File:\t", *otRefFiles)
	}
	fmt.Println("\nLoading target sequences")
	ref := RefLoad(*refFile)
	if *biasHeader != "" {
		ref, err = biasMod(ref, *biasHeader, *biasLvl)
		if err != nil {
			log.Fatal(err)
		}
	}
	fmt.Printf("Getting target sequence kmers\n\n")
	goodKmers := getKmers(ref, *kmerLength)
	fmt.Printf("%s target kmers loaded\n\n", intWithCommas(len(goodKmers)))

	switch {

	case *otRefFiles != "" && *otKmerFile != "":
		fmt.Println("Error: both off-target FASTA files and an off-target kmer file specified.  Please specify only one.")
		errorShutdown()
	case *otRefFiles != "":
		fmt.Printf("Removing off-target kmers\n\n")
		files := strings.Split(*otRefFiles, ",")
		ConcurrentlyProcessSequences(files, goodKmers, *kmerLength, *otKmerLength)
	case *otKmerFile != "":
		fmt.Printf("Removing off-target kmers\n\n")
		removeOffTargetKmersFromGoodKmers(goodKmers, *otKmerFile, *kmerLength)
	}

	fmt.Println("Counting kmers")
	kmerCts := kmerAbun(goodKmers)
	fmt.Println("Finding best construct")
	selConstruct := conBestConstruct(goodKmers, kmerCts, *kmerLength, len(ref), *consLength, *iterations)
	if selConstruct != nil {
		outputResults(goodKmers, kmerLength, selConstruct, ref, *csv)
	} else {
		fmt.Println("Could not identify a dsRNA sense arm sequence.  Check input format and sequence lengths")
	}
}
