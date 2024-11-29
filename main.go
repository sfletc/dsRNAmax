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

var Version = "1.1.12"

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
	log.Printf("dsRNAmax - dsRNA maximizer (Version: %s)\n", Version)

	refFile, otRefFiles, kmerLength, otKmerFile, otKmerLength, consLength, iterations, biasHeader, biasLvl, csv, err := clInput()
	if err != nil {
		log.Fatal(err)
	}

	if *otKmerLength > *kmerLength && *otKmerFile == "" {
		log.Fatalf("Off-target kmer length (%d) must be <= kmer length (%d)", *otKmerLength, *kmerLength)
	}

	if _, err := os.Stat(*refFile); os.IsNotExist(err) {
		log.Fatalf("Target FASTA file does not exist: %s", *refFile)
	}

	log.Printf("Target FASTA File: %s", *refFile)
	if *otRefFiles != "" {
		log.Printf("Off-target FASTA File: %s", *otRefFiles)
	}

	log.Println("Loading target sequences...")
	ref := RefLoad(*refFile)

	if *biasHeader != "" {
		log.Printf("Applying bias modification to sequence '%s' at level %d...", *biasHeader, *biasLvl)
		ref, err = biasMod(ref, *biasHeader, *biasLvl)
		if err != nil {
			log.Fatal(err)
		}
	}

	log.Println("Getting target sequence kmers...")
	goodKmers := getKmers(ref, *kmerLength)
	log.Printf("%s target kmers loaded\n", intWithCommas(len(goodKmers)))

	if *otRefFiles != "" && *otKmerFile != "" {
		log.Fatalln("Error: both off-target FASTA files and an off-target kmer file specified. Please specify only one.")
	}

	if *otRefFiles != "" {
		log.Println("Removing off-target kmers from FASTA files...")
		files := strings.Split(*otRefFiles, ",")
		for _, file := range files {
			if _, err := os.Stat(file); os.IsNotExist(err) {
				log.Fatalf("Off-target FASTA file does not exist: %s", file)
			}
		}
		removeOffTargetKmersFromFasta(files, goodKmers, *kmerLength, *otKmerLength)
	}

	if *otKmerFile != "" {
		log.Println("Removing off-target kmers from kmer file...")
		err := removeOffTargetKmersFromFile(goodKmers, *otKmerFile, *kmerLength)
		if err != nil {
			log.Fatal(err)
		}
	}

	log.Println("Finding best construct...")
	kmerCts := kmerAbun(goodKmers)
	selConstruct := conBestConstruct(goodKmers, kmerCts, *kmerLength, len(ref), *consLength, *iterations)
	if selConstruct != nil {

		outputResults(goodKmers, kmerLength, selConstruct, ref, *csv)
	} else {
		log.Println("Could not identify a dsRNA sense arm sequence. Check input format, increase OT kmer length, and/or try a shorter construct length")
		os.Exit(1)
	}
}

func removeOffTargetKmersFromFasta(files []string, goodKmers map[string][]int, kmerLength int, otKmerLength int) {
	ConcurrentlyProcessSequences(files, goodKmers, kmerLength, otKmerLength)
}

func removeOffTargetKmersFromFile(goodKmers map[string][]int, otKmerFile string, kmerLength int) error {
	return removeOffTargetKmersFromGoodKmers(goodKmers, otKmerFile, kmerLength)
}
