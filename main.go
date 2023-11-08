package main

import (
	"errors"
	"flag"
	"fmt"
	"log"
	"os"
	"strings"
)

var Version = "1.0.3"

func errorShutdown() {
	fmt.Println("\nExiting program")
	os.Exit(1)
}

func clInput() (*string, *string, *int, *int, *int, *int, *string, *int, error) {
	refFile := flag.String("targets", "", "Path to target FASTA file (required)")
	otRefFiles := flag.String("offTargets", "", "Comma-separated list of off-target FASTA file/s")
	kmerLength := flag.Int("kmerLen", 21, "Kmer length")
	otKmerLength := flag.Int("otKmerLen", *kmerLength, "Off-target Kmer length (must be <= kmer length)")
	consLength := flag.Int("constructLen", 300, "dsRNA sense arm length")
	iterations := flag.Int("iterations", 100, "No. of iterations")
	biasHeader := flag.String("biasHeader", "", "Header of target sequence to bias toward")
	biasLvl := flag.Int("biasLvl", 0, "Level of bias to apply")
	flag.Parse()
	if *refFile == "" {
		return refFile, otRefFiles, kmerLength, otKmerLength, consLength, iterations, biasHeader, biasLvl, errors.New("error: no target FASTA file was specificed")
	}
	return refFile, otRefFiles, kmerLength, otKmerLength, consLength, iterations, biasHeader, biasLvl, nil
}

func main() {
	fmt.Println("\ndsRNAmax - dsRNA maximizer")
	fmt.Println("Version:\t", Version)
	fmt.Println("")
	refFile, otRefFiles, kmerLength, otKmerLength, consLength, iterations, biasHeader, biasLvl, err := clInput()
	if err != nil {
		log.Fatal(err)
	}
	//TODO: check if otKmer length <= kmerLength
	fmt.Println("Target FASTA File:\t", *refFile)
	if *otRefFiles != "" {
		fmt.Println("Off-target FASTA File:\t", *otRefFiles)
	}
	fmt.Println("Loading target sequences")
	ref := RefLoad(*refFile)
	if *biasHeader != "" {
		ref, err = biasMod(ref, *biasHeader, *biasLvl)
		if err != nil {
			log.Fatal(err)
		}
	}
	fmt.Println("Getting target sequence kmers")
	goodKmers := getKmers(ref, *kmerLength)
	if *otRefFiles != "" {
		fmt.Println("Loading and removing off-target sequences")
		// if *otKmerLength < *kmerLength {
		// 	goodKmers = otShortRemoval(goodKmers, *otKmerLength, *otRefFiles)
		// } else {
		// 	goodKmers = otRemoval(goodKmers, *otKmerLength, *otRefFiles)
		// }
		files := strings.Split(*otRefFiles, ",")
		ConcurrentlyProcessSequences(files, goodKmers, *otKmerLength)
	}

	fmt.Println("Counting kmers")
	kmerCts := kmerAbun(goodKmers)
	fmt.Println("Finding best construct")
	selConstruct := conBestConstruct(goodKmers, kmerCts, *kmerLength, len(ref), *consLength, *iterations)
	if selConstruct != nil {
		outputResults(goodKmers, kmerLength, selConstruct, ref)
	} else {
		fmt.Println("Could not identify a dsRNA sense arm sequence.  Check input format and sequence lengths")
	}
}
