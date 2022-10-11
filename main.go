package main

import (
	"errors"
	"flag"
	"fmt"
	"log"
	"os"
)

var Version = "1.0.1"


func errorShutdown() {
	fmt.Println("\nExiting program")
	os.Exit(1)
}

func clInput() (*string, *string, *int, *int, *int, *string, *int, error) {
	refFile := flag.String("targets", "", "Path to target FASTA file (required)")
	otRefFile := flag.String("offTargets", "", "Path to off-target FASTA file")
	kmerLength := flag.Int("kmerLen", 21, "Kmer length")
	consLength := flag.Int("constructLen", 300, "dsRNA sense arm length")
	iterations := flag.Int("iterations", 100, "No. of iterations")
	biasHeader := flag.String("biasHeader", "", "Header of target sequence to bias toward")
	biasLvl := flag.Int("biasLvl", 0, "Level of bias to apply")
	flag.Parse()
	if *refFile == "" {
		return refFile, otRefFile, kmerLength, consLength, iterations, biasHeader, biasLvl, errors.New("error: no target FASTA file was specificed")
	}
	return refFile, otRefFile, kmerLength, consLength, iterations, biasHeader, biasLvl, nil
}

func main() {
	fmt.Println("\ndsRNAmax - dsRNA maximizer")
	fmt.Println("Version:\t", Version)
	fmt.Println("")
	refFile, otRefFile, kmerLength, consLength, iterations, biasHeader, biasLvl, err := clInput()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println(*refFile)
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
	if *otRefFile != "" {
		fmt.Println("Loading and removing off-target sequences")
		goodKmers = otRemoval(otRefFile, goodKmers, kmerLength)
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
