package main

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"log"
	"os"
	"strings"
	"sync"
)

// removeOffTargetUint64KmersConcurrent removes off-target kmers from a set of good kmers using concurrent processing.
//
// Parameters:
//   - filename: The path to the file containing off-target kmers in binary format.
//   - goodUint64Kmers: A map of good kmers in uint64 representation.
//   - numWorkers: The number of worker goroutines to use for concurrent processing.
//
// Returns:
//   - A map of removed kmers, where the keys are the kmer sequences.
func removeOffTargetUint64KmersConcurrent(filename string, goodUint64Kmers map[uint64]struct{}, numWorkers int) (map[string]struct{}, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening file: %v", err)
	}
	defer file.Close()

	// Read the first line of the binary file to get the kmer length
	var k64 uint64
	if err := binary.Read(file, binary.LittleEndian, &k64); err != nil {
		return nil, fmt.Errorf("error reading kmer length: %v", err)
	}
	k := int(k64)

	const chunkSize = 1024 * 64 // 64 KB; adjust as needed
	kmerSize := binary.Size(uint64(0))
	kmerChan := make(chan []byte, numWorkers)
	removedKmerChan := make(chan map[string]struct{}, numWorkers)

	var wg sync.WaitGroup

	// Start worker goroutines to process chunks of k-mers
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			localRemovedKmers := make(map[string]struct{})

			for chunk := range kmerChan {
				processChunk(chunk, kmerSize, k, goodUint64Kmers, localRemovedKmers)
			}

			removedKmerChan <- localRemovedKmers
		}()
	}

	// Read large chunks of data and send to worker goroutines
	buf := make([]byte, chunkSize)
	for {
		bytesRead, err := file.Read(buf)
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, fmt.Errorf("error reading file: %v", err)
		}

		if bytesRead > 0 {
			chunk := make([]byte, bytesRead)
			copy(chunk, buf[:bytesRead])
			kmerChan <- chunk
		}
	}

	close(kmerChan)
	wg.Wait()
	close(removedKmerChan)

	// Merge results from workers
	mergedRemovedKmers := make(map[string]struct{})
	for localMap := range removedKmerChan {
		for kmer := range localMap {
			mergedRemovedKmers[kmer] = struct{}{}
		}
	}

	return mergedRemovedKmers, nil
}

func processChunk(chunk []byte, kmerSize, k int, goodUint64Kmers map[uint64]struct{}, localRemovedKmers map[string]struct{}) {
	for j := 0; j < len(chunk); j += kmerSize {
		var kmer uint64
		buf := bytes.NewReader(chunk[j : j+kmerSize])
		if err := binary.Read(buf, binary.LittleEndian, &kmer); err != nil {
			log.Printf("Error reading kmer: %v", err)
			continue
		}

		if _, found := goodUint64Kmers[kmer]; found {
			seq := kmerToSequence(kmer, k)
			localRemovedKmers[seq] = struct{}{}
		}
	}
}

// kmerToSequence converts a uint64 representation of a kmer to its corresponding DNA sequence string.
//
// Args:
//
//	kmer: The uint64 representation of the kmer.
//	k: The length of the kmer.
//
// Returns:
//
//	The DNA sequence string corresponding to the kmer.
func kmerToSequence(kmer uint64, k int) string {
	bases := "ACGT"
	sequence := make([]byte, k)

	for i := k - 1; i >= 0; i-- {
		sequence[i] = bases[kmer&3]
		kmer >>= 2
	}
	return string(sequence)
}

// removeKmersFromGoodKmers removes kmers (and their reverse complements) from the 'goodKmers' map if they are present in the 'removedKmers' map.
//
// Args:
//
//	goodKmers: A map where keys are target kmers and values are presence/absence slices.
//	removedKmers: A map where keys are kmers identified for removal (e.g., off-target kmers).
func removeKmersFromGoodKmers(goodKmers map[string][]int, removedKmers map[string]struct{}) {
	for kmer := range goodKmers {
		_, inRemoved := removedKmers[kmer]
		rc := reverseComplement(kmer)
		_, rcInRemoved := removedKmers[rc]

		// If the k-mer or its reverse complement is in removedKmers, delete it from goodKmers
		if inRemoved || rcInRemoved {
			delete(goodKmers, kmer)
		}
	}
}

// convertGoodKmersToUint64Set converts a map of string-based kmers to a map using their canonical uint64 representations (considering reverse complements).
//
// Args:
//
//	goodKmers: A map where keys are kmers as strings and values are presence/absence slices.
//	k: The length of the kmers.
//
// Returns:
//  1. A map[uint64]struct{} where keys are canonical uint64 representations of kmers.
//  2. An error if any occurs during conversion.
func convertGoodKmersToUint64Set(goodKmers map[string][]int, k int) (map[uint64]struct{}, error) {
	goodUint64Kmers := make(map[uint64]struct{})

	toShift := uint((k - 1) * 2)         // Number of bit positions to shift for reverse complement
	mask := (^uint64(0)) >> uint(64-k*2) // Mask to isolate k-mer bits

	for kmer := range goodKmers {
		var next, nextRC uint64 // next for the k-mer, nextRC for its reverse complement
		for _, b := range []byte(kmer) {
			val := ((b >> 1) ^ ((b & 4) >> 2)) & 3
			next = ((next << 2) | uint64(val)) & mask
			nextRC = ((nextRC >> 2) | ((^uint64(val))<<toShift)&mask)
		}
		// Store the lexicographically smaller of the k-mer and its reverse complement.
		canonicalKmer := min(next, nextRC)
		goodUint64Kmers[canonicalKmer] = struct{}{}
	}

	return goodUint64Kmers, nil
}

// min returns the smaller of two uint64 values (utility function).
func min(a, b uint64) uint64 {
	if a < b {
		return a
	}
	return b
}

// removeOffTargetKmersFromGoodKmers filters out off-target kmers from a set of good kmers (string-based) using a file of off-target kmer references.
// It handles kmer encoding, concurrent reading of off-target kmers, and final removal.
//
// Args:
//
//	goodKmers: A map where keys are target kmers as strings and values are presence/absence slices.
//	offTargetKmersFile: The path to a file containing off-target kmers (likely in uint64 representation).
//	goodKmerLength: The length of the kmers in the 'goodKmers' map.
//
// Returns:
//
//	An error if any occurs during the filtering process
func removeOffTargetKmersFromGoodKmers(goodKmers map[string][]int, offTargetKmersFile string, goodKmerLength int) error {
	file, err := os.Open(offTargetKmersFile)
	if err != nil {
		return fmt.Errorf("error opening file: %v", err)
	}
	defer file.Close()
	//read the first line of the binary file to get the kmer length
	var k64 uint64

	if err := binary.Read(file, binary.LittleEndian, &k64); err != nil {
		return fmt.Errorf("error reading kmer length: %v", err)
	}
	fmt.Println("Kmer length: ", goodKmerLength)
	var OTKmerLen int = int(k64)
	fmt.Println("OT kmer length: ", OTKmerLen)
	switch {
	case OTKmerLen > goodKmerLength:
		return fmt.Errorf("off-target kmer length is greater than target kmer length - it must be equal or lower")
	case OTKmerLen < goodKmerLength:
		removeOffTargetSubKmersFromGoodKmers(goodKmers, offTargetKmersFile, OTKmerLen)
	default:
		// Convert the good k-mers to canonical uint64 representation
		goodUint64Kmers, err := convertGoodKmersToUint64Set(goodKmers, goodKmerLength)
		if err != nil {
			return err
		}

		// Read off-target k-mers and build a map of removed k-mers
		removedKmers, err := removeOffTargetUint64KmersConcurrent(offTargetKmersFile, goodUint64Kmers, 32)
		if err != nil {
			return err
		}

		// Remove the k-mers found in removedKmers from the original goodKmers map
		removeKmersFromGoodKmers(goodKmers, removedKmers)
		fmt.Printf("Total off-target-matching kmers removed: %d\n\n", len(removedKmers))
	}
	return nil
}

// removeOffTargetSubKmersFromGoodKmers removes off-target subkmers from a set of good kmers (string-based) using a file of off-target kmer references.
//
// Parameters:
//   - goodKmers: A map where keys are target kmers as strings and values are presence/absence slices.
//   - offTargetKmersFile: The path to a file containing off-target kmers (likely in uint64 representation).
//   - subKmerLength: The length of the subkmers to consider.
//
// Returns:
//   - An error if any occurs during the filtering process.
func removeOffTargetSubKmersFromGoodKmers(goodKmers map[string][]int, offTargetKmersFile string, subKmerLength int) error {
	// Convert the good k-mers to canonical uint64 representation
	ori_len := len(goodKmers)

	goodSubKmers := generateSubkmers(goodKmers, subKmerLength)

	goodUint64Kmers, err := convertGoodKmersToUint64Set(goodSubKmers, subKmerLength)
	if err != nil {
		return err
	}

	// Read off-target k-mers and build a map of removed k-mers
	removedKmers, err := removeOffTargetUint64KmersConcurrent(offTargetKmersFile, goodUint64Kmers, 32)
	if err != nil {
		return err
	}

	// Remove the k-mers found in removedKmers from the original goodKmers map
	removeSubKmersFromGoodKmers(goodKmers, removedKmers)
	fmt.Printf("Total off-target-matching kmers removed: %d\n\n", ori_len-len(goodKmers))
	return nil
}

// generateSubkmers generates all subkmers of a given length from a set of kmers.
//
// Parameters:
//   - goodKmers: A map where keys are kmers as strings and values are presence/absence slices.
//   - subkmerLength: The length of the subkmers to generate.
//
// Returns:
//   - A map where keys are the generated subkmers and values are nil.
func generateSubkmers(goodKmers map[string][]int, subkmerLength int) map[string][]int {
	subkmers := make(map[string][]int)

	for kmer := range goodKmers {
		for i := 0; i <= len(kmer)-subkmerLength; i++ {
			subkmer := kmer[i : i+subkmerLength]
			subkmers[subkmer] = nil // Explicitly setting value to nil for each subkmer key
		}
	}

	return subkmers
}

// removeSubKmersFromGoodKmers removes entries from goodKmers whose keys contain any subkmer from subkmers.
//
// Parameters:
//   - goodKmers: A map where keys are target kmers and values are presence/absence slices.
//   - subkmers: A map where keys are subkmers to be removed from goodKmers.
func removeSubKmersFromGoodKmers(goodKmers map[string][]int, subkmers map[string]struct{}) {
	for kmer := range goodKmers {
		for subkmer := range subkmers {
			if strings.Contains(kmer, subkmer) || strings.Contains(kmer, reverseComplement(subkmer)) {
				delete(goodKmers, kmer)
				break
			}
		}
	}
}
