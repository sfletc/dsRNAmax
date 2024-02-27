package main

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"sync"
)

func removeOffTargetUint64KmersConcurrent(filename string, goodUint64Kmers map[uint64]struct{}, numWorkers int) (map[string]struct{}, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening file: %v", err)
	}
	defer file.Close()
	//read the first line of the binary file to get the kmer length
	var k64 uint64

	if err := binary.Read(file, binary.LittleEndian, &k64); err != nil {
		return nil, fmt.Errorf("error reading kmer length: %v", err)
	}

	var k int = int(k64)
	const chunkSize = 1024 * 64 // 64 KB; adjust as needed
	kmerSize := binary.Size(uint64(0))
	kmerChan := make(chan []byte, numWorkers)
	removedKmerChan := make(chan map[string]struct{})
	var wg sync.WaitGroup

	// Start worker goroutines to process chunks of k-mers
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			localRemovedKmers := make(map[string]struct{})
			for chunk := range kmerChan {
				for j := 0; j < len(chunk); j += kmerSize {
					var kmer uint64
					buf := bytes.NewReader(chunk[j : j+kmerSize])
					if err := binary.Read(buf, binary.LittleEndian, &kmer); err != nil {
						// Handle error (log or break)
						continue
					}
					if _, found := goodUint64Kmers[kmer]; found {
						seq := kmerToSequence(kmer, k)
						localRemovedKmers[seq] = struct{}{}
					}
				}
			}
			removedKmerChan <- localRemovedKmers
		}()
	}

	// Read large chunks of data and send to worker goroutines
	go func() {
		defer close(kmerChan)
		buf := make([]byte, chunkSize)
		for {
			bytesRead, err := file.Read(buf)
			if err != nil {
				if err == io.EOF {
					break
				}
				// Handle other errors
				break
			}
			if bytesRead > 0 {
				chunk := make([]byte, bytesRead)
				copy(chunk, buf[:bytesRead])
				kmerChan <- chunk
			}
		}
	}()

	// Wait for all workers to finish and close the removedKmerChan
	go func() {
		wg.Wait()
		close(removedKmerChan)
	}()

	// Merge results from workers
	mergedRemovedKmers := make(map[string]struct{})
	for localMap := range removedKmerChan {
		for kmer := range localMap {
			mergedRemovedKmers[kmer] = struct{}{}
		}
	}

	return mergedRemovedKmers, nil
}

func kmerToSequence(kmer uint64, k int) string {
	bases := "ACGT"
	sequence := make([]byte, k)

	for i := k - 1; i >= 0; i-- {
		sequence[i] = bases[kmer&3]
		kmer >>= 2
	}
	return string(sequence)
}

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

// convertGoodKmersToUint64Set converts a set of k-mers from string format to canonical uint64 format.
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

// min returns the smaller of two uint64 values.
func min(a, b uint64) uint64 {
	if a < b {
		return a
	}
	return b
}

func removeOffTargetKmersFromGoodKmers(goodKmers map[string][]int, offTargetKmersFile string, goodKmerLength int) error {
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
	return nil
}
