// BSD 3-Clause License

// Copyright (c) 2023, Stephen Fletcher
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

package main

import (
	"bufio"   // For reading files with a scanner
	"fmt"     // For printing to the console
	"os"      // For opening and closing files
	"strings" // For string manipulation like HasPrefix
	"sync"    // For using WaitGroup and Mutex
)

// getKmers extracts all unique kmers of a specified length from a set of reference sequences.
//
// Args:
//
//	ref: A slice of HeaderRef structures containing reference sequences.
//	kmerLen: The length of kmers to extract.
//
// Returns:
//
//	A map[string][]int where keys are kmers and values are slices indicating presence (1) or absence (0) of the kmer in each input sequence.
func getKmers(ref []*HeaderRef, kmerLen int) map[string][]int {
	refLen := len(ref)
	kmers := make(map[string][]int)

	for i, hr := range ref {
		ref_seq_len := len(hr.Seq)
		for pos := 0; pos <= ref_seq_len-kmerLen; pos++ {
			fwd_seq := hr.Seq[pos : pos+kmerLen]
			kmerIndices, ok := kmers[fwd_seq]
			if !ok {
				// Preallocate a slice with zeroes.
				kmerIndices = make([]int, refLen)
				kmers[fwd_seq] = kmerIndices
			}
			// Update the presence of k-mer for the current sequence only.
			kmerIndices[i] = 1
		}
	}

	return kmers
}

// removeOTKmers removes off-target kmers from the 'goodKmers' map.
//
// Args:
//
//	goodKmers: A map where keys are target kmers and values are presence/absence slices.
//	otKmers: A map where keys are identified off-target kmers.
func removeOTKmers(goodKmers map[string][]int, otKmers map[string]struct{}) {
	for key := range otKmers {
		delete(goodKmers, key)
	}
}

// removeMappedLongOTKmers removes target kmers from the 'goodKmers' map if they contain any of the provided off-target kmers as substrings.
// This helps filter out potential off-target effects even if there's not an exact match.
//
// Args:
//
//	goodKmers: A map where keys are target kmers and values are presence/absence slices.
//	otKmers: A map where keys are identified off-target kmers.
func removeMappedLongOTKmers(goodKmers map[string][]int, otKmers map[string]struct{}) {
	for ok := range otKmers {
		for gk := range goodKmers {
			if strings.Contains(gk, ok) {
				delete(goodKmers, gk)
			}
		}
	}
}

// conGetOTKmers identifies off-target kmers present in a set of off-target sequences, performing the search concurrently for efficiency.
//
// Args:
//
//	kmers: A map where keys are target kmers and values are presence/absence slices.
//	otRef: A slice of HeaderRef structures containing off-target sequences.
//	kmerLen: The length of kmers to search for.
//
// Returns:
//
//	A map[string]struct{} where keys represent off-target kmers found within the provided off-target sequences.
func conGetOTKmers(kmers map[string][]int, otRef []*HeaderRef, kmerLen int) map[string]struct{} {
	wg := &sync.WaitGroup{}
	wg.Add(len(otRef))
	otRefChan := make(chan *HeaderRef, len(otRef))
	for _, header_ref_pair := range otRef {
		otRefChan <- header_ref_pair
	}
	close(otRefChan)
	headerKmerChan := make(chan map[string]struct{})

	for a := 0; a < len(otRef); a++ {
		go workerGo(kmers, otRefChan, kmerLen, headerKmerChan, wg)
	}
	go func(cs chan map[string]struct{}, wg *sync.WaitGroup) {
		wg.Wait()
		close(cs)
	}(headerKmerChan, wg)
	allOTKmers := compileOTKmers(headerKmerChan)
	return allOTKmers
}

// workerGo is a concurrent worker function that processes a single off-target sequence for off-target kmer identification.
//
// Args:
//
//	kmers: A map where keys are target kmers and values are presence/absence slices.
//	otRefChan: A channel receiving HeaderRef structures (off-target sequences).
//	kmerLen: The length of kmers to search for.
//	headerKmerChan: A channel for sending maps of identified off-target kmers.
//	wg: A WaitGroup for synchronization with the main process.
func workerGo(kmers map[string][]int, otRefChan chan *HeaderRef, kmerLen int, headerKmerChan chan map[string]struct{}, wg *sync.WaitGroup) {
	otRef := <-otRefChan
	otKmers := make(map[string]struct{})

	pos := 0
	for pos <= len(otRef.Seq)-kmerLen {
		fseq := otRef.Seq[pos : pos+kmerLen]
		if _, ok := kmers[fseq]; ok {
			otKmers[fseq] = struct{}{}
		}
		rseq := otRef.ReverseSeq[pos : pos+kmerLen]
		if _, ok := kmers[rseq]; ok {
			otKmers[rseq] = struct{}{}
		}
		pos++
	}
	if len(otKmers) > 0 {
		headerKmerChan <- otKmers
	}
	wg.Done()
}

// Kmer maps in the output queue are compiled into a single off-target kmer map
func compileOTKmers(headerKmerChan chan map[string]struct{}) map[string]struct{} {
	allOTKmers := make(map[string]struct{})
	for OTkmers := range headerKmerChan {
		for OTkmer := range OTkmers {
			allOTKmers[OTkmer] = struct{}{}
		}
	}
	return allOTKmers
}

// Kmer abundance (max = no. of input_target_sequences) calculated for each target kmer
func kmerAbun(kmers map[string][]int) map[string]int {
	kmerCts := make(map[string]int)
	for k, v := range kmers {
		tot := 0
		for _, val := range v {
			tot += val
		}
		kmerCts[k] = tot
	}
	return kmerCts
}

func LoadAndSendSeqs(refFile string, seqChan chan<- string, wg *sync.WaitGroup) {
	defer wg.Done()

	f, err := os.Open(refFile)
	if err != nil {
		fmt.Println("Problem opening FASTA reference file", refFile)
		return // or handle error as needed
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	var refSeq strings.Builder

	for scanner.Scan() {
		fastaLine := scanner.Text()
		if strings.HasPrefix(fastaLine, ">") {
			if refSeq.Len() > 0 {
				seqChan <- refSeq.String()
				refSeq.Reset()
			}
		} else if len(fastaLine) > 0 {
			refSeq.WriteString(strings.ToUpper(fastaLine))
		}
	}

	// Send the last sequence if there is one
	if refSeq.Len() > 0 {
		seqChan <- refSeq.String()
	}
}

// TODO: setup for smaller OT kmers
func KmerCheckSeqs(seqChan <-chan string, goodKmers map[string][]int, kmerLen int, wg *sync.WaitGroup, toDeleteChan chan<- map[string]struct{}) {
	defer wg.Done()

	toDelete := make(map[string]struct{}) // Temporary set to store k-mers to delete

	for seq := range seqChan {
		// Compute the reverse complement of the entire sequence once
		rcSeq := reverseComplement(seq)

		// Iterate over the original sequence and the reverse complement
		for _, s := range []string{seq, rcSeq} {
			for pos := 0; pos <= len(s)-kmerLen; pos++ {
				kmer := s[pos : pos+kmerLen]
				// Check if the k-mer is in goodKmers
				if _, exists := goodKmers[kmer]; exists {
					toDelete[kmer] = struct{}{}
				}
			}
		}
	}

	// Send the toDelete map to the channel
	toDeleteChan <- toDelete
}

// smallKmerCheckSeqs processes a sequence, searching for sub-kmers that may indicate potential off-target matches. It sends kmers to be filtered out through a channel.
//
// Args:
//
//	seqChan: A channel receiving DNA sequences.
//	subKmers: A map where keys are sub-kmers and values are lists of their corresponding longer kmers.
//	subKmerLen: The length of the sub-kmers.
//	wg: A WaitGroup for synchronization with the main process.
//	toDeleteChan: A channel for sending maps of sub-kmers (and their associated longer kmers) to be deleted.
func smallKmerCheckSeqs(seqChan <-chan string, subKmers map[string][]string, subKmerLen int, wg *sync.WaitGroup, toDeleteChan chan map[string][]string) {
	defer wg.Done()

	toDelete := make(map[string][]string) // Temporary set to store k-mers to delete

	for seq := range seqChan {
		// Compute the reverse complement of the entire sequence once
		rcSeq := reverseComplement(seq)

		// Iterate over the original sequence and the reverse complement
		for _, s := range []string{seq, rcSeq} {
			for pos := 0; pos <= len(s)-subKmerLen; pos++ {
				kmer := s[pos : pos+subKmerLen]
				// Check if the k-mer is in goodKmers
				if _, exists := subKmers[kmer]; exists {
					toDelete[kmer] = subKmers[kmer]
				}
			}
		}
	}

	// Send the toDelete map to the channel
	toDeleteChan <- toDelete
}

// GenerateSubKmersMap creates a map where keys are shorter substrings (sub-kmers) and values are lists of the original kmers that contain those sub-kmers.
// This is often used for more flexible off-target matching.
//
// Args:
//
//	goodKmers: A map where keys are target kmers and values are presence/absence slices.
//	subKmerLen: The desired length of the sub-kmers to be generated.
//
// Returns:
//
//	A map[string][]string where keys are sub-kmers and values are lists of the original kmers containing them.
func GenerateSubKmersMap(goodKmers map[string][]int, subKmerLen int) map[string][]string {
	subKmers := make(map[string][]string)

	// Iterate through all k-mers in the input map
	for kmer := range goodKmers {
		if len(kmer) < subKmerLen || subKmerLen == 0 { //TODO: prob should throw an error
			return subKmers
		}
		// Generate all possible substrings of the specified length from the kmer
		for i := 0; i <= len(kmer)-subKmerLen; i++ {
			subKmer := kmer[i : i+subKmerLen]

			// Append the original kmer to the list of kmers that contain this subKmer
			subKmers[subKmer] = append(subKmers[subKmer], kmer)
		}
	}

	return subKmers
}

func ConcurrentlyProcessSequences(refFiles []string, goodKmers map[string][]int, kmerLen int, subKmerLen int) {
	ori_len := len(goodKmers)
	seqChan := make(chan string, 100)                         // Buffered channel for better performance
	toDeleteChan := make(chan map[string]struct{}, 20)        // Channel to collect toDelete maps from workers
	toDeleteSubKmerChan := make(chan map[string][]string, 20) // Channel to collect toDelete maps from workers
	var subKmers map[string][]string
	var producerWG sync.WaitGroup // WaitGroup for producers
	var consumerWG sync.WaitGroup // WaitGroup for consumers

	// Set up sequence producers
	for _, refFile := range refFiles {
		producerWG.Add(1)
		go LoadAndSendSeqs(refFile, seqChan, &producerWG)
	}

	// Close the sequence channel once all producers are done
	go func() {
		producerWG.Wait()
		close(seqChan)
	}()
	if subKmerLen < kmerLen {

		subKmers = GenerateSubKmersMap(goodKmers, subKmerLen)
	}
	// Set up sequence consumers
	numConsumers := 20 // Set the number of workers as needed.
	for i := 0; i < numConsumers; i++ {
		consumerWG.Add(1)
		if subKmerLen == kmerLen {
			go KmerCheckSeqs(seqChan, goodKmers, kmerLen, &consumerWG, toDeleteChan)
		} else {
			go smallKmerCheckSeqs(seqChan, subKmers, subKmerLen, &consumerWG, toDeleteSubKmerChan)
		}
	}

	// Wait for all consumers to finish processing
	consumerWG.Wait()

	// Collect toDelete maps, count and delete the kmers from goodKmers after ensuring all consumers are done
	if subKmerLen == kmerLen {
		close(toDeleteChan) // Close the toDelete channel once all consumers are done
		deletedKmerCount := 0
		for toDelete := range toDeleteChan {
			for kmer := range toDelete {
				if _, exists := goodKmers[kmer]; exists {
					deletedKmerCount++      // Increment the counter if the kmer is found in goodKmers
					delete(goodKmers, kmer) // Delete the kmer from goodKmers
				}
			}
		}
		fmt.Printf("Total off-target-matching kmers removed: %d\n\n", ori_len-len(goodKmers))
	} else {
		close(toDeleteSubKmerChan) // Close the toDelete channel once all consumers are done
		deletedKmerCount := 0
		for toDelete := range toDeleteSubKmerChan {
			for _, longKmers := range toDelete {
				for _, longKmer := range longKmers {
					if _, exists := goodKmers[longKmer]; exists {
						deletedKmerCount++          // Increment the counter if the kmer is found in goodKmers
						delete(goodKmers, longKmer) // Delete the kmer from goodKmers
					}
				}
			}
		}
		fmt.Printf("Total off-target-matching kmers removed: %d\n\n", ori_len-len(goodKmers))

	}

}
