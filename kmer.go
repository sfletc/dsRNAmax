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

// Identify all unique kmers in provided sense strand sequences
// TODO: should the reverse complement be included?
// Output is a map with kmer seq as key and a slice of 1s and 0s indicating which
// input target sequence the kmer is present in.

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

func otRemoval(goodKmers map[string][]int, OTkmerLength int, otRefFiles string) map[string][]int {
	return processOTRemoval(goodKmers, OTkmerLength, otRefFiles, removeOTKmers)
}

func otShortRemoval(goodKmers map[string][]int, OTkmerLength int, otRefFiles string) map[string][]int {
	return processOTRemoval(goodKmers, OTkmerLength, otRefFiles, removeMappedLongOTKmers)
}

func processOTRemoval(goodKmers map[string][]int, OTkmerLength int, otRefFiles string, removalFunc func(map[string][]int, map[string]struct{})) map[string][]int {
	files := strings.Split(otRefFiles, ",")
	for _, otRefFile := range files {
		otRef := RefLoad(otRefFile)
		otKmers := conGetOTKmers(goodKmers, otRef, OTkmerLength)
		removalFunc(goodKmers, otKmers)
	}
	return goodKmers
}

func removeOTKmers(goodKmers map[string][]int, otKmers map[string]struct{}) {
	for key := range otKmers {
		delete(goodKmers, key)
	}
}

func removeMappedLongOTKmers(goodKmers map[string][]int, otKmers map[string]struct{}) {
	for ok := range otKmers {
		for gk := range goodKmers {
			if strings.Contains(gk, ok) {
				delete(goodKmers, gk)
			}
		}
	}
}

// Identifies off-target kmers present in the off-target (either orientation) set
// This is done concurrently to improve speed
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

// Concurrent worker function.  Takes an off-target HeaderRef from the queue
// and adds checks if derived kmers are in the target kmer pool.  If so,
// added to a map, when is then added to the output queue upon completion
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
		fmt.Println("Problem opening fasta reference file", refFile)
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

func ConcurrentlyProcessSequences(refFiles []string, goodKmers map[string][]int, kmerLen int) {
	seqChan := make(chan string, 100)                  // Buffered channel for better performance
	toDeleteChan := make(chan map[string]struct{}, 20) // Channel to collect toDelete maps from workers

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

	// Set up sequence consumers
	numConsumers := 20 // Set the number of workers as needed.
	for i := 0; i < numConsumers; i++ {
		consumerWG.Add(1)
		go KmerCheckSeqs(seqChan, goodKmers, kmerLen, &consumerWG, toDeleteChan)
	}

	// Wait for all consumers to finish processing
	consumerWG.Wait()
	close(toDeleteChan) // Close the toDelete channel once all consumers are done

	// Collect toDelete maps, count and delete the kmers from goodKmers after ensuring all consumers are done
	deletedKmerCount := 0
	for toDelete := range toDeleteChan {
		for kmer := range toDelete {
			if _, exists := goodKmers[kmer]; exists {
				deletedKmerCount++      // Increment the counter if the kmer is found in goodKmers
				delete(goodKmers, kmer) // Delete the kmer from goodKmers
			}
		}
	}

	// Print the count of deleted kmers
	fmt.Printf("Total off-target kmers deleted: %d\n", deletedKmerCount)

}
