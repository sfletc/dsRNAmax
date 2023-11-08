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
	"strings"
	"sync"
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
