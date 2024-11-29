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
	"errors"
	"math/rand"
	"sort"
	"sync"
	"time"
)

// construct struct contains kmerHits slice (total present in each input target),
// the geometric mean of the slice, and the sequence of the selected construct
type construct struct {
	kmerHits   []int
	medianHits float64
	seq        string
}

// Concurrent implementation to identify the best construct over multiple iterations
func conBestConstruct(goodKmers map[string][]int, kmerCts map[string]int, kmerLen int, seqLen int, constructLen int, iterations int) *construct {
	wg := &sync.WaitGroup{}
	wg.Add(iterations)
	consSeqsChan := make(chan *construct, iterations)
	for a := 0; a < iterations; a++ {
		go workerBC(goodKmers, kmerCts, kmerLen, seqLen, constructLen, consSeqsChan, wg)
	}
	go func(cs chan *construct, wg *sync.WaitGroup) {
		wg.Wait()
		close(cs)
	}(consSeqsChan, wg)
	selConstruct := compileConsSeqs(consSeqsChan)
	return selConstruct
}

// Worker function for a single iteration of identifying the best construct among the input sequences
// Randomly selects an initKmer, then build forward and backward based on the highest no. of kmer matches
// to each input sequences for a given extension nucleotide until no nucleotides can be added.  kmers are removed
// from the input kmer map upon extension.  The best constrcut with the assembled sequences is selected based on maximising
// the geometric mean of input kmer hits.

func workerBC(goodKmers map[string][]int, kmerCts map[string]int, kmerLen int, seqLen int, constructLen int, consSeqsChan chan *construct, wg *sync.WaitGroup) {
	var kmerSeq []string
	kmerCtsCpy := make(map[string]int)
	for k, v := range kmerCts {
		kmerCtsCpy[k] = v
		kmerSeq = append(kmerSeq, k)
	}
	r := rand.New(rand.NewSource(time.Now().UnixNano()))
	randomIndex := r.Intn(len(kmerSeq))
	initKmer := kmerSeq[randomIndex]
	fcons := buildf(kmerCtsCpy, initKmer, kmerLen)
	bcons := buildr(kmerCtsCpy, fcons, kmerLen)
	construct, _ := bestConstruct(goodKmers, bcons, constructLen, kmerLen, seqLen)
	consSeqsChan <- construct
	wg.Done()
}

// Checks all the generated constrcuts and retains the best (highest geomean)
func compileConsSeqs(consSeqsChan chan *construct) *construct {
	var selConstruct *construct
	best := 0.0
	for eachConstruct := range consSeqsChan {

		if eachConstruct.medianHits > best {
			best = eachConstruct.medianHits
			selConstruct = eachConstruct
		}
	}
	return selConstruct
}

// Build forward from the initial kmer until no nucleotides can be added.  Nucleotide selection
// is random when the kmer abundance for 2 or nucleotides is even.
func buildf(kmerCtsCpy map[string]int, initKmer string, kmerLen int) string {
	consensus := initKmer
	nucs := []string{"A", "C", "G", "T"}
	for {
		nextSub := consensus[len(consensus)-kmerLen+1:]
		bestScore := 0
		bestKmer := ""
		bestNuc := ""
		rand.Seed(time.Now().UnixNano())
		rand.Shuffle(len(nucs), func(i, j int) { nucs[i], nucs[j] = nucs[j], nucs[i] })
		for _, v := range nucs {
			nextKmer := nextSub + v
			if val, ok := kmerCtsCpy[nextKmer]; ok {
				if val > bestScore {
					bestNuc = v
					bestKmer = nextKmer
					bestScore = val
				}
			}
		}
		if bestScore == 0 {
			return consensus
		} else {
			consensus = consensus + bestNuc
			delete(kmerCtsCpy, bestKmer)
		}
	}
}

// Build backward from the completed forward consensus.  Sames rules apply.
func buildr(kmerCtsCpy map[string]int, fcons string, kmerLen int) string {
	consensus := fcons
	nucs := []string{"A", "C", "G", "T"}
	for {
		nextSub := consensus[:kmerLen-1]
		bestScore := 0
		bestKmer := ""
		bestNuc := ""
		rand.Seed(time.Now().UnixNano())
		rand.Shuffle(len(nucs), func(i, j int) { nucs[i], nucs[j] = nucs[j], nucs[i] })
		for _, v := range nucs {
			nextKmer := v + nextSub
			if val, ok := kmerCtsCpy[nextKmer]; ok {
				if val > bestScore {
					bestNuc = v
					bestKmer = nextKmer
					bestScore = val
				}
			}
		}
		if bestScore == 0 {
			return consensus
		} else {
			consensus = bestNuc + consensus
			delete(kmerCtsCpy, bestKmer)
		}
	}
}

// Select the best construct of the specified length from the provided consensus sequence
// by maximising the geometric mean of the number of kmers to match each input target sequence
func bestConstruct(goodKmers map[string][]int, consensus string, constructLen int, kmerLen int, seqLen int) (*construct, error) {
	if len(consensus) < constructLen {
		var bad []int
		return &construct{bad, 0.0, ""}, errors.New("consensus shorter than construct length")
	}
	bestScore := 0.0
	bestPos := 0
	var bestConScores []int
	var allScores [][]int
	for i := 0; i < len(consensus)-kmerLen; i++ {
		s := goodKmers[consensus[i:i+kmerLen]]
		allScores = append(allScores, s)
	}
	for i := 0; i < len(consensus)-constructLen; i++ {
		bestScore, bestPos, bestConScores = bcHelper(seqLen, i, constructLen, kmerLen, allScores, bestScore, bestPos, bestConScores)
	}
	return &construct{bestConScores, bestScore, consensus[bestPos : bestPos+constructLen]}, nil
}

func bcHelper(seqLen int, i int, constructLen int, kmerLen int, allScores [][]int, bestScore float64, bestPos int, bestConScores []int) (float64, int, []int) {
	var conScores []int
	for seq := 0; seq < seqLen; seq++ {
		conScores = append(conScores, 0)
	}

	for j := i; j < i+constructLen-kmerLen+1; j++ {
		for x, y := range allScores[j] {
			conScores[x] += y
		}
	}
	median, err := calculateMedian(conScores)
	if err == nil {
		if median > bestScore {
			bestScore = median
			bestPos = i
			bestConScores = conScores
		}
	}
	return bestScore, bestPos, bestConScores
}

// calculateMedian takes a slice of integers and returns their median as a float64.
// If the slice is empty, it returns an error.
func calculateMedian(numbers []int) (float64, error) {
	// Check if the slice is empty
	if len(numbers) == 0 {
		return 0, errors.New("slice is empty")
	}

	// Sort the slice in ascending order
	sort.Ints(numbers)

	// Calculate the median
	n := len(numbers)
	midIndex := n / 2
	if n%2 == 0 {
		// Even length, take the average of the two middle elements
		return float64(numbers[midIndex-1]+numbers[midIndex]) / 2.0, nil
	} else {
		// Odd length, take the middle element
		return float64(numbers[midIndex]), nil
	}
}
