package main

import (
	"sync"
)

// Identify all unique kmers in provided sense strand sequences
// TODO: should the reverse complement be included?
// Output is a mpa with kmer seq as key and a slice of 1s and 0s indicating which
// input target sequence the kmer is present in.
func getKmers(ref []*HeaderRef, kmerLen int) map[string][]int {
	refLen := len(ref)
	kmers := make(map[string][]int)
	count := 0
	for _, hr := range ref {
		pos := 0
		ref_seq_len := len(hr.Seq)
		for pos <= ref_seq_len-kmerLen {
			fwd_seq := hr.Seq[pos : pos+kmerLen]
			if _, ok := kmers[fwd_seq]; ok {
				kmers[fwd_seq][count] = 1
			} else {
				kmers[fwd_seq] = make([]int, refLen)
				kmers[fwd_seq][count] = 1
			}
			pos++
		}
		count++
	}
	return kmers
}

// Identifies off-target kmers present in the off-target (either orientation) set
// This is done concurrently to improve speed
func conGetOTKmers(kmers map[string][]int, otRef []*HeaderRef, kmerLen int) map[string]bool {
	wg := &sync.WaitGroup{}
	wg.Add(len(otRef))
	otRefChan := make(chan *HeaderRef, len(otRef))
	for _, header_ref_pair := range otRef {
		otRefChan <- header_ref_pair
	}
	close(otRefChan)
	headerKmerChan := make(chan map[string]bool)

	for a := 0; a < len(otRef); a++ {
		go workerGo(kmers, otRefChan, kmerLen, headerKmerChan, wg)
	}
	go func(cs chan map[string]bool, wg *sync.WaitGroup) {
		wg.Wait()
		close(cs)
	}(headerKmerChan, wg)
	allOTKmers := compileOTKmers(headerKmerChan)
	return allOTKmers
}

func otRemoval(otRefFile *string, goodKmers map[string][]int, kmerLength *int) map[string][]int {

	otRef := RefLoad(*otRefFile)
	otKmers := conGetOTKmers(goodKmers, otRef, *kmerLength)
	otRef = nil
	goodKmers = removeOTKmers(goodKmers, otKmers)
	return goodKmers
}

// Concurrent worker function.  Takes an off-target HeaderRef from the queue
// and adds checks if derived kmers are in the target kmer pool.  If so,
// added to a map, when is then added to the output queue upon completion
func workerGo(kmers map[string][]int, otRefChan chan *HeaderRef, kmerLen int, headerKmerChan chan map[string]bool, wg *sync.WaitGroup) {
	otRef := <-otRefChan
	otKmers := make(map[string]bool)

	pos := 0
	for pos <= len(otRef.Seq)-kmerLen {
		fseq := otRef.Seq[pos : pos+kmerLen]
		if _, ok := kmers[fseq]; ok {
			otKmers[fseq] = true
		}
		rseq := otRef.ReverseSeq[pos : pos+kmerLen]
		if _, ok := kmers[rseq]; ok {
			otKmers[rseq] = true
		}
		pos++
	}
	if len(otKmers) > 0 {
		headerKmerChan <- otKmers
	}
	wg.Done()
}

// Kmer maps in the output queue are compiled into a single off-target kmer map
func compileOTKmers(headerKmerChan chan map[string]bool) map[string]bool {
	allOTKmers := make(map[string]bool)
	for OTkmers := range headerKmerChan {
		for OTkmer := range OTkmers {
			allOTKmers[OTkmer] = true
		}
	}
	return allOTKmers
}

// Off-target kmers removed from the target kmer pool
// TODO: combine to the above function
func removeOTKmers(kmers map[string][]int, otKmers map[string]bool) map[string][]int {
	for key := range otKmers {
		delete(kmers, key)
	}
	return kmers
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
