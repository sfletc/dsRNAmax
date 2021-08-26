package main

import (
	"fmt"
	"os"
	"strconv"

	"github.com/adrg/strutil"
	"github.com/adrg/strutil/metrics"
	"github.com/olekukonko/tablewriter"
)

func outputResults(goodKmers map[string][]int, kmerLength *int, selConstruct *construct, ref []*HeaderRef) {
	fmt.Println("\nResults:")
	modKmerHits := outputTable(goodKmers, kmerLength, selConstruct, ref)
	fmt.Println("\nGeometric mean of kmer hits to each target sequence:", geoMean(modKmerHits))
	fmt.Println("\ndsRNA sense-arm sequence:")
	fmt.Println(selConstruct.seq)
	fmt.Println("")
}

func outputTable(goodKmers map[string][]int, kmerLength *int, selConstruct *construct, ref []*HeaderRef) []int {
	kmerLenStr := strconv.Itoa(*kmerLength)
	table := tablewriter.NewWriter(os.Stdout)
	table.SetHeader([]string{"Target sequence header", "Exact " + kmerLenStr + "nt matches", "Smith Waterman Gotoh similarity (%)", "mean GC percentage for matching kmers"})
	modKmerHits := accountForBias(goodKmers, selConstruct, ref, table, *kmerLength)
	table.SetAlignment(tablewriter.ALIGN_LEFT)
	table.Render()
	return modKmerHits
}

//Return GC content of a kmer expressed as a percentage
func gcContent(kmer string) float64 {
	gc := 0.0
	for _, nuc := range kmer {
		if nuc == 'G' || nuc == 'C' {
			gc += 1.0
		}
	}
	return gc * 100 / float64(len(kmer))
}

//Returns the kmers that match each input sequence
//Use only after the best constrcut has been returned
func kmersPerInput(goodKmers map[string][]int, bestConstruct string, kmerLen int, headerNo int) [][]string {
	kmersForInput := make([][]string, headerNo)
	for i := 0; i < len(bestConstruct)-kmerLen+1; i++ {
		kmer := bestConstruct[i : i+kmerLen]
		for index, count := range goodKmers[kmer] {
			if count == 1 {
				kmersForInput[index] = append(kmersForInput[index], kmer)
			}
		}
	}
	return kmersForInput
}

func meanGCforKmers(kmersForInput [][]string) []float64 {
	var meanGC []float64
	gc := 0.0
	for _, kmers := range kmersForInput {
		gc = 0.0
		for _, kmer := range kmers {
			gc += gcContent(kmer)
		}
		meanGC = append(meanGC, gc/float64(len(kmers)))
	}
	return meanGC
}

func accountForBias(goodKmers map[string][]int, selConstruct *construct, ref []*HeaderRef, table *tablewriter.Table, kmerLength int) []int {
	headerMap := make(map[string]bool)
	kmers := kmersPerInput(goodKmers, selConstruct.seq, kmerLength, len(ref))
	meanGC := meanGCforKmers(kmers)
	var modKmerHits []int
	swg := metrics.NewSmithWatermanGotoh()
	swg.GapPenalty = -2
	swg.Substitution = metrics.MatchMismatch{
		Match:    1,
		Mismatch: -2,
	}
	for i, j := range selConstruct.kmerHits {
		if _, ok := headerMap[ref[i].Header]; !ok {
			table.Append([]string{ref[i].Header, strconv.Itoa(j),
				strconv.FormatFloat(strutil.Similarity(selConstruct.seq, ref[i].Seq, swg)*100, 'f', 1, 64),
				strconv.FormatFloat(meanGC[i], 'f', 1, 64)})
			modKmerHits = append(modKmerHits, j)
			headerMap[ref[i].Header] = true
		}
	}
	return modKmerHits
}
