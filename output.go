package main

import (
	"fmt"
	"os"
	"strconv"

	"github.com/adrg/strutil"
	"github.com/adrg/strutil/metrics"
	"github.com/olekukonko/tablewriter"
)

// Output results to commandline for each input sequence and the dsRNA sense arm itself
func outputResults(goodKmers map[string][]int, kmerLength *int, selConstruct *construct, ref []*HeaderRef) {
	fmt.Println("\nResults:")
	modKmerHits := outputTable(goodKmers, kmerLength, selConstruct, ref) // Assuming outputTable returns a []int

	// Calculate and output the geometric mean
	mean, err := geoMean(modKmerHits)
	if err != nil {
		fmt.Println("Error calculating geometric mean:", err)
	} else {
		fmt.Println("\nGeometric mean of kmer hits to each target sequence:", mean)
	}

	// Other output information
	fmt.Println("\ndsRNA sense-arm sequence - " + strconv.FormatFloat(gcContent(selConstruct.seq), 'f', 1, 64) + "% GC content")
	fmt.Println(selConstruct.seq)
	fmt.Println("")
}

// Generate table
func outputTable(goodKmers map[string][]int, kmerLength *int, selConstruct *construct, ref []*HeaderRef) []int {
	kmerLenStr := strconv.Itoa(*kmerLength)
	kmers := kmersPerInput(goodKmers, selConstruct.seq, *kmerLength, len(ref))
	meanGC := meanGCforKmers(kmers)
	kmerStats := kmerStats(kmers)
	table := tablewriter.NewWriter(os.Stdout)
	table.SetHeader([]string{"Target sequence header",
		kmerLenStr + "nt matches",
		"SWG similarity (%)",
		"Kmer mean GC (%)",
		"5'U (%)",
		"5'A (%)",
		"5'C (%)"})
	modKmerHits := generateRowData(kmerStats, meanGC, selConstruct, ref, table, *kmerLength)
	table.SetAlignment(tablewriter.ALIGN_LEFT)
	table.Render()
	return modKmerHits
}

// Returns the kmers that match each input sequence
// Use only after the best constrcut has been returned
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

func kmerStats(kmersPerInput [][]string) [][]float64 {
	var results [][]float64
	// result := []float64{0, 0, 0}
	for _, kmers := range kmersPerInput {
		result := []float64{0, 0, 0}
		for _, kmer := range kmers {
			switch kmer[len(kmer)-1:] {
			case "A":
				result[0] += 1.0 / float64(len(kmers))
			case "T":
				result[1] += 1.0 / float64(len(kmers))
			case "G":
				result[2] += 1.0 / float64(len(kmers))
			}
		}
		results = append(results, result)
	}
	// fmt.Println(results)
	return results
}

func generateRowData(kmerStats [][]float64, meanGC []float64, selConstruct *construct, ref []*HeaderRef, table *tablewriter.Table, kmerLength int) []int {
	headerMap := make(map[string]bool)
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
				strconv.FormatFloat(meanGC[i], 'f', 1, 64),
				strconv.FormatFloat(kmerStats[i][0]*100, 'f', 1, 64),
				strconv.FormatFloat(kmerStats[i][1]*100, 'f', 1, 64),
				strconv.FormatFloat(kmerStats[i][2]*100, 'f', 1, 64),
			})
			modKmerHits = append(modKmerHits, j)
			headerMap[ref[i].Header] = true
		}
	}
	return modKmerHits
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

// Return GC content of a kmer expressed as a percentage
func gcContent(kmer string) float64 {
	gc := 0.0
	for _, nuc := range kmer {
		if nuc == 'G' || nuc == 'C' {
			gc += 1.0
		}
	}
	return gc * 100 / float64(len(kmer))
}
