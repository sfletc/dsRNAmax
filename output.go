package main

import (
	"encoding/csv"
	"fmt"
	"os"
	"strconv"

	"github.com/adrg/strutil"
	"github.com/adrg/strutil/metrics"
	"github.com/olekukonko/tablewriter"
)

// Output results to commandline and a CSV file for each input sequence and the dsRNA sense arm itself
func outputResults(goodKmers map[string][]int, kmerLength *int, selConstruct *construct, ref []*HeaderRef, csvFileName string) {
	fmt.Println("\nResults:")
	modKmerHits, rowData := outputTable(goodKmers, kmerLength, selConstruct, ref) // outputTable will now also return rowData for CSV

	if csvFileName != "" {
		err := writeToCSV(csvFileName, rowData, selConstruct.seq)
		if err != nil {
			fmt.Println("Error writing to CSV:", err)
			return
		}
		fmt.Println("Results written to", csvFileName)
	}
	// Calculate and output median
	median, err := calculateMedian(modKmerHits)
	if err != nil {
		fmt.Println("Error calculating median:", err)
	} else {
		fmt.Println("\nMedian of kmer hits to each target sequence:", median)
	}

	// Other output information
	fmt.Println("\ndsRNA sense-arm sequence - " + strconv.FormatFloat(gcContent(selConstruct.seq), 'f', 1, 64) + "% GC content")
	fmt.Println(selConstruct.seq)
	fmt.Println("")
}

func writeToCSV(filename string, data [][]string, seq string) error {
	// Create a new CSV file
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	// Initialize the CSV writer
	writer := csv.NewWriter(file)
	defer writer.Flush()

	// Write all data to the CSV file
	for _, row := range data {
		if err := writer.Write(row); err != nil {
			return err // returns early if there is an error
		}
	}

	// Write an empty row to create a space after the table
	if err := writer.Write([]string{}); err != nil {
		return err // returns early if there is an error
	}

	// Write the "seq" label in a cell by itself
	if err := writer.Write([]string{"dsRNA sense-arm sequence:"}); err != nil {
		return err // returns early if there is an error
	}

	// Write the sequence on the next line
	if err := writer.Write([]string{seq}); err != nil {
		return err // returns early if there is an error
	}

	return nil // returns nil if everything was written successfully
}

// Generate table and prepare data for CSV
func outputTable(goodKmers map[string][]int, kmerLength *int, selConstruct *construct, ref []*HeaderRef) ([]int, [][]string) {
	kmerLenStr := strconv.Itoa(*kmerLength)
	kmers := kmersPerInput(goodKmers, selConstruct.seq, *kmerLength, len(ref))
	meanGC := meanGCforKmers(kmers)
	kmerStats := kmerStats(kmers)
	table := tablewriter.NewWriter(os.Stdout)
	headers := []string{"Target sequence header",
		kmerLenStr + "nt matches",
		"SWG similarity (%)",
		"Kmer mean GC (%)",
		"5'U (%)",
		"5'A (%)",
		"5'C (%)"}
	table.SetHeader(headers)
	modKmerHits, csvData := generateRowData(kmerStats, meanGC, selConstruct, ref, table, *kmerLength)
	table.SetAlignment(tablewriter.ALIGN_LEFT)
	table.Render()
	return modKmerHits, csvData
}

// Returns the kmers that match each input sequence
// Use only after the best construct has been returned
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
		result := []float64{0, 0, 0, 0}
		for _, kmer := range kmers {
			switch kmer[len(kmer)-1:] {
			case "A":
				result[0] += 1.0 / float64(len(kmers))
			case "T":
				result[1] += 1.0 / float64(len(kmers))
			case "G":
				result[2] += 1.0 / float64(len(kmers))
			}
			result[3] += 1.0
		}
		results = append(results, result)
	}

	return results
}

// generateRowData prepares data for the output table and CSV export
func generateRowData(kmerStats [][]float64, meanGC []float64, selConstruct *construct, ref []*HeaderRef, table *tablewriter.Table, kmerLength int) ([]int, [][]string) {
	headerMap := make(map[string]bool)
	var modKmerHits []int
	var csvData [][]string // Initialize slice to hold CSV data rows

	// Initialize the Smith-Waterman-Gotoh algorithm parameters
	swg := metrics.NewSmithWatermanGotoh()
	swg.GapPenalty = -2
	swg.Substitution = metrics.MatchMismatch{
		Match:    1,
		Mismatch: -2,
	}

	// Prepare headers for the CSV output
	csvHeaders := []string{
		"Target sequence header",
		strconv.Itoa(kmerLength) + "nt matches",
		"SWG similarity (%)",
		"Kmer mean GC (%)",
		"5'U (%)",
		"5'A (%)",
		"5'C (%)",
	}
	csvData = append(csvData, csvHeaders) // Append headers to the CSV data slice

	// Iterate over each target sequence to prepare data for output and CSV
	for i, _ := range selConstruct.kmerHits {
		if _, ok := headerMap[ref[i].Header]; !ok {
			// Prepare row data for the terminal table output
			modKmerHits = append(modKmerHits, int(kmerStats[i][3])) // Collect modified kmer hit counts
			tableRow := []string{
				ref[i].Header,
				strconv.FormatFloat(kmerStats[i][3], 'f', 0, 64),
				strconv.FormatFloat(strutil.Similarity(selConstruct.seq, ref[i].Seq, swg)*100, 'f', 1, 64),
				strconv.FormatFloat(meanGC[i], 'f', 1, 64),
				strconv.FormatFloat(kmerStats[i][0]*100, 'f', 1, 64),
				strconv.FormatFloat(kmerStats[i][1]*100, 'f', 1, 64),
				strconv.FormatFloat(kmerStats[i][2]*100, 'f', 1, 64),
			}
			table.Append(tableRow) // Append row data to the terminal table

			// Prepare row data for CSV output
			csvRow := []string{
				ref[i].Header,
				strconv.FormatFloat(kmerStats[i][3], 'f', 0, 64),
				strconv.FormatFloat(strutil.Similarity(selConstruct.seq, ref[i].Seq, swg)*100, 'f', 1, 64),
				strconv.FormatFloat(meanGC[i], 'f', 1, 64),
				strconv.FormatFloat(kmerStats[i][0]*100, 'f', 1, 64),
				strconv.FormatFloat(kmerStats[i][1]*100, 'f', 1, 64),
				strconv.FormatFloat(kmerStats[i][2]*100, 'f', 1, 64),
			}
			csvData = append(csvData, csvRow) // Append row data to the CSV data slice

			headerMap[ref[i].Header] = true // Mark header as processed
		}
	}

	return modKmerHits, csvData // Return the modified kmer hits and the CSV data for export
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
