package main

import (
	"bufio"
	"bytes"
	"errors"
	"flag"
	"fmt"
	"math"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/adrg/strutil"
	"github.com/adrg/strutil/metrics"
	"github.com/olekukonko/tablewriter"
)

// errorShutdown shuts down the software if there's an error
func errorShutdown() {
	fmt.Println("\nExiting program")
	os.Exit(1)
}

// HeaderRef is a struct comprising a reference sequence header, seques and reverse complement
type HeaderRef struct {
	Header     string
	Seq        string
	ReverseSeq string
}

// RefLoad loads a reference sequence DNA file (FASTA format).
// It returns a slice of HeaderRef structs (individual reference header, sequence and reverse complement).
// Lower case nucleotides are converted to uppercase.
func RefLoad(refFile string) []*HeaderRef {
	var refSlice []*HeaderRef
	var singleHeaderRef *HeaderRef
	var header string
	var refSeq bytes.Buffer
	f, err := os.Open(refFile)
	if err != nil {
		fmt.Println("Problem opening fasta reference file " + refFile)
		errorShutdown()
	}
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		fastaLine := scanner.Text()
		switch {
		case strings.HasPrefix(fastaLine, ">"):
			seq := refSeq.String()
			singleHeaderRef = &HeaderRef{header, seq, reverseComplement(seq)}
			refSlice = append(refSlice, singleHeaderRef)
			header = fastaLine[1:]
			refSeq.Reset()
		case len(fastaLine) != 0:
			refSeq.WriteString(strings.ToUpper(fastaLine))
		}
	}
	seq := refSeq.String()
	singleHeaderRef = &HeaderRef{header, seq, reverseComplement(seq)}
	refSlice = append(refSlice, singleHeaderRef)
	refSlice = refSlice[1:]
	fmt.Println("     --->", len(refSlice), "sequences loaded")
	f.Close()
	return refSlice
}

// Reverse complementary DNA sequence
// Nucleotides must be upper case
// Only complements to ACGT are substituted.  Others remain the same.
func reverseComplement(s string) string {
	var dnaComplement = strings.NewReplacer(
		"A", "T", "T", "A", "G", "C", "C", "G",
	)
	complement := dnaComplement.Replace(s)
	reverseComplement := make([]byte, len(complement))
	for i, j := 0, len(reverseComplement)-1; i < len(reverseComplement); i, j = i+1, j-1 {
		reverseComplement[i] = complement[j]
	}
	return string(reverseComplement)
}

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

func compileOTKmers(headerKmerChan chan map[string]bool) map[string]bool {
	allOTKmers := make(map[string]bool)
	for OTkmers := range headerKmerChan {
		for OTkmer := range OTkmers {
			allOTKmers[OTkmer] = true
		}
	}
	return allOTKmers
}

func removeOTKmers(kmers map[string][]int, otKmers map[string]bool) map[string][]int {
	//Use a pointer so map not copied
	for key := range otKmers {
		delete(kmers, key)
	}
	return kmers
}

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

type construct struct {
	kmerHits []int
	geoMean  float64
	seq      string
}

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

func workerBC(goodKmers map[string][]int, kmerCts map[string]int, kmerLen int, seqLen int, constructLen int, consSeqsChan chan *construct, wg *sync.WaitGroup) {
	var kmerSeq []string
	kmerCtsCpy := make(map[string]int)
	for k, v := range kmerCts {
		kmerCtsCpy[k] = v
		kmerSeq = append(kmerSeq, k)
	}
	rand.Seed(time.Now().UnixNano())
	randomIndex := rand.Intn(len(kmerSeq))
	initKmer := kmerSeq[randomIndex]
	fcons := buildf(kmerCtsCpy, initKmer, kmerLen)
	bcons := buildr(kmerCtsCpy, fcons, kmerLen)
	construct, _ := bestConstruct(goodKmers, bcons, constructLen, kmerLen, seqLen)
	consSeqsChan <- construct
	wg.Done()
}

func compileConsSeqs(consSeqsChan chan *construct) *construct {
	var selConstruct *construct
	best := 0.0
	for eachConstruct := range consSeqsChan {

		if eachConstruct.geoMean > best {
			best = eachConstruct.geoMean
			selConstruct = eachConstruct
		}
	}
	return selConstruct
}

func buildf(kmerCtsCpy map[string]int, initKmer string, kmerLen int) string {
	cons := initKmer
	nucs := []string{"A", "C", "G", "T"}
	for {
		nextSub := cons[len(cons)-kmerLen+1:]
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
			return cons
		} else {
			cons = cons + bestNuc
			delete(kmerCtsCpy, bestKmer)
		}
	}
}

func buildr(kmerCtsCpy map[string]int, fcons string, kmerLen int) string {
	cons := fcons
	nucs := []string{"A", "C", "G", "T"}
	for {
		nextSub := cons[:kmerLen-1]
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
			return cons
		} else {
			cons = bestNuc + cons
			delete(kmerCtsCpy, bestKmer)
		}
	}
}

func bestConstruct(goodKmers map[string][]int, consensus string, constructLen int, kmerLen int, seqLen int) (*construct, error) {

	if len(consensus) < constructLen {
		var bad []int
		return &construct{bad, 0.0, ""}, errors.New("consensus shorter than construct length")
	}
	bestScore := 0.0
	bestPos := 0
	var bestConScores []int
	var allScores [][]int
	for i := 0; i < len(consensus)-kmerLen; i++ { //TODO: check this is correct
		s := goodKmers[consensus[i:i+kmerLen]]
		allScores = append(allScores, s)
	}
	var conScores []int
	for i := 0; i < len(consensus)-constructLen; i++ { //TODO: should this be len(all_scores)
		conScores = nil
		for seq := 0; seq < seqLen; seq++ {
			conScores = append(conScores, 0)
		}
		//A faster method to drop and add end scores could help
		for j := i; j < i+constructLen-kmerLen+1; j++ { //TODO: check this is correct
			for x, y := range allScores[j] {
				conScores[x] += y
			}
		}
		if geoMean(conScores) > bestScore {
			bestScore = geoMean(conScores)
			bestPos = i
			bestConScores = conScores
		}
	}
	return &construct{bestConScores, bestScore, consensus[bestPos : bestPos+constructLen]}, nil
}

func geoMean(input []int) float64 {

	var p float64
	for _, n := range input {
		if n == 0 {
			n = 1
		}
		if p == 0 {
			p = float64(n)
		} else {
			p *= float64(n)
		}
	}
	return math.Pow(p, 1/float64(len(input)))
}

func outputTable(kmerLength *int, selConstruct *construct, ref []*HeaderRef) {
	swg := metrics.NewSmithWatermanGotoh()
	kmerLenStr := strconv.Itoa(*kmerLength)
	table := tablewriter.NewWriter(os.Stdout)
	table.SetHeader([]string{"Target sequence header", "Exact " + kmerLenStr + "nt matches", "Smith Waterman similarity"})
	for i, j := range selConstruct.kmerHits {
		table.Append([]string{ref[i].Header, strconv.Itoa(j), strconv.FormatFloat(strutil.Similarity(selConstruct.seq, ref[i].Seq, swg), 'f', 3, 64)})
	}
	table.SetAlignment(tablewriter.ALIGN_LEFT)
	table.Render()
}

func main() {
	fmt.Println("\ndsRNA construct finder")
	fmt.Println("")
	refFile := flag.String("targets", "", "Path to target FASTA file (required)")
	otRefFile := flag.String("offTargets", "", "Path to off-target FASTA file")
	kmerLength := flag.Int("kmerLen", 21, "Kmer length")
	consLength := flag.Int("constructLen", 300, "dsRNA sense arm length")
	iterations := flag.Int("iterations", 100, "No. of iterations")
	flag.Parse()
	fmt.Println("Loading target sequences")
	ref := RefLoad(*refFile)
	fmt.Println("Getting target sequence kmers")
	goodKmers := getKmers(ref, *kmerLength)
	if *otRefFile != "" {
		fmt.Println("Loading off-target sequences")
		otRef := RefLoad(*otRefFile)
		fmt.Println("Finding and subtracting off-target kmers")
		otKmers := conGetOTKmers(goodKmers, otRef, *kmerLength)
		otRef = nil
		goodKmers = removeOTKmers(goodKmers, otKmers)
	}
	fmt.Println("Counting kmers")
	kmerCts := kmerAbun(goodKmers)
	fmt.Println("Finding best construct")
	selConstruct := conBestConstruct(goodKmers, kmerCts, *kmerLength, len(ref), *consLength, *iterations)
	if selConstruct != nil {
		fmt.Println("\nResults:")
		outputTable(kmerLength, selConstruct, ref)
		fmt.Println("\nGeometric mean of kmer hits to each target sequence:", selConstruct.geoMean)
		fmt.Println("\ndsRNA sense-arm sequence:")
		fmt.Println(selConstruct.seq)
		fmt.Println("")
	} else {
		fmt.Println("Could not identify a dsRNA sense arm sequence.  Check input format and sequence lengths")
	}
}
