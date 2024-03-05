package main

import (
	"fmt"
	"os"
	"reflect"
	"strings"
	"sync"
	"testing"
)

func TestRefLoad(t *testing.T) {
	var should_be []*HeaderRef
	ref1 := &HeaderRef{"ref_1", "AAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTT"}
	ref2 := &HeaderRef{"ref_2", "GGGGGGGGGGGGGGGGGGGGGGGGTAAAAAAAAAAAAAAAAAAAAAAAAG", "CTTTTTTTTTTTTTTTTTTTTTTTTACCCCCCCCCCCCCCCCCCCCCCCC"}
	ref3 := &HeaderRef{"ref_3", "", ""}
	should_be = append(should_be, ref1, ref2, ref3)
	type args struct {
		refFile string
	}
	tests := []struct {
		name string
		args args
		want []*HeaderRef
	}{
		{
			name: "refLoadSuccess",
			args: args{
				refFile: "./testData/testRef.fa",
			},
			want: should_be,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := RefLoad(tt.args.refFile); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("RefLoad() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_reverseComplement(t *testing.T) {
	type args struct {
		seq string
	}
	tests := []struct {
		name string
		args args
		want string
	}{
		{
			name: "reverseComplementSuccess",
			args: args{
				seq: "AACT",
			},
			want: "AGTT",
		},
		{
			name: "badNuc",
			args: args{
				seq: "AACH",
			},
			want: "HGTT",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := reverseComplement(tt.args.seq); got != tt.want {
				t.Errorf("reverseComplement() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_getKmers(t *testing.T) {
	type args struct {
		ref []*HeaderRef
		nt  int
	}
	tests := []struct {
		name string
		args args
		want map[string][]int
	}{
		{
			name: "getKmersSuccess",
			args: args{
				ref: []*HeaderRef{{"test", "ACGTA", "TACGT"}, {"test2", "ACGT", "ACGT"}},
				nt:  4,
			},
			want: map[string][]int{"ACGT": {1, 1}, "CGTA": {1, 0}},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := getKmers(tt.args.ref, tt.args.nt); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("getKmers() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_conGetOTKmers(t *testing.T) {
	type args struct {
		kmers   map[string][]int
		otRef   []*HeaderRef
		kmerLen int
	}
	tests := []struct {
		name string
		args args
		want map[string]struct{}
	}{
		{
			name: "conGetOTKmersSuccess",
			args: args{
				kmers:   map[string][]int{"ACGT": {1, 1}, "CGTA": {1, 0}},
				otRef:   []*HeaderRef{{"test", "ACGTA", "TACGT"}, {"test2", "ACGT", "ACGT"}},
				kmerLen: 4,
			},
			want: map[string]struct{}{"ACGT": {}, "CGTA": {}},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := conGetOTKmers(tt.args.kmers, tt.args.otRef, tt.args.kmerLen); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("conGetOTKmers() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_removeOTKmers(t *testing.T) {
	type args struct {
		kmers   map[string][]int
		otKmers map[string]struct{}
	}
	tests := []struct {
		name string
		args args
		want map[string][]int
	}{
		{
			name: "removeKmerSuccess",
			args: args{
				kmers:   map[string][]int{"ACGT": {1, 1}, "CGTA": {1, 0}, "GTC": {0, 1}, "AATC": {1, 1}},
				otKmers: map[string]struct{}{"ACGT": {}, "CGTA": {}, "GTC": {}},
			},
			want: map[string][]int{"AATC": {1, 1}},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			removeOTKmers(tt.args.kmers, tt.args.otKmers)
			if !reflect.DeepEqual(tt.args.kmers, tt.want) {
				t.Errorf("removeOTKmers() got = %v, want %v", tt.args.kmers, tt.want)
			}
		})
	}
}

func Test_kmerAbun(t *testing.T) {
	type args struct {
		kmers map[string][]int
	}
	tests := []struct {
		name string
		args args
		want map[string]int
	}{
		{
			name: "kmerAbunSuccess",
			args: args{
				kmers: map[string][]int{"ACGT": {1, 1}, "CGTA": {1, 0}, "AATC": {1, 1}},
			},
			want: map[string]int{"ACGT": 2, "CGTA": 1, "AATC": 2},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := kmerAbun(tt.args.kmers); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("kmerAbun() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_conBestConstruct(t *testing.T) {
	type args struct {
		goodKmers    map[string][]int
		kmerCts      map[string]int
		kmerLen      int
		seqLen       int
		constructLen int
		iterations   int
	}
	tests := []struct {
		name string
		args args
		want *construct
	}{
		{
			name: "conBestConstructSuccess",
			args: args{
				goodKmers:    map[string][]int{"ACGT": {1, 1}, "CGTA": {1, 1}, "GTAC": {1, 0}},
				kmerCts:      map[string]int{"ACGT": 2, "CGTA": 2, "GTAC": 1},
				kmerLen:      4,
				seqLen:       2,
				constructLen: 5,
				iterations:   2,
			},
			want: &construct{[]int{2, 2}, 2.0, "ACGTA"},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := conBestConstruct(tt.args.goodKmers, tt.args.kmerCts, tt.args.kmerLen, tt.args.seqLen, tt.args.constructLen, tt.args.iterations); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("conBestConstruct() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_mean(t *testing.T) {
	type args struct {
		input []int
	}
	tests := []struct {
		name    string
		args    args
		want    float64
		wantErr bool
	}{
		{
			name: "medianSuccess",
			args: args{
				input: []int{1, 2, 3},
			},
			want:    2.0,
			wantErr: false,
		},
		{
			name: "medianError",
			args: args{
				input: []int{},
			},
			want:    0,
			wantErr: true,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := calculateMedian(tt.args.input)
			if (err != nil) != tt.wantErr {
				t.Errorf("median() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if got != tt.want {
				t.Errorf("median() = %v, want %v", got, tt.want)
			}
		})
	}
}

// TestRemoveMappedLongOTKmers verifies that the removeMappedLongOTKmers function
// removes keys from goodKmers that contain any key from otKmers as a substring.
func TestRemoveMappedLongOTKmers(t *testing.T) {
	// Define a test case structure
	tests := []struct {
		name      string
		goodKmers map[string][]int
		otKmers   map[string]struct{}
		expected  map[string][]int
	}{
		{
			name: "remove substrings",
			goodKmers: map[string][]int{
				"ACTGG": {1},
				"ACTGA": {1},
				"GGCTC": {1},
				"TGACC": {1},
			},
			otKmers: map[string]struct{}{
				"ACT": {}, // "ACTG" and "ACTGA" should be removed
				"TG":  {}, // "TGAC" should be removed
			},
			expected: map[string][]int{
				"GGCTC": {1}, // Only "GACT" should remain
			},
		},
		{
			name: "no removal when no substrings match",
			goodKmers: map[string][]int{
				"AAAA": {1},
				"CCCC": {1},
			},
			otKmers: map[string]struct{}{
				"GGGG": {},
			},
			expected: map[string][]int{
				"AAAA": {1},
				"CCCC": {1},
			},
		},
		// Add more test cases as needed
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			removeMappedLongOTKmers(tt.goodKmers, tt.otKmers)
			if !reflect.DeepEqual(tt.goodKmers, tt.expected) {
				t.Errorf("removeMappedLongOTKmers() got = %v, want %v", tt.goodKmers, tt.expected)
			}
		})
	}
}

// TestLoadAndSendSeqs tests the LoadAndSendSeqs function to ensure it sends the correct sequences.
func TestLoadAndSendSeqs(t *testing.T) {
	// Create a temporary file to simulate the fasta file input
	tmpFile, err := os.CreateTemp("", "example.*.fasta")
	if err != nil {
		t.Fatalf("Failed to create temp file: %v", err)
	}
	defer os.Remove(tmpFile.Name()) // clean up after the test

	// Write example fasta content to the temp file
	_, err = tmpFile.WriteString(">seq1\nATGC\n>seq2\nCGTA")
	if err != nil {
		t.Fatalf("Failed to write to temp file: %v", err)
	}
	tmpFile.Sync() // Ensure the file is written and available for reading
	tmpFile.Close()

	// Set up the channel and waitgroup
	seqChan := make(chan string) // unbuffered channel
	var wg sync.WaitGroup
	wg.Add(1)

	// Invoke LoadAndSendSeqs in a goroutine
	go LoadAndSendSeqs(tmpFile.Name(), seqChan, &wg)

	// Prepare a set to track received sequences
	receivedSeqs := make(map[string]struct{})

	// Start another goroutine to close the channel once all sequences are processed
	go func() {
		wg.Wait()
		close(seqChan)
	}()

	// Collect sequences from the channel
	for seq := range seqChan {
		receivedSeqs[seq] = struct{}{}
	}

	// Define the expected sequences as a set
	expectedSeqs := map[string]struct{}{
		"ATGC": {},
		"CGTA": {},
	}

	// Check that each expected sequence is in the received set
	for seq := range expectedSeqs {
		if _, ok := receivedSeqs[seq]; !ok {
			t.Errorf("Expected sequence %s not found", seq)
		}
	}

	// Check for any unexpected sequences
	for seq := range receivedSeqs {
		if _, ok := expectedSeqs[seq]; !ok {
			t.Errorf("Unexpected sequence %s received", seq)
		}
	}
}

// mock data for testing

// Mock LoadAndSendSeqs to preload sequences instead of reading from a file
func MockLoadAndSendSeqs(preloadedSequences []string, seqChan chan<- string, wg *sync.WaitGroup) {
	defer wg.Done()
	for _, seq := range preloadedSequences {
		seqChan <- seq
	}
}

// Create a temporary file with the given sequences and return its path
func createTempFastaFile(sequences []string, t *testing.T) string {
	content := strings.Join(sequences, "\n")
	tmpfile, err := os.CreateTemp("", "example.*.fasta")
	if err != nil {
		t.Fatalf("Failed to create temp file: %s", err)
	}
	if _, err := tmpfile.WriteString(content); err != nil {
		t.Fatalf("Failed to write to temp file: %s", err)
	}
	if err := tmpfile.Close(); err != nil {
		t.Fatalf("Failed to close temp file: %s", err)
	}
	return tmpfile.Name()
}

func TestKmerCheckSeqs(t *testing.T) {
	seqChan := make(chan string, 10)                   // Buffered for sending test sequences without blocking.
	toDeleteChan := make(chan map[string]struct{}, 10) // Buffered to receive toDelete maps without blocking.
	wg := &sync.WaitGroup{}

	// Define the goodKmers map with example data.
	goodKmers := map[string][]int{
		"ATGC": {1},
		"GCAT": {1}, // reverse complement of ATGC
		"CGTA": {1},
		"TTTT": {1},
	}

	// Define the input sequences containing the kmers.
	inputSequences := []string{
		"ATGCATAAA", // contains ATGC (which has a reverse complement in the map)
		"TAGCATT",   // does not contain any kmers from the map
		"ACGTAC",    // contains CGTA (which is in the map)
		"GCATGC",    // contains GCAT (which is the reverse complement of ATGC in the map)
	}

	// Start the KmerCheckSeqs in a separate goroutine.
	wg.Add(1)
	go KmerCheckSeqs(seqChan, goodKmers, 4, wg, toDeleteChan)

	// Send the sequences to the channel and close it.
	for _, seq := range inputSequences {
		seqChan <- seq
	}
	close(seqChan)

	// Wait for KmerCheckSeqs to finish processing.
	wg.Wait()
	close(toDeleteChan) // Close the toDeleteChan channel after all workers are done.

	// Process the toDelete maps and update the goodKmers map accordingly.
	for toDelete := range toDeleteChan {
		for kmer := range toDelete {
			delete(goodKmers, kmer)
		}
	}

	// Define the expected results.
	expectedGoodKmers := map[string][]int{
		"TTTT": {1}, // Assuming this is a kmer that should remain.
	}
	fmt.Println(goodKmers)
	// Check if the goodKmers map has been updated correctly.
	if len(goodKmers) != len(expectedGoodKmers) {
		t.Errorf("KmerCheckSeqs failed: expected %d good kmers, got %d", len(expectedGoodKmers), len(goodKmers))
	}

	for kmer, count := range expectedGoodKmers {
		if actualCount, exists := goodKmers[kmer]; !exists || actualCount[0] != count[0] {
			t.Errorf("KmerCheckSeqs failed: kmer '%s' was not processed correctly", kmer)
		}
	}
}

// TestConcurrentlyProcessSequences tests the ConcurrentlyProcessSequences function
func TestConcurrentlyProcessSequences(t *testing.T) {
	// Arrange
	refSequences := []string{
		">seq1",
		"ATCGATCGATCG",
		">seq2",
		"GCTAGCTAGCTA",
	}
	refFile := createTempFastaFile(refSequences, t)
	defer os.Remove(refFile) // clean up

	goodKmers := map[string][]int{
		"ATCG": {1},
		"GCTA": {1},
	}

	kmerLen := 4
	otKmerLen := 4

	// Act
	ConcurrentlyProcessSequences([]string{refFile}, goodKmers, kmerLen, otKmerLen)

	// Assert
	expectedRemainingKmers := map[string][]int{
		// Based on the logic, these k-mers should have been deleted
		// if they were present in the input sequences.
	}

	if len(goodKmers) != len(expectedRemainingKmers) {
		t.Errorf("got %d goodKmers; want %d", len(goodKmers), len(expectedRemainingKmers))
	}

	for kmer := range expectedRemainingKmers {
		if _, exists := goodKmers[kmer]; !exists {
			t.Errorf("expected kmer %s to be present; it was not", kmer)
		}
	}
}

// TestGenerateSubKmersMap provides test cases for the GenerateSubKmersMap function.
func TestGenerateSubKmersMap(t *testing.T) {
	tests := []struct {
		name       string
		goodKmers  map[string][]int // Assuming the type should be map[string]struct{} as per original function
		subKmerLen int
		want       map[string][]string
	}{
		{
			name: "single k-mer",
			goodKmers: map[string][]int{
				"ACGTACGTACGT": {1}}, // Corrected here
			subKmerLen: 10,
			want: map[string][]string{
				"ACGTACGTAC": {"ACGTACGTACGT"},
				"CGTACGTACG": {"ACGTACGTACGT"},
				"GTACGTACGT": {"ACGTACGTACGT"},
			},
		},
		{
			name: "subKmerLen longer than k-mer",
			goodKmers: map[string][]int{
				"ACGT": {1}}, // Assuming this is the intended type
			subKmerLen: 5,
			want:       map[string][]string{},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := GenerateSubKmersMap(tt.goodKmers, tt.subKmerLen); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("GenerateSubKmersMap() = %v, want %v", got, tt.want)
			}
		})
	}
}
