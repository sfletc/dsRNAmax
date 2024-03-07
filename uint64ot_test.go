package main

import (
	"reflect"
	"testing"
)

func TestRemoveOffTargetUint64KmersConcurrent(t *testing.T) {
	// Create a temporary file for testing

	// Define the goodUint64Kmers map
	goodUint64Kmers, _ := convertGoodKmersToUint64Set(map[string][]int{
		"AAAT": {},
		"TTGC": {},
		"AAAC": {},
	}, 4)
	for kmer := range goodUint64Kmers {
		seq := kmerToSequence(kmer, 4)
		t.Logf("Good k-mer: %s", seq)
	}
	// Call the function with test data
	removedKmers, err := removeOffTargetUint64KmersConcurrent("testData/test.kmer", goodUint64Kmers, 2)
	if err != nil {
		t.Fatalf("Unexpected error: %v", err)
	}
	t.Logf("Removed k-mers: %v", removedKmers)
	// Check the expected removed k-mers
	expectedRemovedKmers := map[string]struct{}{
		"AAAC": {},
		"GCAA": {},
	}
	if len(removedKmers) != len(expectedRemovedKmers) {
		t.Errorf("Unexpected number of removed k-mers. Got: %d, Want: %d", len(removedKmers), len(expectedRemovedKmers))
	}
	for kmer := range expectedRemovedKmers {
		if _, found := removedKmers[kmer]; !found {
			t.Errorf("Expected k-mer not found in removed k-mers: %s", kmer)
		}
	}
}

func TestRemoveSubKmersFromGoodKmers(t *testing.T) {
	// Define the goodKmers map
	goodKmers := map[string][]int{
		"AAAT": {},
		"TTGC": {},
		"AAAC": {},
	}
	// Call the function with test data
	removeOffTargetKmersFromGoodKmers(goodKmers, "testData/test_sub.kmer", 4)
	// Check the expected goodKmers
	expectedGoodKmers := map[string][]int{
		"TTGC": {},
	}
	if len(goodKmers) != len(expectedGoodKmers) {
		t.Errorf("Unexpected number of good k-mers. Got: %d, Want: %d", len(goodKmers), len(expectedGoodKmers))
	}
	for kmer := range expectedGoodKmers {
		if _, found := goodKmers[kmer]; !found {
			t.Errorf("Expected k-mer not found in good k-mers: %s", kmer)
		}
	}
}

func Test_removeKmersFromGoodKmers(t *testing.T) {
	tests := []struct {
		name         string
		goodKmers    map[string][]int
		removedKmers map[string]struct{}
		expected     map[string][]int
	}{
		{
			name: "Remove k-mers present in removedKmers",
			goodKmers: map[string][]int{
				"ATCG": {1, 2, 3},
				"CGAT": {4, 5, 6},
				"TGCA": {7, 8, 9},
			},
			removedKmers: map[string]struct{}{
				"ATCG": {},
				"TGCA": {},
			},
			expected: map[string][]int{},
		},
		{
			name: "No k-mers to remove",
			goodKmers: map[string][]int{
				"ATCG": {1, 2, 3},
				"CGAT": {4, 5, 6},
			},
			removedKmers: map[string]struct{}{
				"TGCA": {},
			},
			expected: map[string][]int{
				"ATCG": {1, 2, 3},
				"CGAT": {4, 5, 6},
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			removeKmersFromGoodKmers(tt.goodKmers, tt.removedKmers)
			if !reflect.DeepEqual(tt.goodKmers, tt.expected) {
				t.Errorf("unexpected result\nexpected: %v\ngot: %v", tt.expected, tt.goodKmers)
			}
		})
	}
}
