package main

import (
	"encoding/binary"
	"os"
	"reflect"
	"testing"
)

// Test_bestConstruct_lastWindow is a regression test for the off-by-one bug in
// bestConstruct (construct.go). The two loops used `<` instead of `<=`, so the
// final kmer position and the final construct window were never evaluated.
//
// Here the only matching kmers ("ACGT" at pos 7 and "CGTA" at pos 8) sit in the
// very last construct window (start pos 7 = len-constructLen). The buggy code
// silently skipped that window and returned an inferior construct ("AACGT",
// score 1.0); the fixed code finds the true best ("ACGTA", score 2.0).
func Test_bestConstruct_lastWindow(t *testing.T) {
	goodKmers := map[string][]int{
		"ACGT": {1},
		"CGTA": {1},
	}
	const consensus = "AAAAAAAACGTA" // len 12
	got, err := bestConstruct(goodKmers, consensus, 5 /*constructLen*/, 4 /*kmerLen*/, 1 /*seqLen*/)
	if err != nil {
		t.Fatalf("bestConstruct returned error: %v", err)
	}
	want := &construct{kmerHits: []int{2}, medianHits: 2.0, seq: "ACGTA"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("bestConstruct() = %+v, want %+v", got, want)
	}
}

// Test_conBestConstruct_emptyKmers is a regression test for the panic that
// occurred when no candidate kmers were available (e.g. target sequences shorter
// than the kmer length, or off-target filtering removed everything): each worker
// called rand.Intn(0) and panicked. conBestConstruct must instead return nil so
// the caller can report the failure cleanly.
func Test_conBestConstruct_emptyKmers(t *testing.T) {
	got := conBestConstruct(
		map[string][]int{}, // goodKmers
		map[string]int{},   // kmerCts (empty -> would previously panic)
		21,                 // kmerLen
		2,                  // seqLen
		300,                // constructLen
		100,                // iterations
	)
	if got != nil {
		t.Errorf("conBestConstruct() with no kmers = %+v, want nil", got)
	}
}

// Test_GenerateSubKmersMap_skipsShortKmers is a regression test for the
// early-`return` bug in GenerateSubKmersMap (kmer.go). A kmer shorter than the
// requested sub-kmer length used to abort the whole map (returning whatever had
// accumulated so far, which depended on map iteration order); it must instead
// skip that single kmer and process the rest.
func Test_GenerateSubKmersMap_skipsShortKmers(t *testing.T) {
	goodKmers := map[string][]int{
		"ACGTACGT": {1},
		"AC":       {1}, // too short for a 6-mer sub-kmer; must be skipped, not fatal
	}
	got := GenerateSubKmersMap(goodKmers, 6)
	want := map[string][]string{
		"ACGTAC": {"ACGTACGT"},
		"CGTACG": {"ACGTACGT"},
		"GTACGT": {"ACGTACGT"},
	}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("GenerateSubKmersMap() = %v, want %v", got, want)
	}
}

// Test_removeOffTargetKmersFromGoodKmers_subKmerRoute exercises the sub-kmer
// off-target branch (otKmerLen < kmerLen) of removeOffTargetKmersFromGoodKmers
// (uint64ot.go) end-to-end. This is the branch whose error return was being
// discarded; the test confirms the branch runs, removes the correct target kmer
// via an off-target sub-kmer match, and now reports no error.
func Test_removeOffTargetKmersFromGoodKmers_subKmerRoute(t *testing.T) {
	goodKmers := map[string][]int{
		"ACGT": {1}, // contains sub-kmer "ACG" -> should be removed
		"TTTT": {1}, // unrelated -> should remain
	}

	// Binary off-target kmer file: a uint64 header giving the OT kmer length (3),
	// followed by the canonical uint64 encoding of the 3-mer "ACG" (== 6).
	f, err := os.CreateTemp("", "ot.*.bin")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(f.Name())
	if err := binary.Write(f, binary.LittleEndian, uint64(3)); err != nil {
		t.Fatal(err)
	}
	if err := binary.Write(f, binary.LittleEndian, uint64(6)); err != nil {
		t.Fatal(err)
	}
	f.Close()

	if err := removeOffTargetKmersFromGoodKmers(goodKmers, f.Name(), 4); err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	want := map[string][]int{"TTTT": {1}}
	if !reflect.DeepEqual(goodKmers, want) {
		t.Errorf("goodKmers = %v, want %v", goodKmers, want)
	}
}

// Test_processChunk_unalignedChunk guards against the out-of-range panic that
// occurred when a chunk's length was not a multiple of the 8-byte kmer size
// (e.g. an unaligned/short read). The whole kmers must still be processed and
// the trailing partial bytes ignored without panicking.
func Test_processChunk_unalignedChunk(t *testing.T) {
	const kmerSize = 8
	const k = 3

	// 11 bytes with cap == len, matching a real chunk produced by the reader (a
	// buggy loop bound reslices to chunk[8:16] and panics because cap is only 11):
	// one complete kmer (canonical "ACG" == 6) followed by 3 stray bytes.
	chunk := make([]byte, 11)
	binary.LittleEndian.PutUint64(chunk[:kmerSize], 6)
	chunk[8], chunk[9], chunk[10] = 0x01, 0x02, 0x03 // partial kmer -> must be ignored

	good := map[uint64]struct{}{6: {}}
	removed := map[string]struct{}{}
	processChunk(chunk, kmerSize, k, good, removed) // must not panic

	if _, ok := removed["ACG"]; !ok || len(removed) != 1 {
		t.Errorf("processChunk() removed = %v, want exactly {ACG}", removed)
	}
}
