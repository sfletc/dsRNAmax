package main

import (
	"reflect"
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
			want: map[string]struct{}{"ACGT": struct{}{}, "CGTA": struct{}{}},
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
				otKmers: map[string]struct{}{"ACGT": struct{}{}, "CGTA": struct{}{}, "GTC": struct{}{}},
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

func Test_geoMean(t *testing.T) {
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
			name: "geoMeanSuccess",
			args: args{
				input: []int{1, 2, 4},
			},
			want:    2,
			wantErr: false,
		},
		{
			name: "geoMeanError",
			args: args{
				input: []int{0, 2, 4},
			},
			want:    0,
			wantErr: true,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := geoMean(tt.args.input)
			if (err != nil) != tt.wantErr {
				t.Errorf("geoMean() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if got != tt.want {
				t.Errorf("geoMean() = %v, want %v", got, tt.want)
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
