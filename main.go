// Copyright © 2020 Genome Research Limited
// Author: Sendu Bala <sb10@sanger.ac.uk>.
//
// This file is part of fastq_stats.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"os"

	"compress/gzip"
)

const longestRead = 100000
const numBases = 5
const highestScore = 127

var l = log.New(os.Stderr, "", 0)

func main() {
	// arg handling
	var help = flag.Bool("h", false, "print help text")
	flag.Parse()
	if *help {
		exitHelp("")
	}

	if len(flag.Args()) == 0 {
		exitHelp("ERROR: you must provide fastq input files")
	}
	fastqPaths := flag.Args()

	// prepare arrays to store counts in
	var bases [longestRead][numBases]int
	var quals [longestRead][highestScore]int

	// parse the fastq files to generate our stats
	parseFastqs(fastqPaths, &bases, &quals)

	// print out the results
	printStats(&bases, &quals)
}

// exitHelp prints help text and exits 0, unless a message is passed in which
// case it also prints that and exits 1.
func exitHelp(msg string) {
	if msg != "" {
		l.Println(msg)
		fmt.Printf("\n")
	}

	fmt.Printf(`fastq_stats reports summary stats on fastq files.

You use it to get an overview of what bases were called at what quality at each
base position. Useful for seeing what difference re-calling bases made (compare
the output of this program on original fastqs to the output on re-called
fastqs).

Usage: fastq_stats *.fastq.gz
Options:
  -h          this help text
`)

	if msg != "" {
		os.Exit(1)
	}
	os.Exit(0)
}

func die(err error) {
	l.Printf("ERROR: %s", err.Error())
	os.Exit(1)
}

func parseFastqs(paths []string, bases *[longestRead][numBases]int, quals *[longestRead][highestScore]int) {
	for _, path := range paths {
		f, err := os.Open(path)
		if err != nil {
			die(err)
		}
		gzf, err := gzip.NewReader(f)
		if err != nil {
			die(err)
		}
		calculateStats(gzf, bases, quals)
	}
}

func calculateStats(r io.Reader, bases *[longestRead][numBases]int, quals *[longestRead][highestScore]int) int64 {
	var sum int64
	scanner := bufio.NewScanner(r)
	line := 0
	for scanner.Scan() {
		line++
		if line == 5 {
			line = 1
		}

		b := scanner.Bytes()

		switch line {
		case 1, 3:
			continue
		case 2:
			handleBases(b, bases)
		case 4:
			handleQuals(b, quals)
		}
	}
	return sum
}

func handleBases(b []byte, store *[longestRead][numBases]int) {
	for i, base := range b {
		var j int
		switch base {
		case 'A':
			j = 0
		case 'C':
			j = 1
		case 'G':
			j = 2
		case 'T':
			j = 3
		case 'N':
			j = 4
		default:
			die(fmt.Errorf("unknown base %s\n", string(base)))
		}

		store[i][j]++
	}
}

func handleQuals(b []byte, store *[longestRead][highestScore]int) {
	for i, qual := range b {
		store[i][qual]++
	}
}

func printStats(bases *[longestRead][numBases]int, quals *[longestRead][highestScore]int) {
	fmt.Println("bases (A,C,G,T,N):")
	last := 0
	for i, c := range bases {
		if c[0]+c[1]+c[2]+c[3]+c[4] == 0 {
			last = i
			break
		}
		fmt.Printf("%d,%d,%d,%d,%d,%d\n", i+1, c[0], c[1], c[2], c[3], c[4])
	}

	fmt.Println("quals (0..93):")
	for i, c := range quals {
		if i == last {
			break
		}
		fmt.Printf("%d", i+1)
		for j, n := range c {
			if j < 33 {
				continue
			}
			fmt.Printf(",%d", n)
		}
		fmt.Printf("\n")
	}
}
