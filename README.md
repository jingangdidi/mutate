# mutate
Introduce random mutations to each read in the FASTQ file based on a specified mutation ratio.

# workflow
1. Traverse each read in the FASTQ file, calculating the number of mutation sites based on sequence length and the specified mutation ratio (-r/--ratio).  
2. Randomly select the calculated number of mutation sites within the read sequence (excluding N, only selecting A/T/G/C), and record their positions (-p/--pos).  
3. For each mutation site, randomly choose a different nucleotide from the remaining three bases to substitute the original base (e.g., if the current base is A, select randomly from T, G, C).  
4. After all mutations are applied, save the modified read.

# download releases
[download release](https://github.com/jingangdidi/mutate/releases/tag/v0.1.0)

# usage
```
Usage: mutate.exe -f <fastq> -r <ratio> [-t <thread>] [-p] [-o <outpath>]

mutate fastq

Options:
  -f, --fastq       fastq file (support gz compression), separate with commas
  -r, --ratio       mutation ratio, 0~1
  -t, --thread      max thread number, default: 4
  -p, --pos         add mutation position in header description
  -o, --outpath     output path, default: ./
  --help, help      display usage information
```

# example
1. Randomly select 10% position of each read from FASTQ file for mutation.
```
./mutate -f test.fastq -r 0.1
```
2. Randomly select 10% position of each read from multiple fastq.gz files for mutation, recording the mutation sites within each read's header.
```
./mutate -f test_R1.fastq.gz,test_R2.fastq.gz -r 0.1 -p
```

# compare
**before mutation:**
```
@LH00128:306:22MTVJLT4:2:1101:25772:1056 1:N:0:CAAGGTGA+CGGCTATG
CNTGTATCCAAAGAGGTGGCTCTTGGTAGGCAGGAGGCAAAGGGTGGTATTTGTTCCAGCCAAGAAGGCCAGACCACTTCTCTGTAGCTCATGTTCTTCTGACTCACACCCCTCAGTCCTATGGGCTGAGTAACTNTTCCAGCCTTGTGA
+
I#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIII9IIIIIIIII#IIIIIIIIIIIIII
```
**after mutation (mutation ratio: 10%, mutation number: 150*0.1=15):**
```
@LH00128:306:22MTVJLT4:2:1101:25772:1056 1:N:0:CAAGGTGA+CGGCTATG [18, 21, 23, 39, 40, 41, 45, 52, 58, 77, 117, 124, 128, 138, 145]
CNTGTATCCAAAGAGGTCGCGCATGGTAGGCAGGAGGCCTGGGGAGGTATTAGTTCCCGCCAAGAAGGCCAGACCAGTTCTCTGTAGCTCATGTTCTTCTGACTCACACCCCTCAGGCCTATGTGCTCAGTAACTNTGCCAGCCCTGTGA
+
I#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIII9IIIIIIIII#IIIIIIIIIIIIII
```
**mutation position:**
```
CNTGTATCCAAAGAGGTGGCTCTTGGTAGGCAGGAGGCAAAGGGTGGTATTTGTTCCAGCCAAGAAGGCCAGACCACTTCTCTGTAGCTCATGTTCTTCTGACTCACACCCCTCAGTCCTATGGGCTGAGTAACTNTTCCAGCCTTGTGA
CNTGTATCCAAAGAGGTCGCGCATGGTAGGCAGGAGGCCTGGGGAGGTATTAGTTCCCGCCAAGAAGGCCAGACCAGTTCTCTGTAGCTCATGTTCTTCTGACTCACACCCCTCAGGCCTATGTGCTCAGTAACTNTGCCAGCCCTGTGA
                 *  * *               ***   *      *     *                  *                                       *      *   *         *      *
```

# 
