## ANALYSIS FOR CR9114 DEEP MUTATIONAL SCANNING
This experiment aims to search for CR9114 mutants that can overcome the virus resistance.

### INPUT FILE
* All raw sequencing reads, which can be downloaded from NIH SRA database [PRJNA510700](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA510700), should be placed in fastq/ folder. The filename for read 1 should match those described in [./doc/SampleID.tsv](./doc/SampleID.tsv). The filename for read 2 should be the same as read 1 except "R1" is replaced by "R2".
* [./doc/SampleID.tsv](./doc/SampleID.tsv): Describes the sample identity for each fastq file.

### ANALYSIS PIPELINE
1. [./script/CR9114\_read\_to\_count.py](./script/CR9114\_read\_to\_count.py): Converts raw reads to variant counts.
    - Input files:
      - Raw sequencing reads in fastq/ folder
    - Output files: count/count\_\*.tsv
2. [./script/CR9114\_count\_to\_freq.py](./script/CR9114\_count\_to\_freq.py): Converts counts into frequencies.
    - Input files: count/count\_\*.tsv
    - Output file: [./data/VariantFreqTable.tsv](./data/VariantFreqTable.tsv)

### PLOTTING
1. [./script/CR9114\_plot\_TopClonesFreq.R](./script/CR9114\_plot\_TopClonesFreq.R): Plot the frequencies of top clones in each selection.
    - Input file: [./data/VariantFreqTable.tsv](./data/VariantFreqTable.tsv)
    - Output files: graph/Freq\_YDisplay\_\*.png
