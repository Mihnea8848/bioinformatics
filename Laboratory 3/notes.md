# Laboratory 3

### Exercise 1
we will analyze a specifc virus/bacteria. Design a software app which is able to open a .fasta file containing a viral or bacterial genome and that is able to detect the gene frequencies for that specific codon. Show in output the top10 highest freq for these codons. In order to find the freq for all 64 codons a sliding window of exactly 3 positions in which you add one to a counter for the current codon. Create a dictionnary to associate to count the frequency {"AUG" : 0, ...}