This program models ILUMINA pair end read process on an example sequence. 
The model uses categorical probability distributions derived from the average distribution of nucleotide/quality scores in a given collection of reference files.


The read is assembled by reading the example sequence and copying each nucleotide position. For each position, a random uniform [0,1] variable is drawn. If this variable is below a user input threshold, the nucleotide is considered an error an a nucleotide is drawn from the reference nucleotide distribution. A scaling parameter is used to define the rate at which the error threshold increases as nucleotide position increases. This simulates the increase in probability of read error as sequence length continues. 

Quality scores are calculated independently of nucleotide sequence and is drawn from a categorical probability distribution averaged from reference fastq files.

Arguments

file - (character) a "true" sequence file to be read.
c - (numeric) the number of reads the program should simulate
p_error - (numeric) probability of error from 0-1
dist_scale -  (numeric) a numeric value to scale the error threshold. Error threshold is calculated with the eq. exp(0+n*(dist_scale/100000)) * p_error
pairwise - (boolean) simulate a pairwise read or single read
reference - (boolean) are there RDS files with calculated categorical probabilities? If F then a data/ fold with reference fastq files are necessary.