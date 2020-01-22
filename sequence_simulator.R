library(ShortRead)
library(tidyverse)

#FASTQ simulator using reference sequence

#Read in reference files
dataDir <- "data/"
files <- list.files(dataDir, pattern = ".*.fastq")
temp <- readFastq(paste(dataDir, files[1], sep = ""))

seq_length = 40
#this simulator uses rnorm for 
#Calculate nucleotide probability
##the probability of a given nucleotides is based on a normal distribution given the mean # and sd of 
##that nucleotide in each read of the reference sequence
nuc_BString <- sread(temp)
nuc_freq <- alphabetFrequency(nuc_BString)
nchars <- nchar(nuc_BString)
nuc_prob <- t(t(nuc_freq)/nchars)
nuc_mean <- apply(nuc_prob, MARGIN = 2, FUN = mean)
chars <- sample(size = seq_length,  alphabet(nuc_BString), prob = nuc_mean, replace = T)


#Quality scores is based on the simple probability of a given nucleotide having a specific quality score as
#the reference file
quality_BString <- quality(temp)
shortreads_list <- unlist(strsplit(paste(as.character(nuc_BString), sep = ""), split = ""))
qualityscores_list <- unlist(strsplit(paste(as.character(quality_BString@quality), sep = ""), split = ""))
scores_df <- as.data.frame(cbind(shortreads_list, qualityscores_list)); colnames(scores_df) <- c("shortreads", "qualityscores")
scores_freq <- scores_df %>%
  group_by(shortreads, qualityscores) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(qualityscores = as.character(qualityscores)) %>%
  mutate(freq = n / sum(n)) 


scores <- NULL
for(i in 1:length(chars)){
  nuc <- chars[i] 
  probs <- scores_freq %>% 
    filter(shortreads == nuc)
  scores[i] <- sample(size = 1, probs$qualityscores, probs$freq, replace = T)
}

chars
scores
