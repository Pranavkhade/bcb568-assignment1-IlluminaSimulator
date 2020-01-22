##################################################
## Project: Biostatistics
## Script purpose:  Ilumna Seq simulator
## Date: Tue Jan 21 13:40:48 2020
## Author: Henri C Chung
## Revisions: 
##################################################

## Section A: Set up libraries and load data
##################################################
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ShortRead")){
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install("ShortRead")
}
library(ShortRead)
library(tidyverse)


#take in command line arguments with defaults
args = commandArgs(trailingOnly=TRUE, defaults = c("example.fastq", 1, 0.02, 3, T, F))
true_file = args[1]
c = args[2]
p_error = args[3]
dist_scale = args[4]
pairwise = args[5]
reference = args[6]

#Calculate probability and quality distributions from reference fastq files
#reference = T means a reference RDS file has already been computed and is in the data
{
if(reference == F){ 
  dataDir <- "data/"
  files <- list.files(dataDir, pattern = ".*.fastq")
  temp <- ShortRead::readFastq(paste(dataDir, files[1], sep = ""))

  #Calculate probabilities from reference fasta.
  nuc_mean_list <- list()
  scores_freq_list <- list()

  for(f in 1:length(files)){
    nuc_BString <- ShortRead::sread(temp)
    nuc_freq <- ShortRead::alphabetFrequency(nuc_BString)
    nchars <- nchar(nuc_BString)
    nuc_prob <- t(t(nuc_freq)/nchars)
    nuc_mean <- base::apply(nuc_prob, MARGIN = 2, FUN = mean)
    nuc_mean_list[[f]] <- nuc_mean

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
    scores_freq_list[[f]] <- scores_freq
  }

  nuc_mean_table <- data.frame(do.call(rbind, nuc_mean_list))
  colnames(nuc_mean_table) <- names(nuc_mean_list[[1]])
  nuc_means <- colMeans(nuc_mean_table)

  scores_freq_table<- dplyr::bind_rows(scores_freq_list) %>%
    group_by(shortreads, qualityscores) %>%
    summarise(n = n(), freq = sum(freq)) %>%
    mutate(freq = freq / n) %>%
    ungroup()

  saveRDS(scores_freq_table, paste(dataDir, "scores_freq_table.RDS", sep = "") )
  saveRDS(nuc_means, paste(dataDir, "nuc_means.RDS", sep = "") )
}else{
  scores_freq_table <- readRDS(paste(dataDir, "scores_freq_table.RDS", sep = ""))
  nuc_means <- readRDS(paste(dataDir, "nuc_means.RDS", sep = ""))}
}


if(pairwise == T){c = c * 2} #Because no part of the model relies on nucleotide distance; a pairwise read can be created by 
                            #repeating the first read and inversing it. 

temp <- readFastq(paste(dataDir, true_file, sep = ""))
true_seq <- as.character(ShortRead::sread(temp)[1])
n = nchar(true_seq)

char_list <- list()
score_list <- list()

for(i in 1:c){
temp_char <- NULL
for(j in 1:n){
  err <- runif(1,0,1)
  p_error <- exp(0+n*(dist_scale/100000)) * p_error
  if(err < p_error){
    temp_char[j] <- sample(size = 1, names(nuc_means), prob = nuc_means)
  }else{
    temp_char[j] <- substr(true_seq, j, j)
  }
}

scores <- NULL
for(l in 1:length(temp_char)){
  nuc <- temp_char[l] 
  probs <- scores_freq_table %>% 
    filter(shortreads == nuc)
  scores[l] <- sample(size = 1, probs$qualityscores, probs$freq, replace = T)
}

char_list[[i]] <- base::paste(temp_char, collapse = "")
score_list[[i]] <- base::paste(scores, collapse = "")
}
names(char_list)[seq(1, c, 2)] <- paste("F", 1:length(seq(1,c,2)))
names(char_list)[seq(2, c, 2)] <-  paste("R", 1:length(seq(2,c,2)))
names(score_list)[seq(1, c, 2)] <- paste("F", 1:length(seq(1,c,2)))
names(score_list)[seq(2, c, 2)] <- paste("R", 1:length(seq(2,c,2)))

results <- do.call(rbind, Map(data.frame, read=char_list, quality=score_list))
write.csv(results, "results.csv")
