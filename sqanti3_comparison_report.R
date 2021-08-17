#######################################
#                                     #
#  SQANTI3 output comparison report   #
#             generation              #
#                                     #
#######################################


# Author: Jorge Martinez Tomas & Alejandro Paniagua
# Last modified: 16/08/2021 by Jorge Martinez

# Arguments
# $ sqanti3_comparison_report.R dir_in report_name dir_out
# dir_in: Name of the directory containing the SQANTI3 output files (classification and junction files are required)
# report_name: Output name for the HTML report (without extension)
# dir_out: Output directory for the report and CSV file (working directory as default)


# -------------------- Scan arguments

args <- commandArgs(trailingOnly = TRUE)

directory <- args[1]
if (length(args) == 0) {
  stop("You must supply at least one argument (directory containing classification and junctions files)")
} else if (length(args) == 1) {
  output_name <- "example_html_report"
  output_directory <- "."
} else if (length(args) == 2) {
  output_directory <- "."
} else if (length(args) > 3) {
  stop("You supplied more than the required arguments")
} else {
  output_name <- args[2]
  output_directory <- args[3]
}


# -------------------- Packages

library(DT)
library(entropy)
library(ggvenn)
library(ggplot2)
library(knitr)
library(rmarkdown)
library(tidyverse)
library(UpSetR)


# -------------------- Load data

# Input: Name of a directory with all classification and juntions files from
# the runs to compare
#
# Example:
# $ ls ~/home/dir_in
# sample1_classification.txt sample1_junctions.txt
# sample2_classification.txt sample2_junctions.txt
# sample3_classification.txt sample3_junctions.txt

dir_in <- paste(getwd(), directory, sep="/")

#!!!!!!! Find better implementation (probably with tidyverse)
class_in <-
  list.files(dir_in,
             pattern = "*_classification.txt",
             all.files = FALSE,
             full.names = TRUE)
junct_in <-
  list.files(dir_in,
             pattern = "*_junctions.txt",
             all.files = FALSE,
             full.names = TRUE)

f_in <- list()
for (i in 1:length(class_in)) {
  f <- class_in[[i]]
  start <- stringr::str_locate(f, dir_in)[[2]]
  end <- stringr::str_locate(f, "_classification.txt")[[1]]
  idx <- substring(f, (start+2), (end-1))
  classification <- read.table(class_in[[i]], header = T, sep = "\t")
  junctions <- read.table(junct_in[[i]], header = T, sep = "\t")
  f_in[[idx]] <- list(classification, junctions)
}


# -------------------- Functions

# CREATE TAGS
#TAG structure: Chr_start_end_star_end
isoformTags <- function(junctions_file) {
  df <- junctions_file[, c("isoform", "chrom", "strand")] # df with isoforms in *junctions.txt
  dt <- data.table::data.table(df)
  dt <- dt[,coord:=paste0(junctions_file$genomic_start_coord, "_", junctions_file$genomic_end_coord)]
  dt <-
    dt[, list(tagcoord = paste0(coord, collapse = "_")),
       by = c("isoform", "chrom", "strand")]
  df <- as.data.frame(dt)
  df <- df[order(df$isoform),]
  tag <- paste(df$chrom, df$strand, df$tagcoord, sep = "_")
  return(tag)
}


# DELETE MONOEXONS
filter_monoexon <- function(class_file){
  # Deletes monoexons and transcripts from rare chromosomes
  valid_chrom <- c(paste0("chr", 1:22), paste0("chr", c("X","Y")), paste0("SIRV", 1:7))
  
  filtered_classification <- class_file[class_file$exons > 1 &
                                          class_file$chrom %in% valid_chrom, ]
  filtered_classification[order(filtered_classification$isoform),]
}


# SWAP COORDINATES STRAND
swapcoord <- function(dfclass){
  tss <- dfclass$TTS_genomic_coord[dfclass$strand == "-"]
  tts <- dfclass$TSS_genomic_coord[dfclass$strand == "-"]
  dfclass$TSS_genomic_coord[dfclass$strand == "-"] <- tss
  dfclass$TTS_genomic_coord[dfclass$strand == "-"] <- tts
  return(dfclass)
}

reverseswap <- function(class_file) {
  class_file <- swapcoord(class_file)
  sdtts <- class_file$sdTSS[class_file$strand == "-"]
  sdtss <- class_file$sdTTS[class_file$strand == "-"]
  etts <- class_file$entropyTSS[class_file$strand == "-"]
  etss <- class_file$entropyTTS[class_file$strand == "-"]
  
  class_file$sdTSS[class_file$strand == "-"] <- sdtss
  class_file$sdTTS[class_file$strand == "-"] <- sdtts
  class_file$entropyTSS[class_file$strand == "-"] <- etss
  class_file$entropyTTS[class_file$strand == "-"] <- etts
  return(class_file)
}


# AGGREGATE TAGS | Calculate sd, entropy, min TSS and max TTS
## WARNING: Calculating the entropy is the slowest part
uniquetag <- function(class_file) {
  dt <- data.table::data.table(class_file)
  dt.out <-
    dt[, list(
      minTSS = min(TSS_genomic_coord),
      maxTTS = max(TTS_genomic_coord),
      sdTSS = sd(TSS_genomic_coord),
      sdTTS = sd(TTS_genomic_coord),
      entropyTSS = entropy::entropy(TSS_genomic_coord),
      entropyTTS = entropy::entropy(TTS_genomic_coord)
    ), by = c("tags", "structural_category")]
  dt.out <- as.data.frame(dt.out)
  return(dt.out[order(dt.out$tags, dt.out$structural_category),])
}

addSC <- function(class_file){
  str_cat <- c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog")
  shortSC <- c("FSM", "ISM", "NIC", "NNC")
  
  class_file$SC <- class_file$structural_category
  for (i in 1:length(str_cat)){
    cond <- which(class_file$SC == str_cat[i])
    class_file$SC[cond] <- shortSC[i]
  }
  class_file$tags <- paste0(class_file$SC, "_", class_file$tags)
  class_file$SC <- NULL
  return(class_file)
}


#  COMPARE TAGS
multipleComparison <- function(l_class){
  a <- c(rbind(names(l_class), paste0(names(l_class), "SC")))
  
  l_class %>%
    purrr::map(~ data.frame(col = .$tags, .$tags,.$structural_category, stringsAsFactors = FALSE)) %>%
    purrr::reduce(full_join, by = "col") %>%
    select(-col) %>%
    setNames(a)
}


# FINAL COMPARISON FUNCTION
compareTranscriptomes <- function(l_iso){
  # Given a list of pairs of classification and junction files
  # generates a data.frame with all the common and unique transcripts
  
  # Filter monoexon, add tag column and aggregate
  n <- names(l_iso)
  l_class <- list()
  for ( i in 1:length(l_iso)) {

    
    class.filtered <- filter_monoexon(l_iso[[i]][[1]]) # delete monoexon
    
    tags <- isoformTags(l_iso[[i]][[2]]) # buil tags
    class.filtered[,"tags"] <- tags # add tags
    
    class.swap <- swapcoord(class.filtered) # swap - strand
    class.uniquetags <- as.data.frame(uniquetag(class.swap)) # group by tag
    class.out <- reverseswap(class.uniquetags) # re-swap TSS and TSS - strand
    class.out <- addSC(class.out)
    l_class[[i]] <- class.out # add to list
  }
  
  names(l_class) <- n # add names
  comptags <- multipleComparison(l_class)
  comptags.SC <- cbind(structural_category =
                         do.call(dplyr::coalesce, comptags[,paste0(n,"SC")]),
                       comptags[,n]
  )
  comptags.out <- cbind( TAGS =
                           do.call(dplyr::coalesce, comptags[,n]),
                         comptags.SC
  )
  
  comptags.PA <- comptags.out
  comptags.PA[,3:ncol(comptags.PA)][!is.na(comptags.PA[,3:ncol(comptags.PA)])] <- 1
  comptags.PA[,3:ncol(comptags.PA)][is.na(comptags.PA[,3:ncol(comptags.PA)])] <- 0
  
  output <- list(classifications = l_class, comparison = comptags.out, comparisonPA = comptags.PA)
  
  return(output)
}


# -------------------- Output report

Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc")
rmarkdown::render(
  input = paste(getwd(), "SQANTI3_comparison_report.Rmd", sep = "/"),
  output_dir = output_directory,
  output_file = paste0(output_name, ".html")
)

# -------------------- Output csv

write.csv(res$comparisonPA, paste0(output_directory, "/", output_name, ".csv"))

