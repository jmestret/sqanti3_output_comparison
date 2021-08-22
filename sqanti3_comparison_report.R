#######################################
#                                     #
#  SQANTI3 output comparison report   #
#             generation              #
#                                     #
#######################################


# Author: Jorge Martinez Tomas & Alejandro Paniagua
# Last modified: 18/08/2021 by Jorge Martinez

# Arguments
# $ sqanti3_comparison_report.R dir_in report_name dir_out
# dir_in: Name of the directory containing the SQANTI3 output files (classification and junction files are required)
# report_name: Output name for the HTML report (without extension)
# dir_out: Output directory for the report and CSV file (working directory as default)

# -------------------- Argument parser

library(optparse)

option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = NULL,
              help="directory with input files (classification and junction files)",
              metavar = "DIRIN"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help="Output directory for the report and CSV file [default= %default]",
              metavar = "DIROUT"),
  make_option(c("-n", "--name"), type = "character", default = "comparison_output",
              help="Output name for the HTML report and CSV file (without extension) [default= %default]",
              metavar = "OUTNAME")
)

opt_parser = OptionParser(
  usage = "usage: %prog [-i DIRIN] [-o DIROUT] [-n OUTNAME]",
  option_list=option_list
  )
opt = parse_args(opt_parser)

directory <- opt$dir
output_directory <- opt$outdir
output_name <- opt$name

if (is.null(directory)) {
  stop("You must supply at least the -d argument (directory containing classification and junctions files)")
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

if (length(class_in) != length(junct_in)){
  stop("ERROR: There are different numbers of classification and junction files in the directory")
} else if (length(class_in) == 0){
  stop(paste0("ERROR: No classification and junction files were found in the directory ", dir_in))
}

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
      entropyTTS = entropy::entropy(TTS_genomic_coord),
      isoform = list(isoform)
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

    
    class.filtered <- filter_monoexon(l_iso[[i]][[1]]) # delete monoexons
    
    tags <- isoformTags(l_iso[[i]][[2]]) # build tags
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

# ISOFORM TO GENOME BROWSER

iso2url <- function(id){
  if (id != "NA") {
    id.split <- str_split(id,"_")
    chr <- id.split[[1]][2]
    start <- id.split[[1]][4]
    end <- id.split[[1]][length(id.split[[1]])]
    name <- substr(id, 1, 30)
    url <- paste0(
      "<a href='https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
      chr, "%3A", start, "%2D", end, "&hgsid=1143169919_jAraPbUWtMCdAHfgTHk4sDHQHW7R",
      "'>", name, "...</a>")
    return(url)
  } else { return(NA) }
}


# -------------------- 
# --------------------  RUN COMPARISON
res <- compareTranscriptomes(f_in)


# -------------------- Table and Plot generation

# -------------------- 
# -------------------- 
# TABLE INDEX
# t1: summary table
# t2: presence/ausence table

# -------------------- 
# -------------------- 
# PLOT INDEX
# p1: gene characterization
#   p1.1: isoforms per gene
#   p1.2: exon structure
# p2: structural category distribution
# p3: splice junction distribution for each SC
# p4: distance to TSS
# p5: distance to TTS
# p6: distance to CAGE peak
# p7: bad quality features
# p8: good quality features
# p9: Venn diagrams 
# p10: UpSet plot
# p11: Venn diagrams for SC
# p12: UpSet plot for SC
# p13: TSS standard deviation
# p14: TTS standard deviation
# p15: TSS entropy
# p16: TTS entropy

# -------------------- 
# TABLE 1: summary table

str_cat <- c("antisense", "full-splice_match", "fusion", "genic", "incomplete-splice_match", "intergenic", "novel_in_catalog", "novel_not_in_catalog")
df_summary <- data.frame(ID = names(res[[2]][3:ncol(res[[2]])]))

# Count total isoform from SQANTI3
n <- c()
for (i in f_in) {
  n <- c(n, nrow(i[[1]]))
}
df_summary$total <- n

# Count unique tags
n <- c()
for (i in res[[1]]) {
  n <- c(n, nrow(i))
}
df_summary$uniq_id <- n

# Count unique tags for each category
for (i in str_cat) {
  n <- c()
  for (j in res[[1]]) {
    n <- c(n, nrow(j[j$structural_category == i,]))
  }
  df_summary[,i] <- n
}

colnames(df_summary) <- 
  c("ID","total", "unique_tag","antisense", "FSM", "fusion", "genic", "ISM", "intergenic", "NIC", "NNC")

t1 <- DT::datatable(df_summary,
              extensions = "Buttons",
              options = list(
                order = list(list(5, "asc"),list(1, "desc")),
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'pdf', 'print')
              ),
              escape = FALSE,
              caption = "Table 1. Summary from SQANTI3 output comparison")

# TABLE 2: presence/ausence isoforms

df.PA <- res[[3]]
df.PA$TAGS <- lapply(df.PA$TAGS, iso2url)
t2 <- DT::datatable(df.PA,
              escape = FALSE,
              options = list(
                pageLength = 10,
                autoWidth = TRUE,
                columnDefs = list(list(width = '10px', targets = "_all"))
              ),
              rownames = FALSE,
              caption = "Table 2. Presence/Ausence of all the isoform models")

# -------------------- 
# PLOT 1: gene characterization
# PLOT 1.1: isoforms per gene

countpergene <- c()
exonstructure <- c()
for (data in f_in){
  data.class <- data[[1]]
  
  data.class$novelGene <- "Annotated Genes"
  data.class[grep("novelGene", data.class$associated_gene), "novelGene"] <- "Novel Genes"
  data.class$novelGene = factor(data.class$novelGene,
                                levels = c("Novel Genes","Annotated Genes"),
                                ordered=TRUE)
  
  isoPerGene = aggregate(data.class$isoform,
                         by = list("associatedGene" = data.class$associated_gene,
                                   "novelGene" = data.class$novelGene,
                                   "FSM_class" = data.class$FSM_class),
                         length)
  
  data.class[which(data.class$exons>1), "exonCat"] <- "Multi-Exon"
  data.class[which(data.class$exons==1), "exonCat"] <- "Mono-Exon"
  data.class$exonCat = factor(data.class$exonCat,
                              levels = c("Multi-Exon","Mono-Exon"),
                              ordered=TRUE)
  
  canonical.labels=c("Canonical", "Non-canonical")
  data.class$all_canonical = factor(data.class$all_canonical,
                                    labels=canonical.labels,
                                    levels = c("canonical","non_canonical"),
                                    ordered=TRUE)
  
  countpergene <- c(
    countpergene,
    sum(isoPerGene$x == 1),
    sum(isoPerGene$x == 2 | isoPerGene$x == 3),
    sum(isoPerGene$x == 4 | isoPerGene$x == 5),
    sum(isoPerGene$x >= 6)
  )
  
  exonstructure <- c(
    exonstructure,
    sum(data.class$novelGene == "Novel Genes" & data.class$exonCat == "Mono-Exon"),
    sum(data.class$novelGene == "Novel Genes" & data.class$exonCat == "Multi-Exon"),
    sum(data.class$novelGene == "Annotated Genes" & data.class$exonCat == "Mono-Exon"),
    sum(data.class$novelGene == "Annotated Genes" & data.class$exonCat == "Multi-Exon")
    
  )
  
  
}

sample <- c(rep(names(f_in), each=4))
number <- rep(c("1","2-3","4-5", ">=6"), times=length(f_in))
isoPerGene <- data.frame(sample, number, countpergene)

p1.1 <- ggplot(isoPerGene, aes(fill=number, y=countpergene, x=sample)) +
  geom_bar(position = "fill", stat = "identity")



# PLOT 1.2: exon structure

category <- rep(c("Novel-Mono", "Novel-Multi", "Annotated-Mono", "Annotated-Multi"), times=length(f_in))
exonstructure <- data.frame(sample, category, exonstructure)

p1.2 <- ggplot(exonstructure, aes(fill=category, y=exonstructure, x=sample)) +
  geom_bar(position = "fill", stat = "identity")

# PLOT 2: structural category distribution

df <- df_summary
df$total <- NULL
df$unique_tag <- NULL
df <- df %>% 
  pivot_longer(!"ID", "SC")

p2 <- ggplot(df, aes(fill=SC, y=value, x=ID)) +
  geom_bar(position = "fill", stat = "identity")

# PLOT 3: splice junction distribution for each SC

#!!!!!!!!!!! ADD THESE PLOTS

# PLOT 4: distance to TSS

sample <- c()
dist <- c()
for (i in 1:length(f_in)){
  data.class <- f_in[[i]][[1]]
  cond <- which(data.class$structural_category == "full-splice_match")
  x <- data.class[cond, "diff_to_TSS"]
  sample <- c(
    sample,
    rep(names(f_in)[i], times=length(x))
  )
  dist <- c(
    dist,
    x
  )
}
TSS.dist.FSM <- data.frame(sample, dist)
p4.1 <- ggplot(TSS.dist.FSM, aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample))

sample <- c()
dist <- c()
for (i in 1:length(f_in)){
  data.class <- f_in[[i]][[1]]
  cond <- which(data.class$structural_category == "incomplete-splice_match")
  x <- data.class[cond, "diff_to_TSS"]
  sample <- c(
    sample,
    rep(names(f_in)[i], times=length(x))
  )
  dist <- c(
    dist,
    x
  )
}
TSS.dist.ISM <- data.frame(sample, dist)
p4.2 <- ggplot(TSS.dist.ISM, aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample))
 

# PLOT5: distance to TTS

sample <- c()
dist <- c()
for (i in 1:length(f_in)){
  data.class <- f_in[[i]][[1]]
  cond <- which(data.class$structural_category == "full-splice_match")
  x <- data.class[cond, "diff_to_TTS"]
  sample <- c(
    sample,
    rep(names(f_in)[i], times=length(x))
  )
  dist <- c(
    dist,
    x
  )
}
TSS.dist.FSM <- data.frame(sample, dist)
p5.1 <- ggplot(TSS.dist.FSM, aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample))

sample <- c()
dist <- c()
for (i in 1:length(f_in)){
  data.class <- f_in[[i]][[1]]
  cond <- which(data.class$structural_category == "incomplete-splice_match")
  x <- data.class[cond, "diff_to_TTS"]
  sample <- c(
    sample,
    rep(names(f_in)[i], times=length(x))
  )
  dist <- c(
    dist,
    x
  )
}
TSS.dist.ISM <- data.frame(sample, dist)
p5.2 <- ggplot(TSS.dist.ISM, aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample))

# PLOT 6: distance to CAGE peak

sample <- c()
dist <- c()
for (i in 1:length(f_in)){
  data.class <- f_in[[i]][[1]]
  cond <- which(data.class$structural_category == "full-splice_match")
  x <- data.class[cond, "diff_to_TSS"]
  sample <- c(
    sample,
    rep(names(f_in)[i], times=length(x))
  )
  dist <- c(
    dist,
    x
  )
}
TSS.dist.FSM <- data.frame(sample, dist)
p6.1 <- ggplot(TSS.dist.FSM, aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample))

sample <- c()
dist <- c()
for (i in 1:length(f_in)){
  data.class <- f_in[[i]][[1]]
  cond <- which(data.class$structural_category == "incomplete-splice_match")
  x <- data.class[cond, "dist_to_cage_peak"]
  sample <- c(
    sample,
    rep(names(f_in)[i], times=length(x))
  )
  dist <- c(
    dist,
    x
  )
}
TSS.dist.ISM <- data.frame(sample, dist)
p6.2 <- ggplot(TSS.dist.ISM, aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample))

# PLOT 7: bad quality features
# PLOT 7.1: RT-switching
FSM <- c()
NIC <- c()
NNC <- c()

for (i in 1:length(f_in)){
  data.class <- f_in[[i]][[1]]
  df <- group_by(data.class, structural_category, RTS_stage) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
  FSM.match <- df$count[which(df$structural_category == "full-splice_match")]
  FSM <- c(FSM, ((FSM.match[2]/(FSM.match[1]+FSM.match[2]))*100))
  
  NIC.match <- df$count[which(df$structural_category == "novel_in_catalog")]
  NIC <- c(NIC, ((NIC.match[2]/(NIC.match[1]+NIC.match[2]))*100))
  
  NNC.match <- df$count[which(df$structural_category == "novel_not_in_catalog")]
  NNC <- c(NNC, ((NNC.match[2]/(NNC.match[1]+NNC.match[2]))*100))
}

sample <- names(f_in)
FSM <- data.frame(sample, FSM)
NIC <- data.frame(sample, NIC)
NNC <- data.frame(sample, NNC)

p7.1.1 <- ggplot(FSM, aes(x=sample, y=FSM)) + geom_bar(color="blue",fill=rgb(0.1,0.4,0.5,0.7), stat="identity")
p7.1.2 <- ggplot(NIC, aes(x=sample, y=NIC)) + geom_bar(color="blue", fill=rgb(0.1,0.4,0.5,0.7),stat="identity")
p7.1.3 <- ggplot(NNC, aes(x=sample, y=NNC)) + geom_bar(color="blue", fill=rgb(0.1,0.4,0.5,0.7),stat="identity")

# p8: good quality features

# PLOT 9: Venn diagrams

l <- list()
for (i in 3:ncol(res[[2]])){
  l[[i-2]] <- res[[2]][,i]
}
names(l) <- colnames(res[[2]])[3:ncol(res[[2]])]

p9 <- ggvenn(
  l, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
) + ggtitle("All isoforms") +
  theme(plot.title = element_text(hjust = 0.5))

# PLOT 10: UpSet plot

p10 <-
  UpSetR::upset(
    UpSetR::fromList(l),
    order.by = "freq",
    mainbar.y.label = "Isoform Intersections",
    sets.x.label = "Isoforms Per Sample"
  )

# PLOT 11: Venn diagrams for SC

p11 <- list()
contador <- 1
for (i in 1:length(res[[2]])) {
  a <- res[[2]][res[[2]]$structural_category == str_cat[i], names(res[[2]])[3:length(names(res[[2]]))]]
  l <- list()
  for (j in 1:ncol(a)) {
    l[[j]] <- na.omit(a[,j])
  }
  names(l) <- names(a)
  p11[[contador]] <- ggvenn(
    l, 
    fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 4
  ) + ggtitle(str_cat[contador]) +
    theme(plot.title = element_text(hjust = 0.5))
  contador <- contador + 1
}

# PLOT 12: UpSet plots for SC

p12 <- list()
for (i in 1:length(res[[2]])) {
  a <- res[[2]][res[[2]]$structural_category == str_cat[i], names(res[[2]])[3:length(names(res[[2]]))]]
  l <- list()
  for (j in 1:ncol(a)) {
    l[[j]] <- na.omit(a[,j])
  }
  names(l) <- names(a)
  
  p12[[i]] <-
    UpSetR::upset(
      UpSetR::fromList(l),
      order.by = "freq",
      mainbar.y.label = "Isoform Intersections",
      sets.x.label = "Isoforms Per Sample"
    )
}

# PLOT 13: TSS standard deviation

a <- bind_rows(res$classifications, .id = "pipeline")
x <- colnames(a)[6]
p13 <- ggplot(a, aes(log2(a[,x]))) +
  geom_density(aes(col = pipeline)) + xlab(paste0("log2(",x,")"))

# PLOT 14: TTS standard deviation

x <- colnames(a)[7]
p14 <- ggplot(a, aes(log2(a[,x]))) +
  geom_density(aes(col = pipeline)) + xlab(paste0("log2(",x,")"))

# p15: TSS entropy

x <- colnames(a)[8]
p15 <- ggplot(a, aes(log2(a[,x]))) +
  geom_density(aes(col = pipeline)) + xlab(paste0("log2(",x,")"))

# PLOT 16: TTS entropy

x <- colnames(a)[9]
p16 <- ggplot(a, aes(log2(a[,x]))) +
  geom_density(aes(col = pipeline)) + xlab(paste0("log2(",x,")"))


# -------------------- Output report

Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio/bin/pandoc")
rmarkdown::render(
  input = paste(getwd(), "SQANTI3_comparison_report.Rmd", sep = "/"),
  output_dir = output_directory,
  output_file = paste0(output_name, ".html")
)

# -------------------- Output csv

write.csv(res$comparisonPA, paste0(output_directory, "/", output_name, ".csv"))

