#######################################
#                                     #
#  SQANTI3 output comparison report   #
#             generation              #
#                                     #
#######################################


# Author: Jorge Martinez Tomas & Alejandro Paniagua
# Last modified: 18/08/2021 by Jorge Martinez

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
  stop("\n\nAt least one argument must be supplied.\nThe -d argument is required (directory containing classification and junctions files)")
}
  

# -------------------- Packages

library(DT)
library(gridExtra)
library(knitr)
library(rmarkdown)
library(tidyverse)
library(UpSetR)
library(VennDiagram)


# -------------------- Load data

if (dir.exists(directory)){
  dir_in <- directory
} else {
  dir_in <- paste(getwd(), directory, sep="/")
  if (!dir.exists(dir_in)){
    stop(paste0("\n\nCould not find the input directory (", directory, ").\nPlease enter a valid path"))
  }
}


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
  stop("ERROR: There is a different number of classification and junction files in the directory")
} else if (length(class_in) == 0){
  stop(paste0("ERROR: No classification and junction files were found in the directory: ", dir_in))
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
# Builds isoform ids (tags) with the following structure: Chr_strand_start_end
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
# Deletes monoexons and transcripts from rare chromosomes
filter_monoexon <- function(class_file){
  filtered_classification <- class_file[class_file$exons > 1, ]
  filtered_classification[order(filtered_classification$isoform),]
}

filter_monoexon_sirv <- function(class_file){
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

# REVERSE SWAP
reverseswap <- function(class_file) {
  class_file <- swapcoord(class_file)
  sdtts <- class_file$sdTSS[class_file$strand == "-"]
  sdtss <- class_file$sdTTS[class_file$strand == "-"]

  class_file$sdTSS[class_file$strand == "-"] <- sdtss
  class_file$sdTTS[class_file$strand == "-"] <- sdtts
  return(class_file)
}


# AGGREGATE TAGS | Calculate sd, median TSS and TTS
## WARNING: Calculating the median is the slowest part
uniquetag <- function(class_file) {
  dt <- data.table::data.table(class_file)
  dt.out <-
    dt[, list(
      TSS_genomic_coord=list(TSS_genomic_coord),
      TTS_genomic_coord=list(TTS_genomic_coord),
      medianTSS = round(median(TSS_genomic_coord)),
      medianTTS = round(median(TTS_genomic_coord)),
      sdTSS = sd(TSS_genomic_coord),
      sdTTS = sd(TTS_genomic_coord),
      isoform = list(isoform)
    ), by = c("tags", "structural_category")]
  dt.out <- as.data.frame(dt.out)
  return(dt.out[order(dt.out$tags, dt.out$structural_category),])
}

uniquetag_simple <- function(class_file) {
  dt <- data.table::data.table(class_file)
  dt.out <-
    dt[, list(
      isoform = list(isoform)
    ), by = c("tags", "structural_category")]
  dt.out <- as.data.frame(dt.out)
  return(dt.out[order(dt.out$tags, dt.out$structural_category),])
}

# ADD THE SC TO THE TAG
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


#  COMPARE TAGS | PRESENCE AUSENCE
multipleComparison <- function(l_class){
  a <- c(rbind(names(l_class), paste0(names(l_class), "SC")))
  
  l_class %>%
    purrr::map(~ data.frame(col = .$tags, .$tags,.$structural_category, stringsAsFactors = FALSE)) %>%
    purrr::reduce(full_join, by = "col") %>%
    select(-col) %>%
    setNames(a)
}


# FINAL COMPARISON FUNCTION
# Given a list of pairs of classification and junction files
# generates a data.frame with all the common and unique transcripts
compareTranscriptomes <- function(l_iso){
  # Check for TSS and TTS genomic coords
  TSS_TTS_coord <- TRUE
  for ( i in 1:length(l_iso)){
    if (!("TSS_genomic_coord" %in% colnames(l_iso[[i]][[1]]))){
      TSS_TTS_coord <- FALSE
    }
  }
  
  # Filter monoexon, add tag column and aggregate
  n <- names(l_iso)
  l_class <- list()
  for ( i in 1:length(l_iso)) {
    class.filtered <- filter_monoexon(l_iso[[i]][[1]]) # delete monoexons
    
    tags <- isoformTags(l_iso[[i]][[2]]) # build tags
    if (length(tags) != nrow(class.filtered)){
      class.filtered <- filter_monoexon_sirv(l_iso[[i]][[1]])
    }
    class.filtered[,"tags"] <- tags # add tags
    
    if (TSS_TTS_coord){
      class.swap <- swapcoord(class.filtered) # swap - strand
      class.uniquetags <- as.data.frame(uniquetag(class.swap)) # group by tag
      class.out <- reverseswap(class.uniquetags) # re-swap TSS and TSS - strand
    } else{
      class.out <- uniquetag_simple(class.filtered) # group by tag
    }
    
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


# -------------------- Output functions

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
# --------------------  Unique tag comparison

res <- try({
  compareTranscriptomes(f_in)
}, silent = TRUE)

if (class(res) == "try-error"){
  print("ERROR: An error has ocurred during the comparison. A partial report will be generated. Here's the original error message:")
  print(geterrmessage())
}

# --------------------
# -------------------- Basic comparison with classification files

##*****  Check for TSS and TTS genomic coords

TSS_TTS_coord <- TRUE
for ( i in 1:length(f_in)){
  if (!("TSS_genomic_coord" %in% colnames(f_in[[i]][[1]]))){
    TSS_TTS_coord <- FALSE
  }
}

##*****  Define max number of samples in plots
if (length(f_in) < 6){
  limit <- length(f_in)
} else {limit <- 5}


##*****  Vector of structural categories

str_cat <- c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "antisense", "fusion", "genic", "intergenic")

##***** Generates dataframe with a summary of the SQANTI3 classification files

df_summary.1 <- data.frame(ID = names(f_in))

# Count total isoform from SQANTI3
n <- c()
for (i in f_in) {
  n <- c(n, nrow(i[[1]]))
}
df_summary.1$total <- n

# Count isoforms for each structural category
for (i in str_cat) {
  n <- c()
  for (j in f_in) {
    k <- j[[1]]
    n <- c(n, nrow(k[k$structural_category == i,]))
  }
  df_summary.1[,i] <- n
}

colnames(df_summary.1) <- 
  c("ID","total", "FSM", "ISM", "NIC", "NNC", "antisense", "fusion", "genic", "intergenic")


##***** Generates dataframe with a summary of the unique tag comparison

df_summary.2 <- data.frame(ID = names(res[[2]][3:ncol(res[[2]])]))

# Count unique tags
n <- c()
for (i in res[[1]]) {
  n <- c(n, nrow(i))
}
df_summary.2$uniq_id <- n

# Count unique tags for each category
for (i in str_cat) {
  n <- c()
  for (j in res[[1]]) {
    n <- c(n, nrow(j[j$structural_category == i,]))
  }
  df_summary.2[,i] <- n
}

colnames(df_summary.2) <- 
  c("ID","tags", "FSM", "ISM", "NIC", "NNC", "antisense", "fusion", "genic", "intergenic")


##***** Add GenomeBrowser URL to the P/A table

df.PA <- res[[3]]
df.PA$TAGS <- lapply(df.PA$TAGS, iso2url)


##***** Counts per gene and exon structure

countpergene <- c()
exonstructure <- c()
for (i in 1:limit){
  data <- f_in[[i]]
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

sample <- c(rep(names(f_in[1:limit]), each=4))
number <- rep(c("1","2-3","4-5", ">=6"), times=limit)
isoPerGene <- data.frame(sample, number, countpergene)

category <- rep(c("Novel-Mono", "Novel-Multi", "Annotated-Mono", "Annotated-Multi"), times=limit)
exonstructure <- data.frame(sample, category, exonstructure)


##***** Summary dataframe pivoted

df_SC <- df_summary.1[1:limit,]
df_SC$total <- NULL
df_SC <- df_SC %>% 
  pivot_longer(!"ID", "SC")

##***** Distance to TSS, TTS and CAGE peak

dist.list <- list()
dist.msr <- c("diff_to_TSS", "diff_to_TTS", "dist_to_cage_peak")
dist.SC <- c("full-splice_match", "incomplete-splice_match")
contador <- 1
for (i in dist.msr){
  for (j in dist.SC){
    sample <- c()
    dist <- c()
    for (k in 1:limit){
      data.class <- f_in[[k]][[1]]
      cond <- which(data.class$structural_category == j)
      x <- data.class[cond, i]
      sample <- c(
        sample,
        rep(names(f_in)[k], times=length(x))
      )
      dist <- c(
        dist,
        x
      )
    }
    dist.list[[contador]] <- data.frame(sample, dist)
    names(dist.list)[length(dist.list)] <- paste(i,j,sep="_")
    contador <- contador + 1
  }
}

##***** RT-switching

FSM <- c()
NIC <- c()
NNC <- c()
for (i in 1:limit){
  data.class <- f_in[[i]][[1]]
  df <- group_by(data.class, structural_category, RTS_stage) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
  FSM.match <- df$count[which(df$structural_category == "full-splice_match")]
  FSM <- c(FSM, ((FSM.match[2]/(FSM.match[1]+FSM.match[2]))*100))
  
  NIC.match <- df$count[which(df$structural_category == "novel_in_catalog")]
  NIC <- c(NIC, ((NIC.match[2]/(NIC.match[1]+NIC.match[2]))*100))
  
  NNC.match <- df$count[which(df$structural_category == "novel_not_in_catalog")]
  NNC <- c(NNC, ((NNC.match[2]/(NNC.match[1]+NNC.match[2]))*100))
}

sample <- names(f_in)[1:limit]
FSM.RT <- data.frame(sample, FSM)
NIC.RT <- data.frame(sample, NIC)
NNC.RT <- data.frame(sample, NNC)

##***** List of unique tags for each sample

l <- list()
for (i in 3:ncol(res[[2]])){
  l[[i-2]] <- na.omit(res[[2]][,i])
}
names(l) <- colnames(res[[2]])[3:ncol(res[[2]])]


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


# -------------------- Global plot parameters
# COPY-PASTE FROM SQANTI3 REPORT

myPalette = c("#6BAED6","#FC8D59","#78C679","#EE6A50","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")

mytheme <- theme_classic(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=12) ) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size=12), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=15, hjust = 0.5)) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm"))

# -------------------- 
# TABLE 1: summary table

t1.1 <- DT::datatable(df_summary.1,
              extensions = "Buttons",
              options = list(
                order = list(list(5, "asc"),list(1, "desc")),
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'pdf', 'print')
              ),
              escape = FALSE,
              caption = "Table 1. Summary from SQANTI3 output comparison")

t1.2 <- DT::datatable(df_summary.2,
                      extensions = "Buttons",
                      options = list(
                        order = list(list(5, "asc"),list(1, "desc")),
                        dom = 'Bfrtip',
                        buttons = c('copy', 'csv', 'pdf', 'print')
                      ),
                      escape = FALSE,
                      caption = "Table 2. Summary from isoform id (tag) comparison")

# TABLE 2: presence/ausence isoforms

t2 <- DT::datatable(df.PA,
              escape = FALSE,
              options = list(
                pageLength = 10,
                autoWidth = TRUE,
                columnDefs = list(list(width = '10px', targets = "_all"))
              ),
              rownames = FALSE,
              caption = "Table 3. Presence/Ausence of all the isoform models")


# -------------------- 
# PLOT 1: gene characterization
# PLOT 1.1: isoforms per gene

p1.1 <- ggplot(isoPerGene, aes(fill=number, y=countpergene, x=sample)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = myPalette) + mytheme 

# PLOT 1.2: exon structure

p1.2 <- ggplot(exonstructure, aes(fill=category, y=exonstructure, x=sample)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = myPalette) + mytheme

# PLOT 2: structural category distribution

p2 <- ggplot(df_SC, aes(fill=SC, y=value, x=ID)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = myPalette) + mytheme

# PLOT 3: splice junction distribution for each SC

#!!!!!!!!!!! ADD THESE PLOTS

# PLOT 4: distance to TSS

p4.1 <- ggplot(dist.list[[1]], aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample)) +
  scale_fill_manual(values = myPalette) + mytheme

p4.2 <- ggplot(dist.list[[2]], aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample)) +
  scale_fill_manual(values = myPalette) + mytheme

# PLOT5: distance to TTS

p5.1 <- ggplot(dist.list[[3]], aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample)) +
  scale_fill_manual(values = myPalette) + mytheme

p5.2 <- ggplot(dist.list[[4]], aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample)) +
  scale_fill_manual(values = myPalette) + mytheme

# PLOT 6: distance to CAGE peak

p6.1 <- ggplot(dist.list[[5]], aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample)) +
  scale_fill_manual(values = myPalette) + mytheme

p6.2 <- ggplot(dist.list[[6]], aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample)) +
  scale_fill_manual(values = myPalette) + mytheme

# PLOT 7: bad quality features
# PLOT 7.1: RT-switching

p7.1.1 <- ggplot(FSM.RT, aes(x=sample, y=FSM, fill=sample)) + geom_bar(stat="identity") +
  scale_fill_manual(values = myPalette) + mytheme
p7.1.2 <- ggplot(NIC.RT, aes(x=sample, y=NIC)) + geom_bar(color="blue", fill=rgb(0.1,0.4,0.5,0.7),stat="identity") +
  scale_fill_manual(values = myPalette) + mytheme
p7.1.3 <- ggplot(NNC.RT, aes(x=sample, y=NNC)) + geom_bar(color="blue", fill=rgb(0.1,0.4,0.5,0.7),stat="identity") +
  scale_fill_manual(values = myPalette) + mytheme

# p8: good quality features

# PLOT 9: Venn diagrams

# To not generate the .log files of VennDiagram
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

p9 <- venn.diagram(l, filename = NULL, fill = myPalette[1:length(l)])

# PLOT 10: UpSet plot

p10 <-
  UpSetR::upset(
    UpSetR::fromList(l),
    order.by = "freq",
    mainbar.y.label = "Isoform Intersections",
    sets.x.label = "Isoforms Per Sample",
    nintersects = 20
    #sets.bar.color = myPalette[1:length(l)]
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
  p11[[contador]] <- venn.diagram(l, filename = NULL, fill = myPalette[1:length(l)])
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
      sets.x.label = "Isoforms Per Sample",
      nintersects = 20
      #sets.bar.color = myPalette[1:length(l)]
    )
}

if (TSS_TTS_coord == TRUE) {
  # Calculate UJC SD
  a <- c("tags", rbind(paste0(names(res$classifications), "TSS"), paste0(names(res$classifications), "TTS")))
  TSS_TTS_params <- list()
  for (i in 1:length(res$classifications)){
    TSS_TTS_params[[i]] <- res$classifications[[i]][,c("tags", "TSS_genomic_coord", "TTS_genomic_coord")]
  }
  TSS_TTS_params <- TSS_TTS_params %>% 
    purrr::reduce(full_join, by="tags") %>% 
    setNames(a)
  
  a <- paste0(names(res$classifications), "TSS")
  b <- paste0(names(res$classifications), "TTS")
  allTSS <- TSS_TTS_params[, a]
  allTSS[allTSS == "NULL"] <- NA
  allTTS <- TSS_TTS_params[, b]
  allTTS[allTTS == "NULL"] <- NA
  
  TSS_TTS_df <- data.frame(tags=TSS_TTS_params$tags)
  
  # max and min value
  sapplycolumns <- function(data, func){
    tmp <- list()
    for (i in 1:ncol(data)){
      tmp[[names(data)[i]]] <- sapply(data[,i], func)
    }
    for (i in 1:length(tmp)){
      tmp[[i]][tmp[[i]]=="NULL"] <- NA
    }
    return(as.data.frame(tmp))
  }
  
  minNA <- function(x) ifelse(length(x) > 1, min(x), NA)
  
  
  maxTSS <- sapplycolumns(allTSS, max)
  maxTSS[maxTSS=="-Inf"] <- NA
  
  maxTTS <- sapplycolumns(allTTS, max)
  maxTTS[maxTTS=="-Inf"] <- NA
  
  minTSS <- sapplycolumns(allTSS, minNA)
  minTSS[minTSS=="Inf"] <- NA
  
  minTTS <- sapplycolumns(allTTS, minNA)
  minTTS[minTTS=="Inf"] <- NA
  
  minmaxTSS <- cbind(minTSS, maxTSS)
  minmaxTTS <- cbind(minTTS, maxTTS)
  
  TSS_TTS_df$minmax.SD.TSS <- apply(minmaxTSS, 1, function(x) sd(unlist(x), na.rm = TRUE))

  TSS_TTS_df$minmax.SD.TTS <- apply(minmaxTTS, 1, function(x) sd(unlist(x),na.rm = TRUE))

  # Median value
  
  medianTSS <- sapplycolumns(allTSS, median)
  medianTTS <- sapplycolumns(allTTS, median)
  
  TSS_TTS_df$median.SD.TSS <- apply(medianTSS, 1, function(x) sd(unlist(x),na.rm = TRUE))

  TSS_TTS_df$median.SD.TTS <- apply(medianTTS, 1, function(x) sd(unlist(x),na.rm = TRUE))
  
  # Max SD
  
  TSS_TTS_df$SD.TSS <- apply(TSS_TTS_df[,c("minmax.SD.TSS", "median.SD.TSS")],1,max)
  TSS_TTS_df$SD.TTS <- apply(TSS_TTS_df[,c("minmax.SD.TTS", "median.SD.TTS")],1,max)
  

  
  # PLOT 13: TSS standard deviation per pipeline
  
  a <- bind_rows(res$classifications, .id = "pipeline")
  p13 <- ggplot(a, aes(log2(a[,colnames(a)[8]]))) +
    geom_density(aes(col = pipeline)) + xlab(paste0("log2(",colnames(a)[8],")")) +
    scale_fill_manual(values = myPalette) + mytheme
  
  # PLOT 14: TTS standard deviation per pipeline
  
  p14 <- ggplot(a, aes(log2(a[,colnames(a)[9]]))) +
    geom_density(aes(col = pipeline)) + xlab(paste0("log2(",colnames(a)[9],")")) +
    scale_fill_manual(values = myPalette) + mytheme
  

  # PLOT 15: TSS and TTS SD UJC
    p15.1 <-  ggplot(TSS_TTS_df, aes(log2(TSS_TTS_df[,6]))) +
        geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
        mytheme + ggtitle(names(TSS_TTS_df)[6]) 
    
    p15.2 <-  ggplot(TSS_TTS_df, aes(log2(TSS_TTS_df[,7]))) +
      geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
      mytheme + ggtitle(names(TSS_TTS_df)[7]) 
}
# -------------------- Output report

Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio/bin/pandoc")
rmarkdown::render(
  input = paste(getwd(), "SQANTI3_comparison_report.Rmd", sep = "/"),
  output_dir = output_directory,
  output_file = paste0(output_name, ".html")
)

# -------------------- Output csv

write.csv(res$comparisonPA, paste0(output_directory, "/", output_name, ".csv"))

