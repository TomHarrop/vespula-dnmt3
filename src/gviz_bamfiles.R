#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)
library(Gviz)
library(GenomicFeatures)
library(GenomicRanges)
library(GenomicAlignments)

#############
# FUNCTIONS #
#############

# wrapper for alignment tracks
MakeAt <- function(x, bamfile, st, chromosome, start, end, aln_scheme) {
    my_at <- AlignmentsTrack(range = bamfile,
                             isPaired = FALSE,
                             # referenceSequence = st,
                             name = x,
                             chromosome = chr,
                             start = start,
                             end = end,
                             stacking = "dense",
                             fill = ifelse(grepl("V\\.", x),
                                           set1[2],
                                           set1[3]))
    displayPars(my_at) <- aln_scheme
    return(my_at) }


###########
# GLOBALS #
###########

options(ucscChromosomeNames=FALSE)

fasta_file <- snakemake@input[["fa"]]
gff_file <- snakemake@input[["gff"]]
txdb_file <- snakemake@output[["txdb"]]

# dev
# fasta_file <- "data/Vvulg.assembly.fna"
# gff_file <- "data/Vvulg_final_sorted.gff3"
# txdb_file <- "data/txdb.sqlite"

# bamfiles
# bamfiles <- list.files("output/02_minimap",
#                        pattern = "align.bam$",
#                        recursive = TRUE,
#                        full.names = TRUE)
bamfiles <- snakemake@input[["bamfiles"]]
names(bamfiles) <- gsub("^.*/Vvulg.([[:alpha:]]+)/.*$", "\\1", bamfiles)

spec_order <- c("Pcad" = " \nP. canadensis",
  "Pdomi" = " \nP. dominula",
  "Pdor" = " \nP. dorsalis",
  "Pfus" = " \nP. fuscatus",
  "Pmet" = " \nP. metricus",
  "Vgerm" = " \nV. germanica",
  "Vpens" = " \nV. pensylvanica")

names(bamfiles) <- plyr::revalue(names(bamfiles[names(spec_order)]), spec_order)


##################
# colour schemes #
##################

set1 <- viridis::viridis_pal()(4)

my_scheme <- list(
    shape = "smallArrow",
    background.title = "transparent",
    fontface.title = 1,
    col.title = "black",
    col.frame = "black",
    thinBoxFeature = c("lincRNA"),
    # fontcolor.group = set1[1],
    cex.title = 0.5,
    cex.group = 0.5,
    showTitle=TRUE
)

gat_scheme <- as.list(c(
    my_scheme,
    showTitle = TRUE,
    col = "black",
    fontcolor = "black",
    cex.title = 1,
    col.title = "black",
    cex = 0.5,
    showTitle = TRUE,
    cex.title = 0.5,
    rotation.title = 0,
    fontface.title = "bold"
))

grt_scheme <- as.list(c(
    my_scheme,
    col.frame = set1[1],
    # fontcolor.group = set1[2],
    # fontcolor.title = set1[2],
    size = 1,
    col = set1[1],
    col.line = set1[1],
    fill = set1[1],
    fontface.group = 3,
    fontface.title = "italic",
    transcriptAnnotation = "name",
    cex.title = 0.5,
    rotation.title = 0
))

aln_scheme <- as.list(c(
    my_scheme,
    size = 2,
    col = NA,
    col.line = NA,
    alpha = 1,
    alpha.reads = 1,
    max.height = 1000,
    min.height = 1000,
    stackHeight = 1,
    coverageHeight = 0,
    minCoverageHeight = 0,
    rotation.title = 0,
    fontface.title = "italic",
    # fontcolor.title = set1[3],
    type = c("pileup")))

########
# MAIN #
########

# prepare the txdb
if(!file.exists(txdb_file)) {
    # tidy up gff3
    gff_dt <- fread(cmd = paste('grep -v "^#"', gff_file))
    cruftlines <- gff_dt[, .I[V3 == "exon" & V8 != "."]]
    gff_lesscruft <- unique(gff_dt[!cruftlines])
    
    # rename Dnmt3
    gff_lesscruft[grepl("Vvulg11g01820\\.t1", V9),
                  V9 := gsub("Vvulg11g01820\\.t1",
                             "Vvulg11g01820.t1 (Dnmt3)", V9),
                  by = V9]
    
    # write out cleaned gff
    tmpgff <- tempfile(fileext = ".gff3")
    fwrite(gff_lesscruft, tmpgff, sep = "\t", col.names = FALSE)
    
    # make the database
    txdb <- GenomicFeatures::makeTxDbFromGFF(tmpgff, format = "gff3")
    AnnotationDbi::saveDb(txdb, txdb_file)
} else {
    txdb <- AnnotationDbi::loadDb(txdb_file)
}

# set up transcripts, exons and coding sequences
tx <- transcriptsBy(txdb, "gene")
exons <- exonsBy(txdb, "gene")
cds <- cdsBy(txdb, "gene")

# find dmnt3 in the txdb
dnmt3_transcripts <- tx$Vvulg11g01820
dnmt3_exons <- exons$Vvulg11g01820
dnmt3_cds <- cds$Vvulg11g01820
dnmt3_cds$cds_name <- "Dnmt3"
dnmt3_cds$id <- dnmt3_cds$cds_name

# define ranges to plot
extra_width <- 10e3
chr <- as.character(seqnames(dnmt3_transcripts))[[1]]
start <- min(start(dnmt3_transcripts)) - extra_width
end <- max(end(dnmt3_transcripts)) + extra_width

# the genome
gat <- GenomeAxisTrack(name = chr)
displayPars(gat) <- gat_scheme

# sequence track for reference
fa <- Biostrings::readDNAStringSet(fasta_file)
st <- SequenceTrack(fa, chromosome = chr)

# gene annotation track
my_grt <- GeneRegionTrack(txdb,
                          chromosome = chr,
                          start = start,
                          end = end,
                          name = "V. vulgaris")
displayPars(my_grt) <- grt_scheme

# alignment tracks
ats <- lapply(names(bamfiles), function(x)
    MakeAt(x,
           bamfile = bamfiles[[x]],
           st = st,
           chromosome = chr,
           start = start,
           end = end,
           aln_scheme = aln_scheme))

full_tracklist <- copy(ats)
full_tracklist[[length(full_tracklist) + 1]] <- my_grt
full_tracklist[[length(full_tracklist) + 1]] <- gat

# join tracks and highlight hypervariable region
ht1 <- HighlightTrack(trackList = full_tracklist,
                      start = min(start(dnmt3_transcripts)),
                      end = max(end(dnmt3_transcripts)),
                      chromosome = chr,
                      col = set1[4],
                      fill = set1[4],
                      alpha = 0.5)


wo <- grid::convertUnit(unit(as.integer(snakemake@params[["width"]]), "pt"),
                        "in", valueOnly = TRUE)
ho <- grid::convertUnit(unit(as.integer(snakemake@params[["height"]]), "pt"),
                        "in", valueOnly = TRUE)

cairo_pdf(snakemake@output[["plot"]],
          width = wo,
          height = ho,
          pointsize = 8)

plotTracks(ht1,
           # cex.main = 2,
           # fontface.main = 1,
           add = FALSE,
           from = start,
           to = end,
           title.width = 2.2)

dev.off()

sessionInfo()
