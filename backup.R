library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("Gviz")
#Download the reference genome from refseq


library(Biostrings)
genome<-readDNAStringSet("/Users/ipaz00/Downloads/BIOL4315_R/biol 4315 lab3/BIOL-4315-lab-3-/GCF_000146045.2_R64_genomic.fna")
names(genome)
xome1<- genome[grep("chromosome I",names(genome), ignore.case=TRUE)]
writeXStringSet(xome1, "R64_chr1.fna")
xome1

library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)
library(Gviz)
# BAM file path
bam_file <- "alignment.sorted.bam"
# Open connection to BAM file
bam <- BamFile(bam_file)

# Read all alignments (already restricted to chr1)
alns <- readGAlignments(bam)
options(ucscChromosomeNames = FALSE)

# Convert to GRanges object
gr_alns <- granges(alns)
aln_track <- AnnotationTrack(gr_alns, name = "Contigs", genome = "sacCer3", chromosome = "NC_001133.9")


# Add genome axis for scale
axis_track <- GenomeAxisTrack()
# Plot
plotTracks(list(axis_track, aln_track),
           from = min(start(gr_alns)),
           to = max(end(gr_alns)),
           main = "Contig Alignments to Chromosome 1")

fastaroni<-readDNAStringSet("/Users/ipaz00/Downloads/BIOL4315_R/biol 4315 lab3/BIOL-4315-lab-3-/assembly.fasta")
length<-width(fastaroni)
fileinclusion<-fastaroni[which.max(length)]
writeXStringSet(fileinclusion,"longest_contig.fasta")


library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)
library(Gviz)
# BAM file path
bam_file <- "cds_aln.sorted.bam"
# Open connection to BAM file
bam <- BamFile(bam_file)

# Read all alignments (already restricted to chr1)
alns <- readGAlignments(bam)
options(ucscChromosomeNames = FALSE)

# Convert to GRanges object
gr_alns <- granges(alns)
aln_track <- AnnotationTrack(gr_alns, name = "Contigs", genome = "sacCer3", chromosome = "contig_26")

# Add genome axis for scale
axis_track <- GenomeAxisTrack()
# Plot
plotTracks(list(axis_track, aln_track),
           from = min(start(gr_alns)),
           to = max(end(gr_alns)),
           main = "CDS alignments to Contig 26 ")