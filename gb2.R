library(Gviz)
library(Rsamtools)
library(dplyr)

query.by.symbol <- function(db=NULL, ...) {
  #if (is.na(db)){
  #  stop("Error: need DB for building query")
  #}
  args = list(...)
  result <- data.frame()
  for (i in args){
    result <- rbind(result, filter(db, symbol == i) %>% select(chromosome, start, end, symbol, transcript))
  }
  return(result)
}

generate.query <- function(db=NULL, query=NULL){
  three <- gsub("mir", "miR", query) %>% paste("-3p", sep="")
  five <- gsub("mir", "miR", query) %>% paste("-5p", sep="")
  queries <- query.by.symbol(db=db, query, three, five)
  if (nrow(queries) == 0) {
    stop(sprintf("QueryError: Can not find \"%s\" in db", query))
  }
  return(queries)
}

generate.bam.tracks <- function(bam.dir, gen, chrm, start, end) {
  bam.files <- list.files(bam.dir, "*sorted.bam$")
  track.name <- gsub(pattern=".bam", replacement="_track", x=bam.files)
  c <- 0
  for (file in bam.files){
    c <- c + 1
    assign(x=track.name[c], value=DataTrack(range=paste(bam.dir, bam.files[c], sep=""),
                                            name=file,
                                            genome=gen,
                                            chromosome=chrm,
                                            start=start, end=end,
                                            type="histogram",
                                            col.histogram=NA,
                                            fill.histogram="lightblue",
                                            col.axis="darkgray",
                                            background.title=NA,
                                            col.title="darkgray",
                                            ylim=c(0, 5),
                                            fontsize=15))
  }
  track.list <- list()
  for (i in 1:length(bam.files)) {
    track.list[[i]] <- get(track.name[i])
  }
  return(track.list)
}

mir.db <- read.table("data/miRBase.Gviz", header=T, sep="\t") # TODO: create prepare() func for annotation data
args <- commandArgs(trailingOnly=TRUE)
#args[1] <- "hsa-mir-22"
if (length(args) == 0) {
  stop("Error: need for miRNA name!")
}
args.mir <- args[1]
interest <- generate.query(db=mir.db, query=args.mir)
chrm <- as.character(interest$chromosome)
start <- min(interest$start)
end <- max(interest$end)
gen <- "hg19"

cat("Create Ideogram track...\n")
ig.track <- IdeogramTrack(genome=gen, chromosome=chrm, showBandId=TRUE, bevel=1,
                          col.axis="black", fontsize=15, fontcolor="black")

cat("Create Annotation track...\n")
mir.track <- AnnotationTrack(range=interest, name="miRBase", group=interest$symbol,
                             showId=TRUE,chromosome=chrm, genome=gen,
                             fill="orange", col="orange", shape="smallArrow",
                             col.axis="black", background.title=NA, col.title="black",
                             fontcolor.group="black", fontsize=15)

cat("Create bam tracks...\n")
bam.tracks <- generate.bam.tracks("data/seqs/bams", gen, chrm, start, end) # TODO: parametalized 
all.tracks <- append(c(ig.track, mir.track), bam.tracks)

cat("Draw all tracks...\n")
pdf(file=sprintf("gb.%s.pdf", args.mir), width=10, height=10)

plotTracks(trackList=all.tracks, showId=TRUE, chromosome=chrm, from=start-10, to=end+10,
           transformation=function(x){return(log10(x+1))})
dev.off()

