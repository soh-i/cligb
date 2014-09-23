library(Gviz)
library(Rsamtools)
library(dplyr)


Cligb <- setRefClass(
    Class = "Cligb",

    fields = list(
        genome.ver = "character",
        chromosome = "character",
        start = "numeric",
        end = "numeric",
        interest = "data.frame"
    ),

    methods = list(
        initialize = function(genome.ver=NULL) {
            if (is.null(genome.ver)) {
                warning("Genome version is not specified, set default[hg19]")
                genome.ver <<- "hg19"
            } else {
                genome.ver <<- genome.ver
            }
        },
        
        generate.db = function(file=NULL) {
            if (!is.null(file)){
                df <- read.table(file, sep="\t", header=TRUE)
                return(df)
            }
        },

        query.by.symbol = function(db=NULL, ...) {
            if (is.null(db)) {
                stop("DB is not specified!")
            }
            args = list(...)
            result <- data.frame()
            for (i in args){
                result <- rbind(result, filter(db, symbol == i) %>% select(chromosome, start, end, strand, symbol, transcript))
            }
            result
        },

        find.mirs = function(db=NULL, query=NULL){
            three <- gsub("mir", "miR", query) %>% paste("-3p", sep="")
            five <- gsub("mir", "miR", query) %>% paste("-5p", sep="")
            queries <- query.by.symbol(db=db, query, three, five)
            if (nrow(queries) == 0) {
                stop(sprintf("Query error: Can not find \"%s\" in db", query))
            }
            .set.interest(queries)
            interest <<- queries
        },

        .set.interest = function(q) {
            chromosome <<- as.character(q$chromosome)
            start <<- min(q$start)
            end <<- max(q$end)
        },

        create.ideogram.track = function(...) {
            ig.track <- IdeogramTrack(genome=genome.ver, chromosome=chromosome, showBandId=TRUE, bevel=1,
                                      col.axis="black", fontsize=15, fontcolor="black", ...)
        },

        create.annotation.track = function(...){
            anno.track <- AnnotationTrack(range=interest, name="miRBase", group=interest$symbol,
                                           showId=TRUE,chromosome=chromosome, genome=genome.ver,
                                           fill="orange", col="orange", shape="smallArrow",
                                           col.axis="black", background.title=NA, col.title="black",
                                           fontcolor.group="black", fontsize=15)
        },
        
        create.seq.track = function(...) {
        },

        create.axis.track = function(...) {
        },

        create.data.track = function(file.dir, ...) {
        },

        plot = function(...) {
            plotTracks(list(...), from=start+10, to=end+10)
        }
    )
)

cligb <- Cligb$new(genome.ver="hg19")
mir.db <- cligb$generate.db(file="~/Dropbox/Desktop.osx/ADAR_paper/data/miRBase.Gviz")
cligb$find.mirs(db=mir.db, "hsa-mir-22")
ann.track <- cligb$create.annotation.track()
print(ann.track)
#id.track <- cligb$create.ideogram.track()
cligb$plot(ann.track)

#cligb$plot(id.track)




stop("\n---OK---\n")




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

