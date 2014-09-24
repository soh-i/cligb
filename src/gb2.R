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
        initialize = function(genome.ver = NULL) {
            if (is.null(genome.ver)) {
                warning("Genome version is not specified, set default[hg19]")
                genome.ver <<- "hg19"
            } else {
                genome.ver <<- genome.ver
            }
        },
        
        generateDb = function(file = NULL) {
            if (!is.null(file)){
                df <- read.table(file, sep = "\t", header = TRUE)
                return(df)
            }
        },

        queryBySymbol = function(db = NULL, ...) {
            if (is.null(db)) {
                stop("DB is not specified!")
            }
            args = list(...)
            result <- data.frame()
            for (i in args){
                result <- rbind(result, filter(db, symbol == i) %>%
                                    select(chromosome, start, end, strand, symbol, transcript))
            }
            result
        },

        setInterests = function(db = NULL, query = NULL){
            if (!is.null(db) & !is.null(query)) {
                three <- gsub("mir", "miR", query) %>% paste("-3p", sep = "")
                five <- gsub("mir", "miR", query) %>% paste("-5p", sep = "")
                queries <- queryBySymbol(db = db, query, three, five)
            } else {
                stop("Error: lack a db and query")
            }
            
            if (nrow(queries) == 0) {
                stop(sprintf("Query error: Can not find \"%s\" in db", query))
            } else {
                interest <<- queries
                chromosome <<- as.character(queries$chromosome)
                start <<- min(queries$start)
                end <<- max(queries$end)
            }
        },

        createIdeogramTrack = function(...) {
            ig.track <- IdeogramTrack(genome = genome.ver,
                                      chromosome = chromosome,
                                      showBandId = TRUE,
                                      bevel = 1,
                                      col.axis = "black",
                                      fontsize = 15,
                                      fontcolor = "black",
                                      ...)
        },

        createAnnotationTrack = function(...){
            anno.track <- AnnotationTrack(range = interest,
                                          name = "miRBase",
                                          group = interest$symbol,
                                          showId = TRUE,
                                          chromosome = chromosome,
                                          genome = genome.ver,
                                          fill = "orange",
                                          col = "orange",
                                          shape = "smallArrow",
                                          col.axis = "black",
                                          background.title = NA,
                                          col.title = "black",
                                          fontcolor.group = "black",
                                          fontsize = 15,
                                          ...)
        },
        
        createSeqTrack = function(...) {
            library(BSgenome.Hsapiens.UCSC.hg19)
            seq.track <- SequenceTrack(Hsapiens,
                                       genome = genome.ver,
                                       chromosme = chromosome,
                                       ...)
        },

        createAxisTrack = function(...) {
            axis.track <- GenomeAxisTrack(genome = genome.ver,
                                          chromosome = chromosome,
                                          showId = TRUE,
                                          add53 = TRUE,
                                          exponent = 0,
                                          labelPos = "above",
                                          fontcolor = "black",
                                          ...)
        },

        createDataTrack = function(bam.dir, ...) {
            bam.files <- list.files(bam.dir, "*.bam$")
            
            track.name <- gsub(pattern = ".bam", replacement = "_track", x = bam.files)
            c <- 0
            for (file in bam.files){
                c <- c + 1
                assign(x = track.name[c], value = DataTrack(range = paste(bam.dir, bam.files[c], sep = ""),
                                              name = file,
                                              genome = genome.ver,
                                              chromosome = chromosome,
                                              start = start,
                                              end = end,
                                              type = "histogram",
                                              col.histogram = NA,
                                              fill.histogram = "lightblue",
                                              col.axis = "darkgray",
                                              background.title = NA,
                                              col.title = "darkgray",
                                              ylim = c(0, 5),
                                              fontsize = 15,
                                              ...)
                       )
            }
            track.list <- list()
            for (i in 1:length(bam.files)) {
                track.list[[i]] <- get(track.name[i])
            }
            return(track.list)
        },
        
        plot = function(...) {
            pdf("demo.pdf", height=12, width=10)
            plotTracks(...,
                       from = start - 10,
                       to = end + 10,
                       transformation = function(x){return(log10(x+1))} )
            dev.off()
        }
    )
)


args <- commandArgs(trailingOnly = TRUE)
if (is.na(args[1])) {
  stop("Error: no queries is found")
} else if (is.na(args[2])) {
    stop("Error: no bam directory is found")
}
q <- args[1]
dir <- args[2]

cligb <- Cligb$new(genome.ver = "hg19")

message("Generate annotation track from external data source")
mir.db <- cligb$generateDb(file = "extdata/miRBase.Gviz")
cligb$setInterests(db=mir.db, query=q)

message("Create each tracks...")
ann.track <- cligb$createAnnotationTrack()
id.track <- cligb$createIdeogramTrack()
ax.track <- cligb$createAxisTrack()
data.track <- cligb$createDataTrack(dir)

message("Plotting all tracks...")
all.track <- append(c(id.track, ax.track, ann.track),  data.track)
cligb$plot(all.track)

