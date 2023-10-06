library(Gviz)
library(dplyr)
library(stringr)
library(GenomicFeatures)
library(tidyr)


offs=50000
data <- read.table("data/list2.txt")
colnames(data) <- c("id","chr","start","stop", "strand")
data2 <- read.table("data/transposons_microRNA.txt",
		    sep="\t")

data2 <- data2[,(c(1,2,4,5,6,7))]


data2[,3] <- gsub('C', '-', data2[,3] )

colnames(data2) <- c("id","transposon","strand","chromosome","start","stop")
data2$chromosome <- sub("^", "chr", data2$chromosome )
axisTrack <- GenomeAxisTrack(labelPos="above")

tmp_dict <- read.table("data/transcript_symbols.txt", sep='\t')


for(i in 1:dim(data)[1]){
    ident <- data$id[i]
    chrom <- data$chr[i]
    from <- data$start[i]-offs
    to <- data$stop[i]+offs
    from_orig <- data$start[i]
    to_orig <- data$stop[i]
    strand <- data$strand[i]

    idTrack <- IdeogramTrack(chromosome = chrom, genome = "hg38")
    if (sum(str_detect(strand, "/")) > 0){
    	
	mirTrack <- AnnotationTrack(chrom=c(chrom,chrom),start=c(from_orig-1000, from_orig-1000),end=c(to_orig+1000, to_orig+1000),strand=c("+", "-"),fill="#6E6E14", name="MicroRNAs", background.title="darkgray")
	}
    

    else{
    	mirTrack <- AnnotationTrack(chrom=chrom,start=from_orig-1000,end=to_orig+1000,strand=strand,fill="#6E6E14", name="MicroRNAs", background.title="darkgray")
 	}
    tx = makeTxDbFromGFF("data/short.gff3")
    geneTrack <- GeneRegionTrack(tx, chromosome=chrom, from=from, to=to, name="Genes", min.height=8, background.title="darkgray")
    
    z <- ranges(geneTrack)
   
    z$symbol <- tmp_dict$V2[match(z$symbol, tmp_dict$V1)]

    ranges(geneTrack) <- z

	
    tracklist <- list(idTrack, axisTrack, mirTrack)
    tracklist2 <- list(geneTrack)
    sizes <- list(1,1)
    sizes2 <-list(2)
    if ((any(data2$id==ident))==TRUE){
	
	subdata <- data2[data2$id==ident,]
    	if (sum(str_detect(subdata[,2], "LINE")) > 0){

		subdata_LINE <- subdata[str_detect(subdata[,2], "LINE"),]

		gr_LINE <- GRanges(
		      seqnames = subdata_LINE[,4], range = IRanges(start = subdata_LINE[,5], end = subdata_LINE[,6]),
    		      strand = subdata_LINE[,3]
		      )
		transpLINE <- AnnotationTrack(range=gr_LINE,
					     name="LINEs",
                                             stacking="squish",
					     background.panel = "#E3E5C8",
					     background.title="#8D9064")
		tracklist2 <- append(tracklist2, transpLINE)
		sizes2 <- append(sizes2, 1)

	    }
	if (sum(str_detect(subdata[,2], "SINE")) > 0){
                subdata_SINE <- subdata[str_detect(subdata[,2], "SINE"),]
		gr_SINE <- GRanges(
        	       seqnames = subdata_SINE[,4], range = IRanges(start = subdata_SINE[,5], end = subdata_SINE[,6]),
        	       strand = subdata_SINE[,3]
                     )
		transpSINE <- AnnotationTrack(range=gr_SINE,
					     name="SINEs",
                                             stacking="squish",
					     background.panel = "#E3E5C8",
					     background.title="#8D9064")
		tracklist2 <- append(tracklist2, transpSINE)
		sizes2 <- append(sizes2, 1)
	}
	if (sum(str_detect(subdata[,2], "LTR")) > 0){
                subdata_DNA <- subdata[str_detect(subdata[,2], "LTR"),]
		gr_DNA <- GRanges(
        	       seqnames = subdata_DNA[,4], range = IRanges(start = subdata_DNA[,5], end = subdata_DNA[,6]),
        	       strand = subdata_DNA[,3]
                     )
		transpDNA <- AnnotationTrack(range=gr_DNA,
					    name="LTR",
        	                            stacking="squish",
					    background.panel = "#E3E5C8",
					    background.title="#8D9064")
		tracklist2 <- append(tracklist2, transpDNA)
		sizes2 <- append(sizes2, 1)

       }
	
    }

    highLightTrack <- HighlightTrack(trackList = tracklist2,
                    chromosome = chrom, start = from_orig, end = to_orig,
		    inBackground=FALSE, fill=NA)
    tracklist <- append(tracklist, highLightTrack, after=3)
    pdf(paste0("results/genome_figs/",ident,"_knowngff3_.pdf"), width=10, height=7)
    plotTracks(tracklist,
	       from = from, to = to, 
       	       transcriptAnnotation = "symbol",
               collapseTranscripts = "longest")
    dev.off()
}
