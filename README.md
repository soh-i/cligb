# `cligb`: Command line based genome browser for *-seq data
## Description
_This is currently under development..._

* Quick and simple way to visualize a genomic region of your interest accoss many sequence samples
* Query by gene symbol

## Requirements
* Gviz (Bioconductor package)
 
## Demo
* Clone this repo.
* Downloading a GM12878 (Whole cell/Nucleus) and K562 cell line small RNA-seq data sequenced by ENCODE project from 
 	* `http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlShortRnaSeq/releaseLatest/wgEncodeCshlShortRnaSeqGm12878CellShortAln.bam` (53MB)
 	* `http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlShortRnaSeq/releaseLatest/wgEncodeCshlShortRnaSeqGm12878NucleusShortAln.bam` (52MB)
 	* `http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlShortRnaSeq/releaseLatest/wgEncodeCshlShortRnaSeqK562CellShortAln.bam` (46MB)
* Indexing the seq. data is the following command,
	* `samtools index bam file`
	* move it to `PATH_TO_BAM_DIR`
* Plotting alignment data: 

	```
	$ R --vanilla --slave -f gb2.R --args hsa-mir-99a PATH_TO_BAM_DIR/
	```

![](https://dl.dropboxusercontent.com/u/8677629/gb-demo.png)

## TODO:
- [ ] Packaging for Python command line interface
	* `$ cligb --start 1 --end 999 --bams ... --gtf gene.gtf --plot gb.pdf`
	* `$ cligb --query foxp2 --bams ... --gtf gene.gtf --plot gb.pdf`
- [ ] Support more query types
- [ ] Error handling
- [ ] Support to create AnnotationTrack from data is not provided as  TxDb.* in Bioconductor packages, for example, miRBase or user-generated GTF/GFF/TSV format.
- [ ] More faster and faster to use cache mechanism

