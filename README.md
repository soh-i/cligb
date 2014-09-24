# `cligb`: Command line based genome browser for *-seq data
## Description
_This is currently under development..._

* Quick and simple way to visualize a genomic region of your interest accoss many sequence samples
* Query by gene symbol

 
## Demo
* Downloadiing a small RNA-seq data sequenced by ENCODE project from `http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlShortRnaSeq/releaseLatest/wgEncodeCshlShortRnaSeqGm12878CellShortAln.bam` (53MB).
* Indexing the small RNA-seq data is following command,
	* `samtools index wgEncodeCshlShortRnaSeqGm12878CellShortAln.bam`
	* move it to `PATH_TO_BAM_DIR`
* Plotting alignmenta data: 

	```
	$ R --vanilla --slave -f gb2.R --args hsa-mir-99a PATH_TO_BAM_DIR/
	```

![](https://dl.dropboxusercontent.com/u/8677629/gb.png)

## TODO:
- [ ] Packaging for Python command line interface
	* `$ cligb --start 1 --end 999 --bams ... --gtf gene.gtf --plot gb.pdf`
	* `$ cligb --query foxp2 --bams ... --gtf gene.gtf --plot gb.pdf`
- [ ] Support more query types
- [ ] Error handling
- [ ] Support to create AnnotationTrack from data is not provided as  TxDb.* in Bioconductor packages, for example, miRBase or user-generated GTF/GFF/TSV format.
- [ ] More faster and faster to use cache mechanism

