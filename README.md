# `cligb`: Command line based genome browser for *-seq data
## Description
_This is currently under development..._

* Quick and simple way to visualize a genomic region of your interest accoss many sequence samples
* Query by gene symbol, region[start:end]

 
## Usage:
```
$ R --vanilla --slave -f gb2.R --args hsa-mir-XXX PATH_TO_BAM_DIR/
```

![](https://dl.dropboxusercontent.com/u/8677629/gb.png)

## TODO:
- [] Packaging for Python command line interface
	* `$ cligb --start 1 --end 999 --bams ... --gtf gene.gtf --plot gb.pdf`
	* `$ cligb --query foxp2 --bams ... --gtf gene.gtf --plot gb.pdf`
- [] Support more query types
- [] Error handling
- [] Support to create AnnotationTrack from data is not provided as  TxDb.* in Bioconductor packages, for example, miRBase or user-generated GTF/GFF/TSV format.
- [] More faster and faster to use cache mechanism

