
def to_Gviz(f, label=None):
    print "chromosome\tstart\tend\twidth\tstrand\tfeature\tgene\texon\ttranscript\tsymbol"

    with open(f, "r") as f:
        for line in f:
            data = line.strip().split("\t")
            chrom = data[0]
            start = int(data[3])
            end = int(data[4])
            width = abs(end - start)
            strand = data[6]
            features = data[8].split(";")
            feature = label
            gene = features[3].replace(" Derives_from ", "").replace("\"", "")
            exon = features[0].replace("transcript_id ", "").replace("\"", "")
            transcript = exon
            symbol = features[2].replace(" Name ", "").replace("\"", "")
        
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start, end, width, strand, feature, gene, exon, transcript, symbol)
        

if __name__ == '__main__':
    mature = "mature.gtf" # Download from mirbase.org
    pre = "pre.gtf"

    to_Gviz(pre, label="pre-miRNA")
