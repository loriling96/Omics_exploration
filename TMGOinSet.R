####make a loop to find intersection of Transmembrane GO terms in given hit list

library(argparse)
description = "this is my script to identify proteins with annotated transmembrane domains"
parser= ArgumentParser(description=description)
parser$add_argument("-i", "--infile",required=TRUE, help="full path to csv file of dataset")
parser$add_argument("-o", "--outfile", help="name of csv file of output dataframe containing annotations")
parser$add_argument("-rp", "--resultspath", help= "specify path to folder which will contain results") 
args= parser$parse_args()

infile=args$infile
outfile=args$outfile
resultspath = as.character(args$resultspath)

##read in subject data set
##updated GOterm list##
GOnew=read.delim("~/Documents/Pringlelab/Expt2.05/GO_transmembrane_v2.txt", header = TRUE)
TMGO=as.character(GOnew$GOterm)
print("Subject data set is GO_transmembrane_v2.txt")

#GO_transmembrane = read.delim("~/Documents/Pringlelab/Expt2.05/GO_transmembrane.txt", header=TRUE)
#TMGO=as.character(GO_transmembrane$GOterm)
PNAS= read.csv(file ="~/Documents/Pringlelab/Expt2.05/Aiptasia_annot_PNAS.csv", header = TRUE )


##read in query aipgene list
DF=read.csv(infile, header=TRUE, stringsAsFactors=FALSE)

##merge GO annotations
DF_GO= merge(x = DF, y = PNAS[ , c(1:4,12)], by = "Gene.identifier", all.x=TRUE)


TMinQ=vector(mode= "character", length=0)
TMgenes=vector(mode = "character", length=0)
Query=DF_GO$GO.terms #keep in mind that not every aipgene has annotation hit

count=0

for (i in 1:length(Query)) {
  Query_GO = unlist(strsplit(as.character(Query[i]), ","))
  x = intersect(Query_GO, TMGO)
  
  if (length(x) == 0){
    TMinQ=append(TMinQ,"NA")
  } else {
    TMinQ=append(TMinQ, x, after = i)
    count = count + 1
    TMgenes=append(TMgenes, DF_GO$Gene.identifier[i], after = i)
      }
}
print(paste("For your query", infile))
print(paste("number of transmembrane GOterm in query is ", count))
print(paste("The percent of transmembrane GO annotations is", count/length(Query)))

DF2=subset(DF_GO, DF_GO$Gene.identifier %in% TMgenes)

write.csv(DF2, file=paste(resultspath, outfile, sep = '/'), row.names = FALSE)


