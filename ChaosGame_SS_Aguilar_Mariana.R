
#Modificaciones al script original de Roberto. Le añadí que en los plots ponga de título
#el nombre que viene  en la secuencia fasta y que además guarde el plot en un PDF.
genoma<-readDNAStringSet("Wuhan_seafood_virus.fasta")

x<-c()
y<-c()

omega_mat<-function(v1,v2){
  
  mat<-matrix(c(0.5,0,0,0.5),nrow=2)
  v1<-matrix(v1,nrow=2)
  v2<-matrix(v2,nrow=2)
  return(mat %*% v1 + v2/2)
  
}

chaos_game<-function(genoma){
  x<-c()
  y<-c()
  x[1]<-runif(1)
  y[1]<-runif(1)
  a<-matrix(c(0,0),nrow = 2)
  c<-matrix(c(0,1),nrow=2)
  g<-matrix(c(1,1),nrow=2)
  t<-matrix(c(1,0),nrow=2)
  
  #for(i in 1:length(genome)){
  for(j in 1:length(genoma[[1]])){
    if (genoma[[1]][j]==DNAString("A")){
      x[j+1]<-omega_mat(matrix(c(x[j],y[j]),nrow=2),a)[1,]
      y[j+1]<-omega_mat(matrix(c(x[j],y[j]),nrow=2),a)[2,]
    }else if(genoma[[1]][j]==DNAString("C")){
      x[j+1]<-omega_mat(matrix(c(x[j],y[j]),nrow=2),c)[1,]
      y[j+1]<-omega_mat(matrix(c(x[j],y[j]),nrow=2),c)[2,]
    }else if(genoma[[1]][j]==DNAString("G")){
      x[j+1]<-omega_mat(matrix(c(x[j],y[j]),nrow=2),g)[1,]
      y[j+1]<-omega_mat(matrix(c(x[j],y[j]),nrow=2),g)[2,]
    }else if (genoma[[1]][j]==DNAString("T")){ 
      x[j+1]<-omega_mat(matrix(c(x[j],y[j]),nrow=2),t)[1,]
      y[j+1]<-omega_mat(matrix(c(x[j],y[j]),nrow=2),t)[2,]
    }
    #print(j)
    #plot(x,y)
  }
  plot(x,y,cex=0.1,col="orange",main=names(genoma))
  
  pdf(file="ChaosGame_plot.pdf") # Hay que corregirlo cada que se use una secuencia.
  plot(x,y,cex=0.1,col="orange",main=names(genoma))
  dev.off()
  # }
  
}

chaos_game(genoma)
 
##genomas eucariotas

BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene") #ratón
#Exposes an annotation databases generated from UCSC by exposing these as TxDb objects
#The TxDb class is a container for storing transcript annotations
#https://bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene") #humano

BiocManager::install("GenomicFeatures")
# set of tools and methods for making and manipulating transcript centric annotations.
#With these tools the user can easily download the genomic locations of the transcripts, exons and cds of a given organism, from either the UCSC Genome Browser or a BioMart database 

library("GenomicFeatures")
library("TxDb.Mmusculus.UCSC.mm10.knownGene") 

                ##### objetos TxDb #######

#loading Transcript Data
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene #shorthand (for convenience)
txdb

#Pre-filtering data based on Chromosomes
head(seqlevels(txdb))
seqlevels(txdb) #To determine which chromosomes are currently active
#Will tell you all the chromosomes that are active for the TxDb.Hsapiens.UCSC.hg19.knownGeneTxDb object

seqlevels(txdb) <- "chrX" #So if you ran this, then from this point on in your R session only chromosome 1 would
#be consulted when you call the various retrieval methods.

