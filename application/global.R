
library(shinyjs)
library(shinydashboard)
library(ggplot2)
library(DT)
library(scales)
library(gdata)
library(Biostrings)
library(deconstructSigs)
library(VennDiagram)
library(data.table)
source("loadData.R")

contextColors <- c("#56B4E9", "coral4", "red1", "grey", "olivedrab2", "pink")

mutOrder <- data.frame("type" = c('A>G:T>C','C>T:G>A','C>A:G>T','C>G:G>C','A>T:T>A','A>C:T>G'), 'n' = NA)

contextLabel <- c()
for (t in c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")){
  for (b1 in c("A","C","G","T")){
    for (b2 in c("A","C","G","T")){
      contextLabel <- c(contextLabel ,paste(b1,"[",t,"]",b2,sep=""))
    }
  }
}
contextLabel2 <- c()
for (t in c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")){
  for (b1 in c("A","C","G","T")){
    for (b2 in c("A","C","G","T")){
        dna <-  paste(b1,"[",t,"]",b2,sep="")
        dna1 <- reverseComplement(DNAString(substr(dna,1,1)))
        dna2 <- reverseComplement(DNAString(substr(dna,3,3)))
        dna3 <- reverseComplement(DNAString(substr(dna,5,5)))
        dna4 <- reverseComplement(DNAString(substr(dna,7,7)))
        contextLabel2 <- c(contextLabel2, paste(dna,":",dna1,"[",dna2,">",dna3,"]",dna4, sep=""))
    }
  }
}
contextOrder <- data.frame("type" = contextLabel2)
contextOrder[0:16,"b"] <- c('C>A')
contextOrder[17:32,"b"] <- c('C>G')
contextOrder[33:48,"b"] <- c('C>T')
contextOrder[49:64,"b"] <- c('T>A')
contextOrder[65:80,"b"] <- c('T>C')
contextOrder[81:96,"b"] <- c('T>G')
contextOrder[ , "n"] <- NA
contextOrder


consequence_labels <- c("intergenic_variant",
                        "upstream_gene_variant",
                        "5_prime_UTR_variant",
                        "start_lost",
                        "missense_variant",
                        "synonymous_variant",
                        "stop_gained",
                        "splice_acceptor_variant",
                        "splice_donor_variant",
                        "splice_region_variant",
                        "intron_variant",
                        "stop_lost",
                        "stop_retained_variant",
                        "3_prime_UTR_variant", 
                        "downstream_gene_variant",
                        "NMD_transcript_variant", 
                        "non_coding_transcript_variant", 
                        "non_coding_transcript_exon_variant"
                        )
consOrder <- data.frame("type" = consequence_labels, "n"=NA)



