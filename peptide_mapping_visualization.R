
#Library loading

library(tidyverse)
library(reshape2)
library(scales)
library(splitstackshape)

#Loading peptides.txt from maxquant output
peptides<- as.data.frame(read.delim("peptides.txt", sep="\t",header=TRUE))
names(peptides)

#Select the following columns into new dataframe, "Sequence", "Gene.names", "Start.position", "End.position" 
#and the samples intensities (in this case will be columns between 56:61 )

select_col_peptides<- select(peptides, 1, 39,37,38,56:61)
names(select_col_peptides)

#Clean Up and Melting the dataframe

select_col_peptides[select_col_peptides == 0] <- NA
select_col_peptides[, 5:10] <- log(select_col_peptides[,5:10], 2)
select_col_peptides_melt<- melt(select_col_peptides,id.vars = c("Sequence","Gene.names", "Start.position", "End.position"  ), value.name = "Log2_Intesnsity", variable.name = "Samples" )



#Calculate peptide lenght
select_col_peptides_melt$peptide_length <- (select_col_peptides_melt$End.position -select_col_peptides_melt$Start.position)+1

#remove peptides witgh NA intensities
select_col_peptides_melt<- subset(select_col_peptides_melt, Log2_Intesnsity>0  )

#Since we are only mapping true or false coverage, we remove redundant peptide
select_col_peptides_melt$duplication_mark<-paste(select_col_peptides_melt$Samples,select_col_peptides_melt$Gene.names,select_col_peptides_melt$Sequence)
select_col_peptides_melt<-select_col_peptides_melt[!duplicated(select_col_peptides_melt$duplication_mark), ]
#setting the bar height for visulization
select_col_peptides_melt$Bar_length = 1

#Finally breakdown the peptide into individula aa with corresponding corrdinates  

select_col_peptides_melt= select_col_peptides_melt[complete.cases(select_col_peptides_melt), ]

select_col_peptides_melt_aa = select_col_peptides_melt %>%
  group_by( Sequence,Gene.names,Samples, Start.position  ) %>%
  expandRows("peptide_length") %>%
  mutate(Amino.acid.position = Start.position - 1 + cumsum(Bar_length)) 


#Visuliza the protein coverage for any protein of interest (insert Gene name in Gene.names=="XX" , in the code below)



PABPN1= ggplot( subset(select_col_peptides_melt_aa, Gene.names=="PABPN1" ), aes(x=Amino.acid.position,y= Bar_length, fill=Samples)) +
  geom_bar(stat = "identity", position = "fill")+
  facet_wrap(~Samples,ncol = 1)+
  #below you can set the full lenght of the protein of interest
  xlim(1, 306)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_brewer(palette="Dark2")+
  theme_classic()+
  labs(title = "PABPN1")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        axis.line = element_line(colour = 'black', size = 0.25),
        axis.ticks = element_line(colour = 'black', size = 0.25),
        strip.background  = element_blank())+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black",size=0.25, linetype="solid")+
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, color="black",size=0.25, linetype="solid")+
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color="black",size=0.25, linetype="solid")
ggsave(filename="PABPN1.png", PABPN1, width=5, height=3, units="in")

