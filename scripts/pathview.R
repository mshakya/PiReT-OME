library(pathview) 
options(bitmapType="cairo")

map_data <- read.delim(file = "pathway.stats.txt", stringsAsFactors = FALSE,  header=TRUE)

unique_map<-unique(map_data$Map_ID)
unique_sample<-map_data$GeneID[2]
sample_split<- strsplit(as.character(unique_sample), split="::")
geneid=sapply(sample_split, "[[", 2)
sampleid=sapply(sample_split, "[[", 1)

outsuff=paste (sep = "", sampleid, ".pathview")

data(korg)

for (i in 1:length(unique_map)) 
{
map_split<- strsplit(as.character(unique_map[i]), split=":")
pathwayid=sapply(map_split, "[[", 2)
speciesid=sapply(map_split, "[[", 1)
spidc = match(speciesid, korg[, 1:3])%%nrow(korg)
spnai = is.na(spidc)
LOCUS_IDpos=which(names(map_data)=="LOCUS_ID")
unique_row_set <- subset(map_data, Map_ID == unique_map[i],   select = c((LOCUS_IDpos+1):length(map_data[1,])))
limitscale=max(abs(min(unique_row_set)),  max(unique_row_set))
    if (sum(spnai) == 0 && length(unique_row_set[,1]) >1 )
  {


pathview_map_set <- subset(map_data,  select = c(LOCUS_IDpos:length(map_data[1,])))
de_matrix<-as.matrix(pathview_map_set[,2:length(pathview_map_set[1,])])
rownames(de_matrix)=pathview_map_set[,1]


sample_size<-length(de_matrix[1,])
pv.out<-pathview(new.signature = FALSE, out.suffix = outsuff,limit= list(gene = limitscale, cpd = limitscale), bins = list(gene = 40, cpd = 40), gene.data=de_matrix[,1:sample_size], gene.idtype="kegg", pathway.id=pathwayid, species=speciesid, kegg.native = T, same.layer=T)
  }
}
#unique_map_set <- subset(map_data, Map_ID == unique_map[i],   select = c(LOCUS_ID, reduced_EdgeR_logFC, reduced_Deseq_log2FoldChange))
