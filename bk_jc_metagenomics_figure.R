list.of.packages <- c("DT","tidyverse","pavian","plotly","dplyr","plyr","reshape2","randomcoloR")
lapply(list.of.packages,library,character.only = TRUE)

# Read in clomp_viz files and do some cleaning
setwd("/Users/uwvirongs/Documents/Michelle/bk/bk_clomp/clompviz_no_bad_klebsiella")
path<-getwd()
files<-list.files(path, pattern = '*.tsv',full.names = TRUE)
reports<-pavian::read_reports(files)
filenames<-strsplit(basename(files), '[.]')
reports_merged<-as.data.frame(pavian::merge_reports2(reports,fix_taxnames = TRUE, col_names = filenames))
keep<-which(grepl('clade',colnames(reports_merged)))
only_clade<-reports_merged[,c(1:4, keep)]
new_colnames<-c()
for(i in 5:length(colnames(only_clade))){ 
  #print(colnames(only_clade[i]))
  first<-strsplit(colnames(only_clade)[i],'[..]')[[1]][4]
  new_colnames<-append(new_colnames,first)
}
colnames(only_clade)[keep]<-new_colnames
colnames(only_clade)[1:4]<-c('name','rank','taxID','lineage')

cleaned_df <- only_clade

# Save totals for later for rpm calculations
total <- cleaned_df %>% filter(name=="root")
domains <- cleaned_df %>% filter(rank=="D")

# Grab classifications we want
bk_jc_df <- cleaned_df %>% filter(name == "Bacteria" | name == "Viruses" | name == "Fungi" |
                                    name == "Human polyomavirus 1" | name == "Human polyomavirus 2" | name == "Chordata" | 
                                    name == "Eukaryota" | name == "Viridiplantae" | name == "Archaea" |
                                    name == "Polyomaviridae")
bk_jc_df <- bk_jc_df[-c(2:4)]

# Make sure all values are in numeric format
numeric_cols <- which(sapply(bk_jc_df, is.numeric))
bk_jc_df[is.na(bk_jc_df)] <- 0

# Sum up categories with the same name
bk_jc_df <- ddply(bk_jc_df,"name",numcolwise(sum))

# Rename some categories
bk_jc_df[6,"name"] <- "BK"
bk_jc_df[7,"name"] <- "JC"
bk_jc_df[3,"name"] <- "Human"

bk_jc_df <- bk_jc_df[c(1:10),]

# Convert to usable table and clean up
bk_jc_df_2 <- t(bk_jc_df)
names <- rownames(bk_jc_df_2)
rownames(bk_jc_df_2) <- NULL
bk_jc_df_2 <- janitor::row_to_names(bk_jc_df_2,1)
names <- names[2:46]
names <- gsub("_.*$","",names)
bk_jc_df_2 <- cbind(names,bk_jc_df_2)
#colnames(bk_jc_df_2) <- c("Sample","JC","BK")

# Get rid of samples with NA in both BK and JC
bk_jc_df_2 <- as.data.frame(bk_jc_df_2)
df3 <- bk_jc_df_2[rowSums(is.na(bk_jc_df_2[2:3])) < 2L,]

# Convert everything into numerics again
df3$BK = as.numeric(as.character(df3$BK))
df3$JC = as.numeric(as.character(df3$JC))
df3$Bacteria = as.numeric(as.character(df3$Bacteria))
df3$Viruses = as.numeric(as.character(df3$Viruses))
df3$Fungi = as.numeric(as.character(df3$Fungi))
df3$Eukaryota = as.numeric(as.character(df3$Eukaryota))
df3$Human = as.numeric(as.character(df3$Human))
df3$Viridiplantae = as.numeric(as.character(df3$Viridiplantae))
df3$Archaea = as.numeric(as.character(df3$Archaea))
df3$Polyomaviridae = as.numeric(as.character(df3$Polyomaviridae))

df3 <- as.data.frame(df3)
df3[is.na(df3)] <- 0

colnames(df3)
cols <- c("Archaea","Bacteria","Human","Eukaryota","Fungi","BK","JC","Polyomaviridae","Viridiplantae","Viruses")
rpm_df <- data.frame(matrix(ncol=12,nrow = 0))
colnames <- c("names",cols)
colnames(rpm_df) <- colnames

# Go through and do rpm calculations for each category
for(sample in df3$names) {
  root_reads <- total[,grepl(sample,colnames(total))]
  rpm_sample <- df3 %>% filter (names==sample)
  
  sample_total_reads <- root_reads #+ rpm_sample$unclassified
  rpm_sample$Eukaryota <- rpm_sample$Eukaryota - rpm_sample$Fungi - rpm_sample$Viridiplantae - rpm_sample$Human
  rpm_sample$Polyomaviridae <- rpm_sample$Polyomaviridae - rpm_sample$BK - rpm_sample$JC
  rpm_sample$Viruses <- rpm_sample$Viruses - rpm_sample$BK - rpm_sample$JC - rpm_sample$Polyomaviridae
  
  #rpm_sample$Cellular_organisms <- rpm_sample$Cellular_organisms - rpm_sample$Bacteria - rpm_sample$Human - rpm_sample$Eukaryota - rpm_sample$Fungi - rpm_sample$Viridiplantae - rpm_sample$Archaea
  
  #other_root <- root_reads - rpm_sample$Archaea - rpm_sample$Bacteria - rpm_sample$Human - rpm_sample$Eukaryota - rpm_sample$Fungi - 
  #  rpm_sample$BK - rpm_sample$JC - rpm_sample$Viridiplantae - rpm_sample$Viruses - rpm_sample$Polyomaviridae - rpm_sample$Cellular_organisms
  #rpm_sample$Unclassified <- other_root
  
   rpm_sample[cols] <- lapply(rpm_sample[cols],`/`,sample_total_reads)
   rpm_sample[cols] <- lapply(rpm_sample[cols],`*`,1000000)
   
   other_root <- 1000000 - rpm_sample$Archaea - rpm_sample$Bacteria - rpm_sample$Human - rpm_sample$Eukaryota - rpm_sample$Fungi - 
     rpm_sample$BK - rpm_sample$JC - rpm_sample$Viridiplantae - rpm_sample$Viruses - rpm_sample$Polyomaviridae
   rpm_sample$Unclassified <- other_root
  
  rpm_df <- rbind(rpm_df,rpm_sample)
}

colnames(rpm_df) <- c("Sample","Archaea","Bacteria","Human","Other Eukaryota","Fungi","BK","JC","Polyomaviridae","Viridiplantae","Other Viruses","Unclassified")

df4 <- melt(rpm_df)

df4$value <- as.numeric(as.character(df4$value))

df4_bk <- df4 %>% filter(variable == "BK")
df4_jc <- df4 %>% filter(variable == "JC")

df4_bacteria <- df4%>% filter(variable=="Bacteria")

df4$Sample <- factor(df4$Sample, levels = unique(df4$Sample[rev(order(df4_bk$value,df4_jc$value,df4_bacteria$value))]))
#df4$Sample <- factor(df4$Sample, levels = unique(df4$Sample[(order(df4_bk$value))]))

df4$variable <- factor(df4$variable, levels = rev(c("BK","JC","Polyomaviridae","Other Viruses","Bacteria","Fungi","Viridiplantae","Human","Other Eukaryota",
                                                    "Archaea","Unclassified")))

ggplot(data=df4, aes(x=Sample,y=value,fill=variable)) + geom_col(width = 1) + 
  theme_classic() + ylab("RPM") + labs(fill = "Classification") + 
  #scale_fill_brewer(palette="Spectral", direction=-1, guide=guide_legend(reverse=TRUE)) + 
  #scale_fill_manual(values = rev(c("#F49D6E","#F6D55C","da$rkred","#3CAEA3","#E1E289",
  #                                 "#66a61e","#243E36","#CDE6F5","#A8D8EA","#393E46","#102B3F")), guide=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values = rev(c("#F49D6E","#F6D55C","darkred","orange","#3CAEA3","#E1E289",
                                   "#66a61e","#CDE6F5","darkseagreen","gray","#393E46")), guide=guide_legend(reverse=TRUE)) +
  #scale_fill_manual(values = rev(c("indianred1","firebrick2","darkred","dodgerblue2","darkseagreen","springgreen2",
  #                             "slateblue4","springgreen4","snow3","yellow","orange")), guide=guide_legend(reverse=TRUE)) +
  #theme(legend.position = c(0.45,0.88), legend.direction="horizontal") + 
  theme(legend.position="right") + coord_flip()
 #theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1, size = 7))

ggsave("BK_JC_Barplot.pdf",width=10,height=8)


###### Deprecated bacteria donut figure ######

# df5 <- df3 %>% select(names, BK, JC)
# df5 <- df5[order(df5$BK),]
# 
# df6 <- melt(df5)
# df6 <- df6 %>% group_by(variable) %>% arrange(value)
# 
# df6_bk <- df6 %>% filter(variable == "BK")
# df6$Sample <- factor(df6$Sample, levels=df4_bk$Sample)
# df6$Sample <- factor(df4$Sample, levels = unique(df4$Sample[order(df4_bk$value)]))
# 
# ggplot(data=df4, aes(x=Sample,y=value,fill=variable)) + geom_col(width = 0.8,position=position_dodge()) +
#   theme_classic() + ylab("RPM") + labs(fill = "Virus") + guides(fill = guide_legend(reverse=TRUE)) +
#   theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
# 
# ggsave("BK_JC_Barplot.pdf",width=8.5,height=4)

names <- gsub("_.*$","",colnames(only_clade))
colnames(cleaned_df) <- names

bacteria_rows <- dplyr::filter(cleaned_df, grepl("Bacteria",lineage))

bacteria_totals <- bacteria_rows %>% dplyr::filter(bacteria_rows$rank == "D")
bacteria_rows <- bacteria_rows %>% filter(bacteria_rows$rank == "S")

cols <- grepl("1000",names(bacteria_rows))
bacteria_rows <- bacteria_rows %>% replace(is.na(.),0) %>% mutate(sum=rowSums(.[5:60])) 
bacteria_rows <- bacteria_rows %>% filter(bacteria_rows$sum >= 50)

myColors <- distinctColorPalette(length(bacteria_rows$name))
#names(myColors) <- levels(alldata$Read)
#colScale <- scale_colour_manual(name = NULL, guide = FALSE, values = myColors)
colors <- data.frame(color=myColors)
bacteria_df <- cbind(bacteria_rows,colors)

bacteria_df[cols] <- lapply(bacteria_df[cols],`/`,1000000)
#bacteria_totals[cols] <- lapply(bacteria_totals[cols],`/`,1000000)


i=5
col_names <- colnames(bacteria_df)

sample_tot_bacteria <- bacteria_totals[1,i]
subset_df <- bacteria_df %>% select(1,61,i) %>% top_n(5)
sample_name <- col_names[i]
colnames(subset_df) <- c("name","color","rpm")
subset_df <- subset_df %>% mutate(fraction = rpm / sample_tot_bacteria)

other_rpm <- sample_tot_bacteria - sum(subset_df$rpm)
other_fraction <- 1 - sum(subset_df$fraction)
subset_df <- subset_df %>% add_row(name="Other",color="#D3D3D3",rpm=other_rpm,fraction=other_fraction)

ggplot(subset_df, aes(x = 2, y = rpm, fill = color)) + geom_bar(position = "fill", stat="identity", color="white") + 
  coord_polar(theta="y", start=0) + 
  scale_fill_identity("",guide = "legend",labels=paste(subset_df$name,subset_df$rpm,sep=": "), breaks=subset_df$color) + theme_void() + 
  theme(legend.text = element_text(size=10), legend.position="bottom", legend.margin = margin(t=0,unit="cm")) + guides(fill=guide_legend(ncol=1)) + 
  xlim(0.5,2.5) + #annotate("text",x = 0.5, y=0.2, label=paste(sample_name,"\n","Total RPM: ",sample_tot_bacteria,sep=""))
  annotate("text",x = 0.5, y=0.2, vjust=0,label=paste(sample_name),fontface=2, size = 5) + 
  annotate("text",x=0.5,y=0.2,vjust=2, label=paste("Total RPM: ",sample_tot_bacteria,sep=""))
# + annotate(geom = 'text', x = 0.5, y = 0.2, label=sample_name)

ggsave("Bacteria_Breakdowns.pdf",width=8.5,height=11)