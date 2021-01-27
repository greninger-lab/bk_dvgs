library(plotly)

cleaned_df <- only_clade

total <- cleaned_df %>% filter(name=="root")
domains <- cleaned_df %>% filter(rank=="D")

# Grab common uropathogenic bacteria
bacteria_df <- cleaned_df %>% filter(name == "Citrobacter" | name == "Escherichia coli" |
                                       name == "Enterobacter" | name == "Klebsiella" | name == "Candida" | name == "Acinetobacter baumannii" | 
                                     #| name == "Enterobacter cloacae complex" | 
                                      # name == "Klebsiella pneumoniae" | 
                                       name == "Proteus mirabilis" | name == "Providencia stuartii" | 
                                       name == "speciestaphylococcus aureus" | name == "speciestaphylococcus saprophyticus" | 
                                       name == "Enterococcus faecium" | name == "Enterococcus faecalis" | name == "Pseudomonas aeruginosa group" |
                                       name == "unclassified" )#| name == "Human polyomavirus 1")

# Convert values to numeric
numeric_cols <- which(sapply(bacteria_df, is.numeric))
bacteria_df[is.na(bacteria_df)] <- 0

bk_jc_df <- bacteria_df
bk_jc_df <- bk_jc_df[-c(2:4)]
bk_jc_df_2 <- t(bk_jc_df)
names <- rownames(bk_jc_df_2)
rownames(bk_jc_df_2) <- NULL
bk_jc_df_2 <- janitor::row_to_names(bk_jc_df_2,1)

names <- names[2:46]
names <- gsub("_.*$","",names)
bk_jc_df_2 <- cbind(names,bk_jc_df_2)
#colnames(bk_jc_df_2) <- c("Sample","JC","BK")

bk_jc_df_2 <- as.data.frame(bk_jc_df_2)
df3 <- bk_jc_df_2
df3[] <- cbind(df3[1], lapply(df3[2:15], function(x) as.numeric(as.character(x))))

cols <- colnames(df3)
rpm_df <- data.frame(matrix(ncol=15,nrow = 0))
colnames(rpm_df) <- cols
cols <- cols[2:15]

# Calculate rpms
for(sample in df3$names) {
  root_reads <- total[,grepl(sample,colnames(total))]
  rpm_sample <- df3 %>% filter (names==sample)
  # Actual total is root + unclassified
  sample_total_reads <- root_reads #+ rpm_sample$unclassified
  rpm_sample[cols] <- lapply(rpm_sample[cols],`/`,sample_total_reads)
  rpm_sample[cols] <- lapply(rpm_sample[cols],`*`,1000000)
  rpm_df <- rbind(rpm_df,rpm_sample)
}

rpm_df <- rpm_df %>% select(1,3:15)
colnames(rpm_df)
colnames(rpm_df) <- c("Sample","S.aureus","Enterobacter","Citrobacter","Klebsiella","E.coli","Pseudomonas aeruginosa","Candida","Enterococcus faecalis","Acinetobacter baumannii","S. saprophyticus","Providencia stuartii","Enterococcus faecium","Proteus mirabilis")
col_order <- c("Sample","S.aureus","S. saprophyticus","Enterobacter","Citrobacter","Klebsiella","E.coli","Pseudomonas aeruginosa","Enterococcus faecalis","Enterococcus faecium","Acinetobacter baumannii","Providencia stuartii","Proteus mirabilis","Candida")

df4 <- melt(rpm_df)
df5 <- (rpm_df[,col_order])

df4$value <- as.numeric(as.character(df4$value))
df4_bk <- df4 %>% filter(variable == "BK")


df6 <- df4 %>% filter(variable!="BK")
df6$Sample <- factor(df6$Sample, levels = unique(df6$Sample[rev(order(df4_bk$value))]))

# Count as uropathogenic bacteria present when any bacteria is >10 rpm
df5$UB_present <- apply(df5[, 2:13], MARGIN = 1, function(x) any(x > 10))
#df5$Fungi_present <- apply(df5[,14], MARGIN = 1, function(x) any(x > 10))

df5$Bacteria_Sum <- rowSums(df5[2:13])
df5 <- df5[order(-df5$Bacteria_Sum),]


# Plot for funsies
ggplot(data=df6, aes(x=Sample,y=value,fill=variable)) + geom_col(width = 1) + 
  theme_classic() + ylab("RPM") + labs(fill = "Classification") + 
  scale_fill_brewer(palette="Spectral", direction=-1, guide=guide_legend(reverse=TRUE)) + 
  #theme(legend.position = c(0.45,0.88), legend.direction="horizontal") + 
  theme(legend.position="right") + coord_flip()
#theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1, size = 7))

write.csv(df5,"Uropathogenic_bacteria_fungus_table.csv",row.names=FALSE,quote=FALSE)
ggsave("Uropathogenic_bacteria_BK.pdf",width=10,height=8)
