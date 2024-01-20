# Sincerely grateful to Pat Schloss and his many wonderful tutorials at https://riffomonas.org

library(tidyverse)

d = read_csv("AMR.mobileOG.VFDB.merged.csv")
metadata = read_csv("metadata4.csv")

dj = inner_join(d, metadata, by="Sample") %>%
  mutate(Collection_Type = factor(Collection_Type, levels=c("Autosampler1","Autosampler2","Biofilm","Grab","Moore")))

# first, some statistics out of curiosity
# count contigs that have both AMRfinder and mobileOG result
# dj %>%
#   filter(!is.na(node_id) & !is.na(mobileOG_Result)) %>%
#   add_count(Contig,class,Collection_Type)
# count number of individual AMRfinder, mobileOG, and VFDB results per contig to find the contigs with the most genes
cc = dj %>%
  add_count(Contig)
# don't include contigs without ARG
cc = dj %>%
  filter(!is.na(node_id)) %>%
  add_count(Contig)

# first plot contigs with ARGs without MGEs
# pool counts of contigs within categories and average counts among samples within a type 
subtotals1 = dj %>%
  select(Contig,Sample,Collection_Type,class,`Major mobileOG Category`) %>%
  distinct() %>%
  
  # count ARGs without mobileOG result
  filter(!is.na(class) & is.na(`Major mobileOG Category`)) %>%
  add_count(class,Sample,Collection_Type) %>%
  group_by(Collection_Type,class) %>%   # average counts among samples within collection_type
  mutate(Contigs_with_ARG_without_MGE = mean(n))

# mark rare categories as "other"
others = subtotals1 %>%
  group_by(class) %>%
  summarize(pool = max(Contigs_with_ARG_without_MGE) < 10, 
            max = max(Contigs_with_ARG_without_MGE),
            .groups="drop")

# pool rare categories into Other by summing their counts
subtotals2 = inner_join(subtotals1, others, by="class") %>%
  mutate(class = if_else(pool, "Other", class)) %>%
  add_count(Collection_Type, Sample, class) %>%     
  group_by(Collection_Type, class) %>%
  mutate(Contigs_with_ARG_without_MGE2 = mean(nn)) %>%        # (same as above, this time including Other Category)
  select(Sample,Collection_Type,class,Contigs_with_ARG_without_MGE2) %>%                  # clean up the final table
  distinct() %>%
  ungroup()
  
# find the max count across all samples for each class  
cm = subtotals2 %>% 
  group_by(class) %>%
  summarize(maxC = max(Contigs_with_ARG_without_MGE2), .groups="drop")
final_ARGs_without_MGE = inner_join(subtotals2, cm, by="class")

# how many?
N = sum(others$pool == FALSE)
N

# alternative color scale = pal2
library(RColorBrewer)
pal2 = brewer.pal(N, "Paired")
names(pal2) = others$class[others$pool == FALSE]
pal2["Other"] = "darkgray"

# plot
p = final_ARGs_without_MGE %>%
  mutate(class = factor(class),
         class = fct_reorder(class, maxC)) %>%
  ggplot(aes(x=Collection_Type, y=Contigs_with_ARG_without_MGE2, fill=class)) +
  geom_col() +
  scale_fill_manual(values=pal2) +
  #scale_color_manual(values=c(NA,"black"), guide="none", na.value="transparent") +
  #geom_text(aes(label = l), size=4, position=position_stack(vjust=0.5)) +
  #geom_text(aes(y=CollectionTypeTotalTPM, label=genes, fontface=2), vjust=-0.5) +
  #annotate(geom="text", label="Number of Genes Detected", size=4, fontface=2, x=Inf, y=Inf) +
  #scale_x_discrete(breaks=c("Autosampler1","Autosampler2","Biofilm","Grab","Moore")) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black", face= "bold"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) +
  #theme(plot.margin = margin(t = 20, unit = "pt")) +
  coord_cartesian(clip = "off")
p
ggsave("AMR-mobileOG-contigs-barplot.jpg", p)
  

# plot contigs with both ARGs and MGEs
Both_subtotals1 = dj %>%
  select(Contig,Sample,Collection_Type,class,`Major mobileOG Category`) %>%
  distinct() %>%

  # count ARGs with mobileOG result
  filter(!is.na(class) & !is.na(`Major mobileOG Category`)) %>%
  add_count(class,Sample,Collection_Type) %>%
  group_by(Collection_Type,class) %>%   # average counts among samples within collection_type
  mutate(Contigs_with_ARG_and_MGE = mean(n))

# use same rare categories as for first plot
Both_subtotals2 = inner_join(Both_subtotals1, others, by="class") %>%
  mutate(class = if_else(pool, "Other", class)) %>%
  add_count(Collection_Type, Sample, class) %>%     
  group_by(Collection_Type, class) %>%
  mutate(Contigs_with_ARG_and_MGE2 = mean(nn)) %>%        # (same as above, this time including Other Category)
  select(Sample,Collection_Type,class,Contigs_with_ARG_and_MGE2) %>%                  # clean up the final table
  distinct() %>%
  ungroup()

# use same maxC as first plot to keep same order
final_ARGs_and_MGE = inner_join(Both_subtotals2, cm, by="class")

# plot with same color palette as for first plot
p = final_ARGs_and_MGE %>%
  mutate(class = factor(class),
         class = fct_reorder(class, maxC)) %>%
  ggplot(aes(x=Collection_Type, y=Contigs_with_ARG_and_MGE2, fill=class)) +
  geom_col() +
  scale_fill_manual(values=pal2) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black", face= "bold"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) +
  coord_cartesian(clip = "off")
p
ggsave("AMR-with-mobileOG-contigs-barplot.jpg", p)
   