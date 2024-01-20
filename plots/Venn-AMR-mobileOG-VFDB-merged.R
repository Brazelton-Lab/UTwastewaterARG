# Sincerely grateful to Pat Schloss and his many wonderful tutorials at https://riffomonas.org

library(tidyverse)

d = read_csv("AMR.mobileOG.VFDB.merged.csv")
metadata = read_csv("metadata4.csv")

dj = inner_join(d, metadata, by="Sample") %>%
  mutate(Collection_Type = factor(Collection_Type, levels=c("Autosampler1","Autosampler2","Biofilm","Grab","Moore")))

# count contigs that have both AMRfinder and mobileOG result
c = dj %>%
  filter(!is.na(node_id) & !is.na(mobileOG_Result)) %>%
  add_count(Contig,class,Collection_Type)

# VENN DIAGRAMS
# One Collection_Type at a time, construct a list of contigs in AMRfinder and mobileOG results
v = dj %>%
  select(Collection_Type,Contig,node_id,mobileOG_Result) %>%
  pivot_longer(c("node_id","mobileOG_Result"), values_to="result_ID", names_to="Result" ) %>%
  drop_na()
 
library(ggVennDiagram)
vp1 = v %>%
  filter(Collection_Type == "Autosampler1") %>%
  select(-Collection_Type) %>%
  {split(.$Contig,.$Result)} %>%
  ggVennDiagram(category.names = c("MGEs","ARGs")) + 
  scale_x_continuous(expand = expansion(mult = .2)) +
  labs(title = "Autosampler1")
vp1

vp2 = v %>%
  filter(Collection_Type == "Autosampler2") %>%
  select(-Collection_Type) %>%
  {split(.$Contig,.$Result)} %>%
  ggVennDiagram(category.names = c("MGEs","ARGs")) + 
  scale_x_continuous(expand = expansion(mult = .2)) +
  labs(title = "Autosampler2")
vp2

vp3 = v %>%
  filter(Collection_Type == "Biofilm") %>%
  select(-Collection_Type) %>%
  {split(.$Contig,.$Result)} %>%
  ggVennDiagram(category.names = c("MGEs","ARGs")) + 
  scale_x_continuous(expand = expansion(mult = .2)) +
  labs(title = "Biofilm")
vp3

vp4 = v %>%
  filter(Collection_Type == "Moore") %>%
  select(-Collection_Type) %>%
  {split(.$Contig,.$Result)} %>%
  ggVennDiagram(category.names = c("MGEs","ARGs")) + 
  scale_x_continuous(expand = expansion(mult = .2)) +
  labs(title = "Moore")
vp4

vp5 = v %>%
  filter(Collection_Type == "Grab") %>%
  select(-Collection_Type) %>%
  {split(.$Contig,.$Result)} %>%
  ggVennDiagram(category.names = c("MGEs","ARGs")) + 
  scale_x_continuous(expand = expansion(mult = .2)) +
  labs(title = "Grab")
vp5

# count contigs that have both AMRfinder and VFDB result
c2 = dj %>%
  filter(!is.na(node_id) & !is.na(VFDB_Result)) %>%
  add_count(Contig,class,Collection_Type)



