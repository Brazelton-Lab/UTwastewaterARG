# Sincerely grateful to Pat Schloss and his many wonderful tutorials at https://riffomonas.org

library(tidyverse)

d = read_tsv("amrfinderplus.txt") %>%
  select("Name","Gene symbol","Sequence name","Scope","Element type","Element subtype","Class","Subclass","Method") %>%
  rename(isolate = Name) %>%
  
  # include only ARGs, not stress or virulence genes
  filter(`Element type` == "AMR")

metadata = read_csv("grandeur_summary.csv") %>%
  rename(isolate = sample) %>%
  select(isolate,fastani_reference) %>%
  separate(fastani_reference, sep="_", c("genus","species"), extra="drop") %>%
  separate(genus, sep="/", c(NA,"genus"), extra="drop") %>%
  filter(!isolate == "Undetermined")

# join with metadata and assign factors
dj = inner_join(d, metadata, by="isolate")

# count genes within categories per isolate and average number of genes among isolates within a genus 
t = dj %>%
  group_by(isolate,genus,Class) %>%    
  add_count(name="n_isolate") %>%              # counts genes per class per isolate
  ungroup() %>%
  group_by(genus,Class) %>%
  add_count(name="n_genus") %>%
  ungroup() %>%
  group_by(genus) %>%
  mutate(isolates_per_genus = n_distinct(isolate),
         avg_per_genus = n_genus / isolates_per_genus)

# how many?
N = n_distinct(t$Class)
N

# make palette
library(RColorBrewer)
pal2 = brewer.pal(N, "Paired")
names(pal2) = t %>% group_by(Class) %>% distinct(Class) %>% pull()

# plot
p = t %>%
  mutate(genus = factor(genus)) %>%
  ggplot(aes(x=isolate, y=n_isolate, fill=Class)) +
  geom_col() +
  scale_fill_manual(values=pal2) +
  #scale_x_discrete(breaks=c("Escherichia","Citrobacter","Providencia")) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) +
  #theme(plot.margin = margin(t = 20, unit = "pt")) +
  facet_wrap(~genus, scales="free") +
  ylab("Number of genes")
p
ggsave("isolates-barplot.jpg", p)

# VENN DIAGRAMS
# Construct list of genes for each Collection_Type
vl = t %>%
  select(n_isolate,`Gene symbol`,genus) %>%
  filter(n_isolate> 0) %>%
  {split(.$`Gene symbol`,.$genus)}

library(ggVennDiagram)
# Venn diagram with all Collection Types
vp1 = ggVennDiagram(vl, label="count") + scale_x_continuous(expand = expansion(mult = .2))
vp1

# save Venn diagrams
ggsave("isolates-Venn-all.jpg",vp1)


