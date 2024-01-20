# Sincerely grateful to Pat Schloss and his many wonderful tutorials at https://riffomonas.org

library(tidyverse)

d = read_csv("deeparg.final2.csv")
metadata = read_csv("metadata2.csv")

# join with metadata and assign factors
dj = inner_join(d, metadata, by="Sample") %>%
  mutate(Collection_Type = factor(Collection_Type, levels=c("Autosampler1","Autosampler2","Biofilm","Grab","Moore")))

# VENN DIAGRAMS
# Construct list of genes for each Collection_Type
vl1 = d %>%
  select(`16s-NormalizedReadCount`,`ARG-group`,Sample) %>%
  filter(`16s-NormalizedReadCount`> 0) %>%
  {split(.$`ARG-group`,.$Sample)}

library(ggVennDiagram)
# Venn diagram with all Collection Types
# vp1.1 = ggVennDiagram(vl1, label="count") + scale_x_continuous(expand = expansion(mult = .2))
# vp1.1
# Venn diagram comparing only Autosampler1, Biofilm, and Moore
vp1.2 = ggVennDiagram(vl1[c("2-3.subset5.qtrim.deeparg","12-17.subset2.qtrim.deeparg","4-5.deeparg")], label="count") + scale_x_continuous(expand = expansion(mult = .5))
vp1.2

# Venn diagram comparing Collection_Types BETA-LACTAM only
vl2 = d %>%
  filter(`ARG-category`=="beta-lactam") %>%
  select(`16s-NormalizedReadCount`,`ARG-group`,Sample) %>%
  filter(`16s-NormalizedReadCount`> 0) %>%
  {split(.$`ARG-group`,.$Sample)}

# Venn diagram with all Collection Types
# vp2.1 = ggVennDiagram(vl2, label="count") + scale_x_continuous(expand = expansion(mult = .2))
# vp2.1
# Venn diagram comparing only Autosampler1, Biofilm, and Moore
vp2.2 = ggVennDiagram(vl2[c("2-3.subset5.qtrim.deeparg","12-17.subset2.qtrim.deeparg","4-5.deeparg")], label="count") + scale_x_continuous(expand = expansion(mult = .5))
vp2.2

# get genes unique to one Collection_Type
everything_else = union(vl1$"Autosampler1",vl1$"Autosampler2")
everything_else = union(everything_else,vl1$"Moore")
everything_else = union(everything_else,vl1$"Grab")
u = setdiff(vl1$Biofilm,everything_else)
u

# save Venn diagrams
# ggsave("deeparg-Venn-all.png",vp1.1)
# ggsave("deeparg-Venn-all.jpg",vp1.1)
ggsave("deeparg-Venn-subsampled-3samples.png",vp1.2)
ggsave("deeparg-Venn-subsampled-3samples.jpg",vp1.2)
# ggsave("deeparg-Venn-all-BLA.png",vp2.1)
# ggsave("deeparg-Venn-all-BLA.jpg",vp2.1)
ggsave("deeparg-Venn-subsampled-3samples-BLA.png",vp2.2)
ggsave("deeparg-Venn-subsampled-3samples-BLA.jpg",vp2.2)

# calculate alpha diversity (code from Pat Schloss Code Club) https://www.youtube.com/watch?v=wq1SXGQYgCs
richness <- function(x){
  sum(x>0)
}

simpson <- function(x){
  n <- sum(x)
  sum(x * (x-1) / (n * (n-1)))
  # vegan simpson: 1 - sum((x/n)^2)
}

dp = dj %>%
  group_by(Sample, Collection_Type) %>%
  summarize(num_genes = richness(ReadCount),
            #shannon = shannon(value),
            simpson = simpson(ReadCount),
            invsimpson_diversity = 1/simpson,
            simpsoneven = 1/simpson/num_genes,
            n = sum(ReadCount),
            .groups="drop") %>%
  pivot_longer(cols=c(num_genes, invsimpson_diversity), names_to="metric") %>%
  ggplot(aes(x=n, y=value)) +
  geom_point(aes(color=Collection_Type, size=10)) +
  facet_wrap(~metric, nrow=3, scales="free_y")

ggsave("deeparg-diversity.jpg", dp)

