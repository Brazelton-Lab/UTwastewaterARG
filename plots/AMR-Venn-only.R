# Sincerely grateful to Pat Schloss and his many wonderful tutorials at https://riffomonas.org

library(tidyverse)

d = read_csv("merged.abunds.nodes.counts.csv")  # note that this table contains counts, not TPM, for plotting seq coverage below
metadata = read_csv("metadata.csv")

# traditional transpose matrix
dt = d %>% select(-class,-subclass,-type,-subtype) %>% 
  pivot_longer(cols=-1,) %>% 
  pivot_wider(names_from=node_id,values_from=value) %>%
  rename(Samples = name)

# convert transposed matrix to long format
dl = dt %>% 
  pivot_longer(-Samples, names_to="node_id", values_to="num_reads")

# join metadata
dj = inner_join(dl, metadata, by="Samples") %>%
  mutate(Collection_Type = factor(Collection_Type, levels=c("Autosampler1","Autosampler2","Biofilm","Grab","Moore")))


# VENN DIAGRAMS
# Construct list of genes for each Collection_Type
vl1 = dl %>%
  select(num_reads,node_id,Samples) %>%
  filter(num_reads> 0) %>%
  {split(.$node_id,.$Samples)}

library(ggVennDiagram)
# Venn diagram with all Collection Types
# vp1.1 = ggVennDiagram(vl1, label="count") + scale_x_continuous(expand = expansion(mult = .2))
# vp1.1
# Venn diagram comparing only Autosampler1, Biofilm, and Moore
vp1.2 = ggVennDiagram(vl1[c(12,2,16)], label="count") + scale_x_continuous(expand = expansion(mult = .2))
vp1.2

# Venn diagram comparing Collection_Types BETA-LACTAM only
vl2 = dl %>%
  select(num_reads,node_id,Samples) %>%
  filter(num_reads> 0) %>%
  {split(.$node_id,.$Samples)}

# Venn diagram with all Collection Types
# vp2.1 = ggVennDiagram(vl2, label="count") + scale_x_continuous(expand = expansion(mult = .2))
# vp2.1
# Venn diagram comparing only Autosampler1, Biofilm, and Moore
vp2.2 = ggVennDiagram(vl2[c(12,2,16)], label="count") + scale_x_continuous(expand = expansion(mult = .2))
vp2.2

# get genes unique to one Collection_Type
everything_else = union(vl1$"Autosampler1",vl1$"Autosampler2")
everything_else = union(everything_else,vl1$"Moore")
everything_else = union(everything_else,vl1$"Grab")
u = setdiff(vl1$Biofilm,everything_else)
u

# save Venn diagrams
# ggsave("AMR-Venn-all-node_ids.png",vp1.1)
# ggsave("AMR-Venn-all-node_ids.jpg",vp1.1)
ggsave("AMR-Venn-subsampled-3samples-node_ids.png",vp1.2)
ggsave("AMR-Venn-subsampled-3samples-node_ids.jpg",vp1.2)
# ggsave("AMR-Venn-all-BLA-node_ids.png",vp2.1)
# ggsave("AMR-Venn-all-BLA-node_ids.jpg",vp2.1)
ggsave("AMR-Venn-subsampled-3samples-BLA-node_ids.png",vp2.2)
ggsave("AMR-Venn-subsampled-3samples-BLA-node_ids.jpg",vp2.2)

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
  group_by(Samples, Collection_Type) %>%
  summarize(num_genes = richness(num_reads),
            #shannon = shannon(value),
            simpson = simpson(num_reads),
            invsimpson_diversity = 1/simpson,
            simpsoneven = 1/simpson/num_genes,
            n = sum(num_reads),
            .groups="drop") %>%
  pivot_longer(cols=c(num_genes, invsimpson_diversity), names_to="metric") %>%
  ggplot(aes(x=n, y=value)) +
  geom_point(aes(color=Collection_Type, size=8)) +
  facet_wrap(~metric, nrow=3, scales="free_y")
dp

ggsave("AMR-diversity.jpg", dp)
