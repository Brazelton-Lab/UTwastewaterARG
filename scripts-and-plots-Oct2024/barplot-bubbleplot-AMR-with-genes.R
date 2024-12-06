# Sincerely grateful to Pat Schloss and his many wonderful tutorials at https://riffomonas.org

library(tidyverse)

d = read_csv("amrfinder-contig-cov-merged-all-samples.csv")
metadata = read_csv("metadata.csv")

# include only ARGs, not stress or virulence genes
df = d %>% 
  filter(`Element type` == "AMR") %>%
  select(!c(`Element type`,`Element subtype`))

# pivot to long format
dl = df %>% 
  pivot_longer(c(`A12-17`,`M2-3`,`A22-24`,`B4-5`,`G6-11`, `A2024-take1`, `A2024-take2`, `M2024-take1`, `M2024-take2`, `A2024-hybrid`, `M2024-hybrid`), names_to="Samples", values_to="TPM") %>%
  pivot_longer(c(Class,Subclass), names_to="level", values_to="Category")

# join with metadata
dj = inner_join(dl, metadata, by="Samples")

# subtotal TPM for gene symbol and ARG class
tpm = dj %>%
  filter(level=="Class") %>%
  select(`Gene symbol`,Samples,TPM,Category,Collection_Type) %>%
  group_by(`Gene symbol`,Samples,Category,Collection_Type) %>%
  replace(is.na(.), 0) %>%
  summarize(Gene_TPM = sum(TPM)) 

# write table to file
u = tpm %>% 
  select(!c(Collection_Type)) %>%
  pivot_wider(names_from="Samples", values_from="Gene_TPM")
write.csv(u, "amrfinderplus_comparison.csv", row.names=FALSE)

# mark rare categories as "other"
others = tpm %>%
  group_by(Category) %>%
  summarize(pool = sum(Gene_TPM) < 700, 
            max = sum(Gene_TPM),
            .groups="drop")

# pool rare categories into Other by summing their TPM
tpm = inner_join(tpm, others, by="Category") %>%
  mutate(Category = if_else(pool, "Other", Category)) %>%
  group_by(Collection_Type, Samples, Category) %>%                 
  mutate(total_TPM2 = sum(Gene_TPM)) %>%     
  ungroup()

# find the max total_TPM2 across all samples for each Category  
cm = tpm %>% 
  group_by(Category) %>%
  summarize(maxC = max(total_TPM2), .groups="drop")
tpm = inner_join(tpm, cm, by="Category")

# how many?
N = sum(others$pool == FALSE)
N

# alternative color scale = pal2
library(RColorBrewer)
pal2 = brewer.pal(N, "Paired")
names(pal2) = others$Category[others$pool == FALSE]

# show labels and borders only for the most abundant genes
tpm = tpm %>%
  mutate(l = factor(vector("character", length(tpm$`Gene symbol`)))) %>%
  mutate(l = if_else(Gene_TPM > 100, `Gene symbol`, l)) %>%
  mutate(bc = factor(vector("character", length(tpm$`Gene symbol`)))) %>%
  mutate(bc = if_else(Gene_TPM > 100, "black", bc))

# count number of genes in each Sample
tpm_samples = tpm %>%
  group_by(Samples) %>%
  summarize(SampleGenes = sum(Gene_TPM > 0),
            SampleTotalTPM = sum(Gene_TPM)) 
plotmax = max(tpm_samples$SampleTotalTPM)

# plot
p = tpm %>%
  mutate(Category = factor(Category),
        Category = fct_reorder(Category, maxC)) %>%
  ggplot(aes(x=Samples, y=Gene_TPM, fill=Category, color = bc=="black")) +
  geom_col() +
  scale_fill_manual(values=pal2) +
  scale_color_manual(values=c(NA,"black"), guide="none", na.value="transparent") +
  geom_text(aes(label = l), size=3, position=position_stack(vjust=0.5)) +
  geom_label(inherit.aes = FALSE, data = tpm_samples, aes(label=SampleGenes, x=Samples, y=SampleTotalTPM), nudge_y=100) +
  #annotate(geom="text", label="Number of Genes Detected", size=4, fontface=2, x=Inf, y=Inf) +
  #scale_x_discrete(breaks=c("Autosampler1","Autosampler2","Biofilm","Grab","Moore")) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black", face= "bold", angle=90), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) +
  theme(plot.margin = margin(t = 20, unit = "pt")) +
  coord_cartesian(clip = "off")
p
ggsave("AMRbarplot-allsamples.png", p, width=12)
  
# BETA-LACTAMASE only
bla = tpm %>%
  select(`Gene symbol`,Samples,Category,Collection_Type,Gene_TPM)

# CHOOSE: assign border colors and labels to key genes
bla = bla %>%
  # border color
  mutate(bc = factor(vector("character", length(bla$`Gene symbol`)))) %>%
  # mutate(bc = if_else(`Gene symbol` == "blaKPC-2" & Gene_TPM > 0, "black", bc)) %>%
  # mutate(bc = if_else(`Gene symbol` == "blaKPC" & Gene_TPM > 0, "black", bc)) %>%
  mutate(bc = if_else(`Gene symbol` == "blaVIM-2" & Gene_TPM > 0, "black", bc)) %>%
  mutate(bc = if_else(`Gene symbol` == "blaVIM" & Gene_TPM > 0, "black", bc)) %>%
  # label
  mutate(l = factor(vector("character", length(bla$`Gene symbol`)))) %>%
  # mutate(l = if_else(`Gene symbol` == "blaKPC-2" & Gene_TPM > 0,`Gene symbol`, l)) %>%
  # mutate(l = if_else(`Gene symbol` == "blaKPC" & Gene_TPM > 0, `Gene symbol`, l)) %>%
  mutate(l = if_else(`Gene symbol` == "blaVIM-2" & Gene_TPM > 0, `Gene symbol`, l)) %>%
  mutate(l = if_else(`Gene symbol` == "blaVIM" & Gene_TPM > 0, `Gene symbol`, l))

# OR: show labels for the most abundant genes
bla = bla %>%
  mutate(l = factor(vector("character", length(bla$`Gene symbol`)))) %>%
  mutate(l = if_else(Gene_TPM > 30, `Gene symbol`, l)) %>%
  mutate(bc = factor(vector("character", length(bla$`Gene symbol`)))) %>%
  mutate(bc = if_else(Gene_TPM > 30, "black", bc))

# include only bla genes
bla = bla %>% filter(Category=="BETA-LACTAM")

# count number of genes in each collection type for BETA-LACTAM only
bla_samples = bla %>%
  group_by(Samples) %>%
  summarize(SampleGenes = sum(Gene_TPM > 0),
            SampleTotalTPM = sum(Gene_TPM)) 
blaplotmax = max(bla_samples$SampleTotalTPM)

# BLA only plot
p2 = bla %>%
  ggplot(aes(x=Samples, y=Gene_TPM, fill=Category, color = bc=="black")) +
  geom_col() +
  scale_fill_manual(values=pal2) +
  scale_color_manual(values=c(NA,"black"), guide="none", na.value="gray") +
  geom_text(aes(label = l), size=3, position=position_stack(vjust=0.5)) +
  geom_label(inherit.aes = FALSE, data = bla_samples, aes(label=SampleGenes, x=Samples, y=SampleTotalTPM), nudge_y=100) +
  #annotate(geom="text", label="Number of Genes Detected", size=4, fontface=2, x=Inf, y=Inf) +
  #scale_x_discrete(breaks=c("Autosampler1","Autosampler2","Biofilm","Grab","Moore")) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black", face= "bold", angle=90), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) +
  theme(plot.margin = margin(t = 20, unit = "pt")) +
  coord_cartesian(clip = "off")  
p2
ggsave("BLAbarplot-allsamples-VIM.png", p2, width=12)


# VENN DIAGRAMS
# Construct list of genes for each Collection_Type
vl1 = tf1 %>%
  select(TPM,gene_symbol,Collection_Type) %>%
  filter(TPM> 0) %>%
  {split(.$gene_symbol,.$Collection_Type)}

library(ggVennDiagram)
# Venn diagram with all Collection Types
vp1.1 = ggVennDiagram(vl1, label="count") + scale_x_continuous(expand = expansion(mult = .2))
vp1.1
# Venn diagram comparing only Autosampler1, Biofilm, and Moore
vp1.2 = ggVennDiagram(vl1[c(1,3,5)], label="count") + scale_x_continuous(expand = expansion(mult = .2))
vp1.2

# Venn diagram comparing Collection_Types BETA-LACTAM only
vl2 = tf2 %>%
  select(TPM,gene_symbol,Collection_Type) %>%
  filter(TPM> 0) %>%
  {split(.$gene_symbol,.$Collection_Type)}

# Venn diagram with all Collection Types
vp2.1 = ggVennDiagram(vl2, label="count") + scale_x_continuous(expand = expansion(mult = .2))
vp2.1
# Venn diagram comparing only Autosampler1, Biofilm, and Moore
vp2.2 = ggVennDiagram(vl2[c(1,3,5)], label="count") + scale_x_continuous(expand = expansion(mult = .2))
vp2.2

# get genes unique to one Collection_Type
everything_else = union(vl1$"Autosampler1",vl1$"Autosampler2")
everything_else = union(everything_else,vl1$"Moore")
everything_else = union(everything_else,vl1$"Grab")
u = setdiff(vl1$Biofilm,everything_else)
u

# save Venn diagrams
ggsave("AMR-Venn-all.png",vp1.1)
ggsave("AMR-Venn-all.jpg",vp1.1)
ggsave("AMR-Venn-all-3samples.png",vp1.2)
ggsave("AMR-Venn-all-3samples.jpg",vp1.2)
ggsave("AMR-Venn-all-BLA.png",vp2.1)
ggsave("AMR-Venn-all-BLA.jpg",vp2.1)
ggsave("AMR-Venn-all-3samples-BLA.png",vp2.2)
ggsave("AMR-Venn-all-3samples-BLA.jpg",vp2.2)


# BUBBLE PLOT
b1 = tf1 %>% 
  arrange(desc(TPM)) %>% # arrange(desc()) sorts the bubbles so that smaller bubbles are plotted on top of larger bubbles
  ggplot(aes(x = Collection_Type, y = l)) + # l is the list of abundant genes created above
  geom_point(aes(size = TPM, fill=Category), shape=21, alpha = 0.75) + 
  scale_size_continuous(limits = c(1, 500), range=c(0.1,20), breaks = c(10,100,200,300)) +
  scale_fill_manual(values=pal2) +
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 11, face ="bold", colour ="black"), 
        legend.title = element_text(size = 11, face = "bold"), panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        #plot.margin = margin(t = 50, unit = "pt"),
        #legend.position = "none", 
        panel.grid.major.y = element_line(colour = "grey95")) +
  coord_cartesian(clip = "off") #+
  #scale_y_discrete(limits = rev(levels(tf1$gene_symbol)))
b1

# BUBBLE PLOT BLA only
b2 = tf2 %>% 
  arrange(desc(TPM)) %>% # arrange(desc()) sorts the bubbles so that smaller bubbles are plotted on top of larger bubbles
  ggplot(aes(x = Collection_Type, y = gene_symbol)) + # l is the list of abundant genes created above
  geom_point(aes(size = TPM, fill=Category), shape=21, alpha = 0.75) + 
  scale_size_continuous(limits = c(1, 500), range=c(0.1,20), breaks = c(10,100,200,300)) +
  scale_fill_manual(values=pal2) +
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 11, face = "bold"), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 11, face ="bold", colour ="black"), 
        legend.title = element_text(size = 11, face = "bold"), panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        #plot.margin = margin(t = 50, unit = "pt"),
        #legend.position = "none", 
        panel.grid.major.y = element_line(colour = "grey95")) +
  coord_cartesian(clip = "off") #+
#scale_y_discrete(limits = rev(levels(tf1$gene_symbol)))
b2

# save bubble plots
ggsave("AMR-bubbles.jpg", b1)
ggsave("AMR-bubbles-BLA.jpg", b2)

