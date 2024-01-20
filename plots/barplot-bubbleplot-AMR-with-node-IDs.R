# Sincerely grateful to Pat Schloss and his many wonderful tutorials at https://riffomonas.org

library(tidyverse)

d = read_csv("merged.abunds.tpm.nodes.csv")
metadata = read_csv("metadata3.csv")

types = d %>% select(node_id,class,subclass,type,subtype)

# traditional transpose matrix
dt = d %>% select(-class,-subclass,-type,-subtype) %>% 
  pivot_longer(cols=-1,) %>% 
  pivot_wider(names_from=node_id,values_from=value) %>%
  rename(Samples = name)

# convert transposed matrix to long format
dl = dt %>% 
  pivot_longer(-Samples, names_to="node_id", values_to="TPM")

# join all three data tables (counts, types, metadata)
# and use pivot_longer to make separate rows for types of ARGs
dj = inner_join(dl, metadata, by="Samples") %>%
  inner_join(., types, by="node_id") %>%
  pivot_longer(c("class","subclass","type","subtype"), names_to="level", values_to="Category") %>%
  mutate(Collection_Type = factor(Collection_Type, levels=c("Autosampler1","Autosampler2","Biofilm","Grab","Moore")))

# # keep key genes as separate categories
# dj =  dj %>%
#   # change value in data table dj
#   mutate(Category = ifelse(node_id == "blaKPC", node_id, Category)) %>%
#   mutate(Category = ifelse(node_id == "blaVIM", node_id, Category)) %>%
#   mutate(Category = ifelse(node_id == "blaOXA-48_fam", node_id, Category)) %>%
#   mutate(Category = ifelse(node_id == "blaDIM-SIM-IMP", node_id, Category)) %>%
#   mutate(Category = ifelse(node_id == "blaNDM", node_id, Category))

# pool TPMs within categories and average TPM among samples within a type 
tpm = dj %>%
  filter(level=="class") %>%
  group_by(Collection_Type,Samples,Category) %>%    # delete "Collection_Type" if you want to show individual samples
  mutate(total_TPM = sum(TPM)) %>%                  # sums TPM within a Category while keeping individual node_id TPM
  ungroup() %>%
  group_by(Collection_Type,Category) %>%  
  mutate(avg_TPM = mean(total_TPM)) %>%             # average TPM among samples within a collection type
  ungroup()

# mark rare categories as "other"
others = tpm %>%
  group_by(Category) %>%
  summarize(pool = max(avg_TPM) < 70, 
            max = max(avg_TPM),
            .groups="drop")

# don't pool the key genes into "Other" even though they may be small
# others = others %>%
#   mutate(pool = ifelse(Category == "blaKPC", FALSE, pool)) %>%
#   mutate(pool = ifelse(Category == "blaVIM", FALSE, pool)) %>%
#   mutate(pool = ifelse(Category == "blaOXA-48_fam", FALSE, pool)) %>%
#   mutate(pool = ifelse(Category == "blaDIM-SIM-IMP", FALSE, pool)) %>%
#   mutate(pool = ifelse(Category == "blaNDM", FALSE, pool)) 

# pool rare categories into Other by summing their avg_TPM
tf = inner_join(tpm, others, by="Category") %>%
  mutate(Category = if_else(pool, "Other", Category)) %>%
  group_by(Collection_Type, Samples, Category) %>%                 
  mutate(total_TPM2 = sum(TPM)) %>%     
  ungroup() %>%                                                # sum TPM for all genes within a Category (again, this time including Other Category)
  group_by(Collection_Type, Category) %>%
  mutate(avg_TPM2 = mean(total_TPM2)) %>%                         # average TPM for each Category among all samples within a Collection Type (again, this time including Other Category)
  ungroup()

# find the max avg_TPM across all samples for each Category  
cm = tf %>% 
  group_by(Category) %>%
  summarize(maxC = max(avg_TPM2), .groups="drop")
tf = inner_join(tf, cm, by="Category")

# how many?
N = sum(others$pool == FALSE)
N

# make color palette
# library(hues)
# pal = iwanthue(N)
# names(pal) = others$Category[others$pool == FALSE]
# # assign all bla genes to same color as BETA-LACTAM
# pal["blaKPC"] = pal["BETA-LACTAM"]
# pal["blaVIM"] = pal["BETA-LACTAM"]
# #pal["blaOXA-48_fam"] = pal["BETA-LACTAM"]
# #pal["blaDIM_SIM_IMP"] = pal["BETA-LACTAM"]
# #pal["blaNDM"] = pal["BETA-LACTAM"]
# pal["Other"] = "darkgray"
# save(pal, file="pal_AMR.txt")
# load(file="pal_AMR.txt")

# alternative color scale = pal2
library(RColorBrewer)
pal2 = brewer.pal(N, "Paired")
names(pal2) = others$Category[others$pool == FALSE]
# assign all bla genes to same color as BETA-LACTAM
pal2["blaKPC"] = pal2["BETA-LACTAM"]
pal2["blaVIM"] = pal2["BETA-LACTAM"]
#pal2["blaOXA-48_fam"] = pal2["BETA-LACTAM"]
#pal2["blaDIM_SIM_IMP"] = pal2["BETA-LACTAM"]
#pal2["blaNDM"] = pal2["BETA-LACTAM"]
pal2["Other"] = "darkgray"

# assign border colors to key genes
tf1 = tf %>%
  mutate(bc = factor(vector("character", length(tf$node_id)))) %>%
  mutate(bc = if_else(node_id == "blaKPC", "black", bc)) %>%
  mutate(bc = if_else(node_id == "blaVIM", "black", bc))

# show labels and borders only for the most abundant genes
tf1 = tf1 %>%
  mutate(l = factor(vector("character", length(tf$node_id)))) %>%
  mutate(l = if_else(TPM > 50, node_id, l)) %>%
  mutate(bc = if_else(TPM > 50, "black", bc))

# count number of genes in each collection type
tf1 = tf1 %>%
  #count(node_id, Collection_Type) %>%
  group_by(Collection_Type) %>%
  mutate(genes = sum(TPM > 0),
         CollectionTypeTotalTPM = sum(TPM)) %>%
  ungroup()
plotmax = max(tf1$CollectionTypeTotalTPM)

# plot
p = tf1 %>%
  mutate(Category = factor(Category),
        Category = fct_reorder(Category, maxC)) %>%
        # Category = fct_relevel(Category, "blaVIM", after=Inf),
        # Category = fct_relevel(Category, "blaKPC", after=Inf)) %>%
  ggplot(aes(x=Collection_Type, y=TPM, fill=Category, color = bc=="black")) +
  geom_col() +
  scale_fill_manual(values=pal2) +
  scale_color_manual(values=c(NA,"black"), guide="none", na.value="transparent") +
  geom_text(aes(label = l), size=4, position=position_stack(vjust=0.5)) +
  geom_text(aes(y=CollectionTypeTotalTPM, label=genes, fontface=2), vjust=-0.5) +
  annotate(geom="text", label="Number of Genes Detected", size=4, fontface=2, x=Inf, y=Inf) +
  #scale_x_discrete(breaks=c("Autosampler1","Autosampler2","Biofilm","Grab","Moore")) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black", face= "bold"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) +
  theme(plot.margin = margin(t = 20, unit = "pt")) +
  coord_cartesian(clip = "off")
p
ggsave("AMRbarplot-node_ids.png", p)
ggsave("AMRbarplot-node_ids.jpg", p)
  
# BETA-LACTAMASE only
# assign border colors to key genes
tf2 = tf %>%
  mutate(bc = factor(vector("character", length(tf$node_id)))) %>%
  mutate(bc = if_else(node_id == "blaKPC", "black", bc)) %>%
  mutate(bc = if_else(node_id == "blaVIM", "black", bc))

# also show labels for the most abundant genes
tf2 = tf2 %>%
  mutate(l = factor(vector("character", length(tf2$node_id)))) %>%
  mutate(l = if_else(TPM > 45, node_id, l)) %>%
  mutate(bc = if_else(TPM > 45, "black", bc)) %>%
  filter(Category=="BETA-LACTAM")
         #| Category=="blaKPC" | Category=="blaVIM")

# count number of genes in each collection type for BETA-LACTAM only
tf2 = tf2 %>%
  #count(node_id, Collection_Type) %>%
  group_by(Collection_Type) %>%
  mutate(genes = sum(TPM > 0),
         CollectionTypeTotalTPM = sum(TPM)) %>%
  ungroup()
plotmax = max(tf2$CollectionTypeTotalTPM)

p2 = tf2 %>%
  mutate(Category = fct_reorder(Category, maxC)) %>%
        # Category = fct_relevel(Category, "blaKPC", after=Inf),
        # Category = fct_relevel(Category, "blaVIM", after=Inf)) %>%
  ggplot(aes(x=Collection_Type, y=TPM, fill=Category, color= bc=="black")) +
  geom_col() +
  scale_fill_manual(values=pal2) +
  scale_color_manual(values=c(NA,"black"), guide="none", na.value="gray") +
  geom_text(aes(label=l), position=position_stack(vjust=0.5)) +
  geom_text(aes(y=CollectionTypeTotalTPM, label=genes, colour=TRUE), vjust=-0.5) +
  annotate(geom="text", label="Number of Genes Detected", size=4, fontface=2, x=Inf, y=Inf) +
  #scale_x_discrete(breaks=c("Autosampler1","Autosampler2","Biofilm","Grab","Moore")) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black", face= "bold"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) +
  theme(plot.margin = margin(t = 20, unit = "pt")) +
  coord_cartesian(clip = "off")
p2
ggsave("BLAbarplot-node_ids.png", p2)
ggsave("BLAbarplot-node_ids.jpg", p2)

# VENN DIAGRAMS
# Construct list of genes for each Collection_Type
vl1 = tf1 %>%
  select(TPM,node_id,Collection_Type) %>%
  filter(TPM> 0) %>%
  {split(.$node_id,.$Collection_Type)}

library(ggVennDiagram)
# Venn diagram with all Collection Types
vp1.1 = ggVennDiagram(vl1, label="count") + scale_x_continuous(expand = expansion(mult = .2))
vp1.1
# Venn diagram comparing only Autosampler1, Biofilm, and Moore
vp1.2 = ggVennDiagram(vl1[c(1,3,5)], label="count") + scale_x_continuous(expand = expansion(mult = .2))
vp1.2

# Venn diagram comparing Collection_Types BETA-LACTAM only
vl2 = tf2 %>%
  select(TPM,node_id,Collection_Type) %>%
  filter(TPM> 0) %>%
  {split(.$node_id,.$Collection_Type)}

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
ggsave("AMR-Venn-all-node_ids.png",vp1.1)
ggsave("AMR-Venn-all-node_ids.jpg",vp1.1)
ggsave("AMR-Venn-all-3samples-node_ids.png",vp1.2)
ggsave("AMR-Venn-all-3samples-node_ids.jpg",vp1.2)
ggsave("AMR-Venn-all-BLA-node_ids.png",vp2.1)
ggsave("AMR-Venn-all-BLA-node_ids.jpg",vp2.1)
ggsave("AMR-Venn-all-3samples-BLA-node_ids.png",vp2.2)
ggsave("AMR-Venn-all-3samples-BLA-node_ids.jpg",vp2.2)


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
  #scale_y_discrete(limits = rev(levels(tf1$node_id)))
b1

# BUBBLE PLOT BLA only
b2 = tf2 %>% 
  arrange(desc(TPM)) %>% # arrange(desc()) sorts the bubbles so that smaller bubbles are plotted on top of larger bubbles
  ggplot(aes(x = Collection_Type, y = node_id)) + # l is the list of abundant genes created above
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
#scale_y_discrete(limits = rev(levels(tf1$node_id)))
b2

# save bubble plots
ggsave("AMR-bubbles-node_ids.jpg", b1)
ggsave("AMR-bubbles-BLA-node_ids.jpg", b2)

