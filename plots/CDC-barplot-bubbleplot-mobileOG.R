# Sincerely grateful to Pat Schloss and his many wonderful tutorials at https://riffomonas.org

library(tidyverse)

d = read_csv("merged.abunds.mobileOG.cats.tpm.csv")
metadata = read_csv("metadata4.csv")

types = d %>% select(`mobileOG ID`,`Gene Name`,`Major mobileOG Category`,`Minor mobileOG Category`,`Source Database`,`Evidence Type`)

# traditional transpose matrix
dt = d %>% select(-mobileOG_Result_x,-mobileOG_Result_y,-`Gene Name`,-`Best Hit Accession ID`,-`Major mobileOG Category`,-`Minor mobileOG Category`,-`Source Database`,-`Evidence Type`) %>% 
  pivot_longer(cols=-`mobileOG ID`,) %>% 
  pivot_wider(names_from=`mobileOG ID`,values_from=value) %>%
  rename(Samples = name)

# convert transposed matrix to long format
dl = dt %>% 
  pivot_longer(-Samples, names_to="mobileOG ID", values_to="TPM")

# join transposed and pivoted counts with metadata
dj = inner_join(dl, metadata, by="Samples") %>%
  inner_join(., types, by="mobileOG ID") %>%
  pivot_longer(c("Major mobileOG Category","Minor mobileOG Category","Source Database","Evidence Type"), names_to="level", values_to="Category") %>%
  mutate(Collection_Type = factor(Collection_Type, levels=c("Autosampler1","Autosampler2","Biofilm","Grab","Moore")))

# pool TPMs within categories and average TPM among samples within a type 
tpm = dj %>%
  filter(level=="Major mobileOG Category") %>%
  group_by(Collection_Type,Samples,Category) %>%    # delete "Collection_Type" if you want to show individual samples
  mutate(total_TPM = sum(TPM)) %>%                  # sums TPM within a Category while keeping individual node_id TPM
  ungroup() %>%
  group_by(Collection_Type,Category) %>%  
  mutate(avg_TPM = mean(total_TPM)) %>%             # average TPM among samples within a collection type
  ungroup()

# mark rare categories as "other"
others = tpm %>%
  group_by(Category) %>%
  summarize(pool = max(avg_TPM) < 10, 
            max = max(avg_TPM),
            .groups="drop")

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

# make palette
library(RColorBrewer)
pal2 = brewer.pal(N, "Dark2")
names(pal2) = others$Category[others$pool == FALSE]
pal2["Other"] = "darkgray"

# show labels and borders only for the most abundant genes
tf1 = tf %>%
  mutate(l = factor(vector("character", length(tf$`mobileOG ID`)))) %>%
  mutate(l = if_else(TPM > 500, `Gene Name`, l)) %>%
  mutate(bc = factor(vector("character", length(tf$`mobileOG ID`)))) %>%
  mutate(bc = if_else(TPM > 500, "black", bc))

# count number of genes in each sample type
tf2 = tf1 %>%
  group_by(Collection_Type) %>%
  mutate(genes = sum(TPM > 0),
         CollectionTypeTotalTPM = sum(TPM)) %>%
  ungroup()
plotmax = max(tf2$CollectionTypeTotalTPM)

# plot
p = tf2 %>%
  mutate(Category = factor(Category),
        Category = fct_reorder(Category, maxC)) %>%
  ggplot(aes(x=Collection_Type, y=TPM, fill=Category, color = bc=="black")) +
  geom_col() +
  scale_fill_manual(values=pal2) +
  scale_color_manual(values=c(NA,"black"), guide="none", na.value="transparent") +
  geom_text(aes(label = l), size=4, position=position_stack(vjust=0.5)) +
  geom_text(aes(y=CollectionTypeTotalTPM, label=genes, fontface=2), vjust=-0.5) +
  annotate(geom="text", label="Number of Genes Detected", size=4, fontface=2, x=3, y=Inf) +
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
ggsave("mobileOG-barplot.jpg", p)

# alternative without gene labels
alt = tf %>%
  group_by(Collection_Type) %>%
  mutate(genes = sum(TPM > 0),
         CollectionTypeTotalTPM = sum(TPM)) %>%
  ungroup()
plotmax = max(alt$CollectionTypeTotalTPM)

palt = alt %>%
  mutate(Category = factor(Category),
         Category = fct_reorder(Category, maxC)) %>%
  ggplot(aes(x=Collection_Type, y=avg_TPM2, fill=Category)) +
  geom_col() +
  scale_fill_manual(values=pal2) +
  #scale_color_manual(values=c(NA,"black"), guide="none", na.value="transparent") +
  #geom_text(aes(label = l), size=4, position=position_stack(vjust=0.5)) +
  geom_text(aes(y=CollectionTypeTotalTPM, label=genes, fontface=2), vjust=-0.5) +
  annotate(geom="text", label="Number of Genes Detected", size=4, fontface=2, x=3, y=Inf) +
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
palt
ggsave("mobileOG-barplot-alt.jpg", palt)


# VENN DIAGRAMS
# Construct list of genes for each Collection_Type
vl1 = tf2 %>%
  select(TPM,`mobileOG ID`,Collection_Type) %>%
  filter(TPM> 0) %>%
  {split(.$`mobileOG ID`,.$Collection_Type)}

library(ggVennDiagram)
# Venn diagram with all Sample Types
vp1.1 = ggVennDiagram(vl1, label="count") + scale_x_continuous(expand = expansion(mult = .2))
vp1.1
# Venn diagram comparing only Autosampler1, Biofilm, and Moore
vp1.2 = ggVennDiagram(vl1[c(1,3,5)], label="count") + scale_x_continuous(expand = expansion(mult = .2))
vp1.2

# save Venn diagrams
ggsave("mobileOG-Venn.jpg",vp1.1)
ggsave("mobileOG-Venn-3samples.jpg",vp1.2)

# get genes unique to one Sample_Type
everything_else = union(vl1$"Autosampler1",vl1$"Autosampler2")
everything_else = union(everything_else,vl1$"Moore")
everything_else = union(everything_else,vl1$"Grab")
u = setdiff(vl1$Biofilm,everything_else)
u

# haven't finished this code yet
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

