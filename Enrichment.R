
# install.packages("devtools")
#devtools::install_github("davidsjoberg/ggsankey")
library(tidyverse)
library(dplyr)
library(ggsankey)
library(cols4all)
library(ggplot2)
library(ggforce)
library(randomcoloR)
# Create a list of selected metabolites (PATHWAY_SORTORDER) with its sub-pathway and super pathway based on PATHWAY_SORTORDER (metabolite ID)

#Output from SPLS model
merge_table
load("~/merge_table.Rdata")
load("~/metabolite.info.Rdata")
# Subset metabolite.info where only the selected metabolites are included

selected.metabolite.info <- metabolite.info %>% filter(CHEM_ID%in%merge_table$var)

# Count the number of metabolites invoved in all sub pathway and the n need to be at least >1
subpathway.number <- dplyr::count(metabolite.info, SUB_PATHWAY)
subpathway.number <-filter(subpathway.number,n>1)

metabolite.info<-metabolite.info %>% filter(SUB_PATHWAY %in% subpathway.number$SUB_PATHWAY)
selected.metabolite.info<-selected.metabolite.info %>% filter(SUB_PATHWAY %in% subpathway.number$SUB_PATHWAY)
selected.subpathway.list <- selected.metabolite.info %>% pull(SUB_PATHWAY) %>% unique()

# Count the number of metabolites invovled in the selected sub pathway
selected.subpathway.number <- dplyr::count(selected.metabolite.info, SUB_PATHWAY)
selected.subpathway.number <-filter(selected.subpathway.number,n>1)

#all
# Pathway Enrichment Analysis
enrichment <- selected.subpathway.list %>% data.frame() %>%
  dplyr::rename(SUB_PATHWAY=rlang::eval_tidy(quo("."))) %>%
  left_join(., selected.subpathway.number, by="SUB_PATHWAY") %>%
  dplyr::rename(v1=n) %>%
  left_join(., subpathway.number, by="SUB_PATHWAY") %>%
  dplyr::rename(v2=n) %>%
  mutate(n_N=(selected.metabolite.info %>% nrow())/(metabolite.info %>% nrow()),
         Fold_Enrichment=(v1/v2)/n_N,
         p.value=phyper(q=v1-1, m=v2, n=nrow(metabolite.info)-nrow(selected.metabolite.info),
                        k=nrow(selected.metabolite.info), lower.tail=F),
         p.value.fdr=p.adjust(p.value, "fdr"),
         p.value.bon=p.adjust(p.value, "bon"),
         p.value.size=case_when((p.value>0.05)~"p>0.05",
                                (p.value<0.05 & p.value.fdr>0.05)~"p<0.05",
                                (p.value.fdr<0.05& p.value.bon>0.05)~"p.fdr<0.05",
                                (p.value.bon<0.05)~"p.bon<0.05"),
         p.value.size.num=case_when((p.value>0.05)~1,
                                    (p.value<0.05 & p.value.fdr>0.05)~4,
                                    (p.value.fdr<0.05& p.value.bon>0.05)~6,
                                    (p.value.bon<0.05)~8))
enrichment<-na.omit(enrichment)
g <- enrichment %>%
  mutate(SUB_PATHWAY = fct_reorder(SUB_PATHWAY, desc(Fold_Enrichment))) %>%
  ggplot(aes(x=SUB_PATHWAY, y=log2(Fold_Enrichment), color=-log2(p.value))) +
  geom_point(aes(size=p.value.size)) +
  scale_size_manual(values=sort(unique(enrichment$p.value.size.num), decreasing=T)) +
  scale_colour_gradient2(low="yellow", mid="green", high="red", midpoint=0) +
  geom_hline(yintercept=0, lty=4) +
  coord_flip() +
  labs(x="Fold Enrichment-NDD", color="-log2(p-value)", size="p-value") +
  theme_bw()
g +
  theme_bw(base_size = 15)+
  theme(strip.text.y = element_blank())+
  theme(panel.grid=element_blank())


#FDR multple testing

enrichment <- selected.subpathway.list %>% 
  data.frame() %>%
  dplyr::rename(SUB_PATHWAY = rlang::eval_tidy(quo("."))) %>%
  left_join(., selected.subpathway.number, by = "SUB_PATHWAY") %>%
  dplyr::rename(v1 = n) %>%
  left_join(., subpathway.number, by = "SUB_PATHWAY") %>%
  dplyr::rename(v2 = n) %>%
  mutate(
    n_N = (selected.metabolite.info %>% nrow()) / (metabolite.info %>% nrow()),
    Fold_Enrichment = (v1 / v2) / n_N,
    p.value = phyper(
      q = v1 - 1, 
      m = v2, 
      n = nrow(metabolite.info) - nrow(selected.metabolite.info),
      k = nrow(selected.metabolite.info), 
      lower.tail = F
    ),
    p.value.fdr = p.adjust(p.value, "fdr"),
    p.value.size = case_when(
      (p.value > 0.05) ~ "p>0.05",
      (p.value < 0.05 & p.value.fdr > 0.05) ~ "p<0.05",
      (p.value.fdr < 0.05) ~ "p.fdr<0.05"
    ),
    p.value.size.num = case_when(
      (p.value > 0.05) ~ 1,
      (p.value < 0.05 & p.value.fdr > 0.05) ~ 4,
      (p.value.fdr < 0.05) ~ 6
    )
  )
enrichment<-na.omit(enrichment)
g <- enrichment %>%
  mutate(SUB_PATHWAY = fct_reorder(SUB_PATHWAY, desc(Fold_Enrichment))) %>%
  ggplot(aes(x=SUB_PATHWAY, y=log2(Fold_Enrichment), color=-log2(p.value))) +
  geom_point(aes(size=p.value.size)) +
  scale_size_manual(values=sort(unique(enrichment$p.value.size.num), decreasing=T)) +
  scale_colour_gradient2(low="yellow", mid="green", high="red", midpoint=0) +
  geom_hline(yintercept=0, lty=4) +
  coord_flip() +
  labs(x="Fold Enrichment-NDD", color="-log2(p-value)", size="p-value") +
  theme_bw()
g +
  theme_bw(base_size = 15)+
  theme(strip.text.y = element_blank())+
  theme(panel.grid=element_blank())

#NDD


# Assume df2 is your second data frame
# This creates a new data frame with 'SUB_PATHWAY' and concatenated 'CHEMICAL_NAME'
chemicals_by_pathway <- selected.metabolite.info %>%
  group_by(SUB_PATHWAY) %>%
  summarise(CHEMICALS = paste(CHEMICAL_NAME, collapse = "@")) %>%
  ungroup()
enrichment <- enrichment %>%
  left_join(chemicals_by_pathway, by = "SUB_PATHWAY")
enrichment <- enrichment[order(enrichment$Fold_Enrichment, decreasing = T),]
enrichment$SUB_PATHWAY <- factor(enrichment$SUB_PATHWAY, levels = enrichment$SUB_PATHWAY)

# Assuming enrichment$CHEMICALS contains the metabolites separated by '/'
# Split CHEMICALS into a list of metabolites
enrichment$metabolite_ID <- strsplit(enrichment$CHEMICALS, "@")

# Calculate the number of metabolites per Subpathway
enrichment$num_metabolites <- sapply(enrichment$metabolite_ID, length)

# Create a new sankey_data dataframe
# Calculate the total number of metabolites to set the correct size for sankey_data
# Calculate the total number of metabolites across all subpathways
total_metabolites <- sum(enrichment$num_metabolites)

# Create sankey_data with enough rows to accommodate all metabolites
sankey_data <- data.frame(
  Subpathway = rep(enrichment$SUB_PATHWAY, times=enrichment$num_metabolites),
  Metabolites = vector("character", length = total_metabolites)
)

# Initialize the index for filling sankey_data
n <- 1

# Loop through each Subpathway
for (i in seq_along(enrichment$SUB_PATHWAY)) {
  num_metabolites <- enrichment$num_metabolites[i]
  if (num_metabolites > 0) {
    # Sample metabolites based on the number available
    sampled_metabolites <- sample(enrichment$metabolite_ID[[i]], num_metabolites)
    end_index <- n + num_metabolites - 1
    sankey_data$Metabolites[n:end_index] <- sampled_metabolites
    n <- end_index + 1  # Move index to the next start position
  }
}

metabolites_select <- sankey_data$Metabolites

#Convert to the formate to sankey plot：
sankey_data <- sankey_data %>%
  make_long(Metabolites,Subpathway)

# Mofify the order of the nodes：
sankey_data$node <- factor(sankey_data$node,
                           levels = c(as.character(enrichment$SUB_PATHWAY) %>% unique(),
                                      metabolites_select %>% unique()))

#mycol <- c4a('rainbow_wh_rd', length(unique(sankey_data$node)))
# Custom color palette

custom_colors <- c(
  "#458A74", "#018B38", "#D9A421", "#F5A216", "#57AF37", 
  "#41B9C1", "#008B8B", "#4E5689", "#6A8EC9", "#652884", 
  "#8A7355", "#CC5B45", "#848484", "#E42320", "#B46DA9"
)

# Number of unique nodes in your Sankey diagram
num_unique_nodes <- length(unique(sankey_data$node))

# Generate additional distinct colors if necessary
if (num_unique_nodes > length(custom_colors)) {
  additional_colors_needed <- num_unique_nodes - length(custom_colors)
  additional_colors <- distinctColorPalette(additional_colors_needed)
  mycol <- c(custom_colors, additional_colors)
} else {
  mycol <- custom_colors[1:num_unique_nodes]
}

# Assign the custom colors to your nodes
sankey_data <- sankey_data %>%
  mutate(color = mycol[match(node, unique(node))])

# Create the Sankey diagram with custom colors
p1 <- ggplot(sankey_data, aes(x = x,
                              next_x = next_x,
                              node = node,
                              next_node = next_node,
                              fill = I(color),  # Use the custom color column with I() to interpret colors
                              label = node)) +
  geom_sankey(flow.alpha = 0.5,
              flow.fill = 'grey',
              smooth = 8,
              width = 0.1) +
  geom_sankey_text(size = 5, hjust = 1.2,  # Horizontal alignment: 0=left, 0.5=center, 1=right
                   vjust = 0.5,
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none')

print(p1)


####### bubble plot--------------
go_data <- enrichment %>%
  mutate(ymax = cumsum(num_metabolites)) %>% 
  mutate(ymin = ymax - num_metabolites) %>%
  mutate(label = (ymin + ymax)/2) 
head(go_data)

p2 <- ggplot(go_data) +
  geom_point(aes(x = log2(go_data$Fold_Enrichment),
                 y = label, 
                 size = p.value.size,
                 color = -log2(p.value))) +
  scale_size_continuous(range=c(1, 4)) + 
  scale_size_manual(values=sort(unique(go_data$p.value.size.num), decreasing=T)) +
  scale_colour_distiller(palette = "Reds", direction = 1) + 
  scale_y_continuous(expand = expansion(mult = 0, add = c(0.8, 0.3)))+
  labs(x = "-log2(Fold_Enrichment)",
       y = "") +
  theme_bw()+
  theme(axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_blank(),
        text = element_text(size = 10))

library(cowplot)

# Figure 6
p <- ggdraw() +
  
  draw_plot(p1, 0, 0, .8, 1)+
  draw_plot(p2, 0.6, 0.21, .21, 0.525)
p


