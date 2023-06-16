
library(wTO)
library(tidyverse)

library(igraph)
library(tidygraph)
library(ggraph)




rm(list = ls())

dir <- "~/Documents/ANALYSIS_JUREL/"

ps <- read_rds(paste0(dir, '/inputs_ancom_ps.rds'))

class(ps)

head(TAX <- as(tax_table(ps), "matrix"))

agg_level <- "Family"

ps <- microbiome::aggregate_taxa(ps, level = agg_level)

colnames(TAX)[colnames(TAX) %in% agg_level] <- "unique"

prep_pyseq <- function(ps) {
  
  ps %>% 
    prune_taxa(taxa_sums(.) > 0, .) %>% 
    prune_samples(sample_sums(.) > 0, .) %>%
    transform_sample_counts(., function(x) sqrt(x / sum(x))) %>%
    otu_table() %>%
    as.data.frame()
}


# overlap may be query.taxa observed during ANCOM test

# overlap <- head(TAX[,4]) # Example

from_pyseq_to_wTO <- function(ps, overlap = NULL, n = 100) {
  
  require(wTO)
  require(phyloseq)

  
  metagenomics_abundance <- prep_pyseq(ps)
  
  # Overlap: Set of nodes of interest, where the Overlapping weights will be computed.
  
  if (is.null(overlap)) {
    wTO <- wTO.fast(Data = metagenomics_abundance,
      Overlap = row.names(metagenomics_abundance),
      method = 's', sign = 'sign', n = n, 
      method_resampling = 'Bootstrap')
    
  } else
    
    wTO <- wTO.fast(Data = metagenomics_abundance,
      Overlap = overlap,
      method = 's', sign = 'sign', n = n, 
      method_resampling = 'Bootstrap')
  
  return(wTO)
  
}

# TEST TIME-SERIES ====

table(sample_data(ps)$Time)

dim(Data1 <- ps %>% subset_samples(Time == 0) %>% prep_pyseq())
dim(Data2 <- ps %>% subset_samples(Time == 30) %>% prep_pyseq())
dim(Data3 <- ps %>% subset_samples(Time == 60) %>% prep_pyseq())

Net1 <- wTO.fast(Data = Data1,
  Overlap = rownames(Data1),
  method = 's', sign = 'sign', n = 50,
  method_resampling = 'BlockBootstrap', lag = 2)

Net2 <- wTO.fast(Data = Data2,
  Overlap = rownames(Data2),
  method = 's', sign = 'sign', n = 50,
  method_resampling = 'BlockBootstrap', lag = 2)

Net3 <- wTO.fast(Data = Data3,
  Overlap = rownames(Data3),
  method = 's', sign = 'sign', n = 50,
  method_resampling = 'BlockBootstrap', lag = 2)

# ID = sample_data(ps2)$Time
# method_resampling = 'BlockBootstrap', lag = 2 or ID

# save(Net1, Net2, Net3, file = paste0(dir, "MakeDiffNet.RData"))

load(paste0(dir, "MakeDiffNet.RData"))

library(CoDiNA)

DiffNet <- MakeDiffNet(Data = list(Net1,Net2,Net3), Code = c('Net1', 'Net2', 'Net3') )

table(DiffNet$Group)

# RUN NETWORK ====

# Split by tissue

ps <- ps %>% subset_samples(Time != 60)

table(sample_data(ps)$Tissue)

ps1 <- ps %>% subset_samples(Tissue == "Foregut")

ps2 <- ps %>% subset_samples(Tissue == "Hindgut")

wTO1 <- ps1 %>% from_pyseq_to_wTO(n = 50)

wTO2 <- ps2 %>% from_pyseq_to_wTO(n = 50)

DiffNet <- CoDiNA::MakeDiffNet(Data = list(wTO1,wTO2), Code = c('Foregut', 'Hindgut') )

plot(DiffNet)



wTO1 <- wTO1 %>% mutate(Tissue = "Foregut")
wTO2 <- wTO2 %>% mutate(Tissue = "Hindgut")

WTO <- rbind(wTO1, wTO2 )

# write_csv(WTO, file = paste0(getwd(), "/weighted_topological_overlap_net.csv"))
# WTO <- read_csv(paste0(getwd(), "/weighted_topological_overlap_net.csv"))

# wTOout <- ps %>% subset_samples(Tissue == "Foregut") %>% from_pyseq_to_wTO(n = 50)

head(WTO)

WTO %>%
  mutate(y = -log10(pval.adj-abs(wTO))) %>%
  ggplot(aes(color =  pval.adj)) + 
  facet_grid(~ Tissue) +
  geom_point(aes(x = wTO, y = y), alpha = 0.5) +
  ggsci::scale_color_gsea(reverse = T) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = 'Correlation (Tau)',   
    y = expression(~Log[10]~('p.adj'))) +
  theme_bw(base_size = 16, base_family = "GillSans")

WTO %>%
  filter(pval.adj < 0.05) %>%
  mutate(y = -log10(pval.adj-abs(wTO))) %>%
  drop_na(y) %>% mutate(g = sign(wTO)) %>%
  group_by(Tissue) %>%
  count(g) 

# Most biological networks show a power-law degree distribution, where a few nodes have a very large number of connections, while other nodes have no or few connections

# Positive correlations in co-occurence networks may represent symbiotic or commensal relationships, while negative correlations may represent predator-prey interactions, allelopathy or competition for limited resources.

WTO %>%
  mutate(y = abs(wTO)) %>% 
  ggplot(aes(y)) + 
  geom_histogram() +
  facet_grid(~ Tissue) +
  geom_vline(xintercept = 0.2, linetype = "dashed") +
  theme_bw(base_size = 16, base_family = "GillSans")

WTO %>%
  mutate(y = pval.adj) %>% 
  ggplot(aes(y)) + 
  geom_histogram() +
  facet_grid(~ Tissue) +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  theme_bw(base_size = 16, base_family = "GillSans")

# Refinamiento de la red en base a la significancia estadistica:

WTO <- WTO %>% filter(pval.adj < 0.05)

# .WTO <- WTO %>% filter(abs(wTO) >= 0.2)

Node.1 = as.character(.WTO$Node.1)

Node.2 = as.character(.WTO$Node.2)

length(sort(unique(c(Node.1,Node.2))))

# summary(wTOout$wTO)

# summary(wTOout$pval.adj)

WTO <- WTO %>% filter(Tissue == "Foregut")

wTOcutoff <- function(WTO, cutoff = 0.01, tau = 0.5) {
  
  n <- length(unique(c(WTO$Node.1, WTO$Node.2)))
  
  WTO %>% mutate_if(is.factor, as.character) %>%
    mutate(wTO = ifelse(pval-abs(wTO) < cutoff, wTO, 0 )) %>%
    filter(wTO != 0 ) %>%
    filter(abs(wTO) > tau) %>%
    # as.data.frame() -> out
    as_tibble() -> out
  
  cat("Dimension of df:", nrow(out), "\n")
  cat("Number of significant nodes interacting: ", length(unique(c(out$Node.1, out$Node.2))))
  cat("\nProportion of : ", length(unique(c(out$Node.1, out$Node.2)))/ n, "\n")
  
  return(out)
  
}

exportNet <- function(WTO, infoNodes = TAX, cutoff = 0.05, tau = 0.5) {
  
  WTO %>% wTOcutoff(cutoff, tau) -> WTO
  
  Node.1 = as.character(WTO$Node.1)
  Node.2 = as.character(WTO$Node.2)
  
  nodes <- data.frame(unique = sort(unique(c(Node.1,Node.2))))
  
  infoNodes <- infoNodes %>% as_tibble()
  
  # infoNodes <- infoNodes %>% filter(unique %in% nodes$unique)
  
  nodes %>% mutate_if(is.character, as.factor) -> nodes
  
  # ps %>% tax_table() %>% as.data.frame() -> infoNodes
  
  # ps %>%
  #   transform_sample_counts(., function(x) x / sum(x)) %>%
  #   otu_table() %>% 
  #   as.data.frame() %>% cbind(unique = rownames(.)) %>%
  #   left_join(infoNodes, by = 'unique') %>% 
  #   mutate_if(is.character, as.factor) %>% as_tibble() -> infoNodes
  
  # sum(infoNodes$unique %in% nodes$unique)
  
  nodes %>%
    as_tibble() %>%
    left_join(infoNodes, by = c('unique')) -> nodes
  
  
  WTO %>%
    mutate_if(is.character, as.factor) %>%
    mutate(typeEdge = sign(wTO)) -> edges
  
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = F)
  
  graph %>% activate(nodes) %>%
    mutate(degree = centrality_degree()) %>%
    mutate(value = (degree - min(degree))/(max(degree) - min(degree))) %>%
    mutate(weight = value * 2 + 1) -> graph
  
  gr <- function(x) {ifelse(x > 0, '+', '-')}
  
  graph %>% activate("edges") %>%  mutate(color = gr(wTO)) -> graph
  
  
  return(graph)
  
}

exportNet(WTO, cutoff = 0, tau = 0) -> graph

class(graph)

# igraph::components(graph)

# viz
# igraph centrality measures

graph %>% 
  activate("nodes") %>%
  mutate(betweenness = betweenness(.),
    pageRank = page_rank(.)$vector) %>% 
  as_data_frame(., "vertices") %>%
  mutate(color = ifelse(unique == "Nitrincolaceae", "Nitrincolaceae", "None")) %>%
  select(unique, color, degree, betweenness, pageRank ) %>%
  pivot_longer(cols = c('degree', 'betweenness')) %>%
  ggplot(aes(value, pageRank, color = color)) + 
  geom_point() +
  facet_grid(~name, scales = 'free')


# FILTERING NETWORK =====

graph %>%
  activate("edges") %>% 
  mutate(wTO = ifelse(abs(wTO) < 0.2, 0, wTO)) %>%
  filter(abs(wTO) > 0) -> g
  
g %>% activate("nodes") %>% 
    mutate(degree = degree(.)) %>% 
  filter(degree > 0) %>%
  mutate(betweenness = betweenness(.),
    pageRank = page_rank(.)$vector) %>%
  mutate(unique = ifelse(pageRank > 0.006, unique, NA)) -> g

# g %>% as_data_frame(., "vertices") %>% drop_na(unique)

names <- g %>% activate("nodes") %>% distinct(Phylum) %>% pull()

cols <- ggpubr::get_palette(palette = 'default', length(names))

cols = structure(cols, names = names)

str(cols)

# Hub detection as Layeghifard, M.,et al 2019
# top-ranked network hub taxa:

# Ex. We applied the PageRank algorithm to the microbiome networks to identify key members of CF lung microbiome in each patient. PageRank is a link analysis algorithm with the underlying assumption that hubs are likely to be more connected to other nodes when compared to non-hub nodes. This algorithm assigns a numerical weight to each node of a network as a measure of its relative importance within the network. The numerical weights of the nodes (also called PageRanks of the nodes) were then sorted and the top five or ten taxa with highest PageRank were selected as microbiome network hubs (or key taxa). For the network hub taxa, we recorded both the top five/top ten taxa for further analysis. 
  
# Finding community structure 

g %>% 
    activate("nodes") %>% 
    mutate(betweenness = betweenness(.),
      membership = igraph::cluster_louvain(.)$membership,
      pageRank = page_rank(.)$vector) -> g
  
  
table(igraph::cluster_louvain(g)$membership)

g %>% activate("edges") %>%
  mutate(betweenness = edge.betweenness(.)) -> g

# PREPARE LAYOUT TO DATAVIZ

layout = create_layout(g, layout = 'igraph', algorithm = 'kk')

ggraph(layout) +
  ggforce::geom_mark_hull(
    aes(x, y, group = as.factor(membership), 
      fill = as.factor(membership)), 
    color = NA,
    concavity = 4,
    con.size = 0.3,
    con.linetype = 2,
    expand = unit(2, "mm"),
    alpha = 0.25)  +
  guides(fill = FALSE) -> psaveNet


ggraph(layout) +
  geom_edge_link(aes(color = wTO, alpha = abs(wTO)), # edge_width = betweenness
    arrow = arrow(
      angle = 10,
      length = unit(0.1, "inches"),
      ends = "last",
      type = "closed")) +
  geom_node_point(aes(size = degree, color = Phylum)) +
  geom_node_text(aes(label = unique), repel = T) + 
  coord_fixed() +
  theme_graph(base_family = "GillSans") +
  theme(legend.position = "top") +
  # scale_edge_colour_manual(values = cols, guide = "none") +
  scale_color_manual(values = cols, name = '')

# Compares community structures using various metrics =====
# normalized mutual information
# We used normalized mutual information (NMI, Supplementary Note 2) as a measure of shared information between the detected communities and the true clusters. NMI ranges between 0 to 1, where a higher value indicates higher concordance and 1 indicates an exact match between the detected communities and the true clusters.

exportNet(wTO1, cutoff = 0.05, tau = 0.2) -> g1
exportNet(wTO2, cutoff = 0.05, tau = 0.2) -> g2

# compare( wTO1, wTO2, method = "nmi")

# 2
# Also separete by time-series as independent dataset and then use wTO.complete to search for consensus and difenciated interactions
