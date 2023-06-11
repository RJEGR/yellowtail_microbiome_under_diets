

library(ANCOMBC)
library(phyloseq)
library(microbiome)
library(tidyverse)
library(readxl)

dir <- "/Users/cigom/Documents/ANALYSIS_JUREL/"


hindgut <- read_tsv(paste0(dir, "ANCOMBC_RES_Hindgut.tsv")) %>% mutate(Tissue = 'Hindgut')

foregut <- read_tsv(paste0(dir, "ANCOMBC_RES_foregut.tsv")) %>% mutate(Tissue = 'Foregut')

rbind(hindgut, foregut) %>% # distinct(taxon_id)
  group_by(wrap, group, Tissue) %>%
  tally() %>% view()

library(ggupset)

# treatL <- c("CONTROL", "1","2", "3")

rbind(hindgut, foregut) %>% 
  group_by(taxon_id, Tissue, group) %>%
  summarise(across(wrap, .fns = list), n = n()) %>% #view()
  ggplot(aes(x = wrap, fill = group)) +
  geom_bar(position = position_dodge(width = 1), width = 0.8) +
  geom_text(stat='count', aes(label = after_stat(count)), 
    position = position_dodge(width = 1), vjust = -0.5, family = "GillSans", size = 2) +
  scale_x_upset(order_by = "degree", reverse = T) +
  theme_classic(base_family = "GillSans") +
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
    combmatrix.panel.line.size = 0, base_family = "GillSans") +
  labs(x = '', y = 'Number of Families') +
  ggsci::scale_fill_aaas() +
  ggsci::scale_color_aaas() +
  # scale_fill_manual(name = "", values = c('#1f78b4', '#a6cee3')) +
  # scale_color_manual(name = "", values = c('#1f78b4', '#a6cee3')) +
  theme(legend.position = 'top') -> psave

# psave + facet_grid( Tissue ~ .)

# ggsave(psave, filename = 'DA_upset.png', path = dir,
#   width = 5,height = 3.7)



rbind(hindgut, foregut) %>% distinct(taxon_id) %>% pull() -> keepTaxa

obj %>% filter(Family %in% keepTaxa) %>% drop_na(Family) -> fam_top

hclust_df <- data.frame(fam_top[-1], row.names = fam_top$Family)
hclust_df <- hclust(dist(hclust_df), "complete")

tax_hclust <- hclust_df$labels[hclust_df$order]

lineage %>%
  select(Class, taxon_id) %>% distinct() %>%
  right_join(fam_top , by = c('taxon_id'='Family')) -> fam_top

raf <- function(x) x/sum(x) * 100

colNames <- fam_top %>% select_if(is.double) %>% names()

fam_top %>%
  mutate_at(colNames, raf) %>% select() %>% 
  pivot_longer(cols = colNames, values_to = 'ra', names_to = 'ID') %>%
  mutate(taxon_id = factor(taxon_id, levels = tax_hclust),
    ra = ifelse(ra == 0, NA, ra)) %>%
  left_join(sam) %>%
  mutate(Tissue = recode_factor(Tissue, IP = "Hindgut", 
    IA = "Foregut")) -> dataviz

dataviz %>%
  ggplot(aes(y = taxon_id, x = ID, fill = ra)) +
  geom_tile() +
  geom_tile(color = "black") +
  ggsci::scale_fill_material("light-green",  
    name="Relative\nAbundance\n", 
    na.value = 'white',
    limits = c(0,100),
    labels = scales::percent_format(scale = 1)) +
  ggh4x::facet_nested(Class ~ Tissue+Treatment,
    scales = "free", space = "free" , nest_line = F) +
  labs(x = NULL, y = NULL) + 
  theme_classic(base_family = "GillSans", base_size = 14) +
  guides(fill = guide_colorbar(barheight = unit(5, "in"), 
    barwidth = unit(0.5, "in"),
    ticks.colour = "black", 
    frame.colour = "black",
    label.theme = element_text(size = 14)) ) +
  theme(
    title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12),
    strip.text.y = element_text(
      angle = 0, 
      size = 14),
    strip.background = element_rect(
      colour = "black", 
      fill = "transparent",
      size = 0.4),
    panel.spacing = unit(0.007, "lines")) -> heatPlot

ggsave(heatPlot, filename = "ANCOMBC_heatmap.png", path = dir, 
  width = 16, height = 10, dpi = 300)

bar_df <- fam_top %>%
  pivot_longer(cols = colNames, values_to = 'ab', names_to = 'ID') %>%
  left_join(sam) %>% 
  group_by(Tissue, taxon_id) %>%
  mutate(ra = ab/sum(ab))  

# sanity check
# bar_df %>% group_by(Treatment, taxon_id) %>% summarise(sum(ra)) 


# Sanity check
fam_top %>% filter(grepl('Microtrichaceae', taxon_id)) %>% 
  pivot_longer(cols = colNames, values_to = 'ab', names_to = 'ID') %>%
  left_join(sam) 

library(NatParksPalettes)

bar_df %>%
  arrange(taxon_id) %>%
  ggplot() +
  geom_col(aes(x = reorder(taxon_id, desc(taxon_id)), 
    y = ra, fill = Treatment),
    position = position_fill(reverse = T)) +
  # facet_grid(Phylum ~ ., scales = "free", space = "free") +
  ggh4x::facet_nested(.  ~ Tissue, scales = "free", space = "free") +
  labs(x = NULL, y = NULL) +
  coord_flip() +
  theme_classic(base_size = 16, base_family = "GillSans") +
  theme(strip.background.y = element_blank(),
    axis.text.y.left = element_text(angle = 0, 
      hjust = 1, vjust = 0, size = 16),
    axis.text.y.right = element_text(angle = 0, 
      hjust = 1, vjust = 0, size = 6)) -> p2


values <- natparks.pals("Yosemite", 4, direction = 1) # "Yellowstone"

p2 + scale_fill_manual("", values=values)  -> p2

ggsave(p2, filename = "ANCOMBC_barplot.png", path = dir, 
  width = 8, height = 8)

# Step 2 (omit): adjust the log observed abundances by subtracting the estimated  sampling fraction from log observed abundances of each sample.

adj_ab <- function(out_ancombc, phyloseq) {
  samp_frac = out_ancombc$samp_frac
  # Replace NA with 0
  samp_frac[is.na(samp_frac)] = 0 
  
  # Add pesudo-count (1) to avoid taking the log of 0
  log_obs_abn = log(abundances(phyloseq) + 1) 
  # Adjust the log observed abundances
  log_obs_abn_adj = t(t(log_obs_abn) - samp_frac)
  
  log_obs_abn_adj %>% as_tibble(rownames = "Family") %>% pivot_longer(-Family, names_to = "Index")
}

adj_ab(out1,ancombc_data1) %>%
  right_join(adj_ab(out2, ancombc_data2), by = c("Family", "Index"))

# re-doing heatmap

ancombc_res <- rbind(hindgut, foregut) %>% arrange(taxon_id)

diff_ab_list <- ancombc_res %>% pull(taxon_id) %>% unique()

# intersected_taxa
phyloseq <- read_rds(paste0(dir, '/inputs_ancom_ps.rds'))

phyloseq %>%
  aggregate_taxa(., "Family") %>%
  # transform_sample_counts(function(x) sqrt(x / sum(x))) %>%
  subset_taxa(Family %in% diff_ab_list) %>%
  psmelt() %>%
  rename("taxon_id" = OTU, "ab" = Abundance) %>% 
  as_tibble() -> dataHeat

# sanity check

dataHeat %>%
  group_by(taxon_id, Tissue) %>%
  summarise(Prev = sum(ab > 0), Tab = sum(ab)) %>%
  # mutate(Tab = log(Tab + 1)) %>%
  # left_join(mtd) %>%
  left_join(lineage) %>%
  arrange(desc(Tab)) %>%
  group_by(Tissue) %>%
  mutate(Rank = rank(Tab)) -> dataViz

# labels <- dataViz %>% pull(Phylum) %>% unique() %>%
dataViz %>% group_by(Phylum) %>% summarise(t = sum(Tab)) %>% arrange(desc(t)) %>% pull(Phylum) -> labels 

colourCount = length(labels)

library(ggsci)

if(colourCount > 7) {
  getPalette <- colorRampPalette(pal_locuszoom(7))(colourCount)
} else
  getPalette <- pal_locuszoom()(colourCount)  



dataViz %>%
  ggplot(aes(Tab, Prev)) +
  geom_point(aes(color = Tissue)) +
  scale_x_log10() +
  geom_smooth(se = F) 

dataViz %>%
  arrange(taxon_id) %>%
  group_by(taxon_id) %>%
  mutate(ra = Tab / sum(Tab)) %>%
  mutate(Phylum = factor(Phylum,  levels = labels)) %>%
  mutate(taxon_id = forcats::fct_reorder(taxon_id, Rank)) %>%
  ggplot() +
  geom_col(aes(x = taxon_id ,
    y = ra, fill = Tissue)) +
  ggsci::scale_fill_rickandmorty() +
  # facet_grid(Phylum ~ ., scales = "free", space = "free") +
  ggh4x::facet_nested(Phylum  ~., scales = "free", space = "free", switch = "y") +
  labs(x = NULL, y = NULL) +
  coord_flip() +
  theme_classic(base_size = 16, base_family = "GillSans") +
  theme(strip.background.y = element_blank(),
    # strip.text = element_text(margin = margin(1, 1, 1, 1)),
    axis.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0, size = 16),
    axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 6)) -> p2


ggsave(p2, filename = "bar_ANCOM_list.png", path = dir, 
  width = 8, height = 8)
