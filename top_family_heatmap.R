
# Barplots were described  At the class level (figure 2 paper)

rm(list = ls())

dir <- "/Users/cigom/Documents/ANALYSIS_JUREL"

files <- list.files(path = dir, 
  pattern = "Tabla abundancias familias tesis Rocio Valenzuela", full.names = TRUE)

tax_f <- list.files(path = dir, pattern = "Taxonomy", full.names = TRUE)

library(readxl)

library(tidyverse)

library(NatParksPalettes)

obj <- read_xlsx(files, sheet = 1, skip = 0, col_names = T)

mtd <- read_xlsx(files, sheet = 2, skip = 0, col_names = T) %>% mutate(ID = as.character(ID))

# 
write_rds(tax, file = paste0(dir, '/lineage.rds'))

# no contamos con todos los linajes, o los linajes estan salteados y rellenados en familia

sum(unique(obj$Family) %in% unique(tax$Family))

unique(obj$Family)[!grepl('eae', unique(obj$Family))]

sum(!grepl('eae', unique(obj$Family)))


# prepare data to select top

obj %>% select_if(is.double) %>% names() -> colNames


agglom_lev <- "Family"

which_vars <- c(agglom_lev, colNames)
# obj %>% inner_join(tax) %>% distinct(Family, .keep_all = T) -> obj

obj %>%
  select_at(vars(all_of(which_vars))) %>%
  rename("Level" = agglom_lev) %>%
  # filter(grepl('ceae', Level)) %>% # keep all 
  group_by(Level) %>%
  summarise_at(vars(colNames), sum) %>%
  ungroup() -> agg_wide

# test prevalence/ab

prevelancedf = apply(X = agg_wide[-1],
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})

df <- data.frame(Prevalence = prevelancedf, 
                TotalAbundance = rowSums(agg_wide[-1]),
                Level = agg_wide$Level)

names(df)[3] <- agglom_lev

df %>%
  # parse taxonomy
  left_join(tax) %>% 
      distinct_at(agglom_lev, .keep_all = TRUE) %>%
  as_tibble() %>%
  rename("Level" = agglom_lev) %>%
  arrange(desc(TotalAbundance)) -> df


# intersected_taxa
df %>% arrange(Prevalence) 


pd <- position_dodge(0.1)

df %>% group_by(Phylum) %>% summarise(n = sum(TotalAbundance), 
                                      mean = mean(Prevalence),
                                      N = sum(!is.na(Phylum)),
                                      sd = sd(Prevalence, na.rm = T),
                                      se = sd  / sqrt(N)) %>% 
  arrange(desc(n)) -> summary_prevalence

summary_prevalence %>%
  # group_by(Phylum) %>%
  mutate(pct = N/ sum(N) * 100) %>%
  ggplot(aes(n, mean)) +
  scale_x_log10() +
  geom_point(aes(size = pct)) +
  ggrepel::geom_text_repel(aes(label = Phylum), size = 4) +
  # labs(y = "Prevalence (mean)", x = "Total Abundance") +
  labs(x = "x", y = 'y') +
  theme_bw(base_size = 16) +
  theme(legend.position = "top",
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank())

keepPhyla <- summary_prevalence$Phylum[summary_prevalence$n/sum(summary_prevalence$n) >= 0.001]

df %>%
  mutate(col = ifelse(Prevalence == 1, 'exclusive', 'Intersected')) %>%
  # filter(Phylum %in% keepPhyla) %>%
  mutate(Level = ifelse(is.na(Level), "Incomplete", "Complete")) %>%
  mutate(Phylum= ifelse(Phylum %in% keepPhyla, Phylum, "Low Taxa")) %>%
  mutate(Phylum = factor(Phylum, levels = c(keepPhyla, "Low Taxa"))) %>%
  # mutate(Top = ifelse(Level %in% fam_top$Level, TRUE, FALSE)) %>%
  mutate(Prevalence = Prevalence/length(colNames)) %>%
  ggplot(aes(TotalAbundance, Prevalence)) + 
  geom_point(size = 2, alpha = 0.7, aes(color = col)) + 
  geom_hline(yintercept = 0.1, alpha = 0.5, linetype = 2) +
  scale_x_log10() +
  labs(y = "Prevalence (Frac. Samples)", 
    x = "Total Abundance (log10)", 
    color = "") +
  scale_color_manual(values = c('red', 'blue')) +
  theme_bw(base_size = 17) +
  theme(
    legend.position = "top",
    title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12)) -> saveP


df %>% 
  # drop_na(Phylum) %>% 
  # group_by(Phylum) %>%
  rstatix::cor_test(Prevalence,TotalAbundance, method = 'spearman') %>%
  rstatix::add_significance()

saveP + theme(panel.border = element_blank()) -> saveP

ggsave(saveP, filename = "Family_prevalence.png", path = dir, 
       width = 4, height = 4)

save(df, keepPhyla, file = paste0(dir, "ANCOM_BC_inputs.RData"))

pick_top <- function(x, y, top = 10) {
  
  # x <- vector of abundance
  # y <- vector of taxonomic name
  
  ordered <- order(x, decreasing = T)
  topPos <- head(ordered, top)
  taxPos <- y[topPos]
  
  return(taxPos)
}




obj %>%
  inner_join(tax) %>% distinct(Family, .keep_all = T) %>%
  select_at(vars(ranks)) %>%
  distinct_at(agglom_lev, .keep_all = T) %>%
  rename("Level" = agglom_lev) %>%
  left_join(agg_wide) %>%
  filter(Phylum %in% keepPhyla) %>%
  select_at(vars(c("Level",colNames))) %>%
  mutate(Level = ifelse(is.na(Level), "Incomplete", Level)) -> agg_wide_filtered


apply(agg_wide_filtered[-1], 2, pick_top, top = 10, 
      y = agg_wide_filtered$Level) %>%
  as_tibble() %>%
  pivot_longer(all_of(colNames), names_to = 'Index', 
               values_to = "Level") %>%
  distinct(Level) %>%
  inner_join(agg_wide_filtered) -> fam_top

obj %>% filter(Family %in% keepTaxa) -> fam_top


# make tax clustering 

library(superheat)

m <- data.frame(fam_top[-1])
rownames(m) <- fam_top$Level

raf <- function(x) x/sum(x) * 100

superheat(apply(m, 2, raf),
          scale = F,
          row.dendrogram = T,
          col.dendrogram = T,
          clustering.method = 'hierarchical',
          dist.method = 'euclidean',
          print.plot = F) -> sh

tax_hclust <- sh$order.rows
tax_hclust <- rownames(m)[tax_hclust]

obj %>%
  inner_join(tax) %>% distinct(Family, .keep_all = T) %>%
  mutate(Family = ifelse(is.na(Family), "Incomplete", Family)) %>%
  select_at(vars(ranks)) %>%
  distinct_at(agglom_lev, .keep_all = T) %>%
  rename("Level" = agglom_lev) %>%
  inner_join(fam_top, by = "Level") %>%
  mutate_at(colNames, raf) %>%
  pivot_longer(cols = colNames, 
               names_to = "ID", values_to = 'ra') %>%
  filter(ra > 0) %>% inner_join(mtd) -> dataHeat

# sanity check ----



# set left-panel (Phylum) ordering 
dataHeat %>% group_by(Phylum) %>% summarise(t = sum(ra)) %>% arrange(desc(t)) %>% pull(Phylum) -> PhylumLevel


dataHeat %>%
  # mutate(Time = str_replace_all(Time, c("Day" = ""))) %>%
  # mutate(Time = factor(Time, levels = TimeLev)) %>%
  mutate(Level = factor(Level, levels = tax_hclust),
         Phylum = factor(Phylum, levels = PhylumLevel),
         ra = ifelse(ra == 0, NA, ra)) %>%
  ggplot(aes(y = Level, x = ID, fill = ra)) +
  geom_tile() +
  # facet_grid( ~ as.factor(Tissue) , scales = "free", space = "free", switch = "x") +
  geom_tile(color = "black") +
  ggsci::scale_fill_material("blue-grey",  
                             name="Relative\nAbundance\n", 
                             na.value = 'white',
                             limits = c(0,100),
                             labels = scales::percent_format(scale = 1)) +
  # ggh4x::facet_nested(~ Tissue + Time, scales = "free", space = "free") +
  ggh4x::facet_nested(Phylum + Class ~ Tissue + Time,
                      scales = "free", space = "free" ) +
  labs(x = NULL, y = NULL) + 
  theme_classic(base_family = "GillSans", base_size = 14) +
  guides(fill = guide_colorbar(barheight = unit(9, "in"), 
                               barwidth = unit(0.5, "in"),
                               ticks.colour = "black", 
                               frame.colour = "black",
                               label.theme = element_text(size = 14))
  ) +
  
  theme(
    title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12),
    strip.text.y = element_text(
      angle = 0, 
      size = 14),
    strip.background = element_rect(colour = "black", 
                                    fill = "transparent",
                                    size = 0.4),
    panel.spacing = unit(0.007, "lines")) -> heatPlot

ggsave(heatPlot, filename = "heatmap.png", path = dir, 
       width = 22, height = 16)

save(obj, mtd, colNames, keepPhyla,
     tax_hclust, ranks, file = paste0(dir, "/objects.RData"))

