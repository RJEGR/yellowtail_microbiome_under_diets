
library(tidyverse)

# load(paste0(dir, "/objects.RData"))

# pre-process tax 

tax_f <- list.files(path = dir, pattern = "Taxonomy", full.names = TRUE)


tax <- read_xlsx(tax_f, sheet = 1, skip = 0, col_names = T) %>% distinct(Taxonomy)

ranks <- c("Phylum", "Class", "Order", "Family", "Genus")

tax %>% separate(col = "Taxonomy", into = ranks, sep = "[|]") -> tax

tax[tax == ""] <- NA  

# clean taxonomy
# 

tax %>%
  mutate_at(ranks, funs(str_replace_all(., c("_[1-9]" = "", ".+_Incertae_Sedis"=NA_character_)))) -> tax

tax %>% distinct() -> tax


# obj %>% select_at(all_of(ranks)) %>%
#   data.frame(row.names = obj$`Feature ID`) %>%
#   mutate_all(., funs(str_replace_all(., c(".+_unclassified"=NA_character_, "_"=" ", "Unclassified"=NA_character_, "NA"=NA_character_, "sp."=NA_character_, "Undetermined"=NA_character_, "Unknown"=NA_character_)))) -> tax

# Resolution ----
max.rank <-  ncol(tax) 

tail(!is.na(tax))

# apply(!is.na(tax), 1, function(x) which(x == T))

tail(res <- rowSums(!is.na(tax))) # max.rank - 

features <- c(nrow(tax) - table(res))
pct <- c(features / nrow(tax))

data.frame(ranks, features, pct) %>%
  mutate(ranks = factor(ranks, levels = ranks)) %>%
  ggplot(aes(x = ranks, y = features)) +
  geom_path(size = 1.5, alpha=0.6, group = 1) +
  geom_point(size = 3, alpha=0.6) +
  ylim(250, 360) + 
  xlim(ranks[-5]) +
  geom_text(aes(label = paste0(round(pct*100, digits = 2), "%")), 
    size = 4, vjust = -1) +
  labs(y = "Taxonomic groups", x = '') +
  theme_bw(base_size = 20) -> ps 

ps + theme(panel.border = element_blank()) -> ps

ggsave(ps, filename = "resolution.png", path = dir, 
  width = 5.5, height = 5)

# fill ---
tax %>%
  select_at(ranks) %>% 
  mutate(id = 1:nrow(tax)) %>%
  pivot_longer(cols = all_of(ranks)) %>% fill(value) %>% # con esto rellenas espacios vacios o NA !!!
  pivot_wider(names_from = name) %>%
  select(-id) -> tax


# la esperanza es que todos los asvs/otus lleguen a la clasificacion mas profunda (species), evaluemos cuanto fue asi:

# divergence
length(unique(df$Phylum))
df %>%
  drop_na(Phylum) %>%
  group_by(Phylum) %>%
  summarise(
    n = n(),
    nR = length(unique(Level)),
    divergence = 1-(1/length(unique(Level))),
    nTotal = sum(TotalAbundance)) %>%
  mutate(pct = nTotal/ sum(nTotal) * 100) %>%
  arrange(n) %>%
  mutate(cs = cumsum(n), csR = cumsum(nR)) -> dataV 


dataV %>% view()

dataV %>%
  filter(divergence > 0.5) %>% # pull(Phylum)
  summarise(cor(nR, n, method = "pearson"),
    cor(divergence, nR, method = "pearson"),
    cor(divergence, nTotal))

dataV %>%
  mutate(facet = ifelse(csR > 15, "A", "B")) %>%
  ggplot(aes(divergence, csR)) +
  geom_point(aes(size = pct)) +
  scale_size(name = "Abundance (%)") +
  # ggsci::scale_color_gsea() +
  ggrepel::geom_text_repel(aes(label = Phylum), size = 4) +
  facet_grid(facet~., scales = "free_y",space = "fixed") +
  # scale_color_manual(values = c('black', 'red')) +
  labs(x = "x", y = 'y') +
  theme_bw(base_size = 16) +
  theme(legend.position = "top",
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()) -> p1


ggsave(p1, filename = "divergence.png", path = dir, 
  width = 5, height = 6)

# include the prevalence and abundance of data

obj %>% select_at(all_of(colNames)) %>% 
  data.frame(row.names = obj$`Feature ID`)-> ab

prevelancedf = apply(X = ab,
  MARGIN = 1,
  FUN = function(x){sum(x > 0)})

identical(rownames(tax), rownames(ab)) # sanity check

df = data.frame(Prevalence = prevelancedf, 
  TotalAbundance = rowSums(ab),
  tax) %>% as_tibble(rownames = "id")

#

# Now lets investigate low prevelance/abundance phylum and subset them out.

df %>%
  group_by(Phylum) %>%
  summarise(cor = cor(Prevalence, TotalAbundance, method = "pearson"),
    mean_prevalence = mean(Prevalence), 
    total_abundance = sum(TotalAbundance)) %>%
  drop_na(Phylum) %>%
  arrange(desc(total_abundance)) -> summary_prevalence

# OrderPhylumLevel <- summary_prevalence$Phylum 
# Using the table above, determine the phyla to filter based on the 0.001 threshold
threshold <- 0.001

sum(summary_prevalence$total_abundance)*threshold

table(summary_prevalence$total_abundance/sum(summary_prevalence$total_abundance) >= threshold)

keepPhyla <- summary_prevalence$Phylum[summary_prevalence$total_abundance/sum(summary_prevalence$total_abundance) >= threshold]

summary_prevalence <- summary_prevalence[summary_prevalence$Phylum %in% keepPhyla,]

summary_prevalence 

# Individual taxa filtering
# Subset to the remaining phyla by prevelance.

df %>%
  filter(TotalAbundance > 1) %>%
  filter(Phylum %in% keepPhyla) %>%
  group_by(Phylum) %>%
  summarise(cor = cor(Prevalence, TotalAbundance))

df %>%
  filter(TotalAbundance > 1) %>%
  # mutate(F_Resolution = ifelse(Resolution >= 4, TRUE, FALSE)) %>%
  # left_join(r) %>%
  # mutate(rank = factor(rank, levels = ranks)) %>%
  # filter(Phylum %in% keepPhyla) %>%
  mutate(Phylum= ifelse(Phylum %in% keepPhyla, Phylum, "Low Taxa")) %>%
  mutate(Phylum = factor(Phylum, levels = c(keepPhyla, "Low Taxa"))) %>%
  mutate(Prevalence = Prevalence/ncol(ab)) %>%
  ggplot(aes(TotalAbundance, Prevalence)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_hline(yintercept = 0.10, alpha = 0.5, linetype = 2) +
  # geom_vline(xintercept = 0.1, alpha = 0.5, linetype = 2) +
  scale_x_log10() +
  # geom_smooth(se = F, color = "red") +
  labs(y = "Prevalence (Frac. Samples)", x = "Total Abundance (log10)", color = "Family\nResolution") +
  facet_wrap(~Phylum) +
  theme_bw(base_size = 17) +
  theme(
    title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12)) -> p2

p2 + theme(panel.border = element_blank()) 

ggsave(p2, filename = "Prevalence.png", path = dir, 
  width = 8, height = 8)

# exit


load(paste0(dir, "lulu_curated_result.Rdata"))

otu_map <- curated_result$otu_map %>% as_tibble(rownames = "id")

otu_map %>%
  select(id, curated) %>%
  left_join(df) -> df

df %>% view()
df%>% 
  summarise(cor(Prevalence, TotalAbundance))

transform_sample_counts(ps, function(x) sqrt(x / sum(x)))

df %>%
  filter(Resolution > 0) %>%
  # group_by(Family)
  mutate(normAb = sqrt(TotalAbundance / sum(TotalAbundance))) %>%
  ggplot() +
  geom_point(aes(x = Prevalence, y = normAb, color = curated)) +
  facet_wrap(~ as.factor(Resolution), scales = "free")

r <- c(0:5)
names(r) <- ranks
data.frame(Resolution = r) %>% as_tibble(rownames = "rank") -> r

sum(df$TotalAbundance)
df %>%
  group_by(Resolution) %>% 
  summarise(t = sum(TotalAbundance)) %>%
  arrange(desc(t)) %>%
  mutate(ct = cumsum(t)) %>%
  mutate(reads_pct = ct/sum(t)) %>%
  left_join(r) %>%
  select(reads_pct,rank) -> pct_reads

df %>%
  group_by(Resolution) %>% # curated, 
  tally(sort = T) %>%
  # summarise(n = n(), t = sum(TotalAbundance)) %>%
  mutate(cs = cumsum(n)) %>%
  mutate(pct = cs/nrow(df)) %>%
  left_join(r) %>%
  left_join(pct_reads) %>%
  mutate(rank = factor(rank, levels = ranks)) -> pct_asvs

cor(pct_asvs$pct, pct_asvs$reads_pct)

pct_asvs %>%
  ggplot(aes(x = rank, y = cs, color = reads_pct, fill = reads_pct)) +
  # ggplot(aes(x = rank, y = n, color = curated, fill = curated, group = curated)) +
  geom_path(size=1, alpha=0.6, group = 1) +
  geom_point(size=2, alpha=0.6) +
  geom_text(aes(label = paste0(round(pct*100, digits = 2), "%")), size=4, vjust=-1) +
  labs(y = "# feature (Cumulative)", x = '') +
  theme_bw(base_size = 17) 


ggsave(p2, filename = "resolution.png", path = dir, 
  width = 6, height = 6)

# geom_col(aes(x = rank, y = n, color = curated), position = position_dodge(width = 0.4))

# Gracias papa
