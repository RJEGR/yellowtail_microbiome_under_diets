

rm(list = ls())

# 
#           Treatment
# Tissue     1  2  3 CONTROL
# Hindgut   10 10 10       9
# Foregut   10 10 10      15

# ... to identify the significant taxa between the different sections of the intestine and the treatment.

# Analysis composition of microbiomes w/ ANCOM-BC ----

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ANCOMBC")

# http://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
#
library(ANCOMBC)
library(phyloseq)
library(microbiome)
library(tidyverse)
library(readxl)

dir <- "/Users/cigom/Documents/ANALYSIS_JUREL/"


phyloseq <- read_rds(paste0(dir, '/inputs_ancom_ps.rds'))


lineage <- read_rds(paste0(dir, '/lineage.rds')) %>% 
  rename('taxon_id' = Family) %>% distinct(taxon_id, .keep_all = T)

out_df <- function(out_ancombc) {
  
  res = out_ancombc$res
  
  df_fig1 = data.frame(res$beta * res$diff_abn, check.names = FALSE) %>% 
    rownames_to_column("taxon_id") %>%
    pivot_longer(-taxon_id, names_to = "group", values_to = "logFC")
  
  df_fig2 = data.frame(res$se * res$diff_abn, check.names = FALSE) %>% 
    rownames_to_column("taxon_id") %>% 
    pivot_longer(-taxon_id, names_to = "group", values_to = "SE")
  
  df_fig3 = data.frame(res$q_val * res$diff_abn, check.names = FALSE) %>% 
    rownames_to_column("taxon_id") %>%
    pivot_longer(-taxon_id, names_to = "group", values_to = "q_val") %>%
    mutate(star = ifelse(q_val <.001, "***", 
      ifelse(q_val <.01, "**",
        ifelse(q_val <.05, "*", ""))))
  
  # table(df_fig3$star)
  # colnames(df_fig2)[-1] = paste0(colnames(df_fig2)[-1], "SD")
  
  df_fig1 %>% left_join(df_fig2) %>% left_join(df_fig3)
}


# prep dat (omit) ----

files <- list.files(path = dir, 
  pattern = "Tabla abundancias familias tesis Rocio Valenzuela", full.names = TRUE)

obj <- read_xlsx(files, sheet = 1, skip = 0, col_names = T)
mtd <- read_xlsx(files, sheet = 2, skip = 0, col_names = T) %>% mutate(ID = as.character(ID))

obj %>% select_if(is.double) %>% names() -> colNames

sam <- data.frame(mtd, row.names = mtd$ID) %>% 
  arrange(match(ID, colNames)) %>% mutate_if(is.character, as.factor)

sam %>% mutate(Tissue = recode_factor(Tissue, IP = "Hindgut", 
  IA = "Foregut")) -> sam

sam %>% mutate(Treatment = recode_factor(Treatment, `1` = "CD+4%", 
  `2` = "CD+10%", `3` = "CD+25%", CONTROL = 'CD')) -> sam

TreatmentL <- c("CD","CD+4%", "CD+10%", "CD+25%")

sam %>% mutate(Treatment = factor(Treatment, levels = TreatmentL)) -> sam

sam %>% with(., table(Tissue, Treatment, Time)) 


obj %>% select(Family) %>% data.frame(row.names = obj$Family) -> tax

tax <- tax %>% 
  left_join(lineage, by = c("Family"="taxon_id")) %>% 
  select(-Genus)

rownames(tax) <- tax$Family

dat <- obj %>% select_at(colNames) %>% data.frame(row.names = obj$Family)

names(dat) <- colNames

identical(names(dat), rownames(sam))

identical(rownames(dat), rownames(tax))

# and parse

phyloseq = phyloseq(otu_table(dat, taxa_are_rows = TRUE), 
                    tax_table(as(tax, 'matrix')), 
                    sample_data(sam)) 
  # transform_sample_counts(function(x) sqrt(x / sum(x)))

write_rds(phyloseq, file = paste0(dir, '/inputs_ancom_ps.rds'))

# continue ----

# Choose only intersected taxa to this test

contrasts(factor(sam$Treatment), sparse = T)

prevelancedf = apply(X = otu_table(phyloseq),
  MARGIN = 1,
  FUN = function(x){sum(x > 0)})

df <- data.frame(Prevalence = prevelancedf, 
  TotalAbundance = rowSums(otu_table(phyloseq)),
  Level = names(prevelancedf)) %>% as_tibble()

rstatix::cor_test(df, Prevalence, TotalAbundance, method = 'spearman')

# Filtering by abundance

threshold <- 0.001

keepTaxa <- df$Level[df$TotalAbundance/sum(df$TotalAbundance) >= threshold]

# sum(keepTaxa %in% intersected_taxa)

sum(df$TotalAbundance) # 71,309

df %>% filter(Level %in% keepTaxa) %>% arrange(Prevalence) -> df

sum(df$TotalAbundance)/71309

# exclusive taxa there?

df %>% 
  arrange(Prevalence) %>% 
  filter(Prevalence == 1)  # not!

length(keepTaxa) # 70 families 

phyloseq <- phyloseq %>%
  subset_taxa(Family %in% keepTaxa) %>%
  prune_taxa(taxa_sums(.) > 0, .)

# multiple contrast by Foregut ----

ancombc_data1 <-  phyloseq %>%
  subset_samples(Tissue=="Foregut") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  aggregate_taxa(., "Family")

# by default Control vs * ----

formula <- "Treatment" 

out1 <- ANCOMBC::ancombc(phyloseq = ancombc_data1, formula = formula,
                        p_adj_method = "holm", zero_cut = 0.98, lib_cut = 100,
                        group =  formula, struc_zero = TRUE, neg_lb = TRUE,
                        tol = 1e-5, max_iter = 1000, conserve = TRUE,
                        alpha = 0.05, global = FALSE)


out_df(out1) %>%
  # filter(abs(logFC) > 0.05) %>% 
  filter(abs(logFC) != 0) %>% 
  arrange(desc(logFC)) %>%
  mutate(group = str_replace_all(group, c("Treatment" = ""))) %>% 
  mutate(wrap = paste0("CD : ", group)) %>%
  mutate(group = ifelse(logFC > 0, "g1", "g2")) -> df_fig_1


# CD+4%GP vs * ----

ancombc_data2 <-  phyloseq %>%
  subset_samples(Tissue=="Foregut" & Treatment != "CD") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  aggregate_taxa(., "Family")

out2 <- ANCOMBC::ancombc(phyloseq = ancombc_data2, formula = formula,
                        p_adj_method = "holm", zero_cut = 0.98, lib_cut = 100,
                        group =  formula, struc_zero = TRUE, neg_lb = TRUE,
                        tol = 1e-5, max_iter = 1000, conserve = TRUE,
                        alpha = 0.05, global = FALSE)

out_df(out2) %>%
  filter(logFC != 0) %>% arrange(desc(logFC)) %>% 
  mutate(group = str_replace_all(group, c("Treatment" = ""))) %>% 
  mutate(wrap = paste0("CD+4% : ", group)) %>%
  mutate(group = ifelse(logFC > 0, "g1", "g2")) -> df_fig_2

# CD+10%GP vs * ----

ancombc_data3 <-  phyloseq %>%
  subset_samples(Tissue=="Foregut" & !Treatment %in% c("CD", "CD+4%")) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  aggregate_taxa(., "Family")

out3 <- ANCOMBC::ancombc(phyloseq = ancombc_data3, formula = formula,
  p_adj_method = "holm", zero_cut = 0.98, lib_cut = 100,
  group =  formula, struc_zero = TRUE, neg_lb = TRUE,
  tol = 1e-5, max_iter = 1000, conserve = TRUE,
  alpha = 0.05, global = FALSE)

out_df(out3) %>%
  filter(logFC != 0) %>% arrange(desc(logFC)) %>% 
  mutate(group = str_replace_all(group, c("Treatment" = ""))) %>% 
  mutate(wrap = paste0("CD+10% : ", group)) %>%
  mutate(group = ifelse(logFC > 0, "g1", "g2")) -> df_fig_3



table(df_fig_1$wrap)
table(df_fig_2$wrap)
table(df_fig_3$wrap)


# Data are represented by effect size (log fold change) and 95% confidence interval bars (two-sided; Bonferroni adjusted) derived from the ANCOM-BC model.\nAll effect sizes with adjusted p<0.05 are indicated:

write_tsv(rbind(df_fig_1, df_fig_2, df_fig_3), 
  file = paste0(dir, "ANCOMBC_RES_Foregut.tsv"))


# data viz ----
unique(df_fig_3$wrap)

wrapL <- c("CD : CD+4%", "CD : CD+10%" , "CD : CD+25%",
  "CD+4% : CD+10%", "CD+4% : CD+25%", "CD+10% : CD+25%")

rbind(df_fig_1, df_fig_2, df_fig_3) %>% 
  mutate(wrap = factor(wrap, levels = wrapL)) %>%
  # inner_join(lineage) %>% distinct(.keep_all = T) %>%
  mutate(
    y_star = logFC + (0.2+SE)*sign(logFC),
    ymin = (abs(logFC) - SE) * sign(logFC),
    ymax = (abs(logFC) + SE) * sign(logFC)) %>% 
  # mutate(Phylum = factor(Phylum,  levels = labels)) %>%
  ggplot(data = ., 
    aes(x = reorder(taxon_id, desc(taxon_id)), y = logFC, fill = group)) + 
  geom_bar(stat = "identity", width = 0.5, 
    position = position_dodge(width = 0.4)) +
  ggh4x::facet_nested(  ~ wrap, scales = "free", space = "free", switch = "y") +
  coord_flip() +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2,
    position = position_dodge(0.05), color = "black") + 
  geom_text(aes(y = y_star, label=star), 
    vjust=.7, color="black", position=position_dodge(width = .5)) +
  ggsci::scale_fill_aaas() +
  guides(color = "none") +
  labs(x = NULL, y = "Log fold change") +
  theme_classic(base_size = 16, base_family = "GillSans") + 
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    strip.background.y = element_blank(),
    axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 6),
    axis.text.x = element_text(angle = 0, hjust = 1)) -> p


ggsave(p, filename = "ANCOMBC_foregut_multiple_contrast.png", path = dir, 
  width = 12, height = 8, dpi = 300)

# multiple contrast by Hindgut ----

ancombc_data1 <-  phyloseq %>%
  subset_samples(Tissue=="Hindgut") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  aggregate_taxa(., "Family")

# by default Control vs * ----

formula <- "Treatment" 

out1 <- ANCOMBC::ancombc(phyloseq = ancombc_data1, formula = formula,
  p_adj_method = "holm", zero_cut = 0.98, lib_cut = 100,
  group =  formula, struc_zero = TRUE, neg_lb = TRUE,
  tol = 1e-5, max_iter = 1000, conserve = TRUE,
  alpha = 0.05, global = FALSE)



# plot(abs(out_df(out1)$logFC), out_df(out1)$SE)

out_df(out1) %>%
  # filter(abs(logFC) > 0.05) %>% 
  filter(abs(logFC) != 0) %>% 
  arrange(desc(logFC)) %>%
  mutate(group = str_replace_all(group, c("Treatment" = ""))) %>% 
  mutate(wrap = paste0("CD : ", group)) %>%
  mutate(group = ifelse(logFC > 0, "g1", "g2")) -> df_fig_1


# CD+4%GP vs * ----

ancombc_data2 <-  phyloseq %>%
  subset_samples(Tissue=="Hindgut" & Treatment != "CD") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  aggregate_taxa(., "Family")

out2 <- ANCOMBC::ancombc(phyloseq = ancombc_data2, formula = formula,
  p_adj_method = "holm", zero_cut = 0.98, lib_cut = 100,
  group =  formula, struc_zero = TRUE, neg_lb = TRUE,
  tol = 1e-5, max_iter = 1000, conserve = TRUE,
  alpha = 0.05, global = FALSE)

out_df(out2) %>%
  filter(logFC != 0) %>% arrange(desc(logFC)) %>% 
  mutate(group = str_replace_all(group, c("Treatment" = ""))) %>% 
  mutate(wrap = paste0("CD+4% : ", group)) %>%
  mutate(group = ifelse(logFC > 0, "g1", "g2")) -> df_fig_2

# CD+10%GP vs * ----

ancombc_data3 <-  phyloseq %>%
  subset_samples(Tissue=="Hindgut" & !Treatment %in% c("CD", "CD+4%")) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  aggregate_taxa(., "Family")

out3 <- ANCOMBC::ancombc(phyloseq = ancombc_data3, formula = formula,
  p_adj_method = "holm", zero_cut = 0.98, lib_cut = 100,
  group =  formula, struc_zero = TRUE, neg_lb = TRUE,
  tol = 1e-5, max_iter = 1000, conserve = TRUE,
  alpha = 0.05, global = FALSE)

out_df(out3) %>%
  filter(logFC != 0) %>% arrange(desc(logFC)) %>% 
  mutate(group = str_replace_all(group, c("Treatment" = ""))) %>% 
  mutate(wrap = paste0("CD+10% : ", group)) %>%
  mutate(group = ifelse(logFC > 0, "g1", "g2")) -> df_fig_3



table(df_fig_1$wrap)
table(df_fig_2$wrap)
table(df_fig_3$wrap)


# Data are represented by effect size (log fold change) and 95% confidence interval bars (two-sided; Bonferroni adjusted) derived from the ANCOM-BC model.\nAll effect sizes with adjusted p<0.05 are indicated:

write_tsv(rbind(df_fig_1, df_fig_2, df_fig_3), 
  file = paste0(dir, "ANCOMBC_RES_Hindgut.tsv"))


# data viz ----

unique(df_fig_3$wrap)
wrapL <- c("CD : CD+4%", "CD : CD+10%" , "CD : CD+25%",
  "CD+4% : CD+10%", "CD+4% : CD+25%", "CD+10% : CD+25%")

rbind(df_fig_1, df_fig_2, df_fig_3) %>% 
  mutate(wrap = factor(wrap, levels = wrapL)) %>%
  # inner_join(lineage) %>% distinct(.keep_all = T) %>%
  mutate(
    y_star = logFC + (0.2+SE)*sign(logFC),
    ymin = (abs(logFC) - SE) * sign(logFC),
    ymax = (abs(logFC) + SE) * sign(logFC)) %>% 
  # mutate(Phylum = factor(Phylum,  levels = labels)) %>%
  ggplot(data = ., 
    aes(x = reorder(taxon_id, desc(taxon_id)), y = logFC, fill = group)) + 
  geom_bar(stat = "identity", width = 0.5, 
    position = position_dodge(width = 0.4)) +
  ggh4x::facet_nested(  ~ wrap, scales = "free", space = "free", switch = "y") +
  coord_flip() +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2,
    position = position_dodge(0.05), color = "black") + 
  geom_text(aes(y = y_star, label=star), 
    vjust=.7, color="black", position=position_dodge(width = .5)) +
  ggsci::scale_fill_aaas() +
  guides(color = "none") +
  labs(x = NULL, y = "Log fold change") +
  theme_classic(base_size = 16, base_family = "GillSans") + 
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    strip.background.y = element_blank(),
    axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0, size = 6),
    axis.text.x = element_text(angle = 0, hjust = 1)) -> p


ggsave(p, filename = "ANCOMBC_Hindgut_multiple_contrast.png", path = dir, 
  width = 12, height = 8, dpi = 300)

# exit ----
