
dir <- "/Users/cigom/Documents/ANALYSIS_JUREL/"

files <- list.files(path = dir, 
  pattern = "Tabla abundancias familias tesis Rocio Valenzuela", full.names = TRUE)

# prep dat (omit) ----

obj <- read_xlsx(files, sheet = 1, skip = 0, col_names = T)
mtd <- read_xlsx(files, sheet = 2, skip = 0, col_names = T) %>% mutate(ID = as.character(ID))

obj %>% select_if(is.double) %>% names() -> colNames

# group by tissue ---

library(ggupset)

treatL <- c("CONTROL", "1","2", "3")

obj %>% 
  pivot_longer(-Family, values_to = 'ab', names_to = "ID") %>%
  left_join(mtd) %>%
  mutate(Tissue = recode_factor(Tissue, IP = "Hindgut", 
    IA = "Foregut")) %>%
  filter(ab > 0) %>%
  group_by(Family, Tissue, Treatment) %>%
  # summarise by factors
  summarise(ab = sum(ab)) %>% 
  mutate(Treatment = factor(Treatment, levels = treatL)) %>%
  summarise(across(Treatment, .fns = list), n = n()) %>% #view()
  ggplot(aes(x = Treatment, fill = Tissue)) +
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

ggsave(psave, filename = 'upset.png', path = dir,
  width = 5,height = 3.7)


obj %>% 
  pivot_longer(-Family, values_to = 'ab', names_to = "ID") %>%
  left_join(mtd) %>%
  mutate(Tissue = recode_factor(Tissue, IP = "Hindgut", 
    IA = "Foregut")) %>%
  filter(ab > 0) %>%
  group_by(Family, Tissue, Treatment) %>%
  # summarise by factors
  summarise(ab = sum(ab)) %>% 
  mutate(Treatment = factor(Treatment, levels = treatL)) %>%
  group_by(Family, Treatment) %>%
  summarise(across(Tissue, .fns = list), n = n()) %>% #view()
  ggplot(aes(x = Tissue, group = Treatment)) +
  geom_bar(position = position_dodge(width = 1), width = 0.8, fill = 'black') +
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

psave <- psave + facet_grid(~ Treatment)

ggsave(psave, filename = 'upset_2.png', path = dir,
  width = 7,height = 3)


# 3 -----

obj %>% 
  pivot_longer(-Family, values_to = 'ab', names_to = "ID") %>%
  left_join(mtd) %>%
  mutate(Tissue = recode_factor(Tissue, IP = "Hindgut", 
    IA = "Foregut")) %>%
  filter(ab > 0) %>%
  group_by(Family, Tissue) %>%
  # summarise by factors
  summarise(ab = sum(ab)) %>% 
  # mutate(Treatment = factor(Treatment, levels = treatL)) %>%
  group_by(Family) %>%
  summarise(across(Tissue, .fns = list), n = n()) %>% #view()
  ggplot(aes(x = Tissue)) +
  geom_bar(position = position_dodge(width = 1), width = 0.8, fill = 'black') +
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


ggsave(psave, filename = 'upset_3.png', path = dir,
  width = 3,height = 3)
