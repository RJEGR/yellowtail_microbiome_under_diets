rm(list = ls())

dir <- "/Users/cigom/Documents/ANALYSIS_JUREL"

load(paste0(dir, "objects.RData"))

# Check redundant abreviation nomenclature ----
# raw datavis

# sum(table(obj$`Feature ID`) > 1)

agglom_lev <- "Family"

obj %>% 
  filter(!grepl('ceae', Family)) %>%
  pivot_longer(cols = colNames, 
               names_to = "ID", values_to = 'ab') %>%
  rename("Level" = agglom_lev) %>%
  filter(ab > 0) %>%
  arrange(desc(Level)) %>%
  # mutate(Level = ifelse(!grepl('ceae', Level), NA) %>%
  group_by(Level) %>%
  select(-ID) -> df1

color_by = "Phylum"
labels <- df1 %>% pull(color_by) %>% unique() %>% sort()
colourCount = length(labels)

library(ggsci)

if(colourCount > 7) {
  getPalette <- colorRampPalette(pal_locuszoom()(7))(colourCount)
} else
  getPalette <- pal_locuszoom()(colourCount)


df1 %>% 
  summarise(Freq = sum(ab > 0), ab = sum(ab)) %>%
  left_join(., df1 %>% 
              distinct(Level, .keep_all = T) %>% 
              select_if(is.character), by = "Level") %>% 
  arrange(desc(Freq)) %>%
  mutate(wrap = ifelse(ab > 60, "A", "B")) %>%
  ggplot() +
  geom_point(aes(Freq, ab, color = Phylum)) +
  facet_grid(wrap ~., scales = "free") +
  ggrepel::geom_label_repel(aes(Freq, ab, 
                                label = Level,
                                color = Phylum)) +
  theme_classic(base_size = 16) +
  theme(
    # legend.position = "none",
    strip.text.y = element_blank(), 
    strip.background = element_blank()) +
  labs(x = "Features Frequency (# ASVs)", 
       y = "Abundance (Reads)") +
  scale_color_manual(labels = labels, values = getPalette) -> p1

df1 %>% 
  ggplot() +
  geom_freqpoly(aes(ab, color = Phylum), 
                binwidth = 10) +
  # facet_grid(Phylum ~., scales = "free_y") +
  coord_flip() + labs(x = "", 
                      y = "") +
  scale_color_manual(labels = labels, values = getPalette) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none") -> p2

# devtools::install_github("thomasp85/patchwork")
library(patchwork)

p2 + p1  + 
  plot_layout(widths = c(1, 2), 
              heights = c(1,NULL)) -> p



ggsave(p, filename = "redundant_abreviations.png", path = dir, 
       width = 10, height = 5)

