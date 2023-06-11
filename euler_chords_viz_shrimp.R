dir <- '~/Documents/Shrimp_Estefany/'
load(paste0(dir, "euler_outputs.RData"))
load(paste0(dir, "objects.RData"))
# phyloseq <- readRDS(paste0(dir, 'phyloseq_ancom.rds'))

library(tidyverse)

mtd %>% mutate(Tissue = recode_factor(Tissue, Intestine = "Hindgut", Hepatopancreas = "Midgut", Stomach = "Foregut")) -> mtd 

obj %>% select_at(all_of(ranks)) %>% 
  mutate(id = 1:nrow(obj)) %>%
  pivot_longer(cols = ranks) %>% fill(value) %>%
  pivot_wider(names_from = name) %>%
  select(-id) %>%
  data.frame(row.names = obj$`Feature ID`) -> tax

     
euler %>% as_tibble(rownames = "Family") %>%
  left_join(tax %>% select_at(ranks[2:5]) %>% distinct(Family, .keep_all = T)) %>%
  pivot_longer(cols = names(euler), names_to = "Tissue", values_to = "pa") -> df

obj %>% 
  filter(Phylum %in% keepPhyla) %>%
  pivot_longer(cols = colNames, 
               names_to = "Index", values_to = 'ab') %>%
  filter(ab > 0) %>%
  # select_at(ranks) %>% 
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = ranks) %>% fill(value) %>%
  pivot_wider(names_from = name) %>%
  select(-id) %>%
  left_join(mtd) %>%
  group_by(Tissue,Family) %>%
  summarise(ab = sum(ab), pa = sum(ab>0)) %>% arrange(desc(ab)) %>%
  ungroup() %>%
  filter(Family %in% exclusive_taxa) %>%
  left_join(tax %>% select_at(ranks[2:5]) %>% distinct(Family, .keep_all = T)) -> dff

# colSums(euler[rownames(euler) %in% exclusive_taxa, ])
# table(dff$Tissue)

# no sirve!
# circos.clear()
# dff %>% select(Phylum, Tissue, pa) %>% chordDiagram(link.arr.type = "big.arrow") 
# circos.track(ylim = c(-1, 1), panel.fun = function(x, y) {
#      value = dff$ab
#      circos.barplot(value, 1:length(value), col = NA, bar_width = 0.3)
#    })

df %>% 
  filter(pa > 0) %>%
  filter(Family %in% exclusive_taxa) %>% 
  with(., table(Phylum, Tissue)) %>% t()
df %>% 
  filter(pa > 0) %>%
  filter(Family %in% intersected_taxa) %>%
  with(., table(Phylum, Tissue)) %>% t()

# df %>% filter(Family %in% exclusive_taxa) -> df
df %>% filter(Family %in% intersected_taxa) -> df

library(circlize)

set.seed(260220)

n <- length(unique(df$Tissue))
grid.col <- ggsci::pal_rickandmorty(alpha = 0.8)(n)

names(grid.col) <- c("Hindgut", "Midgut", "Foregut")# sort(unique(df$Tissue),decreasing = T)

# establecer paleta de colores de filos

df %>% pull(Phylum) %>% sort() %>% unique() -> labels 

colourCount = length(labels)

library(ggsci)

if(colourCount > 7) {
  getPalette <- colorRampPalette(pal_locuszoom()(7))(colourCount)
} else
  getPalette <- pal_locuszoom(alpha = 1)(colourCount)  

names(getPalette) <- labels

# filename = paste0(dir, "/intersected_taxa_chord.png")
# intercalate the filename
filename = paste0(dir, "/exclusive_taxa_chord.png")
width= 1500;height=1500;res = 150

png(filename, width = width, height = height, res = res)

circos.clear()
circos.par(start.degree = 0, gap.degree = 4, 
           track.margin = c(-0.01, 0.01), 
           points.overflow.warning = FALSE)
df %>% 
  filter(pa > 0) %>%
  select(Tissue, Phylum) %>%
  with(., table(Phylum, Tissue)) %>%
  chordDiagram(grid.col = c(grid.col, getPalette),
               directional = -1,
               diffHeight = mm_h(5), target.prop.height = mm_h(4),
               # annotationTrack = "grid", 
               # direction.type = c("arrows", "diffHeight"),
               preAllocateTracks = 1,
               small.gap = 10, big.gap = 15)

dev.off()
# circos.track(track.index = 1, panel.fun = function(x, y) {
#   circos.text(CELL_META$xcenter, CELL_META$ylim[1],
#               CELL_META$sector.index,
#               facing = "clockwise",
#               niceFacing = TRUE, adj = c(0, 0.5))
# }, bg.border = NA)

# abline(h = 0.05, lty = 2, col = "#00000080")

# later test some of this https://www.nature.com/articles/s41522-018-0077-y#Sec7
