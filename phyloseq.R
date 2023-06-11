
rm(list = ls())

# lets EDA 

library(tidyverse)
library(phyloseq)

dir <- "/Users/cigom/Documents/ANALYSIS_JUREL"

# attach(paste0(dir, "dna-sequences-fitGTR.RData"))

tree <- get("fitGTR", pos=2)$tree

load(paste0(dir, "objects.RData"))


obj %>% select_at(all_of(ranks)) %>%
  data.frame(row.names = obj$`Feature ID`) -> tax

obj %>% select_at(all_of(colNames)) %>% 
  data.frame(row.names = obj$`Feature ID`)-> ab

rownames(mtd) <- mtd$Index

# identical(names(ab),rownames(mtd))
#identical(rownames(ab), rownames(tax))

phyloseq = phyloseq(otu_table(ab, taxa_are_rows = TRUE), 
                    tax_table(as(tax, 'matrix')), 
                    sample_data(mtd),
                    phy_tree(tree))

saveRDS(phyloseq, paste0(dir, "phyloseq.rds"))

# Lets generate a prevelance table (number of samples each taxa occurs in) for each taxa.
# https://ucdavis-bioinformatics-training.github.io/2019_September_UCD_Microbial_Community_Analysis_Workshop/MCA_Workshop_R/phyloseq

# phyloseq <- ps

prevelancedf = apply(X = otu_table(phyloseq),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame

# # instead of top 10 of taxa, lets take the prevalence feature in account! ----

prevelancedf = data.frame(Prevalence = prevelancedf,
                          TotalAbundance = taxa_sums(phyloseq),
                          tax_table(phyloseq))
colnames(prevelancedf) <- c("Prevalence", 
                            "TotalAbundance", 
                            colnames(tax_table(phyloseq)))
# plot(prevelancedf$Prevalence, log10(prevelancedf$TotalAbundance))
# prevelancedf %>% view()

prevelancedf %>%
  ggplot() +
  geom_point(aes(x = Prevalence, y = TotalAbundance)) +
  facet_wrap(~Phylum)

# Now lets investigate low prevelance/abundance phylum and subset them out.

summary_prevalence <- plyr::ddply(prevelancedf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})

# summary_prevalence %>% arrange(desc(total_abundance)) %>% view()

# Using the table above, determine the phyla to filter based on the 0.001 threshold

sum(summary_prevalence$total_abundance)*0.001
table(summary_prevalence$total_abundance/sum(summary_prevalence$total_abundance) >= 0.001)

keepPhyla <- summary_prevalence$Phylum[summary_prevalence$total_abundance/sum(summary_prevalence$total_abundance) >= 0.001]

physeq = subset_taxa(phyloseq, Phylum %in% keepPhyla)

summary_prevalence <- summary_prevalence[summary_prevalence$Phylum %in% keepPhyla,]

summary_prevalence

# Individual taxa filtering
# Subset to the remaining phyla by prevelance.

prevelancedf1 = subset(prevelancedf, Phylum %in% get_taxa_unique(physeq, taxonomic.rank = "Phylum"))
ggplot(prevelancedf1, aes(TotalAbundance,Prevalence / nsamples(physeq),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.10, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# here we can see a break above the 0.25 % with exception in proteobacteria

#  Define prevalence threshold as 10% of total samples ~ set of replicates
prevalenceThreshold = 0.10 * nsamples(physeq)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevelancedf1)[(prevelancedf1$Prevalence >= prevalenceThreshold)]
length(keepTaxa)

physeq2 = prune_taxa(keepTaxa, physeq)
physeq2


# Agglomerate taxa at the Genus level (combine all with the same name) keeping all taxa without genus level assignment

length(get_taxa_unique(physeq2, taxonomic.rank = "Family"))
physeq2_glom = tax_glom(physeq2, "Genus", NArm = FALSE)
physeq2_glom

sum(colSums(otu_table(physeq2_glom))) /sum(colSums(otu_table(phyloseq)))

# Now lets filter out samples (outliers and low performing samples)
# Do some simple ordination looking for outlier samples, first we variance stabilize the data with a log transform, the perform PCoA using bray’s distances



logt  = transform_sample_counts(physeq2_glom, function(x) log(1 + x) )
out.pcoa.logt <- ordinate(logt, method = "MDS", distance = "bray")
evals <- out.pcoa.logt$values$Eigenvalues
plot_ordination(logt, out.pcoa.logt, type = "samples",
                color = "Tissue", shape = "Time") + 
  labs(col = "Tissue") +
  coord_fixed(sqrt(evals[2] / evals[1]))

# Look for low perfroming samples

qplot(colSums(otu_table(physeq2_glom)),bins=30) + 
  xlab("Logged counts-per-sample")

physeq3 <- prune_samples(sample_sums(physeq2_glom)>=4000, physeq2_glom)
physeq3

# investigate transformations.

## for 
plot_abundance = function(physeq, meta, title = "",
                          Facet = "Phylum", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  # p1f = subset_taxa(physeq, Phylum %in% c("Proteobacteria"))
  p1f = physeq
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = meta,y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

plot_abundance(ps, "Tissue", title="original")
