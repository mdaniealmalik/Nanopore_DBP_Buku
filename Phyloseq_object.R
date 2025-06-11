# Tentukan "working directory" file kerja selama menjalankan R ini
setwd("C:/Users/danie/Downloads/Compressed/test_nanopore/Nanopore_DBP_Buku-main")

# Aktifkan package 
library(tidyverse)
library(phyloseq)

# Panggil OTU Table file
otu_table <- read.delim("otu_table.tsv", header=TRUE, row.names=1)
colnames(otu_table) <- c("Sample_1", "Sample_2", "Sample_3")
View(otu_table)

# Panggil Taxon Table file
taxonomy <- read.csv("taxon_table.csv", sep=";", row.names=1)
View(taxonomy)
taxonomy <- as.matrix(taxonomy)

# Panggil Metadata/Map file
metadata <- read.csv("map.csv", row.names=1, sep=";")

# Assign Phyloseq Object berdasarkan 3 data diatas (OTU Table, Taxon Table, dan Map/Metadata)
OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy)
META = sample_data(metadata)

physeq1 = phyloseq(OTU, TAX, META)

# Hilangkan singleton data (Optional)
physeq1_without_singleton = filter_taxa(physeq1, function(x) sum(x) > 1, TRUE)
physeq1_without_singleton

# Telusuri Taxa berdasarkan level taxon
get_taxa_unique(physeq1_without_singleton, "Family")
get_taxa_unique(physeq1_without_singleton, "Genus")

### Membuat Stacked Barplot (Family), Script (line 36 - 50) bisa disesuaikan berdasarkan level taxon (ex: Genus atau Species atau lainnya)
physeq_family <- physeq1_without_singleton %>%
  tax_glom(taxrank = "Family") %>%                        
  transform_sample_counts(function(x) {x/sum(x)} ) %>%    
  psmelt() %>%                                                                
  arrange(Family)

unique(physeq_family$Family)

Class_colors <- c("coral1", "yellow2", "cadetblue", "burlywood2", 
                  "brown", "green3", "pink", "darkturquoise", 
                  "darkorchid2", "darkorange2", "red", "blue")

ggplot(physeq_family, aes(x = Sample,  y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values =Class_colors) 

## Melihat Alpha diversity (Shannon dan Simpson)
p <- plot_richness(physeq1_without_singleton, x="Sample", measures = c("Shannon", "Simpson")) +  
  geom_point(size=7, alpha=0.5) + xlab("")+ ylab("") + theme_bw()

estimate_richness(physeq1_without_singleton, measures = c("Shannon", "Simpson"))

