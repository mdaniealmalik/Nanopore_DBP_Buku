# Tentukan "working directory" file kerja selama menjalankan R ini
setwd("D:/DMBP_Undip/eDNA_project_analysis/Tutorial_Nanopore_Book_DBP")

library(tidyverse)

# Panggil file hasil blastn
result_blastn <- read.delim("result_blastn.txt", header=FALSE) %>%
  select(V1, V2, V3) %>%
  rename(qseqid = V1, sseqid = V2, pident = V3)

# panggil taxonomy data dan gabungkan dengan file hasil blastn 
taxonomy_names <- read.delim("database/database.txt", header=FALSE) %>%
  rename(sseqid = V1, taxon = V2)

edna_blastn <- left_join(result_blastn, taxonomy_names, by = "sseqid")

# Filter hanya menampilkan yang memiliki pident lebih dari 90
pident_90 <- edna_blastn %>% 
  arrange(qseqid) %>%
  group_by(qseqid, pident, taxon) %>% 
  slice_head(n=1) %>% 
  filter(pident >= 90) %>%
  ungroup() %>%
  group_by(qseqid) %>%
  slice_max(order_by = pident, n =1)

# Kolom yang ingin diperiksa
column_to_check <- "qseqid"

# Hitung frekuensi masing-masing nilai di kolom yang ingin Anda periksa
value_counts <- table(pident_90[[column_to_check]])

value_morethan_1 <- as.data.frame(value_counts) %>%
  filter(Freq > 1)

# edit nama dengan menghapus nama spesies dengan OTU yang teridentifikasi 2 species dengan score ident sama
pident_90_replicate <- pident_90 %>% filter((qseqid %in% value_morethan_1$Var1))

# Function to remove the last word before the last semicolon
remove_last_word <- function(s) {
  parts <- strsplit(s, ";")[[1]]
  if (length(parts) > 1) {
    parts <- parts[-length(parts)]  
  }
  modified_string <- paste(parts, collapse = ";")
  modified_string <- paste0(modified_string, ";")
  return(modified_string)
}

# Apply the function to the second column
pident_90_replicate$taxon <- sapply(pident_90_replicate$taxon, remove_last_word)

# delete row yang memiliki OTU replicate
pident_non_replicate_otu <- anti_join(pident_90, pident_90_replicate, by = "qseqid")

p_ident_dereplicate <- full_join(pident_non_replicate_otu, pident_90_replicate)

taxon_table_clean <- p_ident_dereplicate %>% select(qseqid, taxon) %>% data.frame()

taxon_table_clean_edt <-taxon_table_clean %>% separate(taxon, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep = ';')

unique(taxon_table_clean_edt$Species)

colnames(taxon_table_clean)[1] <- "OTU"

# combine with otu table to assign non-identify OTU as unassigned
otu_table <- read.delim("otu_table.tsv")

taxon_table_90_ident <- full_join(taxon_table_clean, otu_table, by ="OTU") %>% 
  select(OTU, taxon)

taxon_table_90_ident_edt <-taxon_table_90_ident %>% separate(taxon, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep = ';')

# Simpan file dengan nama taxon_table.csv
write.csv2(taxon_table_90_ident_edt, "taxon_table.csv", row.names = F)
