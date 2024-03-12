#!/usr/bin/env Rscript
library("radEmu")
library("phyloseq")
library("reshape2")

# Read in abundance information
print("Reading input counts")
read_counts <- read.table(
  "counts.csv",
  sep = ",",
  header = TRUE,
  row.names = 1
)
print(head(t(head(t(read_counts)))))
print(paste("Rows:", nrow(read_counts), "Columns:", ncol(read_counts)))

# Read in metadata annotations
print("Reading metadata")
metadata <- read.table(
  "metadata.csv",
  sep = ",",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)
print(head(t(head(t(metadata)))))
print(paste("Rows:", nrow(metadata), "Columns:", ncol(metadata)))

if(length("${params.filter}") > 0){
  print("Filtering samples by ${params.filter}")
  metadata <- dplyr::filter(metadata, ${params.filter})
  print(paste("Rows:", nrow(metadata), "Columns:", ncol(metadata)))
}

# Read in the fake taxonomy table
print("Reading taxonomy")
taxonomy <- read.table(
  "taxonomy.csv",
  sep = ",",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)
print(head(taxonomy))

# Construct the phyloseq object
phy <- phyloseq(
  otu_table(as.matrix(read_counts), taxa_are_rows = FALSE),
  tax_table(as.matrix(taxonomy)),
  sample_data(metadata)
)

# Run radEmu
print("Running radEmu")
res <- emuFit(
  formula = ~ ${params.formula},
  Y = phy
)
print(res)

print("Writing results")
write.table(
  data.frame(res\$coef),
  file = "radEmu_results.csv",
  quote = FALSE,
  sep = ",",
  row.names = FALSE
)
print("Done")
file.copy(".command.log", "${task.process}.log")
