#!/usr/bin/env Rscript
library("corncob")
library("phyloseq")
library("reshape2")

# Read in abundance information
print("Reading input counts")
read_counts <- read.table(
  "/share/data/merged.csv",
  sep = ",",
  header = TRUE,
  row.names = 1
)
print(head(read_counts))

# Read in metadata annotations
print("Reading metadata")
metadata <- read.table(
  "/share/data/metadata.csv",
  sep = ",",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)
print(head(metadata))

# Read in the fake taxonomy table
print("Reading taxonomy")
taxonomy <- read.table(
  "/share/data/taxonomy.csv",
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

# Run corncob
print("Running corncob")
dv_analysis <- differentialTest(
  formula = ~ CRC + PRJEB10878 + PRJEB6070 + PRJEB7774,
  phi.formula = ~ CRC + PRJEB10878 + PRJEB6070 + PRJEB7774,
  formula_null = ~ 1,
  phi.formula_null = ~ 1,
  data = phy,
  test = "LRT",
  boot = FALSE,
  fdr_cutoff = 0.05
)

reformat_dv <- function(m, n){
  if(class(m) == "summary.bbdml"){
    res <- m$coefficients %>% melt
    colnames(res) <- c("parameter", "variable", "value")
    res$id <- n
  }else{
    print(paste("No results found for", n))
    res <- data.frame(parameter=c(), variable=c(), value=c(), id=c())
  }
  res
}

res <- do.call(
  rbind,
  lapply(
    c(seq_len(ncol(read_counts))),
    function(i){
      reformat_dv(dv_analysis$all_models[[i]], colnames(read_counts)[i])
    }
  )
)
print(head(res))

print("Writing results")
write.table(
  res,
  file = "/share/data/corncob_results.csv",
  quote = FALSE,
  sep = ",",
  row.names = FALSE
)
print("Done")