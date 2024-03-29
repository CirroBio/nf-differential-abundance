#!/usr/bin/env Rscript
library("corncob")
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

# Run corncob
print("Running corncob")
dv_analysis <- differentialTest(
  formula = ~ ${params.formula},
  phi.formula = ~ ${params.formula},
  formula_null = ~ 1,
  phi.formula_null = ~ 1,
  data = phy,
  test = "Wald",
  boot = FALSE,
  fdr_cutoff = ${params.fdr_cutoff},
  taxa_are_rows=TRUE,
  allow_noninteger=TRUE
)

reformat_dv <- function(m, n){
  if(class(m) == "summary.bbdml"){
    res <- m\$coefficients %>% melt
    colnames(res) <- c("parameter", "variable", "value")
    res\$id <- n
  }else{
    print(paste("No results found for", n))
    res <- data.frame(parameter=c(), variable=c(), value=c(), id=c())
  }
  res
}

res <- do.call(
  rbind,
  lapply(
    c(seq_len(length(dv_analysis\$all_models))),
    function(i){
      reformat_dv(dv_analysis\$all_models[[i]], names(dv_analysis\$p)[i])
    }
  )
)
print(head(res))

print("Writing results")
write.table(
  res,
  file = "corncob_results.csv",
  quote = FALSE,
  sep = ",",
  row.names = FALSE
)
print("Done")
file.copy(".command.log", "${task.process}.log")
