# nf-differential-abundance
Nextflow workflow for differential abundance analysis

## Inputs

The format of input data can vary widely. Because the parameters
used to specify the inputs may vary by input type, this repository
will use different `main-script` entrypoints for the different
input types.

### Shared Parameters

Parameters which are used in common across all input types.

- `formula`: String like `"colA + colB"` used for statistical modeling
- `method`: Statistical method used (options: `corncob` (default), `mannwhitneyu`)
- `data_output`: Output folder used for human-readable data
- `web_output`: Output folder used for web-readable data

### Microbial 16S (nf-core/ampliseq)

Script: `ampliseq.nf`

For the analysis of microbial 16S sequencing data, input files are expected to be
formatted following the pattern used in the [nf-core/ampliseq](https://nf-co.re/ampliseq)
workflow.

**Counts Table (e.g. qiime2/abundance_tables/feature_table.tsv)**:

Tab-delimited file with the number of counts for each ASV in each sample.

```
# Constructed from biom file
#OTU ID	Sample1	Sample2	Sample3
8e121f6e82f82d6e7e2f831b9ce9e75c	0.0	0.0	71.0
22f4ee9a41a4d73580bf7ade8e9e017a	400.0	168.0	237.0
2555ab896130cd0373de325088a94110	4.0	5.0	77.0
b4ed8ec9d3313b1ae0b929812511444b	10.0	19.0	5.0
70d55baf78e9ac4d0babeac5dcbae5c2	37.0	62.0	53.0
2eedf896e95f02090d22dd3e68742ce7	64.0	30.0	147.0
bc44453105cf321ef72daa44e95a16ab	0.0	14.0	0.0
...
```

**Taxonomy Table (e.g. qiime2/taxonomy/taxonomy.tsv)**:

Tab-delimited file with the taxonomic assignment for each ASV.

```
Feature ID	Taxon	Confidence
c728ad6f5d183cb36fa06b6a3a47758b	d__Bacteria; p__Firmicutes; c__Clostridia; o__Oscillospirales; f__Ruminococcaceae; g__Faecalibacterium	0.9999993830869841
bbae6ed124f4d6b48435a964a95c8418	d__Bacteria; p__Firmicutes; c__Clostridia; o__Oscillospirales; f__Ruminococcaceae; g__Faecalibacterium	0.9999984489550084
22f4ee9a41a4d73580bf7ade8e9e017a	d__Bacteria; p__Firmicutes; c__Clostridia; o__Oscillospirales; f__Ruminococcaceae; g__Faecalibacterium	0.9999973729608546
ec6c9233b72090035eb87c6ff39a84a7	d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_coprocola	0.9999852779388605
8d478879084ed7c88660c8b0c4b24923	d__Bacteria; p__Firmicutes; c__Clostridia; o__Oscillospirales; f__Ruminococcaceae	0.999599801457001
b264ac8ff5f9aff74f0b9aa084d9a9f0	d__Bacteria; p__Firmicutes; c__Clostridia; o__Lachnospirales; f__Lachnospiraceae; g__Agathobacter	0.9911833494102088
b6a1cfc87c3ffa7f96227d4a372d8dd7	d__Bacteria; p__Firmicutes; c__Clostridia; o__Oscillospirales; f__Ruminococcaceae; g__Subdoligranulum; s__uncultured_bacterium	0.7169473785430964
...
```

**Metadata Table**:

Comma-delimited file with the metadata information assigned to each sample.
This information is used for the `formula` parameter and statistical modeling.

```
sample,batch
Sample1,1
Sample2,1
Sample3,2
Sample4,2
```

Parameters:

- `counts`: URI to counts table
- `taxonomy`: URI to taxonomy table
- `metadata`: URI to metadata table

### MetaPhlAn

Script: `metaphlan.nf`

Each sample has input taxonomic abundance data in MetaPhlAn output format
which resembles:

```
#mpa_vJan21_CHOCOPhlAnSGB_202103
#/usr/local/bin/metaphlan --nproc 16 --input_type bowtie2out --biom SRR27087601.biom -t rel_ab_w_read_stats --unclassified_estimation --bowtie2db db --index mpa_vJan21_CHOCOPhlAnSGB_202103 --sample_id_key SRR27087601 --sample_id SRR27087601 SRR27087601.bowtie2.bz2 SRR27087601.metaphlan
#2948649 reads processed
#SRR27087601	SRR27087601
#estimated_reads_mapped_to_known_clades:1788459
#clade_name	clade_taxid	relative_abundance	coverage	estimated_number_of_reads_from_the_clade
UNCLASSIFIED	-1	38.89555	-	1160190
k__Bacteria	2	61.10445	0.29824	1788459
k__Bacteria|p__Firmicutes	2|1239	61.10445	0.29824	1788459
k__Bacteria|p__Firmicutes|c__Bacilli	2|1239|91061	61.10445	0.29824	1788459
k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales	2|1239|91061|1385	61.10445	0.29824	1788459
k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillaceae	2|1239|91061|1385|186817	61.10445	0.29824	1788459
k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillaceae|g__Bacillus	2|1239|91061|1385|186817|1386	61.10445	0.29824	1788459
k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillaceae|g__Bacillus|s__Bacillus_cereus	2|1239|91061|1385|186817|1386|1396	61.10445	0.29824	1788459
k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillaceae|g__Bacillus|s__Bacillus_cereus|t__SGB7697_group	2|1239|91061|1385|186817|1386|1396|	61.10445	0.29824	1788459
```

> Note: Data for `estimated_number_of_reads_from_the_clade` is required

Parameters:

- `samplesheet`: URI to CSV with columns `sample` and `file`, and any additional metadata columns used for the `formula`
- `tax_level`: `'species'` by default

### Sourmash

Script: `sourmash.nf`

Each sample has input taxonomic abundance data in Sourmash CSV summary output format
which resembles:

```
query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank,query_ani_at_rank,total_weighted_hashes
SRR26891235,superkingdom,0.5497580475489164,d__Bacteria,5ad71395,inputs/SRR26891235_2.fastq.gz,0.9729654193376661,5226000,0.9719126143642576,225415
SRR26891235,superkingdom,0.45024195245108356,unclassified,5ad71395,inputs/SRR26891235_2.fastq.gz,0.027034580662333885,4280000,,225415
SRR26891235,phylum,0.5497580475489164,d__Bacteria;p__Bacillota,5ad71395,inputs/SRR26891235_2.fastq.gz,0.9729654193376661,5226000,0.9719126143642576,225415
SRR26891235,phylum,0.45024195245108356,unclassified,5ad71395,inputs/SRR26891235_2.fastq.gz,0.027034580662333885,4280000,,225415
SRR26891235,class,0.5497580475489164,d__Bacteria;p__Bacillota;c__Bacilli,5ad71395,inputs/SRR26891235_2.fastq.gz,0.9729654193376661,5226000,0.9719126143642576,225415
SRR26891235,class,0.45024195245108356,unclassified,5ad71395,inputs/SRR26891235_2.fastq.gz,0.027034580662333885,4280000,,225415
SRR26891235,order,0.5497580475489164,d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales,5ad71395,inputs/SRR26891235_2.fastq.gz,0.9729654193376661,5226000,0.9719126143642576,225415
SRR26891235,order,0.45024195245108356,unclassified,5ad71395,inputs/SRR26891235_2.fastq.gz,0.027034580662333885,4280000,,225415
SRR26891235,family,0.5497580475489164,d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae_G,5ad71395,inputs/SRR26891235_2.fastq.gz,0.9729654193376661,5226000,0.9719126143642576,225415
SRR26891235,family,0.45024195245108356,unclassified,5ad71395,inputs/SRR26891235_2.fastq.gz,0.027034580662333885,4280000,,225415
SRR26891235,genus,0.5497580475489164,d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A,5ad71395,inputs/SRR26891235_2.fastq.gz,0.9729654193376661,5226000,0.9719126143642576,225415
SRR26891235,genus,0.45024195245108356,unclassified,5ad71395,inputs/SRR26891235_2.fastq.gz,0.027034580662333885,4280000,,225415
SRR26891235,species,0.5442878182200716,d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__Bacillus_A cereus,5ad71395,inputs/SRR26891235_2.fastq.gz,0.9664796042854291,5174000,0.9714499051552751,225415
SRR26891235,species,0.00547022932884494,d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__Bacillus_A tropicus,5ad71395,inputs/SRR26891235_2.fastq.gz,0.006485815052236985,52000,0.7803437885603269,225415
SRR26891235,species,0.45024195245108356,unclassified,5ad71395,inputs/SRR26891235_2.fastq.gz,0.027034580662333885,4280000,,225415
```

Parameters:

- `samplesheet`: URI to CSV with columns `sample` and `file`, and any additional metadata columns used for the `formula`
- `tax_level`: `'species'` by default

### curatedMetagenomicData

Script: `curatedMetagenomicData.nf`

Input data is the `"*.relative_abundance` table output by the
[`curatedMetagenomicDataTerminal` tool](https://github.com/waldronlab/curatedMetagenomicDataTerminal).

> Note: It is strongly recommended to use the `--counts` flag so that the
integer values can be used in the beta-binomial model used by corncob
(which take into account sequencing depth information which is otherwise
not available).

Using the curatedMetagenomicData tool, input data is combined into a single
table which includes both organism abundances and sample metadata.

The workflow will split this data up before analysis, taking advantage of
the pattern that organism abundances are formatted using the MetaPhlAn
nomenclature (e.g. `s__Escherichia_coli`).

Parameters:

- `input`: URI to CSV with output from the curatedMetagenomicDataTerminal` tool

### Relative Abundance Table

Script: `table.nf`

Arbitrary inputs can be processed from any tabular data file.
The relative abundance information must be in a single CSV file in which:

- Rectangular table including a header row and index column
- Each column contains a sample (observation)
- Each row contains a feature (species, metabolite, gene, etc.)
- Each cell contains a numeric value (integer or float)

| index    | sampleA | sampleB | sampleC |
| -------- | ------- | ------- | ------- |
| species1 | 0.1     | 0.08    | 0.15    |
| species2 | 0.85    | 0.72    | 0.68    |
| species3 | 0.05    | 0.20    | 0.17    |

> Note: While the example table shows proportional values which add up to 1
for each sample, there is no such requirement or expectation for this data type.

In addition, a metadata file must be provided which contains
the columns of data which are used in the `formula` input.
The values in the first column of the metadata table must match the
sample identifiers listed in the first column of the relative abundance table.

Parameters:

- `input`: Path to relative abundance table
- `metadata`: Path to metadata table
- `read_kwargs`: Any additional keyword arguments which should be provided
to `pandas.read_csv` as JSON (default: `'{"sep": ",", "index_col": 0}'`) for reading the
`input` file

> Note that while a CSV input is expected, TSV inputs can be parsed by
supplying `read_kwargs`: `'{"sep": "\t,", "index_col": 0}'`


### gig-map/align_reads.nf

script: `gig_map_align_reads.nf`

Gene-level read counts can be processed from the output of the
[gig-map](https://github.com/FredHutch/gig-map) workfow, specifically
using the `align_reads.nf` workflow.
The output of that workflow is the `read_alignments.csv.gz` file, which can
be used as the `input` to this workflow.

Samples are identified in the `specimen` column in that file.
However, the expected `metadata` file will use the `sample` column header
to align with the behavior of the other entrypoints in this repository.

Parameters:

- `input`: Path to `read_alignment.csv.gz` produced by `gig-map/align_reads.nf`
- `metadata`: Path to metadata file (with a `sample` column)
- `min_coverage`: Minimum gene coverage to use for inclusion (range: 0-1, default: 0.75)
- `use_unaligned`: If `true`, calculate the proportional abundance of each gene as the fraction
of reads aligning to that gene divided by the total number of reads collected from the specimen.
If `false` (by default), use the total number of reads which align to any gene in this analysis
as the denominator instead.
