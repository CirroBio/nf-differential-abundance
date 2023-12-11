# nf-differential-abundance
Nextflow workflow for differential abundance analysis

## Inputs

The format of input data can vary widely. Because the parameters
used to specify the inputs may vary by input type, this repository
will use different `main-script` entrypoints for the different
input types.

### Shared Parameters

Parameters which are used in common across all input types.

- `metadata`: CSV with sample metadata (required column name: `sample`)
- `formula`: String like `"colA + colB"` used for statistical modeling
- `method`: Statistical method used (options: `corncob` (default), `wilcoxon`)
- `data_output`: Output folder used for human-readable data
- `web_output`: Output folder used for web-readable data

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

- `samplesheet`: URI to CSV with columns `sample` and `file`
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

- `samplesheet`: URI to CSV with columns `sample` and `file`
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
