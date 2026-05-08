
## How to Run

### 1. Modify the input parameters

At the beginning of the script, update the following two parameters to match your local paths:

```r
GTF_FILE <- "/path/to/your_long_read_isoforms.gtf"
OUT_DIR  <- "/path/to/output/longread_scRBP_regions"


#for example
GTF_FILE <- "/mnt/isilon/project/longread/merged_isoforms.gtf"
OUT_DIR  <- "/mnt/isilon/project/longread/scRBP_regions"
```
### 2. Run the script

```bash
Rscript Gene_coordinate_longread_isoform_scRBP.R
```
#### It is recommended to save the running log:
```bash
Rscript Gene_coordinate_longread_isoform_scRBP.R > Gene_coordinate_longread_isoform_scRBP.log 2>&1
```
#### To run the script in the background:
```bash
nohup Rscript Gene_coordinate_longread_isoform_scRBP.R > Gene_coordinate_longread_isoform_scRBP.log 2>&1 &
```
####  To monitor the progress:
```bash
tail -f Gene_coordinate_longread_isoform_scRBP.log
```
## Output Files
> For each transcript domain, the script generates multiple output files. Using CDS as an example, the following files will be produced:

```bash
CDS_hg38.scan_regions.bed
CDS_hg38.region_to_transcript.tsv
CDS_hg38.transcript_segments.bed
CDS_hg38.transcript_segments.metadata.tsv
```

### 1. BED file for motif scanning
```bash
CDS_hg38.scan_regions.bed
```
> This is the recommended BED file for downstream sequence extraction using bedtools getfasta and motif scanning tools such as Cluster-Buster.

#### The BED file has the following format:
```bash
chrom    start    end    region_id    .    strand
```

#### for example
```bash
chr1    10000    10200    CDS_region_000000001    .    +
```
> This file is coordinate-deduplicated. In other words, if multiple isoforms share the same genomic region, that region is included only once for motif scanning.

### 2. Region-to-transcript/gene mapping table

```bash
CDS_hg38.region_to_transcript.tsv
```
>This file maps each scan region back to its corresponding transcript and gene annotation.
>It contains the following columns:
```bash
region_id
domain
transcript_id
gene_id
gene_name
part_index
chrom
start
end
strand
width
transcript_segment_id
```
#### Downstream motif hits can be mapped as follows:
```bash
motif hit
  ↓
region_id
  ↓
transcript_id / gene_id
  ↓
isoform-level ranking or gene-level ranking
```
### 3. Transcript-level BED file

```bash
CDS_hg38.transcript_segments.bed
```
> This file preserves transcript-specific fragments.
> For example, if the same genomic region is shared by five isoforms, this file will contain five rows corresponding to those isoform-specific segments.
> This file is generally not recommended for large-scale motif scanning because it may repeatedly scan the same genomic sequence. However, it is useful for checking isoform-level annotation structures.


## Downstream FASTA Extraction
> Assume the reference genome FASTA file is:
```bash
/mnt/isilon/reference/hg38.fa
```
#### extract strand-aware FASTA sequences using:
```bash
GENOME=/mnt/isilon/reference/hg38.fa
BEDDIR=/mnt/isilon/project/longread/scRBP_regions

for domain in fiveUTR threeUTR CDS intron exon
do
  bedtools getfasta \
    -fi ${GENOME} \
    -bed ${BEDDIR}/${domain}_hg38.scan_regions.bed \
    -s \
    -nameOnly \
    -fo ${domain}_hg38.sense.fa
done
```
> The -s option is important because it extracts strand-aware sequences according to the strand column in the BED file.

#### If your version of bedtools does not support -nameOnly, use -name instead:
```bash
bedtools getfasta \
  -fi ${GENOME} \
  -bed ${BEDDIR}/CDS_hg38.scan_regions.bed \
  -s \
  -name \
  -fo CDS_hg38.sense.fa
```











































