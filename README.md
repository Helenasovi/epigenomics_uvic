
## Section 4 – ATAC-seq downstream analyses

```bash
cd epigenomics_uvic/ATAC-seq
mkdir -p analyses data/bigBed.files data/bed.files analyses/peaks.analysis

../bin/download.metadata.sh "PASTE_METADATA_URL_HERE"

grep -F "bigBed_narrowPeak" metadata.tsv | \
grep -F "pseudoreplicated_peaks" | \
grep -F "GRCh38" | \
awk 'BEGIN{FS=OFS="\t"}{print $1,$11,$23}' | \
sort -k2,2 -k1,1r | sort -k2,2 -u \
> analyses/ATAC.bigBed.peaks.ids.txt

column -t analyses/ATAC.bigBed.peaks.ids.txt

cut -f1 analyses/ATAC.bigBed.peaks.ids.txt | \
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done

../bin/selectRows.sh <(cut -f1 analyses/ATAC.bigBed.peaks.ids.txt) metadata.tsv | \
cut -f1,46 > data/bigBed.files/md5sum.txt

cat data/bigBed.files/md5sum.txt | \
while read filename original_md5; do
  md5sum data/bigBed.files/"$filename".bigBed | \
  awk -v f="$filename" -v o="$original_md5" 'BEGIN{OFS="\t"}{print f,o,$1}'
done > tmp && mv tmp data/bigBed.files/md5sum.txt

# mismatches (should be empty)
awk '$2!=$3' data/bigBed.files/md5sum.txt

cut -f1 analyses/ATAC.bigBed.peaks.ids.txt | \
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done

# total peaks (so results are reproducible)
cut -f-2 analyses/ATAC.bigBed.peaks.ids.txt | \
while read filename tissue; do
  echo -n "$tissue total peaks: "
  wc -l data/bed.files/"$filename".bed
done

# peaks intersecting promoters
cut -f-2 analyses/ATAC.bigBed.peaks.ids.txt | \
while read filename tissue; do
  echo -n "$tissue peaks in promoters: "
  bedtools intersect -a data/bed.files/"$filename".bed \
    -b ../annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -u | wc -l
done

# peaks outside gene bodies
cut -f-2 analyses/ATAC.bigBed.peaks.ids.txt | \
while read filename tissue; do
  echo -n "$tissue peaks outside genes: "
  bedtools intersect -a data/bed.files/"$filename".bed \
    -b ../annotation/gencode.v24.protein.coding.gene.body.bed -v | wc -l
done
```



**Results**

- **Stomach**
  - Total peaks: 103,609
  - Peaks in promoters: 49,061
  - Peaks outside genes: 34,537

- **Sigmoid colon**
  - Total peaks: 110,999
  - Peaks in promoters: 52,229
  - Peaks outside genes: 37,035




## Section 5 – Distal regulatory activity

```bash
cd epigenomics_uvic
mkdir -p regulatory_elements
cd regulatory_elements

for tissue in stomach sigmoid_colon; do
  bedtools intersect -a ../ATAC-seq/analyses/peaks.analysis/outside_gene_bodies/${tissue}.outside_genes.bed \
    -b ../ChIP-seq/<REAL_H3K27ac_${tissue}_peaks>.bed -u | \
  bedtools intersect -a - \
    -b ../ChIP-seq/<REAL_H3K4me1_${tissue}_peaks>.bed -u \
  > ${tissue}.candidate.distal.regulatory.bed

  wc -l ${tissue}.candidate.distal.regulatory.bed
done

awk 'BEGIN{FS=OFS="\t"}$1=="chr1"{print $4,$2}' stomach.candidate.distal.regulatory.bed \
> regulatory.elements.starts.tsv

awk 'BEGIN{FS=OFS="\t"}$1=="chr1"{print}' ../annotation/gencode.v24.protein.coding.gene.body.bed \
> chr1.genes.bed

awk 'BEGIN{FS=OFS="\t"}{if ($6=="+"){start=$2} else {start=$3}; print $4,start}' chr1.genes.bed \
> gene.starts.tsv

python ../bin/get.distance.py --input gene.starts.tsv --start 980000

cat regulatory.elements.starts.tsv | while read element start; do
  python ../bin/get.distance.py --input gene.starts.tsv --start "$start" | \
  awk -v e="$element" -v s="$start" 'BEGIN{OFS="\t"}{print e,s,$1,$2,$3}'
done > regulatoryElements.genes.distances.tsv

Rscript -e 'd<-read.table("regulatoryElements.genes.distances.tsv",sep="\t"); x<-d[,5]; cat("mean",mean(x),"\n"); cat("median",median(x),"\n")'
```

**Results**

**Candidate distal elements**
- Stomach: 8,022
- Sigmoid colon: 14,215

**Distance to closest gene (chr1, stomach)**
- Mean distance: 45,227 bp
- Median distance: 27,735 bp

