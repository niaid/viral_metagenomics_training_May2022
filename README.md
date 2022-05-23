---
title: "Introduction to identifying viral sequences in bulk metagenomic data"
author: "David B. Stern, Ph.D."
---

[Bioinformatics and Computational Biosciences Branch](https://bioinformatics.niaid.nih.gov/)


## Learning Objectives:

* Highlight Locus resources for viral metagenomic analyses
* Discuss approaches and bioinformatic tools for viral sequence identification in metagenomic assemblies
* Implement a simplified workflow to identify phage sequences in metagenomic data including:
    * Classification of metagenomic contigs as viral / non-viral
    * Genome quality assessment
    * Gene annotation
    * Taxonomic classification (if time allows)

Much of this workflow was adapted from the [Sullivan Lab](https://u.osu.edu/viruslab/) [SOP](https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3?version_warning=no)


#### Pipelines
VirusSeeker, PAIPline, LAZYPIPE, ViroMatch, VIP, virMine

# checkv db: checkv-db-v1.0
#checkv download_database .
# DRAMv: db-dramv
#DRAM-setup.py prepare_databases --skip_uniref --output_dir db-dramv

## Start interactive session and set up working directory

[Link to temporary account information](https://nih-my.sharepoint.com/:x:/g/personal/sterndb_nih_gov/EcAbW7ESV9JOsFFtR1Q9fJsBiwMBSiqjJgAdd4A9aYKNSQ?e=hxL0rq&wdLOR=c7C88D67F-EBFA-E94B-9009-9C8DA514C2FB)

```bash
ssh <username>@ai-submit1.niaid.nih.gov
qrsh -pe threaded 8 -l h_vmem=2G
```

The workflow is designed to take place after metagenome assembly (as in Metagenomics I). For this tutorial, we are using a very small subset (5 sequences) that were output from the previous session (i.e. sample 6 from the CAMI dataset)

Copy data to working directory

```bash
# if using training account ->
    cp /classhome/classroom/test.fa .
# if using personal account
    cp /hpcdata/scratch/viral_metagenomics_training/test.fa .
```

Take a quick look at the data


```bash
less test.fa
```

## Run VirSorter2
[](figs/virsorter2_fig1.png)

[VirSorter](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00990-y) is a multi-classifier method that uses genomic features to assign sequences a 'viralness' score. This applied both to the whole sequence and to sliding windows to identify partial viral sequences (e.g. proviruses). Importantly, different random-forest classifiers were trained for five different groups of viruses with different genomic characteristics, biology, evolutionary origins, etc. (dsDNA phages, NCLDV, RNA, ssDNA, Laviviruses)

### Three steps performed automatically:  
1. Preprocess sequences and identify circular contigs
2. Extract features from input sequences
    - CDS identification with Prodigal (ref) and annotated with HMMER3 against Pfam and custom viral database
        - Manually annotated "viral hallmark genes"
        - Genes assigned to virus, bacteria, archaea, eukaryotes, mixed
    - Circular sequences
    - Gene size, density, and overlapping frequency
    - Strand switching frequency
    - Start codon usage
    - GC content
    - Ribosomal binding site motifs
3. Score each sequence independently using a set of classifiers customized for different viral groups
    - Random forest classifiers trained on high-quality reference genomes
    - Aggregate scores into a single prediction

VirSorter2 is a written as a SnakeMake workflow, and each step of the workflow is printed to the screen when running.

The main command is `virsorter run` and has several options that can be viewed by typing `virsorter run -h`

```bash
module load virsorter/2.2.3-Python-3.8.10
# The reference database is already downloaded and configured on Locus here: /hpcdata/bio_data/virsorter/db-2.2.3/
# If it were not, the command to download and configure the database is:
# virsorter setup -d /path/to/database/db-vs2 -j 8

virsorter run -i test.fa -w vs2-pass1 --keep-original-seq --include-groups dsDNAphage,ssDNA --min-length 1000 --min-score 0.5 -j 8 all

```

Flags:
- `-i test.fa`: specify name and path to input file  
- `-w vs2-pass1`: specify name and path to output directory  
- `--keep-original-seq`: partial viral sequences are not trimmed from the whole contig. Instead we will use CheckV to remove cellular sequence and identify proviruses.
- `--min-score 0.5`: VirSorter2 suggests that a score above 0.9 is strong evidence that the sequence is viral. However, we will be using a low cutoff to capture more sequences of putative viral origin and use checkV to collect additional information.


This produces several output files in vs2-pass1 directory:
- final-viral-score.tsv: contig information and scores

  > This table can be used for further screening of results. It includes the following columns:
  >   - sequence name
  >   - score of each viral sequences across groups (multiple columns)
  >   - max score across groups
  >   - max score group
  >   - contig length
  >   - hallmark gene count
  >   - viral gene %
  >   - nonviral gene %

- final-viral-combined.fa: all viral sequences in fasta format

  > identified viral sequences, including three types:
  > - full sequences identified as viral (identified with suffix `||full`);
  > - partial sequences identified as viral (identified with suffix `||{i}_partial`); here `{i}` can be numbers starting from 0 to max number of viral fragments found in that contig;
  > - short (less than two genes) sequences with hallmark genes identified as viral (identified with suffix `||lt2gene`);

- final-viral-boundary.tsv: table with ORF coordinates and information

 > only some of the columns in this file might be useful:
  >   - seqname: original sequence name
  >   - trim\_orf\_index\_start, trim\_orf\_index\_end:  start and end ORF index on orignal sequence of identified viral sequence
  >   - trim\_bp\_start, trim\_bp\_end:  start and end position on orignal sequence of identified viral sequence
  >   - trim\_pr: score of final trimmed viral sequence
  >   - partial:  full sequence as viral or partial sequence as viral; this is defined when a full sequence has score > score cutoff, it is full (0), or else any viral sequence extracted within it is partial (1)
  >   - pr\_full:  score of the original sequence
  >   - hallmark\_cnt:  hallmark gene count
  >   - group: the classifier of viral group that gives high score; this should **NOT** be used as reliable classification

Sequence names are appended with `||full` or `{i}_partial`. `||full` means that the entire contig has strong viral signal, while `{i}_partial` sequences have some viral and some cellular signal.


# run checkv to qc virsorter results and trim host regions left at the end of proviruses
**insert checkv figure**

CheckV estimates completeness of viral genomes assembled from metagenomic data and also removes host (bacterial) contamination. Alternative methods include VIBRANT (ref.) which assesses completeness based on viral hallmark genes and viralComplete (ref.) which compares sequence length to related sequences in NCBI RefSeq.

For completeness, CheckV uses two methods

```bash
module load checkv/0.8.1
```
CheckV has several modules, each with their own set of options that can be viewed with `checkv -h`. CheckV is a pipeline that consists of several steps


Running `checkv end_to_end` will run the entire pipeline of estimating completeness, contamination, and identify closed genomes.



```bash
checkv end_to_end vs2-pass1/final-viral-combined.fa checkv -t 16 -d /hpcdata/bio_data/checkv/checkv-db-v1.0

# concatenate viral sequences extracted from mostly bacterial sequence (proviruses) and contigs that are mostly viral
cat checkv/proviruses.fna checkv/viruses.fna > checkv/combined.fna

module unload checkv
```
Let's take a look at the main results in `less checkv/quality_summary.tsv`.


## Let's run DRAMv annotate the sequences and help filter some false-positive
**insert dramv figure**

DRAM is a annotation tool for bacterial and viral (DRAM-v) sequences. It


```bash
module load dram/1.3.4-Python-3.10.4
#python -m pip install click
time virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off \
        --prep-for-dramv -i checkv/combined.fna \
        -w vs2-pass2 --include-groups dsDNAphage,ssDNA --min-length 1000 --min-score 0.5 -j 16 all

DRAM-setup.py print_config
```
**no need to run `DRAM-v.py annotate` now, it will take too long**
We can copy a pre-prepared output `cp xxx`

```bash
# run DRAMv to annotate identified sequences
# step 1 annotate
DRAM-v.py annotate -i vs2-pass2/for-dramv/final-viral-combined-for-dramv.fa \
                    -v vs2-pass2/for-dramv/viral-affi-contigs-for-dramv.tab \
                    -o dramv-annotate --skip_trnascan --threads 16 --min_contig_size 1000
## took 21m..
#step 2 summarize anntotations
DRAM-v.py distill -i dramv-annotate/annotations.tsv -o dramv-distill
# took 4s
```

# Compile and filter results
Now we have output from virsorter2, checkv, and dramv. Let's combine these results into a single table using R and filter to
a set of high-confidence viral sequences.

This could be run locally in RStudio. For now, we can run R in Locus in the terminal.

```bash
module load R/4.1.0  
R
```

```R
library(dplyr)
library(stringi)
library(stringr)
#` read in data
vs2_res <- read.delim("vs2-pass1/final-viral-score.tsv",h=T)

checkv_res <- read.delim("checkv/quality_summary.tsv",h=T)
colnames(checkv_res)[1] <- "seqname"

dramv_vMAG_stats <- read.delim("dramv-distill/vMAG_stats.tsv",h=T)
dramv_vMAG_stats$seqname <- dramv_vMAG_stats$X %>%
                    str_replace('-.+','') %>%
                    str_replace('full_\\d','full') %>%
                    str_replace('partial_\\d','partial') %>%
                    stri_replace_last_fixed('__','||')

#` merge tables
res_tmp <- merge(vs2_res,checkv_res,by="seqname")
res <- dramv_vMAG_stats %>%
        select(seqname,potential.AMG.count,Viral.genes.with.host.benefits) %>%
        merge(res_tmp,by="seqname")
```


```R
#` filter by criteria suggested by Sullivan lab
keep1 <- res %>%
        filter(viral_genes > 0)
keep2 <- res %>%
        filter(viral_genes == 0 & (host_genes == 0 | max_score > 0.95 | hallmark > 2))


```

Some genes are common in both viruses and hosts, and can cause false positives in the Keep2 category. We can screen the DRAMv results for those suspicious genes and then check those sequences manually.


```R
dramv_annotations <- read.delim("dramv-annotate/annotations.tsv",h=T)

suspicious_genes <- c('carbohydrate kinase',
'carbohydrate-kinase',
'glycosyltransferase',
'glycosyl transferase',
'glycosyl transferaseendonuclease',
'nucleotide sugar epimerase',
'nucleotide sugar-epimerase',
'nucleotide-sugar epimerase',
'nucleotide-sugar-epimerase',
'nucleotidyltransferase',
'nucleotidyl transferase',
'nucleotidyl-transferase',
'plasmid stability',
'endonuclease')

to_rm <- data.frame()
for (gene in suspicious_genes){
    out <- filter(dramv_annotations, grepl(gene,pfam_hits,ignore.case = TRUE))
    to_rm <- rbind(to_rm,out)
}

sus_contigs <- unique(to_rm$scaffold) %>%
    stri_replace_last_fixed("__","||") %>%
    str_replace("-cat.*","")

#` filter the keep2 list
keep2_good <- filter(keep2, !seqname %in% sus_contigs)

final <- rbind(keep1,keep2_good)

write.table(final, 'good_viral_contigs.txt',quote=F,row.names=F,sep='\t')
```

Let's collect the viral sequences in fasta format
```bash
module load bbmap

cut -d $'\t' -f 1 good_viral_contigs.txt > keep_seqlist.txt

filterbyname.sh in=checkv/combined.fna out=good_viral_contigs.fa include=t substring=name names=keep_seqlist.txt
```


## Blast some of these

## Assigning taxonomy
Viruses do no have a single, universal marker gene (like 16S, ITS, COI) that can be used for phylogenetic or similarity-based taxonomic assignment (although group specific markers do exist, e.g. RDP for RNA viruses)
Bacteriophage taxonomy is xxx. Recent papers have proposed several new groups and classification systems, e.g.:
1.
2.
3.

Two recent methods for assigning viral genomes to taxonomic groups are:
1. vConTACT is xxx
2. GRAViTy

vConTACT is available on Locus, but takes a long time to run and process the output. For reference, one can load the module with `module load xxx`

GRAViTy

## VConTACT2

Predict genes using prodigal
```bash
module load prodigal/2.6.3
module load vcontact2/0.9.19-Python-3.7.10

prodigal -p meta -i good_viral_contigs.fa -a good_viral_contigs.prodigal_out.fa -o good_viral_contigs.prodigal_out.txt
grep '>' good_viral_contigs.prodigal_out.fa > good_viral_contigs.prodigal_names.txt

sed 's/>//g' good_viral_contigs.prodigal_names.txt | \
awk -F' ' '{print $1 "," $1 "REF," $9}' - | \
sed 's/_[0-9]\+REF/REF/g' > g2g.tmp.csv

echo 'protein_id,contig_id,keywords' | cat - g2g.tmp.csv > g2g.csv

awk -F '\t' '{print $1 "," $3 "," $25}' dramv-annotate/annotations.tsv | tail -n +2 - > g2g.dramv.tmp.csv
echo 'protein_id,contig_id,keywords' | cat - g2g.dramv.tmp.csv > g2g.dramv.csv  
```

Ideally, we would also combine our sequences with additional reference sequences with associated taxonomy, as the vConTACT reference library is somewhat limited.

vConTACT takes a long time to run, even on a small dataset. Here is the command to run it, but no need to do that now.

```bash
time vcontact2 --raw-proteins good_viral_contigs.prodigal_out.fa \
            --rel-mode 'Diamond' \
            --proteins-fp g2g.csv \
            --db 'ProkaryoticViralRefSeq201-Merged' \
            --pcs-mode MCL --vcs-mode MCL \
            --output-dir vcontact_out

vcontact2 --raw-proteins good_viral_contigs.prodigal_out.fa \
            --rel-mode 'Diamond' \
            --proteins-fp g2g.csv \
            --db 'ProkaryoticViralRefSeq201-Merged' \
            --pcs-mode MCL --vcs-mode ClusterONE --c1-bin ~/local/bin/cluster_one-1.0.jar \
            --output-dir vcontact_out

vcontact2 --raw-proteins dramv-annotate/genes.faa \
            --rel-mode 'Diamond' \
            --proteins-fp g2g.dramv.csv \
            --db 'ProkaryoticViralRefSeq201-Merged' \
            --pcs-mode MCL --vcs-mode ClusterONE --c1-bin ~/local/bin/cluster_one-1.0.jar \
            --output-dir vcontact_out_dram
```

Let's take a look at the output. Two files are informative:
