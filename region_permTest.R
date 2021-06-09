## Function to measure enrichment analysis of regions with gene-sets
## Written by Eduardo Maury (eduardo_maury@hms.harvard.edu)

library(regioneR)
library(AnnotationHub)
library(biomaRt)
library(GenomeInfoDb)
library(GenomicRanges)

run_pt <- function(gene_set, regions, niter){
  ## input: 
  ## - gene_set: vector of gene symbols
  ## - regions: GenomeRanges object with regions to test
  ## Note: currently only set for autosomal regions, but can be adapted to others. 
  ## Note: currently adapted to the grch37 but can be adapted based on biomart host. 
  ## output:
  ## - permtest object from the regioneR. 
  ensembl <- useMart("ensembl")
  ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  x = getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand','hgnc_symbol'),
        filters=c('hgnc_symbol'),
        values=gene_set,
        mart=ensembl)
  gene_ranges = toGRanges(x)
  gene_ranges = gene_ranges[(gene_ranges@seqnames %in% paste(1:22))]
  seqlevels(gene_ranges) = paste(1:22)
  newStyle <- mapSeqlevels(seqlevels(gene_ranges), "UCSC")
  gene_ranges <- renameSeqlevels(gene_ranges, newStyle)
  human.genome <- getGenomeAndMask("hg19")
  
  set.seed(777)
  perm_test = overlapPermTest(A=regions, ntimes=niter, B=gene_ranges, 
    genome = filterChromosomes(human.genome$genome, organism="hg",chr.type="autosomal"),
    mask=human.genome$mask, per.chromosome=FALSE)
  return(perm_test)
}


## Example:
# regions = toGRanges(data.frame(chr=paste0("chr",df$CHROM), start=df$BEG_GRCh37, end=df$END_GRCh37))

# ## Synaptic genes 
# synaptome_genes = read_table("data/metadata/synaptome.txt",col_names="genes")
# head(synaptome_genes)
# pt_synaptic <- run_pt(synaptome_genes, regions, niter=1000)
# pt_synaptic