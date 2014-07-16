# a function that returns genomic sequence by given chromosome and its start and end coordinates
get_seq_by_genomic_coord <- function(genome, chr, start, end) {
	# genome is a genome library (has to be from bioconductor)
	# e.g. for human: "BSgenome.Hsapiens.UCSC.hg19"
	do.call(library, list(genome))
	# get the sequence
	# chr should have "chr1" format
	if(grep("Hsapiens", genome)) {
		seq <- getSeq(Hsapiens, chr, start, end)
	}
	return(toString(seq))
}

# a function that returns all GO terms (with a short description) of a given set of genes
get_genes_go_terms <- function(genes) {
	# genes is a vector of Ensembl gene names
	library(biomaRt)
	# define biomaRt data base
	mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
	# retrieve go_id and name_1006 for each gene
	res <- getBM(attributes = c("hgnc_symbol", "go_id", "name_1006"),
	    filters = c("hgnc_symbol"),
	    values = list(genes, TRUE), mart = mart)
	return(res)
}

# a function that returns all GO terms (with a short description) of all genes
get_all_genes_go_terms <- function() {
	library(biomaRt)
	# define biomaRt data base
	mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
	gene.data <- getBM(attributes=c('ensembl_gene_id', 'go_id', "name_1006"),
	    mart = mart)
}
