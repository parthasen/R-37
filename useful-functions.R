get_seq_by_genomic_coord <- function(genome, chr, start, end) {
	# import the genome library (has to be from bioconductor)
	# e.g. for human: "BSgenome.Hsapiens.UCSC.hg19"
	do.call(library, list(genome))
	# get the sequence
	# chr should have "chr1" format
	if(grep("Hsapiens", genome)) {
		seq <- getSeq(Hsapiens, chr, start, end)
	}
	return(toString(seq))
}

get_genes_go_terms <- function(genes) {
	# genes is a vector of genes
	library(biomaRt)
	# define biomaRt data base
	mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
	# retrieve go_id and name_1006 for each gene
	res <- getBM(attributes = c("hgnc_symbol", "go_id", "name_1006"),
	    filters = c("hgnc_symbol"),
	    values = list(genes, TRUE), mart = mart)
	return(res)
}
