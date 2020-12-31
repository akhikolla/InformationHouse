#' Gene set enrichment test
#'
#' @param S Ranks of gene set
#' @param N Integer value. Only required if \code{r} is not specified.
#' @param r Rank/correlation scores. If \code{S} is \code{character}, then \code{r} must be named by gene or be a character vector of the gene names in rank order (necessarily containing \code{S}).
#' @param p Weighting of ranking/correlations, see Subramanian et. al 2005.
#' @param min_its Minimum number of null permutations to compare.
#' @param max_its Maximum number of null permutations to compare.
#' @param significance_threshold Maximum p-value of significant result.
#' @param log_dismiss Threshold log probability of returning a significant result, below which function returns current p-value.
#' @param raw_score Logical value determining whether to return the raw value of the gene set enrichment score.
#' @return Numeric value - p-value of enrichment.
#' @examples 
#' gset(S=1:5 * 2, N=1000)
#' gset(S=letters[1:3], r=letters)
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom stats setNames
#' @useDynLib gsEasy
gset <- function(
	S,
	N=NULL,
	r=NULL,
	p=1,
	min_its=2e2,
	max_its=1e5,
	significance_threshold=0.05,
	log_dismiss=-10,
	raw_score=FALSE
) {
	stopifnot(is.numeric(p))
	stopifnot(length(p)==1)
	if (is.null(N) & is.null(r))
		stop("Must specify either N or r!")

	if (!((is.character(S) & (!is.null(names(r)) | is.character(r))) | (is.numeric(S) & is.numeric(r)) | (is.numeric(S) & !is.null(N))))
		stop("Either: S must be 'character' and r 'character' or named; S 'numeric' and N given; S and r 'numeric'")

	if (is.null(N))
		N <- length(r)

	if (is.null(r))
		r <- (N:1)/N

	r_order <- if (is.character(r)) seq(length(r)) else order(decreasing=TRUE, r)

	r_sorted <- abs(if (is.character(r)) setNames(nm=r, (length(r):1/length(r))) else r[r_order])

	stopifnot(is.vector(S) & (all(S %in% names(r_sorted)) | (is.numeric(S) & all(S > 0 & S <= length(r)))))

	S_cpp <- sort((if (is.character(S)) match(S, names(r_sorted)) else match(S, r_order))-1)

	if (!raw_score) {
		gset_raw(
			N,
			S_cpp,
			r_sorted^p,
			min_its,
			max_its,
			significance_threshold,
			log_dismiss
		)
	} else {
		es_raw(
			S_cpp,
			r_sorted^p
		)
	}
}

#' Create list of gene sets defined by ontological annotation 
#'
#' @param ontology \code{ontology_index} object.
#' @param gene Character vector of genes.
#' @param term Character vector of term IDs annotated to corresponding genes.
#' @param min_genes Minimum number of genes in gene sets.
#' @param max_genes Maximum number of genes in gene sets.
#' @return List of character vectors of term IDs.
#' @export
#' @importFrom ontologyIndex get_ancestors
get_ontological_gene_sets <- function(
	ontology,
	gene,
	term,
	min_genes=1,
	max_genes=500
) {
	gene.anno <- lapply(split(term, gene), get_ancestors, ontology=ontology)
	genes.by.term <- lapply(FUN=as.character, X=split(unlist(mapply(SIMPLIFY=FALSE, FUN=rep, names(gene.anno), sapply(gene.anno, length))), unlist(gene.anno)))
	Filter(x=genes.by.term, f=function(x) length(x) <= max_genes & length(x) >= min_genes)
}

#' Create list of gene sets defined by GO term annotation
#'
#' Note, this function takes several minutes to execute.
#'
#' @param GO_annotation_file File path of annotation file, which should contain a column of genes and a column of terms. Can be downloaded from at http://geneontology.org/gene-associations/gene_association.goa_human.gz.
#' @param GO_file File path of gene ontology.
#' @param min_genes Minimum number of genes in gene sets.
#' @param max_genes Maximum number of genes in gene sets.
#' @param verbose Print progress.
#' @return List of character vectors of term IDs.
#' @export
#' @importFrom ontologyIndex get_ontology
#' @importFrom utils read.table
get_GO_gene_sets <- function(
	GO_annotation_file,
	GO_file="http://purl.obolibrary.org/obo/go.obo",
	min_genes=15,
	max_genes=500,
	verbose=TRUE
) {
	if (verbose) cat("reading ontology file...\n")
	go <- get_ontology(GO_file)

	if (verbose) cat("reading annotation file...\n")
	anno.df <- read.table(GO_annotation_file, sep="\t", quote="", stringsAsFactors=FALSE, comment.char="!", skipNul=TRUE)

	if (verbose) cat("creating gene set list...\n")
	get_ontological_gene_sets(ontology=go, term=anno.df[,5], gene=anno.df[,3], min_genes=min_genes, max_genes=max_genes)
}

#' List of gene sets annotated by each GO term
#'
#' Based on gene-GO term annotations downloaded from \code{geneontology.org}. Only contains gene sets for terms with up to 500 genes.
#'
#' @name GO_gene_sets
#' @title GO term gene sets 
#' @docType data
#' @format List of character vectors of genes per GO term, and named by term ID.
NULL


