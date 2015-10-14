main <- function() {
	args <- commandArgs(trailingOnly = TRUE)
	# tsv file with following columns: chrom, pos, ref_allele, alt_allele, 
	# locus, name, ref_count, alt_count
	infile = args[1]
	# Output file for locus level information
	locus_outfile = args[2]
	# Output file for SNV level information
	snv_outfile = args[3]
	# Sample name
	sample = args[4]
	# Whether the input is phased (all ref alleles assumed to be same phase)
	isPhased = as.logical(args[5])
	# Number of simulation (1e6 in manual)
	numSim = as.numeric(args[6])
	threads = as.numeric(args[7])

	library(MBASED)
	set.seed(03252015)
	register(MulticoreParam(threads), default=TRUE)

	t = read.table(infile, header=TRUE)

	ranges = IRanges(start=t$pos, width=1)

	snvs = GRanges(
	    seqnames=as.character(t$chrom),
	    ranges=ranges,
	    aseID=as.character(t$locus),
	    allele1=as.character(t$ref_allele),
	    allele2=as.character(t$alt_allele)
	)
	
	names(snvs) = as.character(t$name)
	
	se <- SummarizedExperiment(
	    assays=list(
	        lociAllele1Counts=matrix(
	            t$ref_count,
	            ncol=1,
	            dimnames=list(
	            names(snvs),
	            sample
	            )
	        ),
	        lociAllele2Counts=matrix(
	            t$alt_count,
	            ncol=1,
	            dimnames=list(
	            names(snvs),
	            sample
	            )
	        )
	    ),
	rowData=snvs
	)
	
	results = runMBASED(se, isPhased=isPhased, numSim=numSim,
			    BPPARAM=bpparam())

	locus_df <- data.frame(
	    major_allele_freq=assays(results)$majorAlleleFrequency[,1],
	    p_val_ase=assays(results)$pValueASE[,1],
	    p_val_het=assays(results)$pValueHeterogeneity[,1]
	)

	write.table(locus_df, file=locus_outfile, sep="\t", quote=FALSE,
		    col.names=NA)

	snv_gr <- rowData(exptData(results)$locusSpecificResults)
	snv_gr$ref_is_major <- 
		assays(exptData(results)$locusSpecificResults)$allele1IsMajor[,1]
	snv_gr$maf <- 
		assays(exptData(results)$locusSpecificResults)$MAF[,1]
	snv_list <- 
		split(snv_gr, factor(snv_gr$aseID,
				     levels=unique(snv_gr$aseID)))
	snv_df = as.data.frame(snv_list)
	snv_df = snv_df[c("seqnames", "start", "aseID", "allele1", "allele2",
			  "ref_is_major", "maf")]
	colnames(snv_df)[1] = "chrom"
	colnames(snv_df)[2] = "position"
	colnames(snv_df)[3] = "locus"
	colnames(snv_df)[4] = "ref"
	colnames(snv_df)[5] = "alt"

	write.table(snv_df, file=snv_outfile, sep="\t", quote=FALSE,
		    col.names=NA)

	sessionInfo()
}

main()
