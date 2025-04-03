#' Prepare copy number data for downstream analysis
#'
#'
#' @param X a matrix or a data frame of copy number data. The rows and columns
#'   of X correspond to genes and subjects, respectively. X must have gene
#'   symbols as its row names.
#' @param Y an optional matrix or data frame of copy number data to compare with
#'   X. As is the case for X, the rows and columns of Y correspond to genes and 
#'   subjects, respectively. Y must have gene symbols as its row names.
#' @param species a value specifying species for ensembl database to be used;
#'   Default is "human".
#'
#' @return a list of processed X and Y with a data frame containing gene
#'   annotation information.
#'
#' @details This helper function is designed to prepare input matrices or data frames
#'   X and Y containing DNA copy number data for analysis. The downstream functions
#'   that X and Y are indexed by a common set of genes that appear in genomic order.
#'   In addition, the peelingOne and peelingTwo functions require information about
#'   the cytoband for each gene, and the resultsProcess function uses information
#'   about gene position. The dataPrep function queries biomaRt so users are not
#'   required to provide this information.
#'
#'
#' @examples
#' \donttest{
#' #This runtime for this code slightly exceeds the limits imposed by CRAN.
#' data(DiNAMIC.Duo)
#' output = dataPrep(X=luadSubset,Y=NULL)
#' }
#'
#' @import biomaRt
#' @import curl
#' @export
dataPrep = function(X, Y=NULL, species=c("human","mouse"), ensemblHost="http://aug2020.archive.ensembl.org")
	{
	#Check connection to host
 	statusCode = curl::curl_fetch_memory(ensemblHost)$status_code
	if (statusCode != 200)
		{
		stop(paste(ensemblHost, "currently not available"))
		}

  	species = match.arg(species,c("human","mouse"))

  	if (species == "human") 
		{
    		ensembl = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl",
          		host = ensemblHost)
    		attribute = "hgnc_symbol"
  		} else if (species=="mouse") 
			{
    			ensembl = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl",
              		host = ensemblHost)
    			attribute = "mgi_symbol"
  			}


  	#Re-order the rows of X (and Y)
  	if (is.null(rownames(X)))
		{
    		stop("Row names of X are missing")
  		}
	
  	# Check if genes are repeated. If repeated, remove.
  	X = X[(! rownames(X) %in% names(which(table(rownames(X))>1))),]
  	X = X[order(rownames(X)),]
  	if (is.null(Y))
		{
    		geneX = rownames(X)
  		} else {
    			if (is.null(rownames(Y))) 
				{
      			stop("Row names of Y are missing")
    				}
    			#Check if genes are repeated. If repeated, remove.
    			Y = Y[(! rownames(Y) %in% names(which(table(rownames(Y))>1))),]
    			Y = Y[order(rownames(Y)),]
    			geneX = rownames(X)[(rownames(X) %in% rownames(Y))] # common genes in X and Y
    			X = X[(rownames(X) %in% geneX),]
    			Y = Y[(rownames(Y) %in% geneX),]
  			}

  	#Get the chromosomes and genomic positions for the genes in tcga using getLDS().  Note that
  	#the result contain chromosomes other than 1-22, X, and Y.  Thus the results need to be
  	#filtered.
  	posDT = data.frame(trimws(as.matrix(getLDS(attributes = attribute, values = geneX, mart = ensembl,
                                             attributesL = c("chromosome_name", "start_position", "end_position", "band"), martL = ensembl, uniqueRows = TRUE))))
  	colnames(posDT) = c("symbol","chr","start","end","cytoband")
  	posDT = posDT[(posDT$chr %in% as.character(c(1:22))),]
  	# Remove rows of posDT corresponding to repeated gene names
  	posDT = posDT[(! posDT$symbol %in% names(which(table(posDT$symbol) > 1))),]

  	#Compute the mean of the start and end position for each gene, then add it to posDT
  	mean.pos = floor(0.5*(as.numeric(posDT$start) + as.numeric(posDT$end)))
 	posDT = cbind(posDT,mean.pos)
  	colnames(posDT)[ncol(posDT)] = "mean_pos"
  	rownames(posDT) = posDT$symbol
  	posDT = posDT[,c("chr","start","end","mean_pos","cytoband")]

  	# Select common genes and reorder the genes according to their genomic position
  	genes = sort(intersect(geneX, rownames(posDT)))
  	posDT = posDT[(rownames(posDT) %in% genes),]
  	posDT = posDT[order(rownames(posDT)),]
  	pos.perm = order(as.numeric(posDT$chr),posDT$mean_pos)
  	posDT = posDT[pos.perm,]

  	X = X[(rownames(X) %in% genes),]
  	X = X[pos.perm,]
  	if (! is.null(Y))
		{
    		Y = Y[(rownames(Y) %in% genes),]
    		Y = Y[pos.perm,]
  		}

  	return(list(X=X,Y=Y,posDT=posDT))
	}
