#' A Function for Plotting Mean Copy Number Values and Differences Across the Genome
#'
#' This function plots mean copy number values from one or two cohorts at a common set of markers across the genome.
#'
#' @param inputList A list produced by dataPrep.
#'
#' @param lwdVec A vector of line widths.  Default = rep(1, 3).  See \code{\link{par}}.
#'
#' @param ltyVec A vector of line types.  Default = c(1:3).  See \code{\link{par}}.
#'
#' @param lineColorVec  A vector of line colors.  Default = c("red", "blue", "black").  See \code{\link{par}}.
#'
#' @param ylimLow  The lower limit of the y-values in the plot.  Default = -1.  See \code{\link{plot}}.
#'
#' @param ylimHigh  The upper limit of the y-values in the plot.  Default = 1.  See \code{\link{plot}}.
#'
#' @param chrLabel  Binary value determining whether or not chromosomes are labeled.  Default = TRUE.
#'
#' @param xaxisLabel  Label for the x-axis in the plot.  Default = "Chromosome".  See \code{\link{plot}}.
#'
#' @param yaxisLabel  Label for the y-axis in the plot.  Default = NULL.  See \code{\link{plot}}.
#'
#' @param mainLabel  Main label in the plot.  Default = NULL.  See \code{\link{plot}}.
#'
#' @param rectColors  Background colors for different chromosomes.  Default = c("light gray", "gray").
#'
#' @param axisCex  Point size for the scale on the axis.  Default = 1.  See \code{\link{par}}.
#'
#' @param labelCex  Point size for the axis label.  Default = 1.  See \code{\link{par}}.
#'
#' @param xaxisLine  Numerical value used to specify the location (line) of the x-axis label.  Default = 2.5.  See \code{\link{mtext}}.
#'
#' @param yaxisLine  Numerical value used to specify the location (line) of the y-axis label.  Default = 2.5.  See \code{\link{mtext}}.
#'
#' @param mainLine  Numerical value used to specify the location (line) of the main.label.  Default = 0.  See \code{\link{mtext}}.
#'
#' @param marginVec Numerical vector specifying margin sizes.  Default = c(4, 4, 3, 3).  See \code{\link{par}}.
#'
#' @param legendText Character vector used to legend.  Only shown if showLegend = TRUE.  Default = NULL.  See \code{\link{legend}}.
#'
#' @param highThreshold Numerical value representing the position of the upper horizontal line, e.g., a threshold for
#'   assessing statistical significance.  Default = NULL.
#'
#' @param lowThreshold Numerical value representing the position of the lower horizontal line, e.g., a threshold for
#'   assessing statistical significance.  Default = NULL.
#'
#' @param showLegend Binary value determining whether or not the legend is shown.  Default = FALSE.  See \code{\link{legend}}.
#'
#' @param legendXQuantile Quantile to specify the "x" location of the legend.  Only relevant if showLegend = TRUE  Default = 0.55.  See \code{\link{legend}}.
#'
#' @param legendYCoord Numerical value to specify the "y"location of the legend.  Only relevant if showLegend = TRUE.  Default = 1.  See \code{\link{legend}}.
#'
#' @importFrom graphics abline axis legend lines mtext par rect segments
#'
#' @importFrom stats quantile
#'
#' @return Creates a genomewide plot of mean copy number values and differences.
#'
#' @details This function is used to visualize copy number values and copy number alterations
#'   across the genome.  If Y = NULL in the input list, then the plot shows a single line
#'   corresponding to the mean DNA copy number values based on the entries in X.
#'   If both X and Y are specified, the plot shows three lines corresponding to the mean
#'   DNA copy number values in X, the mean DNA copy number values in Y, and the difference
#'   of the mean DNA copy number values.
#'
#'
#' @examples genomeChrPlot(inputList = pD, ylimLow = -1.4, ylimHigh = 1.4)
#'
#' @export

genomePlot = function(
	inputList,
	lwdVec = rep(1, 3),
	ltyVec = c(1:3),
	lineColorVec = c("red", "blue", "black"),
	ylimLow = -1,
	ylimHigh = 1,
	chrLabel = TRUE,
	xaxisLabel = "Chromosome",
	yaxisLabel = NULL,
	mainLabel = NULL,
	rectColors = c("light gray", "gray"),
	axisCex = 1,
	labelCex = 1,
	xaxisLine = 2.5,
	yaxisLine = 2.5,
	mainLine = 0,
	marginVec = c(4, 4, 3, 3),
	legendText = NULL,
	highThreshold = NULL,
	lowThreshold = NULL,
	showLegend = FALSE,
	legendXQuantile = 0.55,
	legendYCoord = 1
	)
	{
	#Define terms for plotting
	chr = as.numeric(inputList[["posDT"]]$chr)
	pos = as.numeric(inputList[["posDT"]]$mean_pos)
	chrPosPerm = order(chr, pos)
	m = length(chrPosPerm)
	dists = pos[chrPosPerm][2:m] - pos[chrPosPerm][1:m-1]
	dists[dists < 0] = 0
	cumDists = cumsum(dists)
	cumDists = c(0, cumDists)/1e6
	
	#Reset to default par() after exiting since par is redefined below
	oldPar = par(no.readonly = TRUE)
	on.exit(par(oldPar))

	par(mar = marginVec)
	plot(cumDists, 
		rep(0, length(cumDists)), 
		xlim = range(cumDists), 
		ylim = c(ylimLow, ylimHigh), 
		axes = FALSE, 
		ylab = "", 
		xlab = "", 
		type = "n")

	for (i in unique(chr))
		{
 		rect(min(cumDists[which(chr == i)]),
 			ylimLow,
			max(cumDists[which(chr == i)]),
   			ylimHigh,
			col = rectColors[1 + i%%2], 
			border = rectColors[1 + i%%2])
					
		if (chrLabel == T)
			{
			if ((i %% 2) == 1)
				{
				segments(min(cumDists[which(chr == i)]),
 					ylimLow[1],
					max(cumDists[which(chr == i)]),
   					ylimLow[1],
					lwd = 2)

				mtext(i,
					side = 1,
					line = -.5,
	 				at = .5 *(min(cumDists[which(chr == i)]) +
	      				max(cumDists[which(chr == i)])), 
					adj = .5)
				}
			}
		}

	plotVals1 = rowMeans(inputList[["X"]])
	if(!is.null(inputList[["Y"]]))
		{
		plotVals2 = rowMeans(inputList[["Y"]])
		} else plotVals2 = NULL
	if(!is.null(inputList[["Y"]]))
		{
		plotVals3 = plotVals1 - plotVals2
		} else plotVals3 = NULL

	lines(cumDists, rep(0, length(cumDists)), lty = 2, col = "gray50", lwd = lwdVec[1])	
	lines(cumDists, plotVals1, lwd = lwdVec[1], lty = ltyVec[1], col = lineColorVec[1])
	lines(cumDists, plotVals2, lwd = lwdVec[2], lty = ltyVec[2], col = lineColorVec[2])
	lines(cumDists, plotVals3, lwd = lwdVec[3], lty = ltyVec[3], col = lineColorVec[3])
	if (!is.null(highThreshold))
		{
		abline(h = highThreshold, lwd = 2, lty = 2, col = "gray50")
		}
	if (!is.null(lowThreshold))
		{
		abline(h = lowThreshold, lwd = 2, lty = 2, col = "gray50")
		}
	axis(2, cex.axis = axisCex)
	mtext(xaxisLabel, side = 1, cex = labelCex, line = xaxisLine)
	mtext(yaxisLabel, side = 2, cex = labelCex, line = yaxisLine)
	mtext(mainLabel, side = 3, cex = labelCex, line = mainLine)

	if (showLegend)
		{
		legend(x = quantile(cumDists, legendXQuantile), y = legendYCoord, lwd = lwdVec, lty = ltyVec, 
			col = lineColorVec, legend = legendText, bty = "n")
		}
	}



