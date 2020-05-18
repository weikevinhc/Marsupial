## Allele frequency inference and fitting

#' LOESS fit over allele frequency
#'
#' \code{AFfitloess} Takes the position and allele frequency at diagnostic SNP sites,
#' and returns a list that contains two loess fit objects for sites to the left and right of
#' the selected locus.
#' 
#' 
#' @param sites Numeric vector. The positions of diagnostic sites (in Mb).
#' @param freq Numeric vector. Allele frequency at each diagnostic sites; must be the same length as sites.
#' @param depth Numeric vector. The depths at each diagnostic sites; must be the same length as sites.
#' @param locus Numeric. Position of the selected locus (in Mb).
#' @param locus_freq Numeric. Expected allele frequency of the selected locus.
#' @param span Numeric. The loess smoothing factor. Defaults to 0.2. 
#' Lower span increases resolution at the cost of increased noise.
#' @return Returns a list with five elements (left, right, AF, AF_pos,locus, locus_freq):
#'   $left is the fit object for positions to the left of locus
#'   $right is the fit object for position to the right of locus
#'   $AF is the predicted allele frequency at at $AF_pos
#'   $AF_pos is the ordered positions of the diagnostic sites with the addition of the selected locus
#'   $locus is the position of the selected locus
#'   $locus_freq is the expected frequency of the selected locus 
#' 
#' @export
AFfitloess <- function(sites, freq, depth, locus, locus_freq, span = 0.2, degree = 2) {
	bolleft <- sites < locus
	bolright <- sites > locus
	spleft <- span/(sum(bolleft)/length(sites))
	spright <- span/(sum(bolright)/length(sites))
	locusweight <- 1000000
	loesleft <- loess(c(freq[bolleft], locus_freq) ~ c(sites[bolleft], locus), span = spleft, weights = c(depth[bolleft], locusweight), degree = degree)
	loesright <- loess(c(locus_freq, freq[bolright]) ~ c(locus, sites[bolright]), span = spright, weights = c(locusweight, depth[bolright]), degree = degree)
	return(list("left" = loesleft,
							"right" = loesright,
							"AF" = c(predict(loesleft, sites[bolleft]), locus_freq, predict(loesright, sites[bolright])),
							"AF_pos" = c(sites[bolleft], locus, sites[bolright]),
							"locus" = locus,
							"locus_freq" = locus_freq))
}


#' Cubic splines fit over allele frequency
#'
#' \code{AFfitsplines} Takes the position and allele frequency at diagnostic SNP sites,
#' and returns a list that contains two cubic spline fit objects for sites to the left and right of
#' the selected locus.
#' 
#' 
#' @param sites Numeric vector. The positions of diagnostic sites (in Mb).
#' @param freq Numeric vector. Allele frequency at each diagnostic sites; must be the same length as sites.
#' @param depth Numeric vector. The depths at each diagnostic sites; must be the same length as sites.
#' @param locus Numeric. Position of the selected locus (in Mb).
#' @param locus_freq Numeric. Expected allele frequency of the selected locus.
#' @param df Numeric. Degrees of freedom. Defaults to value determined by crossvalidation. 
#' Lower df increases resolution at the cost of increased noise. 
#' @return Returns a list with five elements (left, right, AF, AF_pos,locus, locus_freq):
#'   $left is the fit object for positions to the left of locus
#'   $right is the fit object for position to the right of locus
#'   $AF is the predicted allele frequency at at $AF_pos
#'   $AF_pos is the ordered positions of the diagnostic sites with the addition of the selected locus
#'   $locus is the position of the selected locus
#'   $locus_freq is the expected frequency of the selected locus 
#' 
#' @export
AFfitsplines <- function(sites, freq, depth, locus, locus_freq, df = "cv") {
	bolleft <- sites < locus
	bolright <- sites > locus
	locusweight <- 1000000
	if (df == "cv") {
		splineleft <- smooth.spline(c(sites[bolleft], locus), c(freq[bolleft], locus_freq),  w = c(depth[bolleft], locusweight))
		splineright <- smooth.spline(c(locus, sites[bolright]) , c(locus_freq, freq[bolright]),  w = c(locusweight, depth[bolright]))
		return(list("left" = splineleft,
								"right" = splineright,
								"AF" = c(predict(splineleft, sites[bolleft]), locus_freq, predict(splineright, sites[bolright])),
								"AF_pos" = c(sites[bolleft], locus, sites[bolright]),
								"locus" = locus,
								"locus_freq" = locus_freq))
	} else {
		kleft <- round(sum(bolleft)/length(sites)*df)
		kright <- round(sum(bolright)/length(sites)*df)
		splineleft <- smooth.spline(c(sites[bolleft], locus), c(freq[bolleft], locus_freq),  w = c(depth[bolleft], locusweight), df = kleft)
		splineright <- smooth.spline(c(locus, sites[bolright]) , c(locus_freq, freq[bolright]),  w = c(locusweight, depth[bolright]), df = kright)
		return(list("left" = splineleft,
								"right" = splineright,
								"AF" = c(predict(splineleft, sites[bolleft]), locus_freq, predict(splineright, sites[bolright])),
								"AF_pos" = c(sites[bolleft], locus, sites[bolright]),
								"locus" = locus,
								"locus_freq" = locus_freq))
	}
	
#	plot(sites, freq, pch = 20, cex = 0.1)
#	points(splineleft, type = "l", col = "red")
#	points(splineright, type = "l", col = "red")
	
#	fit.mgcvleft <- gam(c(freq[bolleft], locus_freq)~s(c(sites[bolleft], locus), bs="ps"), weights = c(depth[bolleft], locusweight), )
#	fit.mgcvright <- gam(c(locus_freq, freq[bolright])~s(c(locus, sites[bolright]), bs="ps"), weights = c(locusweight, depth[bolright]))
#	plot(sites, freq, pch = 20, cex = 0.1)
#	points(c(sites[bolleft], locus), predict(fit.mgcvleft), type = "l", col = "red")
#	points(c(locus, sites[bolright]), predict(fit.mgcvright), type = "l", col = "red")
	
}


#' Linear regressions over allele frequency in windows
#'
#' \code{AFwinlm} Takes the position and allele frequency at diagnostic SNP sites,
#' and returns a list that contains two cubic spline fit objects for sites to the left and right of
#' the selected locus.
#' 
#' 
#' @param sites Numeric vector. The positions of diagnostic sites (in Mb).
#' @param freq Numeric vector. Allele frequency at each diagnostic sites; must be the same length as sites.
#' @param size Numeric. Size of the window (in Mb).
#' @param slide Numeric. Number of base pairs in which windows move (in Mb).
#' @return Returns a list with five elements (left, right, AF, AF_pos,locus, locus_freq):
#'   $left is the fit object for positions to the left of locus
#'   $right is the fit object for position to the right of locus
#'   $AF is the predicted allele frequency at at $AF_pos
#'   $AF_pos is the ordered positions of the diagnostic sites with the addition of the selected locus
#'   $locus is the position of the selected locus
#'   $locus_freq is the expected frequency of the selected locus 
#' 
#' @export
AFwinlm <- function(sites, freq, window = 0.5, slide = 0.1) {
	last <- sort(sites, decreasing = T)[1]
	size <- window/2
	v <- sapply(seq(0, last, by = slide), function(x){
		startend <- c(ifelse(x-size < 0, 0, x-size), ifelse((x+size-0.000001) > last, last, x+size-0.000001))
		winlm <- lm(formula = freq~pos, data = data.frame("pos" = sites[sites >= startend[1] & sites < startend[2]], 
																											"freq" = AF[sites >= startend[1] & sites < startend[2]]))
		return(unlist(c(x, predict(winlm, data.frame(pos = c(x))), as.numeric(winlm$coefficients[2]))))
	})
	return(list(pos = v[seq(1, length(v), by = 3)],
							freq = v[seq(2, length(v), by = 3)],
							slope = v[seq(3, length(v), by = 3)]))
}



#' Predict allele frequency at any position
#'
#' \code{predictAF} Predict allele frequency given AFfitloess or
#'   AFfitsplines for any position.
#' 
#' 
#' @param AFfit List. List object must be produced by AFfitloess
#' @param sites Numeric vector. Positions to be predicted.
#' @return Returns a numeric vector the same length as sites of the predicted allele frequency.
#' 
#' @export
predictAF <- function(AFfit, sites) {
	sitebol <- sites <= AFfit$locus
	if (is.list(predict(AFfit$left))) {
		##predict splines AF
		
		return(c(predict(AFfit$left, sites[sitebol])$y,
						 predict(AFfit$right, sites[!sitebol])$y))
	} else {
		##predict loess AF
		return(c(predict(AFfit$left, sites[sitebol]),
						 predict(AFfit$right, sites[!sitebol])))
	}
	
}

#' Predict the slope of allele frequency at any position
#'
#' \code{predictAF} Predict the slope of allele frequency given
#'   AFfitsplines for any position.
#' 
#' 
#' @param AFfit List. List object must be produced by AFfitsplines
#' @param sites Numeric vector. Positions to be predicted.
#' @return Returns a numeric vector the same length as sites of the predicted allele frequency.
#' 
#' @export
predictm <- function(AFfit, sites) {
	sitebol <- sites <= AFfit$locus
	return(c(predict(AFfit$left, sites[sitebol], deriv = 1)$y,
					 predict(AFfit$right, sites[!sitebol], deriv = 1)$y))
}



#' Estimate recombination rate from slope
#'
#' \code{d2rslope} Estimates recombination rate from slope of 
#'   genetic distance.
#' 
#' 
#' @param d Numeric vector. Vector of genetic distances
#' @param sites Numeric vector. Positions of the genetic distances.
#' @param locus Numeric. Position of the selected locus
#' @return Returns a numeric vector of length(sites) - 1 of the 
#'   slopes of distances.
#' 
#' @export
d2rslope <- function(d, sites, locus) {
	win <- sites[-1] - sites[-length(sites)]
	dif <- (d[-1] - d[-length(d)])/win
	dif[which(sites[-length(sites)] <= locus)] <- -dif[which(sites[-length(sites)] <= locus)]
	return(dif)
}