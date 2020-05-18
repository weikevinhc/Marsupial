##Mathematical properties of r, d, D and AF
##recomb rate (r) and genetic distance (d) conversion given a locus

#' Recombination rate to genetic distance from locus
#'
#' \code{r2d.locus} returns a list the chromosome-wide windows and the recombination rate
#' and genetic distance with recombination rate in cM/Mb as a function of chromosome position.
#' 
#' 
#' @param rate Function or dataframe. A R function object for recombination rate as as function of 
#'   chromosome position, or a bedfile of rates in different intervals. Object must be in the form of e.g. function(x) x^2.
#'   If recombination rate is static (e.g. 4cM/Mb across the chromosome) use function(x) x*0 + 4. 
#'   Bed file must have four columns with chr, start, end, and rate.
#' @param l Numeric. The position of the locus under selection in Mb.
#' @param size Numeric. The size of the chromosome in Mb.
#' @param pos Numeric or numeric vector. If a single number in Mb, the number is used as the non-overlapping window size to 
#'   construct a postional vector from zero to size of chromosome. If a numeric vector, numbers in Mb are used as positions. 
#'   Defaults to 0.01.
#' @param start Numeric. Start position on the chromosome in which recombination rate is measured.
#'   Positions less than start will have a rate of 0. Defaults to 0.
#' @param end Numeric. End position on the chromosome in which recombination rate is measured.  
#'   Positions greater than end will have a rate of 0. Defaults to chr_size.
#' @param bed Dataframe. Bed file in dataframe format with recombination rate (cM/Mb) in nonoverlapping windows.
#'   If a position is not covered, it will have a rate of 0. Coloumns should be chr, start, end, rate. 
#' @return Returns a list with three elements (r, pos, and d).
#'   r and d is the recombination rate and genetic distance from locus (cM) at each window in pos (Mb).
#' @examples 
#' example <- r2d.locus(function(x) -0.03*x^2 + 0.6*x + 1.15, 3, 22.422827, start = 1.22, end = 21.21)
#' plot(example$pos, example$d, type = "l")
#' @export

r2d.locus <- function(rate, l, size, pos = 0.01, start = 0, end = 0) {
	x <- c()
	winsize <- c()
	if (length(pos) > 1) {
		x <- pos
		winsize <- c(0, pos[-1] - pos[-length(pos)])
	} else {
		x <- seq(0, size, by = pos)
		winsize <- rep(pos, length(x))
	}
		
	if (is.matrix(rate) | is.data.frame(rate)) {
		
		r <- sapply(x, function(y){ 
			if (any(as.integer(y*1000000) >= rate[,2] & as.integer(y*1000000) <= rate[,3])) {
				return(rate[which(rate[,2] <= as.integer(y*1000000) & rate[,3] >= as.integer(y*1000000)),4])
#				return(1)
			} else {
				return(0)
			}
		})
		ldex <- tail(which(x < l), n = 1) + 1
		locusd <- rev(sapply((ldex-1):1, function(dex){sum(r[ldex:dex]*winsize[ldex:dex])}))
		locusd <- c(locusd, sapply(ldex:length(x), function(dex){sum(r[ldex:dex]*winsize[ldex:dex])}))
		return(list("pos" = x, "r" = r, "d" = locusd, "l" = l, "range" = c(x[1], size), "winsize" = winsize))
		
	} else if (is.function(rate)){
		if (end == 0) {
			end <- size
		}
		
		r <- rate(x)
		r[r < 0] <- 0
		r[x < start | x > end] <- 0
		locusd <- rep(0, length(r))
		dvec <- c(sapply(x[x >= start & x <= l], function (y) {integrate(rate, lower = y, upper = l)$value}),
							sapply(x[x > l & x <= end], function (y) {integrate(rate, lower = l, upper = y)$value}))
		locusd[!(x < start | x > end)] <- dvec
		locusd[x < start] <- dvec[1]
		locusd[x > end] <- tail(dvec, n = 1)
		return(list("pos" = x, "r" = r, "d" = locusd, "l" = l, "range" = c(start, end), "winsize" = winsize))
	} else {
		stop("rate must be a function or a bed file dataframe")
	}
}

#' Genetic distance to recombinant fraction conversion
#'
#' \code{d2D} returns a numeric vector of recombinant fraction.
#' 
#' 
#' @param d Numeric vector. Genetic distance (d)
#' @param method String. Name of mapping function. Either "haldane" or "kosambi".
#' @return Returns a numeric vector of recombinant fraction given genetic distance
#' @export
d2D <- function(d, method) {
	if (method == "haldane") {
		return((1 - exp(1)^(-2*d))/2)
	} else if (method == "kosambi") {
		return(0.5*(exp(1)^(4*d)-1)/(exp(1)^(4*d)+1))
	} else {
		stop("method must be 'haldane' or 'kosambi'")
	}
} 

#' Recombinant fraction to genetic distance conversion
#'
#' \code{D2d} returns a numeric vector of recombinant fraction.
#' 
#' 
#' @param D Numeric vector. Recombinant fraction
#' @param method String. Name of mapping function. Either "haldane" or "kosambi".
#' @return Returns a numeric vector of genetic distance given recombinant fraction.
#' @export
D2d <- function(D, method) {
	if (method == "haldane") {
		return(-log(1-2*D)/2)
	} else if (method == "kosambi") {
		return(1/4*log((1+2*D)/(1-2*D)))
	} else {
		stop("method must be 'haldane' or 'kosambi'")
	}
} 


#' Recombinant fraction to Allele Frequency
#'
#' \code{D2AF()} returns in a vector of the allele frequency given a vector of the recombinant fraction.
#' 
#' 
#' @param D Numeric vector. Recombinant fraction
#' @param fitness Numeric. Fitness differential at the selected locus, must be between -1 and 1.
#' @param ploidy String. The ploidy of the chromosome "haploid" or "diploid". 
#'   X-linked male selection should be haploid; otherwide diploid.
#' @param mendel_rate Numeric. The expected mendelian ratio if no selection was done.
#'   If should be one of 0.25, 5, and 0.75. If haploid, defaults to 0.5. This rate informs the program
#'   how to include the genetic contribution of the father or inbred parent.
#' @return Returns a numeric vector of the allele frequency given each recombinant fraction
#' @examples
#' D2AF(c(0, 0.2,0.23,0.24, 0.25), 1, ploidy = "diploid", mendel_rate = 0.25)
#' @export
D2AF <- function(D, fitness, ploidy = "haploid", mendel_rate = 0.5) {
	s <- fitness ## fitness differential between alleles at locus of selection, between -1 and 1.
	m <- mendel_rate ## mendelian rate should be 0.5 for haploid, and either 0.25 or 0.75 for diploid
	if (s < -1 | s > 1) {
		stop("fitness_diff must be between -1 and 1")
	}
	if (!(ploidy %in% c("haploid", "diploid"))) {
		stop("ploidy must be 'haploid' or 'diploid'")
	} else {
		if (ploidy == "haploid") {
			m <- 0.5
		}
	}
	if (mendel_rate %in% c(0.25, 0.5, 0.75)) {
	} else {
		stop("mendel_rate must be one of 0.25, 0.5, 0.75")
	}
	ql <- 0.5 + 0.5*s
	q <- (1-D)*ql + D*(1-ql)
	if (mendel_rate == 0.75) {
		q <- q/2 +0.5
	}
	if (mendel_rate == 0.25) {
		q <- q/2
	}
	return(q)
}

#' Convert allele frequency to recombinant fraction
#'
#' \code{AF2D} Covert allele frequency to recombinant fraction 
#' 
#' @param AF Numeric vector. Allele frequency
#' @param fitness Numeric. Fitness differential at the selected locus, must be between -1 and 1.
#' @param ploidy String. The ploidy of the chromosome "haploid" or "diploid". 
#'   X-linked male selection should be haploid; otherwide diploid.
#' @param mendel_rate Numeric. The expected mendelian ratio if no selection was done.
#'   If should be one of 0.25, 5, and 0.75. If haploid, defaults to 0.5. This rate informs the program
#'   how to include the genetic contribution of the father or inbred parent.
#' @return Returns a numeric vector of the recombinant fraction given genetic distance
#' @examples
#' AF2D(c(0.500, 0.400, 0.385, 0.380, 0.375), 1, ploidy = "diploid", mendel_rate = 0.25)
#' @export
AF2D <- function(AF, fitness, ploidy = "haploid", mendel_rate = 0.5) {
	s <- fitness ## fitness differential between alleles at locus of selection, between -1 and 1.
	m <- mendel_rate ## mendelian rate should be 0.5 for haploid, and either 0.25 or 0.75 for diploid
	if (s < -1 | s > 1) {
		stop("fitness_diff must be between -1 and 1")
	}
	if (!(ploidy %in% c("haploid", "diploid"))) {
		stop("ploidy must be 'haploid' or 'diploid'")
	} else {
		if (ploidy == "haploid") {
			m <- 0.5
		}
	}
	ql <- 0.5 + 0.5*s
	AFadj <- AF
	if (mendel_rate == 0.5) {
		
	} else if (mendel_rate == 0.75) {
		AFadj <- (AFadj - 0.5)*2
	} else if (mendel_rate == 0.25) {
		AFadj <- (AFadj)*2
	} else{
		stop("mendel_rate must be one of 0.25, 0.5, 0.75")
	}
	D <- (AFadj-ql)/(1-2*ql)
	return(D)
}

#' Recombination rate to genetic distance for two loci
#'
#' \code{r2d.3pt} returns a list the chromosome-wide windows and the 
#' genetic distance at two positions give recombination rate in cM/Mb 
#' as a function of chromosome position. this function is to create object
#' for two loci selection that create a effectively a 3point cross.
#' 
#' 
#' @param ratefunc A R function object for recombination rate as as function of 
#'   chromosome position. Object must be in the form of e.g. function(x) x^2.
#'   If recombination rate is static (e.g. 4cM/Mb) use function(x) x*0 + 4.
#' @param l1 Numeric. The position of the first/main locus under selection.
#' @param l2 Numeric. The position of the second locus under selection.
#' @param size Numeric. The size of the chromosome.
#' @param winsize Numeric. The non-overlapping window size. Defaults to 0.01Mb
#' @param start Numeric. Start position on the chromosome in which recombination rate is measured.
#'   Positions less than start will have a rate of 0. Defaults to 0.
#' @param end Numeric. End position on the chromosome in which recombination rate is measured.  
#'   Positions greater than end will have a rate of 0. Defaults to chr_size.
#' @return Returns a list with six elements (pos, d1, d2, d12, l1, and l2).
#'   $d1 and $d2 are the genetic distance from l1 and l2 (cM) at each window in pos (Mb).
#'   $d12 is the genetic distance between l1 and l2
#' @examples 
#' example <- r2d.locus(function(x) -0.03*x^2 + 0.6*x + 1.15, 3, 10, 22.422827, start = 1.22, end = 21.21)
#' plot(example$pos, example$d1, type = "l")
#' points(example$pos, example$d2, type = "l")
#' @export
r2d.3pt <- function(ratefunc, l1, l2, size, winsize = 0.01, start = 0, end = 0) {
	if (end == 0) {
		end <- size
	}
	locusd <- rep(0, length(r))
	dvec <- c(sapply(seq(start, l1, by = winsize), function (y) {integrate(ratefunc, lower = y, upper = l1)$value}),
						sapply(seq(l1 + winsize, end, by = winsize), function (y) {integrate(ratefunc, lower = l1, upper = y)$value}))
	locusd[!(x < start | x > end)] <- dvec
	locusd[x < start] <- dvec[1]
	locusd[x > end] <- tail(dvec, n = 1)
	d1 <- locusd
	locusd <- rep(0, length(r))
	dvec <- c(sapply(seq(start, l2, by = winsize), function (y) {integrate(ratefunc, lower = y, upper = l2)$value}),
						sapply(seq(l2 + winsize, end, by = winsize), function (y) {integrate(ratefunc, lower = l2, upper = y)$value}))
	locusd[!(x < start | x > end)] <- dvec
	locusd[x < start] <- dvec[1]
	locusd[x > end] <- tail(dvec, n = 1)
	d2 <- locusd
	d12 <- integrate(ratefunc, lower = l1, upper = l2)$value
	
	return(list("pos" = x, "d1" = d1, "d2" = d2, "d12" = d12, "l1" = l1, "l2" = l2 ))
}


#' Recombinant fraction to allele frequency in two loci selection
#'
#' \code{d2AF.3pt()} returns a numeric vector of alelle frequency give a r2D.3pt object.
#' This function calculates the chromosome-wide allele frequency when two loci are under selection.
#' 
#' 
#' @param list3pt List. Object created by r2D.3pt
#' @param fitness1 Numeric. Fitness differential for locus 1. Between -1 and 1. Defaults to 1.
#' @param fitness2 Numeric. Fitness differential for locus 2. Between -1 and 1.
#' @param method String. Name of mapping function. Either "haldane" or "kosambi".
#' @return Returns a numeric vector of allele frequency of gametes.
#' @export
d2AF.3pt= function(list3pt, fitness1 = 1, fitness2, method = "haldane") {
	pos <- list3pt$pos
	Dil <- d2D(list3pt$d1/100, method = method)
	Dio <- d2D(list3pt$d2/100, method = method)
	Dlo <- d2D(list3pt$d12/100, method = method)
	locusl <- list3pt$l1
	locus2 <- list3pt$l2
	sl <- fitness1 ## strength of selection -1 to 1
	so <- fitness2
	ql <- 0.5 + 0.5*sl ## proportion of gametes
	qo <- 0.5 + 0.5*so ## proportion of gametes
	fqilo <- c()
	fqlio <- c()
	fqloi <- c()
	if (method == "haldane") {
		fqilo <- ((1-Dil)*(1-Dlo)*ql*qo + (1-Dil)*(Dlo)*ql*(1-qo) +
								(Dil)*(1-Dlo)*(1-ql)*(1-qo) + (Dil)*(Dlo)*(1-ql)*qo)/
			((1-Dil)*(1-Dlo)*ql*qo + (1-Dil)*(Dlo)*ql*(1-qo) +
			 	(Dil)*(1-Dlo)*(1-ql)*(1-qo) + (Dil)*(Dlo)*(1-ql)*qo +
			 	(1-Dil)*(1-Dlo)*(1-ql)*(1-qo) + (1-Dil)*(Dlo)*(1-ql)*(qo) +
			 	(Dil)*(1-Dlo)*(ql)*(qo) + (Dil)*(Dlo)*(ql)*(1-qo))
		fqlio <- ((1-Dil)*(1-Dio)*ql*qo + (Dil)*(1-Dio)*(1-ql)*qo +
								(1-Dil)*(Dio)*ql*(1-qo) + (Dil)*(Dio)*(1-ql)*(1-qo))/
			((1-Dil)*(1-Dio)*ql*qo + (Dil)*(1-Dio)*(1-ql)*qo +
			 	(1-Dil)*(Dio)*ql*(1-qo) + (Dil)*(Dio)*(1-ql)*(1-qo) +
			 	(1-Dil)*(1-Dio)*(1-ql)*(1-qo) + (Dil)*(1-Dio)*(ql)*(1-qo) +
			 	(1-Dil)*(Dio)*(1-ql)*(qo) + (Dil)*(Dio)*(ql)*(qo))
		fqloi <- ((1-Dio)*(1-Dlo)*qo*ql + (1-Dio)*(Dlo)*qo*(1-ql) +
								(Dio)*(1-Dlo)*(1-qo)*(1-ql) + (Dio)*(Dlo)*(1-qo)*ql)/
			((1-Dio)*(1-Dlo)*qo*ql + (1-Dio)*(Dlo)*qo*(1-ql) +
			 	(Dio)*(1-Dlo)*(1-qo)*(1-ql) + (Dio)*(Dlo)*(1-qo)*ql +
			 	(1-Dio)*(1-Dlo)*(1-qo)*(1-ql) + (1-Dio)*(Dlo)*(1-qo)*(ql) +
			 	(Dio)*(1-Dlo)*(qo)*(ql) + (Dio)*(Dlo)*(qo)*(1-ql))
		return(c(fqilo[pos <= locusl], fqlio[pos > locusl & pos <= locus2], fqloi[pos > locus2]))
	} else if (method == "kosambi") {
		Da <- Dil[pos <= locusl]
		Db <- Dlo
		Dab <- Dio[pos <= locusl]
		dcab <- (Da + Db - Dab)
		Dab_t <- (Da + Db)/(1+4*Da*Db)
		Cc <- 2*Dab
		Ca <- (1-Db*Cc)/(1-Db)
		Cb <- (1-Da*Cc)/(1-Da)
		fqilo <- ((1-Da)*(1-Db)*ql*qo + (1-Da)*(Db)*ql*(1-qo)*Cb +
								(Da)*(1-Db)*(1-ql)*(1-qo)*Ca + Da*Db*(1-ql)*qo*Cc)/
			((1-Da)*(1-Db)*ql*qo + (1-Da)*(Db)*ql*(1-qo)*Cb +
			 	(Da)*(1-Db)*(1-ql)*(1-qo)*Ca + Da*Db*(1-ql)*qo*Cc +
			 	(1-Da)*(1-Db)*(1-ql)*(1-qo) + (1-Da)*(Db)*(1-ql)*(qo)*Cb +
			 	(Da)*(1-Db)*(ql)*(qo)*Ca + Da*Db*(ql)*(1-qo)*Cc)
		
		Da <- Dil[pos > locusl & pos <= locus2]
		Db <- Dio[pos > locusl & pos <= locus2]
		Dab <- Dlo
		dcab <- (Da + Db - Dab)
		Dab_t <- (Da + Db)/(1+4*Da*Db)
		Cc <- 2*Dab
		Ca <- (1-Db*Cc)/(1-Db)
		Cb <- (1-Da*Cc)/(1-Da)
		fqlio <- ((1-Da)*(1-Db)*ql*qo + (Da)*(1-Db)*(1-ql)*qo*Ca +
								(1-Da)*(Db)*ql*(1-qo)*Cb + Da*Db*(1-ql)*(1-qo)*Cc)/
			((1-Da)*(1-Db)*ql*qo + (Da)*(1-Db)*(1-ql)*qo*Ca +
			 	(1-Da)*(Db)*ql*(1-qo)*Cb + Da*Db*(1-ql)*(1-qo)*Cc +
			 	(1-Da)*(1-Db)*(1-ql)*(1-qo) + (Da)*(1-Db)*(ql)*(1-qo)*Ca +
			 	(1-Da)*(Db)*(1-ql)*(qo)*Cb + Da*Db*(ql)*(qo)*Cc)
		
		Da <- Dlo
		Db <- Dio[pos > locus2]
		Dab <- Dil[pos > locus2]
		dcab <- (Da + Db - Dab)
		Dab_t <- (Da + Db)/(1+4*Da*Db)
		Cc <- 2*Dab
		Ca <- (1-Db*Cc)/(1-Db)
		Cb <- (1-Da*Cc)/(1-Da)
		fqloi <- ((1-Da)*(1-Db)*ql*qo + (Da)*(1-Db)*(1-ql)*qo*Ca +
								(1-Da)*(Db)*(1-ql)*(1-qo)*Cb + (Da)*(Db)*ql*(1-qo)*Cc)/
			((1-Da)*(1-Db)*ql*qo + (Da)*(1-Db)*(1-ql)*qo*Ca +
			 	(1-Da)*(Db)*(1-ql)*(1-qo)*Cb + (Da)*(Db)*ql*(1-qo)*Cc +
			 	(1-Da)*(1-Db)*(1-ql)*(1-qo) + (Da)*(1-Db)*(ql)*(1-qo)*Ca +
			 	(1-Da)*(Db)*(ql)*(qo)*Cb + (Da)*(Db)*(1-ql)*(qo)*Cc)
		return(c(fqilo, fqlio, fqloi))
	} else {
		stop("method must be 'haldane' or 'kosambi'")
	}
	
}
