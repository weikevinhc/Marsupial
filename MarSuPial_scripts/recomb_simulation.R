## simulations
#' Crossover simulation generating allele frequencies
#'
#' \code{COsim} Simulates crossovers based on recombination rate in a r2d.locus object and
#' different number of individuals (n) and number of trials after marker selection. Assumes
#' no crossover interference.
#'
#' 
#' @param r2d List. Output from r2d.locus.
#' @param fitness Numeric. Fitness differential. Must be between -1 and 1.
#' @param n Numeric. Expected allele frequency of the selected locus.
#' @param trials Numeric. Number of trials. For a table of AF in individuals, n = 1, and trials = number of individuals.
#' @return Returns a bolean matrix. Each row is the averaged allele frequency of a trial with n number individuals.
#'   each column is a position in the $pos element of the r2d object. True and False cells represent different alleles
#' 
#' @examples 
#' example <- r2d.locus(function(x) -0.03*x^2 + 0.6*x + 1.15, 3, 22.422827, start = 1.22, end = 21.21)
#' sim.matrix <- COsim(example, 1, 1000, 10)
#' plot(example$pos, sim.matrix[1,], type = "l")
#' sapply(2:nrow(sim.matrix), function(x){points(example$pos, sim.matrix[x,], type = "l")})
#'
#' @export
COsim <- function(r2d, fitness, n = 1, trials = 100) {
	if (fitness < -1 | fitness > 1) {
		stop("fitness must be between -1 and 1")
	}
	locdex <- sum(r2d$pos <= r2d$l)
	winprob <- r2d$r/100*r2d$winsize
	wtot <- length(winprob)
	simm <- matrix(nrow = trials, ncol = wtot)
	ql <- 0.5 + 0.5*fitness
	wp <- matrix(winprob, nrow = n, ncol = wtot, byrow = T)
	for (j in 1:trials) {
		hits <- rbinom(n, 1, prob = ql) ## random draws of individuals with selected allele
		p <- matrix(runif(n * wtot), nrow = n, ncol = wtot) <= wp
		
		for (i in 1:n) {
			#	which(ind == 1)
			if (hits[i] == 1) {
				if (sum(p[i,]) > 0 ) {
					codex <- as.vector(unlist(sapply(seq(1,length(which(p[i,])), by = 2),
																					 function(x){seq(which(p[i,])[x], by = 1, length.out = c(which(p[i,]),wtot)[x+1] - which(p[i,])[x] + 1)})))
					if (locdex %in% codex) {
						p[i,] <- F
						p[i,codex] <- T
					} else {
						p[i,] <- T
						p[i,codex] <- F
					}
				} else {
					p[i,] <- T
				}
			} else {
				if (sum(p[i,]) > 0) {
					codex <- as.vector(unlist(sapply(seq(1,length(which(p[i,])), by = 2),
																					 function(x){seq(which(p[i,])[x], by = 1, length.out = c(which(p[i,]),wtot)[x+1] - which(p[i,])[x] + 1)})))
					if (locdex %in% codex) {
						p[i,] <- T
						p[i,codex] <- F
					} else {
						p[i,] <- F
						p[i,codex] <- T
					}
				} else {
					p[i,] <- F
				}
			}
		}
		#	points(afsim, type = "l", col = "red")
		simm[j,] <- sapply(1:ncol(p), function(x){sum(p[,x])})/n
	}
	return(simm)
}

#' Crossover simulatin for individual genotypes
#'
#' \code{COsim.ind} Simulates crossovers based on recombination rate in a r2d.locus object and
#' different number of individuals (n). Assumes no crossover interference.
#'
#' 
#' @param r2d List. Output from r2d.locus.
#' @param fitness Numeric. Fitness differential. Must be between -1 and 1.
#' @param n Numeric. Expected allele frequency of the selected locus.
#' @return Returns a matrix. Each row is the genotype of an individuals with True and False as the two allele states.
#'   each column is a position in the $pos element of the r2d object. 
#' 
#' @examples 
#' example <- r2d.locus(function(x) -0.03*x^2 + 0.6*x + 1.15, 3, 22.422827, start = 1.22, end = 21.21)
#' genotype.matrix <- COsim.ind(example, 1, 100)
#' plot(example$pos, colSums(genotype.matrix), type = "l")
#' 
#' @export
COsim.ind <- function(r2d, fitness, n = 1) {
	if (fitness < -1 | fitness > 1) {
		stop("fitness_diff must be between -1 and 1")
	}
	locdex <- sum(r2d$pos <= r2d$l)
	winprob <- r2d$r/100*r2d$winsize
	wtot <- length(winprob)
	ql <- 0.5 + 0.5*fitness
	wp <- matrix(winprob, nrow = n, ncol = wtot, byrow = T)
	hits <- rbinom(n, 1, prob = ql) ## random draws of individuals with selected allele
	p <- matrix(runif(n * wtot), nrow = n, ncol = wtot) <= wp
	
	for (i in 1:n) {
		#	which(ind == 1)
		if (hits[i] == 1) {
			if (sum(p[i,]) > 0 ) {
				codex <- as.vector(unlist(sapply(seq(1,length(which(p[i,])), by = 2),
																				 function(x){seq(which(p[i,])[x], by = 1, length.out = c(which(p[i,]),wtot)[x+1] - which(p[i,])[x] + 1)})))
				if (locdex %in% codex) {
					p[i,] <- F
					p[i,codex] <- T
				} else {
					p[i,] <- T
					p[i,codex] <- F
				}
			} else {
				p[i,] <- T
			}
		} else {
			if (sum(p[i,]) > 0) {
				codex <- as.vector(unlist(sapply(seq(1,length(which(p[i,])), by = 2),
																				 function(x){seq(which(p[i,])[x], by = 1, length.out = c(which(p[i,]),wtot)[x+1] - which(p[i,])[x] + 1)})))
				if (locdex %in% codex) {
					p[i,] <- T
					p[i,codex] <- F
				} else {
					p[i,] <- F
					p[i,codex] <- T
				}
			} else {
				p[i,] <- F
			}
		}
	}
	return(p)
}


#' Read counts and AF simulation based on recombinant fraction
#'
#' \code{D2WGSsim} Simulates the allele frequeny based on the expected 
#' recombinant fraction. This simulation is much faster that simulating
#' individual crossovers and can accomodate values of D with crossover interference.
#' But it assumes (unrealistically) that allele frequency at different sites are 
#' indepndent, thus underestimating noise that results from genotype sampling.
#' 
#' @param D Numeric vector. Vector of recombinant fraction.
#' @param fitness Numeric. Fitness differential. Must be between -1 and 1.
#' @param mendel_rate Numeric. The expected mendelian ratio if no selection was done.
#'   If should be one of 0.25, 5, and 0.75. If haploid, defaults to 0.5. This rate informs the program
#'   how to include the genetic contribution of the father or inbred parent.
#' @param n Numeric. Expected allele frequency of the selected locus.
#' @return Returns a list with two elements "counts" and "dp", corresponding to vectors of the number
#'   read counts of the selected allele and the depth for each element in D. 
#' 
#' @export
D2WGSsim <- function(D, fitness, n, dp, mendel_rate = 0.5) {
	s <- fitness ## fitness differential between alleles at locus of selection, between -1 and 1.
	if (s < -1 | s > 1) {
		stop("fitness_diff must be between -1 and 1")
	}
	ql <- 0.5 + 0.5*s
	q <- (1-D)*ql + D*(1-ql)
	qdraw <- rbinom(prob = q, size = n, n = length(q))/n ## binomial sample of genotypes given allele frequency
	if (mendel_rate == 0.5) {
		
	} else if (mendel_rate == 0.75) {
		qdraw <- qdraw/2 + 0.5
	} else if (mendel_rate == 0.25) {
		qdraw <- qdraw/2
	} else {
		stop("mendel_rate must be one of 0.25, 0.5, 0.75")
	}
	poisize <- rpois(n = length(q), lambda = dp) ## poisson draw for read counts at each position
	counts <- rbinom(prob = qdraw, size = poisize, n = length(q))
	return(list("counts" = counts, "dp" = poisize))
}


#' Read count simulation based on Allele frequency
#'
#' \code{AF2WGSsim} Simulates read counts given allele frequency. Read counts are
#'  simulated with binomial draws.
#' 
#' @param AF Numeric vector. Vector of allele frequencies.
#' @param dp Numeric vector. If only one element, it is used as average depth of coverage of a poisson
#'   distribution. If length of dp is the same as length of AF, the simulation assumes each AF will have 
#'   the corresponding dp. 
#' @return Returns a list with two elements "counts" and "dp", corresponding to vectors of the number
#'   read counts of the selected allele and the depth for each element in D. 
#' 
#' @export
AF2WGSsim <- function(AF, dp) {
	
	if (length(dp) == 1) {
		dp <- rpois(n = length(AF), lambda = dp) ## poisson draw for read counts at each position
	} else if (length(dp) == length(AF)) {
	} else {
		stop("dp must have length of 1 or the same length as AF")
	}
	counts <- rbinom(prob = AF, size = dp, n = length(AF))
	return(list("counts" = counts, "dp" = dp))
}