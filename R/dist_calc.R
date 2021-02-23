#' Title
#'
#' @param data the numeric vector. Ideally data from a two component mixture model.
#' @param  gridsize An integer. Number of equally spaced points (between 0 and 1) to evaluate the distance function. Larger values are more computationally intensive but also lead to more accurate estimates, default is 600.
#'
#' @importFrom Iso pava
#' @return An object of class `dist.calc`, a list including the elements:
#' \item{distance}{A vector of the distance computed at a grid of alpha as created by `gridsize`.}
#' \item{grid.pts}{A vector of equally spaced points between 0 and 1 of length `gridsize`}
#' @export
#'
dist.calc <- function(data,gridsize=200){
	q <- 0.6792
	n <- length(data) ## Length of the data set
	data <- sort(data) ## Sorts the data set
	data.1 <- unique(data) ## Finds the unique data points
	Fn <- ecdf(data) ## Computes the empirical DF of the data
	Fn.1 <- Fn(data.1) ## Empirical DF of the data at the data points
	## Calculate the known F_b at the data points
	## Note: for Uniform(0,1) F_b(x) = x
	## Usually would need to CHANGE this
	Fb <- data.1
	## Compute the weights (= frequency/n) of the unique data values, i.e., dF_n
	Freq <- diff(c(0,Fn.1))
	distance <- rep(0,floor(gridsize*.12))
	distance[0]<- sqrt(t((Fn.1-Fb)^2)%*%Freq)

	grid.pts <- {1:floor(gridsize)} /gridsize
	for(ii in 1:length(grid.pts)){
		# print(i)
		aa <- grid.pts[ii]
		F.hat <- (Fn.1-(1-aa)*Fb)/aa ## Computes the naive estimator of F_s
		F.is=Iso::pava(F.hat,Freq,decreasing=FALSE) ## Computes the Isotonic Estimator of F_s
		F.is[which(F.is<=0)]=0
		F.is[which(F.is>=1)]=1
		distance[ii] <- aa*sqrt(t((F.hat-F.is)^2)%*%Freq);
	}
	# Lower.Cfd.Bound <- sum(distance>q/sqrt(n))/gridsize
	ret <- list( distance = distance,
				gridsize = gridsize,
				grid.pts = grid.pts,
				F.n = Fn.1,
			 	F.n.x = data.1,
				Freq = Freq,
				F.b = Fb,
				n = n)
	class(ret) <- "dist.calc"
	ret$call <- match.call()

	return(ret)
}
#'@rdname dist.calc
#' @export
print.dist.calc <- function(x,...){
	cat("Call:\n")
	print(x$call)
	t.mat<- cbind(x$grid.pts , x$distance)
	colnames(t.mat) = c("gamma", "distance")
	print(t(t.mat))
}
#'@rdname dist.calc
#' @export
plot.dist.calc <- function(x,...){
	plot(x$grid.pts, x$distance, ylab= "distance", xlab= "gamma" , type ="l")
}
