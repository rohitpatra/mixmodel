#' A function to estimate  mixing proportion and F.s (the cdf for the unknown component) given a data vector from a two component mixture model. This also provides inference for the mixing proportion.
#'
#'
#' @param data the numeric vector. Ideally data from a two component mixture model.
#' @param  method  A  string indicating which method to use and what is the goal. (a) `fixed` : it will compute the estimate based on the value of c.n; (b) `cv' : use cross-validation for choosing c.n (tuning parameter); (c) `lwr.bnd` : computes 95\% lower confidence bound for the mixing proportion.
#' @param  c.n  a positive number, default value is `.1 log(log(n)`, where `n` is the length of the vector.
#' @param  folds Number of folds used for cross-validation, default is 10.
#' @param  reps Number of replications for cross-validation, default is 1.
#' @param  cn.s  A vector of numbers to be used for cross-validation, a vector of values. default is equally spaced grid of 100 values between `.001 * log(log(n))` to `.2* log(log(n))`, where n is length(data).
#' @param  cn.length A number of equally spaced tuning parameter (between `.001 * log(log(n))` and `.2* log(log(n))`) values to search from, default is 100.
#' @param  gridsize An integer. Number of equally spaced points (between 0 and 1) to evaluate the distance function. Larger values are more computationally intensive but also lead to more accurate estimates, default is 600.
#' @return An object of class `mix.model`, a list including the elements:
#' \item{alp.hat}{estimate of alpha. NULL if `method' is `lwr.bnd'.}
#' \item{alp.Lwr}{95\% lower confidence bound of alpha. NULL if `method` is `cv` or `fixed`.}
#' \item{Fs.hat}{A list containing with elements `x` and `y` values for the function estimate of `F.s`}
#' \item{dist.out}{An object of the class `dist.calc` using the complete data.gen}
#' \item{c.n}{ Value of the tuning parameter used to compute the final estimate.}
#' \item{cv.out}{An object of class `cv.mix.model`. The object is NULL if method is `fixed`.}
#' @importFrom graphics abline
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom stats approx
#' @importFrom stats ecdf
#' @export
#' @examples
#' data.gen<- function( n, alpha){
#'   ind<- rbinom(n, 1, alpha)
#'   data<- ind* rbeta(n, 1,5) + (1-ind)* runif(n,0,1)
#'   return(data)
#' }
#' data.1<- data.gen(n = 5000, alpha = .1)
#' est.lwr.bnd <- mix.model(data.1, method = "lwr.bnd", gridsize = 600)

mix.model <-function(data, method=c("lwr.bnd", "fixed", "cv"),   c.n = NULL, folds =10, reps = 1, cn.s =NULL, cn.length =100, gridsize= 600){
	if(!is.vector(data)) stop(" 'data' has to be a vector.")

	n<- length(data)
	if(is.null(method)){
		stop("'method' can not be NULL")
	}
	if(method == 'lwr.bnd'){
		dist.out <- dist.calc(data, gridsize = gridsize)
		q<- 0.6792
		alp.Lwr <- sum(dist.out$distance > q/ sqrt(n))/gridsize
		alp.hat <- NULL
		Fs.hat.fun <- NULL
		c.n <- NULL
	} else{
		alp.Lwr<- NULL
	}
	if(method == "fixed"){
		if(is.null(c.n)){
			warning("'c.n' is not given. Fixing it to be '0.1*log(log(n))")
			c.n <- 0.1 *log(log(n))
		}
		dist.out <- dist.calc(data, gridsize = gridsize)
		alp.hat <- sum(dist.out$distance > c.n/ sqrt(n))/gridsize
		if(alp.hat >0){
			F.hat <- (dist.out$F.n-(1-alp.hat)*dist.out$F.b)/alp.hat ## Computes the naive estimator of F_s
			Fs.hat=Iso::pava(F.hat,dist.out$Freq,decreasing=FALSE) ## Computes the Isotonic Estimator of F_s
			Fs.hat[which(Fs.hat<=0)]=0
			Fs.hat[which(Fs.hat>=1)]=1
		} else {
			Fs.hat = dist.out$F.b
		}

		Fs.hat.fun<- NULL
		Fs.hat.fun$y = Fs.hat
		Fs.hat.fun$x =  dist.out$F.n.x
	} else if(method == "cv"){
		out.cv <- cv.mix.model(data, folds = folds, reps = reps, cn.s = cn.s, cn.length =cn.length , gridsize = gridsize)
		alp.hat <- out.cv$alp.hat
		Fs.hat.fun <- out.cv$Fs.hat
		dist.out <- out.cv$dist.out
		c.n <- out.cv$cn.cv
	}

	ret <- list(alp.hat = alp.hat,
				Fs.hat = Fs.hat.fun,
				dist.out = dist.out,
				c.n = c.n,
				alp.Lwr =alp.Lwr,
				n = n)
	if(method == "cv"){
		ret$cv.out <-  out.cv
	} else{
		ret$cv.out <- NULL
	}
	ret$method = method
	ret$call <- match.call()

	class(ret) <- "mix.model"
	return(ret)
}
#'@rdname mix.model
#' @export
print.mix.model <- function(x, ...){
	cat("Call:")
	print(x$call)
	if(x$method != "lwr.bnd"){
		print(paste("Estimate of alp is" , x$alp.hat))
		print(paste(" The chosen value c_n is", x$c.n))
		if( !is.null(x$cv.out)){
			par(mfrow=c(1,2))
			plot(x$cv.out)
		}
		plot(x$dist.out)
	} else if(x$method == 'lwr.bnd'){
		plot(x$dist.out)
		print (paste("The  '95%' lower confidence for alp_0 is ", x$alp.Lwr))
	}
}
#'@rdname mix.model
#' @export
plot.mix.model <- function(x, ...){
	if(x$method != "lwr.bnd"){
		plot(x$dist.out)
		abline(h= {x$c.n /sqrt(x$n)}, col="red", lwd= 1.3)
		abline(v= x$alp.hat, lty=2, lwd =1, col="black")
	} else {
		plot(x$dist.out)
		abline(h = 0.6792/sqrt(x$n), col="red", lwd= 1.3)
		abline(v= x$alp.Lwr, lty=2, lwd =1, col="black")
	}
}

