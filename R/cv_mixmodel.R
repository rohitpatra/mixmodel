

#' #' A function to estimate alp and F.s given a data vector via Cross validation
#'
#' @param data the numeric vector. Ideally data from a two component mixture model.
#' @param  folds Number of folds used for cross-validation, default is 10.
#' @param  reps Number of replications for cross-validation, default is 1.
#' @param  cn.s  A vector of numbers to be used for cross-validation, a vector of values. default is equally spaced grid of 100 values between `.001 * log(log(n))` to `.2* log(log(n))`, where n is length(data).
#' @param  cn.length A number of equally spaced tuning parameter (between `.001 * log(log(n))` and `.2* log(log(n))`) values to search from, default is 100.
#' @param  gridsize An integer. Number of equally spaced points (between 0 and 1) to evaluate the distance function. Larger values are more computationally intensive but also lead to more accurate estimates, default is 600.
#'
#' @return An object of class `cv.mix.model`, a list including the elements:
#' \item{alp.hat}{estimate of alpha.}
#' #' \item{cn.cv}{value of the chosen tuning parameter via cross validation for the final estimate of alpha. }
#' \item{Fs.hat}{A list containing with elements `x` and `y` values for the function estimate of `F.s`}
#' \item{dist.out}{An object of the class `dist.fun` using the complete data.gen}
#' \item{score}{ vector of cross validation scores at the input values of `cn.s`.}
#' @export
#'
#' @examples
#' data.gen<- function( n, alpha){
#'   ind<- rbinom(n, 1, alpha)
#'   data<- ind* rbeta(n, 1,5) + (1-ind)* runif(n,0,1)
#'   return(data)
#' }
#' data.1<- data.gen(n = 200, alpha = .1)
#' est.cv <- cv.mix.model(data.1, folds = 10, reps = 1, cn.length =20 , gridsize = 200)
cv.mix.model <- function(data, folds = 10, reps = 1, cn.s = NULL, cn.length = NULL, gridsize = 200){
	if(!is.vector(data)) stop(" 'data' has to be a vector.")
	warning("Make sure that data is transformed such that F_b is Uniform(0,1)")
	n<- length(data)

	# cvFolds(length(data), K =folds, R= reps, type= "random")

	if (is.null(cn.s)){
		if( is.null(cn.length)) stop("Both 'cn.s' and 'cn.length' can not be null")
		cn.s <- seq(.001 * log(log(n)), .2* log(log(n)), length.out =cn.length)
	} else if(length(cn.s) != cn.length){
		stop("Length of 'cn.s' is different  from cn.length")
	}
	score <-  rep(0, length(cn.s))
	for( rr in 1:reps){
		cv.ind <- rep(1:folds, times = ceiling(n/folds) )
		cv.ind <- cv.ind[sample.int(length(cv.ind))]
		cv.ind <- cv.ind[1:n]
		for( kk in 1:folds){
			t.1 <- Sys.time()
			t.data <- dist.calc(data[cv.ind != kk], gridsize = gridsize)
			for (cc in 1: length(cn.s)){
				score[cc] <- score[cc] +  cv.score(t.data, test.data = data[cv.ind == kk], c.n = cn.s[cc] )
			}
			t.2<- Sys.time()
			if(kk==1) {
				print("Expected time for completion")
				print( folds*reps*(t.2-t.1))
			}
		}
	}
	cn.cv <-  cn.s[which.min(score)]
	tot.out <- dist.calc(data, gridsize = gridsize)
	alp.hat.cv <- sum(tot.out$distance>cn.cv/sqrt(tot.out$n))/tot.out$gridsize
	if(alp.hat.cv >0){
		F.hat <- (tot.out$F.n-(1-alp.hat.cv)*tot.out$F.b)/alp.hat.cv ## Computes the naive estimator of F_s
		Fs.hat=Iso::pava(F.hat,tot.out$Freq,decreasing=FALSE) ## Computes the Isotonic Estimator of F_s
		Fs.hat[which(Fs.hat<=0)]=0
		Fs.hat[which(Fs.hat>=1)]=1
	} else {
		Fs.hat = tot.out$F.b
	}
	Fs.hat.fun<- NULL
	Fs.hat.fun$y = Fs.hat
	Fs.hat.fun$x = tot.out$F.n.x

	ret <-  list(cn.cv = cn.cv,
				 score = score,
				 cn.s = cn.s,
				 Fs.hat = Fs.hat.fun,
				 alp.hat = alp.hat.cv,
				 folds = folds,
				 rep = rep,
				 data =data,
				 dist.out = tot.out )
	ret$call <- match.call()
	class(ret)  <- "cv.mix.model"

	return(ret)
}
#'@rdname cv.mix.model
#' @export
print.cv.mix.model <- function(x,...){
	cat("Call:")
	print(x$call)
	print(paste("Cross validated estimate of alp is" , x$alp.hat))
	print(paste(" The cross-validated choice of c_n is", x$cn.cv))
}
#'@rdname cv.mix.model
#' @export
plot.cv.mix.model <- function(x,...){
	plot(x$cn.s, x$score, ylab= "cross validation score", xlab= "c_n" )
}




#' Computes the cross validataion score at c.n. Function used inside `cv.mix.model`
#'
#' @param tr.data Training data vector
#' @param test.data Testing data vector
#' @param c.n a positive number, default value is `.1 log(log(n)`, where `n` is the length of the vector.
#'
#' @return Numeric value of the cross validation score at the input
#' @export
#'
cv.score<- function(tr.data, test.data, c.n){
	if(class(tr.data) != "dist.calc") stop("'tr.data' is of the wrong type.")

	alp.hat <- sum(tr.data$distance>c.n/sqrt(tr.data$n))/tr.data$gridsize
	if(alp.hat >0){
		F.hat <- (tr.data$F.n-(1-alp.hat)*tr.data$F.b)/alp.hat ## Computes the naive estimator of F_s
		F.is=Iso::pava(F.hat,tr.data$Freq,decreasing=FALSE) ## Computes the Isotonic Estimator of F_s
		F.is[which(F.is<=0)]=0
		F.is[which(F.is>=1)]=1
		# Fhat.train = alp.hat * F.is + (1-alp.hat) * tr.data$F.b
	} else {
		F.is = tr.data$F.b
	}

	test.data <- sort(test.data) ## Sorts the data set
	test.data.1 <- unique(test.data) ## Finds the unique data points
	test.Fn <- ecdf(test.data) ## Computes the empirical DF of the data
	test.Fn.1 <- test.Fn(test.data.1)
	test.Fb <- test.data.1
	test.Freq <- diff(c(0,test.Fn.1))

	F.is.interpolate <- approx(x = tr.data$F.b, y = F.is, xout = test.Fb, yleft = 0, yright = 1)
	loss <- sum( {alp.hat* F.is.interpolate$y + (1- alp.hat) *test.Fb- test.Fn.1}^2* test.Freq)
	return(loss)
}