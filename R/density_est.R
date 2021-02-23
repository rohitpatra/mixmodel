#' A function to compute the estimate of f.s (the density corresponding to F.s) when f.s is known to be either decreasing or increasing.
#'
#' @param input An element of class `cv.mix.model` or `mix.model`
#' @param dec.density A logical variable. If FALSE then density assumed to be increasing. Default is TRUE.
#'
#' @return A piecewise function as a list.
#' \item{x}{data points where the piecewise density has changes.}
#' \item{y}{value of density at the above points}
#' @importFrom fdrtool gcmlcm
#'
#' @export
#' @examples
#' data.gen<- function( n, alpha){
#' 	ind<- rbinom(n, 1, alpha)
#' 	data<- ind* rbeta(n, 1,5) + (1-ind)* runif(n,0,1)
#' 	return(data)
#'  }
#'  data.1<- data.gen(n = 1000, alp = .2)
#' est.fixed <- mix.model(data.1, method = "fixed"
#' , c.n = .05*log(log(length(data.1))), gridsize = 600)
#' out <- den.mix.model(est.fixed)
#' plot(est.fixed$Fs.hat, type= "l", xlab= "x", ylab= "F_s")
#' plot(out, type="l", xlab="x", ylab= "f_s")

den.mix.model<- function(input, dec.density =TRUE){
	if(class(input) == "cv.mix.model" || class(input) == "mix.model"){
		Fs.hat <-  input$Fs.hat
	} else {
			stop("This function only works on objects of class 'cv.mix.model' or 'mix.model'. See functions 'mix.model' or 'cv.mix.model'.")
	}
	if(dec.density== TRUE){
		ll <- fdrtool::gcmlcm(Fs.hat$x,Fs.hat$y, type="lcm")
	} else if(dec.density== FALSE){
		ll <- fdrtool::gcmlcm(Fs.hat$x,Fs.hat$y, type="gcm")
	}
	fs.hat <- NULL
	fs.hat$x=rep(ll$x.knots,each=2) #data points for density
	fs.hat$y=c(0,rep(ll$slope.knots,each=2),0) #value of density
	return(fs.hat)
}



