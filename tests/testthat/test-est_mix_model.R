test_that("Computes lower confidence bound for a simple 2 component mixture model", {
  data.gen<- function( n, alpha){
   ind<- rbinom(n, 1, alpha)
   data<- ind* rbeta(n, 1,5) + (1-ind)* runif(n,0,1)
   return(data)
  }
  set.seed(1)
 data.1<- data.gen(n = 500, alpha = .1)
 est.lwr.bnd <- mix.model(data.1, method = "lwr.bnd", gridsize = 50)
 expect_identical(est.lwr.bnd$alp.Lwr,0.06)
})
