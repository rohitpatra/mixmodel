test_that("Computes density estimate in a simple 2 component example", {
  data.gen<- function( n, alpha){
   ind<- rbinom(n, 1, alpha)
   data<- ind* rbeta(n, 1,5) + (1-ind)* runif(n,0,1)
   return(data)
  }
  set.seed(1)
 data.1<- data.gen(n = 500, alpha = .1)
 est.fixed <- mix.model(data.1,method = "fixed", c.n = .05*log(log(length(data.1))), gridsize = 50)
 out <- den.mix.model(est.fixed)
 expect_identical(floor(sum(out$y)),659)
})
