print.rnewton <- function(x, ...){
  x$info <- switch(as.character(x$info), "1" = "1: gradient tolerance reached." , "2"= "2: maximum number of unsuccesful iterations reached.", "3"="3: insufficient improvement in objective function (ran into tol.obj).", "4"="4: reached regularisation limit (ran into tol.mu2).", "5"="5: reached maximum interations.", "6"="6: step tolerance reached.")
  cat(" Code", x$info, "\n")
  
  sts <- data.frame(c("Objective:", "Iterations:", "Evaluations:","max(|grad|):"), c(format(x$objective, nsmall=2,digits=2),x$iterations,x$evalg, signif(x$maxgr, 2)))
  colnames(sts) <- NULL
  print(sts, row.names=FALSE, right = FALSE, digits = 2)
  
  if(!x$convergence)warning("Optimisation unsuccesful.")
}

