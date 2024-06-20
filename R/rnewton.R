#'@aliases rnewt rnewton
#' @export
#' @importFrom Rcpp evalCpp
#' 
rnewton <-function(x0, fn, gr, he = NULL, 
                  quasi = TRUE, method = "LBFGS", verbose = FALSE, return.hess = FALSE, control = list(maxit = 1e3, m = 5, sigma1 = 0.5, sigma2 = 4, c1 = 1e-3, c2 = 0.9, pmin = 1e-3, tol.g = 1e-8, tol.gamma = 1e-5, tol.obj = 1e-8, tol.mu = 1e-4, tol.mu2 = 1e15, tol.c = 1e-8, report.iter = 10, grad.reject = FALSE, max.reject = 50, mu0 = 5), ...){
  
  if(!quasi && is.null(he))stop("Function for hessian must be provided with 'quasi = FALSE'")
  
  fill_control = function(x){
  if (!("maxit" %in% names(x))) 
    x$maxit = 1e3
  if (!("m" %in% names(x))) 
    x$m = 5
  if (!("sigma1" %in% names(x))) 
    x$sigma1 = 0.5
  if (!("sigma2"%in% names(x)))
    x$sigma2 = 4
  if (!("c1"%in% names(x)))
    x$c1 = 1e-3
  if (!("pmin"%in% names(x)))
    x$pmin = 1e-3
  if (!("c2"%in% names(x)))
    x$c2 = 0.9
  if (!("tol.g"%in% names(x)))
    x$tol.g = 1e-8
  if (!("tol.gamma"%in% names(x)))
    x$tol.gamma = 1e-5
  if (!("tol.obj"%in% names(x)))
    x$tol.obj = 1e-8
  if (!("tol.mu"%in% names(x)))
    x$tol.mu = 1e-4
  if (!("tol.mu2"%in% names(x)))
    x$tol.mu2 = 1e15
  if (!("tol.c"%in% names(x)))
    x$tol.c = 1e-8
  if (!("report.iter"%in% names(x)))
    x$report.iter = 10
  if (!("max.reject"%in% names(x)))
    x$max.reject = 50
  if (!("grad.reject"%in% names(x)))
    x$grad.reject = FALSE
  if (!("mu0"%in% names(x)))
    x$mu0 = 5
  x
}

control <- fill_control(control)

maxit = control$maxit;m=control$m;sigma1=control$sigma1;sigma2=control$sigma2;c1=control$c1;pmin=control$pmin;c2=control$c2;tol.g=control$tol.g;tol.gamma=control$tol.gamma;tol.obj=control$tol.obj;tol.mu=control$tol.mu;tol.mu2=control$tol.mu2;tol.c=control$tol.c;report.iter=control$report.iter;max.reject=control$max.reject;grad.reject=control$grad.reject;mu0=control$mu0

# Defensive programming for tolerances
if(pmin>1|pmin<0)stop("pmin must be between 0 and 1.")
if(c1>1|c1<0)stop("c1 must be between 0 and 1.")
if(c2<c1|c2>1)stop("c2 must be between c1 and 1.")
if(sigma1<0|sigma1>1)stop("sigma must be between 0 and 1.")
if(sigma2<1 && sigma2!=0)stop("sigma2 must be larger than 1.")

method = match.arg(method, c("LBFGS", "LSR1", "LPSB", "BFGS", "SR1", "PSB"))
method = switch(method, "LBFGS" = 1, "LSR1" = 2, "LPSB" = 3, "BFGS" = 4, "SR1" = 5, "PSB" = 6)

grdold = matrix(gr(x0), nrow = 1)
d = -t(grdold / norm(grdold,"f"))

lerr <- try({
res <- cvsrch(fn, gr, length(x0), x0, fn(x0), grdold, d, 1, 1e-4, 0.9, 1e-16, 1e-12, 1e20, 20);

d = res$x-x0;
pars = res$x;
}, silent = TRUE)
# linesearch failed for some reason, try to use provided start

if(inherits(lerr,"try-error")){
  pars = x0
  grdold = matrix(gr(x0-d), ncol = 1)
}

if(quasi)he <- function()matrix(0)
  
system.time(opt <- rnewt(x0=pars, fn=fn, gr=gr, he=he, gr0=t(grdold), d0=d, quasi=quasi, method = method, maxit = maxit, m = m, mu0 = mu0, sigma1 = sigma1, sigma2 = sigma2, c1 = c1, c2 = c2, pmin = pmin, tolg = tol.g, tolgamma = tol.gamma, tolobj = tol.obj, tolmu = tol.mu, tolmu2 = tol.mu2, tolc = tol.c , verbose = verbose, riter = report.iter, maxreject = max.reject, grdre = grad.reject, returnhess = return.hess))

opt$par <- c(opt$par)
if(!is.null(names(x0))) names(opt$par) <- names(x0)

opt$convergence = FALSE
# report succesful convergence if: reached gradient tolerance (1) or reached objective tolernace (3)
if(opt$info %in% c(0, 1, 3))opt$convergence = TRUE 

class(opt) <- "rnewton"
return(opt)
}