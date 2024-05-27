#' @title Regularized quasi-Newton optimization
#' @description Performs regularized (quasi-)Newton optimisation with BFGS, SR1, or PSB updates.
#'
#' @param x0 Initial values for the parameters.
#' @param fn A function to be minimized.
#' @param gr A function that returns the gradient.
#' @param he A function that returns the hessian (only used when quasi = FALSE).
#' @param quasi logical. Defaults to TRUE. If FALSE implements regularised Newton optimization.
#' @param method The method to be used when \code{quasi = TRUE}. Defaults to "LBFGS", alternatives are "LPSB", "LSR1" and ther full-memory alternatives "BFGS", "SR1", "PSB". The latter three options should probably not be used in practice (see details).
#' @param verbose logical. Defaults to FALSE. If TRUE prints reports on each iteration.
#' @param return.hess logical. Defaults to FALSE. If TRUE prints returns the (approximated for quasi-Newton methods) Hessian.
#' @param control a \code{"\link{list}"} of control constants.
#' \itemize{
#'  \item{\emph{maxit}: }{ The maximum number of iterations. Defaults to 1000.}
#'  \item{\emph{m}: }{ The number of gradients to remember from previous optimisation steps. Defaults to 5.}  
#'  \item{\emph{sigma1}: }{ Step decrement factor. Defaults to 0.5. Must be smaller than 1 but larger than 0.}
#'  \item{\emph{sigma2}: }{ Step increment factor. Defaults to 4. Must be larger than 1.}
#'  \item{\emph{c1}: }{ First constant for determining step success. Defaults to 1e-3. See details.}
#'  \item{\emph{c2}: }{ Second constant for determining step success. Defaults to 0.9. See details.}
#'  \item{\emph{pmin}: }{ Third constant for determining (lack of) step success. Defaults to 1e-3.}
#'  \item{\emph{tol.g}: }{ Convergence tolerance for gradient. Defaults to 1e-8.}
#'  \item{\emph{tol.gamma}: }{ Threshold for gamma paramter. Defaults to 1e-5.}
#'  \item{\emph{tol.obj}: }{ Convergence tolerance for reduction in the objective. Defaults to 1e-8.}
#'  \item{\emph{tol.mu}: }{ Minimum threshold for the regularisation parameter. Defaults to 1e-4.}
#'  \item{\emph{tol.mu2}: }{ Maximum threshold for the regularisation parameter. Defaults to 1e15.}
#'  \item{\emph{tol.c}: }{ Tolerance for cautious updating. Defaults to 1e-8.}
#'  \item{\emph{report.iter}: }{ If 'verbose = TRUE', how often should a report be printed? Defaults to every 10 iterations.}
#'  \item{\emph{max.reject}: }{ Maximum number of consecutive rejections before algorithm terminates.}
#'  \item{\emph{grad.reject}: }{ Logical. If TRUE the gradient is evaluated at every iteration and information of rejected steps is incorporated in limited-memory methods. Defaults to TRUE.}
#' }
#' @param ... Not used.
#'
#' @details
#' This function implements some of the regularised (quasi-)Newton optimisation methods presented in Kanzow and Steck (2023). The full-memory options that are implemented rely on explicitly inverting the approximated Hessian and regularisation penalty, are thus slow, and should probably not be used in practice.
#' 
#' The function start with a single More-Thuente line search along the normalized negative gradient direction. 
#' The code for this was originally written in matlab by Burdakov et al. (2017), translated to python by Kanzow and Steck (2023), and separately translated to R code for this package. 
#' 
#' A step is considered somewhat successful for \eqn{c_1<\rho\leq c_2}, where \eqn{\rho} is the proportion of achieved and predicted reduction in the objective function. Note the requirement \eqn{c_1 \in (0,1)} and \eqn{c_2 \in (c_1,1)}.
#' A step is considered highly successful for \eqn{c_2 < \rho}, where rho is the proportion of achieved and predicted reduction in the objective function.
#' 
#' The \eqn{\sigma_1} constant controls the decrement of the regularisation parameter \eqn{\mu} on a highly succesful step.
#' The \eqn{\sigma_2} constant controls the increment of the regularisation parameter \eqn{\mu} on a unsuccesful step. A step is defned as unsuccesful if 1) the predicted reduction less than pmin times the product of the l2 norm for step direction and gradient, or 2) if \eqn{\rho \leq c_1}.
#'
#' @return An object of class "rnewton" including the following components:
#' \itemize{
#'  \item{objective: }{ The value of fn corresponding to par.}
#'  \item{iterations: }{ Number of completed iterations.}
#'  \item{evalg: }{ Number of calls to gr.}
#'  \item{par: }{ The best set of parameters found.}
#'  \item{info: }{ Convergence code. See \code{"\link"{print.rnewton}} for details.}
#'  \item{maxgr: }{ Maximum absolute gradient component.}
#'  \item{convergence: }{ logical, TRUE indicating succesful convergence (reached tol.obj or tol.g).}
#' }
#' 
#' @author Bert van der Veen
#' @references
#' Burdakov, O., Gong, L., Zikrin, S., & Yuan, Y. X. (2017). On efficiently combining limited-memory and trust-region techniques. Mathematical Programming Computation, 9, 101-134.
#' 
#' Kanzow, C., & Steck, D. (2023). Regularization of limited memory quasi-Newton methods for large-scale nonconvex minimization. Mathematical Programming Computation, 15(3), 417-444.
#'
#'@aliases rnewt rnewton
#' @useDynLib rnewton
#' @export
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' 
rnewton <-function(x0, fn, gr, he = NULL, 
                  quasi = TRUE, method = "LBFGS", verbose = FALSE, return.hess = FALSE, control = list(maxit = 1e3, m = 5, sigma1 = 0.5, sigma2 = 4, c1 = 1e-3, c2 = 0.9, pmin = 1e-3, tol.g = 1e-8, tol.gamma = 1e-5, tol.obj = 1e-8, tol.mu = 1e-4, tol.mu2 = 1e15, tol.c = 1e-8, report.iter = 10, grad.reject = TRUE, max.reject = 50, mu0 = 5), ...){
  
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
    x$grad.reject = TRUE
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