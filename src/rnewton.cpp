// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// simple example of creating two matrices and
// returning the result of an operation on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
//
// [[Rcpp::export]]
Rcpp::List rnewton(const Eigen::VectorXd &x0, 
                   Rcpp::Function fn,
                   Rcpp::Function gr,
                   Rcpp::Function he, 
                   const Eigen::MatrixXd &gr0,
                   const Eigen::MatrixXd &d0,
                   const bool &regularize,
                   const bool &quasi,
                   const int &method,
                   const int &maxit,
                   const int &m,
                   double &mu0,
                   const double &sigma1,
                   const double &sigma2,
                   const double &c1,
                   const double &c2,
                   const double &pmin,
                   const double &tolg,
                   const double &tolgamma,
                   const double &tolobj,
                   const double &tolmu,
                   const double &tolmu2,
                   const double &tolc,
                   const bool &verbose,
                   const int &riter
) {
  // define objects
  MatrixXd x(x0.rows(), 1);x.setZero(); //pars
  x.col(0) = x0;
  MatrixXd d(x0.rows(), 1);d.setZero(); //step direction
  // MatrixXd p(x0.rows(), 1);p.setZero();
  MatrixXd grd(x0.rows(), 1);grd.setZero();
  MatrixXd grdold(x0.rows(), 1);grdold.setZero();
  MatrixXd hess;
  if(!quasi){
    hess = Rcpp::as<MatrixXd>(he(x0));
  }
  int griter = 0; //nr grd. evals
  bool pass = true;
  bool firstpass = false;
  int info = 0;
  double ared = 0; //achieved improvement
  double pred = 0; //predicted improvement
  double objnew = 0; //best so far
  double obj = Rcpp::as<double>(fn(x0)); //obj of step
  double gamma = 1;double newgamma = 0;
  int s = 0;
  double rho = 1;
  if(method == 1)s = 2*m;
  if(method == 2)s = m;
  if(method == 3)s = 2*m;
  int si = 1;
  MatrixXd y(x0.rows(), 1);y.setZero();
  double sy = 0;
  MatrixXd Y(x0.rows(), s);Y.setZero(); //diff in grd. of updates
  MatrixXd S(x0.rows(), s);S.setZero(); //diff in pars. of updates
  MatrixXd Q(1, 1);MatrixXd SY(1, 1);MatrixXd StS(1, 1);MatrixXd L(1, 1);
  MatrixXd A(x.rows(), 1);
  
  grd.col(0) = Rcpp::as<VectorXd>(gr(x0));
  Y.col(0) = grd - gr0;
  S.col(0) = d0;
  int i = 1;
  while(i<maxit){
    if(i % 10 == 0){
      Rcpp::checkUserInterrupt();
    }
    if(pass){
      grdold = grd;
      grd.col(0) = Rcpp::as<VectorXd>(gr(x));
      griter++;
      if(grd.col(0).lpNorm<Eigen::Infinity>()<=tolg){ // infinity norm, above eq 37
        info = 1;
        break;
      }
      
      if(quasi){
        //tol2: cautious updating(see eq23)
        y = grd-grdold;
        sy = (y.transpose()*d).value();
        if(sy>=(tolc*d.col(0).squaredNorm()) & firstpass){
          if(si == s){
          // max memory: matrices are full
          Y.leftCols(s-1) = Y.rightCols(s-1); //first s columns
          Y.col(s-1) = y; //last column
          S.leftCols(s-1) = S.rightCols(s-1);
          S.col(s-1) = d;
          }else{
          Y.col(si) = y;
          S.col(si) = d;
          si++;
          }
          newgamma = (y.transpose()*y).value()/sy;
          gamma = std::min(1/tolgamma, newgamma);
        }

        if(!std::isnan(rho)){
          if(si<=s && SY.cols()<s){
            SY.resize(si,si);
          }
          if(si<=s && StS.cols()<s){
            StS.resize(si,si);
          }
          if(si<=s && L.cols()<s){
            L.resize(si,si);
          }
          SY = S.leftCols(si).transpose()*Y.leftCols(si);
          StS = S.leftCols(si).transpose()*S.leftCols(si);
          L = SY;
          //diagonal and entries above to zero
          for(int l=0; l<(si-1); l++){
            L(l,l) = 0; // diagonal
            for(int k=l+1; k<(si-1); k++){
              L(l,k) = 0; //entries above diagonal
            }
          }
          Eigen::DiagonalMatrix<double,  Eigen::Dynamic> D = SY.diagonal().asDiagonal();
          
          if(method == 1){
            if(si<=s && A.cols()<(2*s)){
              A.resize(Eigen::NoChange, 2*si);
            }
            A.leftCols(si) = S.leftCols(si);
            A.rightCols(si) = Y.leftCols(si);
            if(si<=s  && Q.cols()<(2*s)){
              Q.resize(2*si, 2*si);
            }
            Q.block(0,0,si,si) = -pow(gamma,-1)*StS;
            Q.block(si,0,si,si) = -pow(gamma, -1)*L.transpose();
            Q.block(0,si,si,si) = -pow(gamma, -1)*L;
            Q.block(si,si,si,si) = D;
          } else if(method == 2){ // probably just runs into non pd updated
            if(si<=s && A.cols()<s){
              A.resize(Eigen::NoChange, si);
            }
            if(si<=s && A.cols()<s){
              Q.resize(si, si);
            }
            A = Y.leftCols(si)-gamma*S.leftCols(si);
            Q = D.toDenseMatrix()+L+L.transpose()-gamma*StS;
          }else if(method == 3){
            if(si<=s && A.cols()<(2*s)){
              A.resize(Eigen::NoChange, 2*si);
            }
            if(si<=s && Q.cols()<(2*s)){
              Q.resize(2*si, 2*si);
            }
            MatrixXd U = StS; // nonstrict upper triangular
            //diagonal and entries above to zero
            for(int l=0; l<(si-1); l++){
              for(int k=l+1; k<(si-1); k++){
                U(k,l) = 0; //entries under diagonal
              }
            }
            A.leftCols(si) = S.leftCols(si);
            A.rightCols(si) = Y.leftCols(si);
            Q.block(0,0,si,si) = MatrixXd::Zero(U.rows(), U.cols());
            Q.block(si,0,si,si) = U.transpose();
            Q.block(0,si,si,si) = U;
            Q.block(si,si,si,si) = D.toDenseMatrix()+gamma*StS.diagonal().asDiagonal().toDenseMatrix()+L+L.transpose();
            }
      }
    }
    }
    
    if(quasi){
      // p = (Q+pow(gamma+mu,-1)*A.transpose()*A).lu().solve(A.transpose()*grd);
      // d = -pow(gamma+mu, -1)*grd + pow(gamma+mu, -2)*A*p; // this inversion can still be improved I think following the paper
      d = -pow(gamma+mu0, -1)*grd + pow(gamma+mu0, -2)*A*(Q+pow(gamma+mu0,-1)*A.transpose()*A).lu().solve(A.transpose()*grd);
    }else{ // do newton update is quasi is starting
    if(pass)hess = Rcpp::as<MatrixXd>(he(x)); //only recalculate hessian if parameters were updated
    d = -(hess + mu0*MatrixXd::Identity(hess.rows(),hess.cols())).llt().solve(grd);
    //d = -(hess + mu0*MatrixXd::Identity(hess.rows(),hess.cols())).llt().solve(grd);
    }
    pred = 0.5*mu0*(d.transpose()*d).value()-0.5*(d.transpose()*grd).value();
    if(pred<=(pmin*d.norm()*grd.norm())){ //sum(abs(grd))*sum(abs(d)))){#absolute norm
      mu0 = std::max(sigma2*mu0, tolmu);
      }else{
        objnew = Rcpp::as<double>(fn(x+d));
        ared = obj-objnew; // actual improvement
        
        if(ared<=tolobj && ared>0){
          info = 3;
          break;
        }
        rho = ared/pred;
        if(std::isnan(rho) || rho<=c1){ // negative difference always goes here, unsucessful
          pass = false;
          mu0 = std::max(sigma2*mu0, tolmu);
        }else  if(c1<rho & rho <= c2){ // somewhat succesful
          pass = true;
          firstpass = true;
          x = x+d;
          obj = objnew;
        }else if(c2<rho){ // highly succesful
          pass = true;
          firstpass = true;
          x = x+d;mu0 = std::max(sigma1*mu0, tolmu);
          obj = objnew;
        }
        if(mu0>=tolmu2){
          info = 4;
          break;
        }

        if(pass && verbose && (i % riter == 0)){
          Rcpp::Rcout << "Iter: " << i << " Objective: " << std::fixed << std::setprecision(2) << obj << std::setprecision(4) << std::scientific << " Achieved reduction: " << ared << " mu: " << mu0  << std::scientific << std::endl;
        }else if(verbose && (i % riter == 0)){
          Rcpp::Rcout << "Iter: " << i << std::endl;
        }
      }
      i++;
      }
  double grmax = grd.cwiseAbs().maxCoeff();
  if(i==maxit)info = 5;
  return Rcpp::List::create(Rcpp::Named("objective")=Rcpp::wrap(obj),
                            Rcpp::Named("iterations")=Rcpp::wrap(i),
                            Rcpp::Named("evalg")=Rcpp::wrap(griter),
                            Rcpp::Named("par")=Rcpp::wrap(x),
                            Rcpp::Named("info")=Rcpp::wrap(info),
                            Rcpp::Named("maxgr")=Rcpp::wrap(grmax)
                              );
  
}