// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in
#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]
// via the exports attribute we tell Rcpp to make this function
// available from R
//
//
// [[Rcpp::export]]
Rcpp::List rnewt(const Eigen::VectorXd &x0, 
                   Rcpp::Function fn,
                   Rcpp::Function gr,
                   Rcpp::Function he, 
                   const Eigen::MatrixXd &gr0,
                   const Eigen::MatrixXd &d0,
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
                   const double &tolstep,
                   const double &tolmu,
                   const double &tolmu2,
                   const double &tolc,
                   const int &maxreject,
                   const bool &grdre,
                   const bool &verbose,
                   const int &riter,
                   const bool &returnhess
) {
  // define objects
  MatrixXd x(x0.rows(), 1);x.setZero(); //pars
  x.col(0) = x0;
  MatrixXd d(x0.rows(), 1);d.setZero(); //step direction
  d = d0;
  MatrixXd p(1,1);
  // MatrixXd p(x0.rows(), 1);p.setZero();
  MatrixXd grd(x0.rows(), 1);grd.setZero();
  MatrixXd grdold(x0.rows(), 1);
  grdold = gr0;
  MatrixXd hess = MatrixXd::Identity(x.rows(),x.rows());
  // if(method == 4 | method == 5 | method == 6){
  //   // BFGS
  //   hess = MatrixXd::Identity(x.rows(),x.rows());
  // } 
  if(!quasi){
    hess = Rcpp::as<MatrixXd>(he(x0));
  }
  int griter = 0; //nr grd. evals
  bool pass = true;
  bool firstpass = false;
  int info = 0;
  int reject = 0; //rejection counter
  double ared = 0; //achieved improvement
  double pred = 0; //predicted improvement
  double arediter = 0;//cumulative improvement for printing messages
  double objnew = 0; //best so far
  double obj = Rcpp::as<double>(fn(x0)); //obj of step
  double gamma = 1;double newgamma = 0;
  int s = 0;
  double rho = 1;
  if((method == 1) || (method == 3))s = 2*m;
  if(method == 2)s = m;
  int si = 0;
  MatrixXd Y(x0.rows(), s);Y.setZero(); //diff in grd. of updates
  MatrixXd S(x0.rows(), s);S.setZero(); //diff in pars. of updates
  MatrixXd Q(1, 1);MatrixXd SY(1, 1);MatrixXd L(1, 1);
  MatrixXd A(x.rows(), 1);
  Eigen::LLT<Eigen::MatrixXd> lltH;
  grd.col(0) = gr0;
  MatrixXd y(x0.rows(), 1);
  
  double sy;

  int i = 0;
  while(i<maxit){
    i++;
    
    if(i % 5 == 0){ //check every 5 iterations for interruption
      Rcpp::checkUserInterrupt();
    }
    
    //if(pass){
    if(quasi && grdre || pass){
      //incorporate rejection information for limited memory
      grdold = grd;
      grd = Rcpp::as<VectorXd>(gr(x));
      griter++;
      if(grd.lpNorm<Eigen::Infinity>()<=tolg){ // infinity norm, above eq 37
        info = 1;
        break;
      }
    }
      if(quasi){
        //tol2: cautious updating(see eq23)
        y = grd-grdold;
        sy = (y.transpose()*d).value();
        
        // newgamma = (y.transpose()*y).value()/sy;
        // gamma = std::min(1/tolgamma, newgamma);
        // newgamma = sy/y.squaredNorm();//doesnt quite work, need to check how this echos into Q
        // gamma = std::min(1/tolgamma, newgamma/(1+newgamma));
        
        if(((method == 1)| (method == 2)| (method ==3))){
          
          if(!firstpass || ((method == 2) | (method == 3)) | ((method == 1) && (sy>=(tolc*d.squaredNorm())))){
          //limited memory matrices
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
        }
          if(sy>=(tolc*d.squaredNorm())){
          newgamma = (y.transpose()*y).value()/sy;
          gamma = std::min(1/tolgamma, newgamma);
        }
        }
      }
      
      if(quasi && pass){
        //not convinced i should be cautiously updating gamma.

        if(!std::isnan(rho)){
        
        if((method == 1)| (method == 2)| (method ==3)){
          //limited memory updates
        
          if(si<=s && SY.cols()<s){
            SY.resize(si,si);
          }
          if(si<=s && L.cols()<s){
            L.resize(si,si);
          }
          SY = S.leftCols(si).transpose()*Y.leftCols(si);
          L = SY.triangularView<Eigen::StrictlyLower>();
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
            Q.topLeftCorner(si,si) = -pow(gamma,-1)*S.leftCols(si).transpose()*S.leftCols(si);
            Q.bottomLeftCorner(si,si) = -pow(gamma, -1)*L.transpose();
            Q.topRightCorner(si,si) = -pow(gamma, -1)*L;
            Q.bottomRightCorner(si,si) = D;
          } else if(method == 2){
            if(si<=s && A.cols()<s){
              A.resize(Eigen::NoChange, si);
            }
            if(si<=s && A.cols()<s){
              Q.resize(si, si);
            }
            A = Y.leftCols(si)-gamma*S.leftCols(si);
            Q = D.toDenseMatrix()+L+L.transpose()-gamma*S.leftCols(si).transpose()*S.leftCols(si);
            //This, SR1, still needs a rejection rule for ill-conditioned updates..
          }else if(method == 3){
            if(si<=s && A.cols()<(2*s)){
              A.resize(Eigen::NoChange, 2*si);
            }
            if(si<=s && Q.cols()<(2*s)){
              Q.resize(2*si, 2*si);
            }
            MatrixXd StS = S.leftCols(si).transpose()*S.leftCols(si);
            MatrixXd U = StS.triangularView<Eigen::Upper>(); // nonstrict upper triangular
            // diagonal and entries above to zero
            // for(int l=0; l<(si-1); l++){
            //   for(int k=l+1; k<(si-1); k++){
            //     U(k,l) = 0; //entries under diagonal
            //   }
            // }
            A.leftCols(si) = S.leftCols(si);
            A.rightCols(si) = Y.leftCols(si);
            Q.topLeftCorner(si,si) = MatrixXd::Zero(U.rows(), U.cols());
            Q.bottomLeftCorner(si,si) = U.transpose();
            Q.topRightCorner(si,si) = U;
            Q.bottomRightCorner(si,si) = D.toDenseMatrix()+gamma*StS.diagonal().asDiagonal().toDenseMatrix()+L+L.transpose();
           }
        }else if((method == 4) | (method == 5) | (method == 6)){
        //full hessian approximation
         if(i==1)hess.diagonal().array() *= newgamma;
         if(method == 4){
          //dhessd = -d*grd
          //due to secant condition dhess = -grad
          hess += ((y*y.transpose()).array()/sy - ((hess*d)*(d.transpose()*hess.transpose())).array()/(d.transpose()*hess*d).sum()).matrix(); //update hessian
        }else if(method == 5){
          MatrixXd yBs = (y-hess*d);
          if((yBs.transpose()*d).sum()>0||i==1)hess += ((yBs*yBs.transpose()).array()/(yBs.transpose()*d).sum()).matrix();
        }else if(method == 6){
          MatrixXd yBs = (y-hess*d);
          hess += ((yBs*d.transpose()+d*yBs.transpose()).array()/(d.transpose()*d).sum()- (yBs.transpose()*d).sum()/pow((d.transpose()*d).sum(),2)*(d*d.transpose()).array()).matrix();
        }
        
    }
    }
    }
    
    if(quasi && method <4){
      // p = (Q+pow(gamma+mu0,-1)*A.transpose()*A).lu().solve(A.transpose()*grd);
      d = -pow(gamma+mu0, -1)*grd + pow(gamma+mu0, -2)*A*(Q+pow(gamma+mu0,-1)*A.transpose()*A).lu().solve(A.transpose()*grd);  
    }else if(quasi){
      // lltH.compute(hess);
      d = -(hess + mu0*MatrixXd::Identity(hess.rows(),hess.cols())).lu().solve(grd);
    }else{
      if(pass)hess = Rcpp::as<MatrixXd>(he(x)); //only recalculate hessian if parameters were updated
      lltH.compute(hess + mu0*MatrixXd::Identity(hess.rows(),hess.cols()));
      d = -lltH.solve(grd);
    }
    pred = 0.5*mu0*(d.transpose()*d).value()-0.5*(d.transpose()*grd).value();
    if(pred<=(pmin*d.norm()*grd.norm())){ //sum(abs(grd))*sum(abs(d)))){#absolute norm
      mu0 = std::max(sigma2*mu0, tolmu);
      pass = false;
      reject ++;
      }else{
        objnew = Rcpp::as<double>(fn(x+d));
        ared = obj-objnew; // actual improvement
        
        // if((ared)/std::max(1., abs(obj))<=tolobj && ared>0){
        //   info = 3;
        //   break;
        // }
        if(abs(ared)<=tolobj){
          info = 3;
          break;
        }

        rho = ared/pred; // trust region ratio
        if((std::isnan(rho) || rho<=c1)){ // reject, no improvement
          reject++;
          if(reject>maxreject){
            info = 2;
            break;
          }
          //might want to consider a linesearch here if rho isnan to prevent getting stuck indefinitely
          //although that will also cause a break with info 4 in a few iterations
          pass = false;
          mu0 = std::max(sigma2*mu0, tolmu);
        }else if((c1<rho) & (rho <= c2)){ // accept, somewhat succesful updates, might also want to go here if sufficien reduction in gradient (see TMB's newton)
          //might want to check if we can take a longer step here to prevent excess iterations.
          pass = true;reject=0;
          firstpass = true;
          x = x+d;
          obj = objnew;
        }else if(c2<rho){ // accept, highly succesful, decrease regularisation
          reject = 0;
          pass = true;
          firstpass = true;
          x = x+d;mu0 = std::max(sigma1*mu0, tolmu);
          obj = objnew;
        }
        
        if(d.squaredNorm()<=tolstep){
          info = 6;
          break;
        }
        
        if(mu0>=tolmu2){
          info = 4;
          break;
          // could try something radical and set mu=m0, if we dont get out on the second try we can still break
          
        }

        if(verbose && pass){
          arediter += ared;
        }
        if(pass && verbose && (i % riter == 0)){
          Rcpp::Rcout << "Iter: " << i << " Objective: " << std::fixed << std::setprecision(2) << obj << std::setprecision(4) << std::scientific << " Achieved reduction: " << arediter << " mu: " << mu0  << std::scientific << std::endl;
          arediter = 0;
        }
      }
      
      if(!pass && verbose && (i % riter == 0)){
        Rcpp::Rcout << "Iter: " << i << " No improvement" << std::endl;
      }
      }
  double grmax = grd.cwiseAbs().maxCoeff();
  if(i==maxit)info = 5;
  
  if(!returnhess){
    return Rcpp::List::create(Rcpp::Named("par")=Rcpp::wrap(x),
                              Rcpp::Named("objective")=Rcpp::wrap(obj),
                              Rcpp::Named("iterations")=Rcpp::wrap(i),
                              Rcpp::Named("evalg")=Rcpp::wrap(griter),
                              Rcpp::Named("info")=Rcpp::wrap(info),
                              Rcpp::Named("maxgr")=Rcpp::wrap(grmax)
    );
  }else{
    if(quasi && ((method == 1) | (method == 2) | (method == 3))){
    if(verbose)Rcpp::Rcout << "Calculating hessian..." << std::endl;
    hess = gamma*MatrixXd::Identity(x.rows(),x.rows()) + A*Q.lu().solve(A.transpose());
    }
    return Rcpp::List::create(Rcpp::Named("par")=Rcpp::wrap(x),
                              Rcpp::Named("objective")=Rcpp::wrap(obj),
                              Rcpp::Named("iterations")=Rcpp::wrap(i),
                              Rcpp::Named("evalg")=Rcpp::wrap(griter),
                              Rcpp::Named("info")=Rcpp::wrap(info),
                              Rcpp::Named("maxgr")=Rcpp::wrap(grmax),
                              Rcpp::Named("hess")=Rcpp::wrap(hess)
    );
  }
  
  
}