//Testing to see that the covariance from the sum of two vectors with different covariances
//is the sum of the covariances
//In general it seems to be: 
//cov(x_i-y_i,x_j-y_j) = cov(x_i,x_j)+cov(y_i,y_j)-cov(x_i,y_j)-cov(x_j,y_i)

void testCovSum()
{
 double x1,x2,y1,y2,sigX1,sigX2,sigY1,sigY2;
 x1 = 1.0; sigX1 = 0.2;
 y1 = -2.0; sigY1 = 0.3;
 x2 = 3.0; sigX2 = 0.1;
 y2 = 2.0; sigY2 = 0.2;

 const int nEle = 2; //number of elements in vector
 const int nVec = 10000;
 double mean1[nEle] = {x1,y1};
 double mean2[nEle] = {x2,y2};
 double sigma1[nEle] = {sigX1,sigY1};
 double sigma2[nEle] = {sigX2,sigY2};
 TVectorD v1(nEle);
 TVectorD v2(nEle);

 TMatrixD corr1(nEle,nEle);
 TMatrixD corr2(nEle,nEle);

 corr1(0,0) = 1.0;
 corr1(0,1) = 0.1;
 corr1(1,0) = 0.1;
 corr1(1,1) = 1.0;

 corr2(0,0) = 1.0;
 corr2(0,1) = 0.3;
 corr2(1,0) = 0.3;
 corr2(1,1) = 1.0;

 TDecompChol chol1(corr1);
 if(!chol1.Decompose()) cout << "Can't Decompose corr1!" << endl;
 TDecompChol chol2(corr2);
 if(!chol2.Decompose()) cout << "Can't Decompose corr2!" << endl;

 TMatrixD cholUp1 = chol1.GetU();
 cholUp1.T(); 
 TMatrixD cholUp2 = chol2.GetU();
 cholUp2.T(); 
 
 TRandom3 gen;
 vector<TVectorD> allCorr1;
 vector<TVectorD> allCorr2;
 double vec[nVec][nEle],vec1[nVec][nEle],vec2[nVec][nEle];
 for(int i=0;i<nVec;i++){
  TVectorD vUn1(nEle);
  TVectorD vUn2(nEle);
  
  //Make uncorrelated vectors
  for(int j=0;j<nEle;j++){
   vUn1(j) = gen.Gaus(0,1);
   vUn2(j) = gen.Gaus(0,1);
  }

  //Make correlated vectors
  TVectorD vCorrUnit1 = cholUp1*vUn1;
  TVectorD vCorrUnit2 = cholUp2*vUn2;
  
  //Impose selected means and uncertainties
  TVectorD vCorr1(nEle);
  TVectorD vCorr2(nEle);
  for(int j=0;j<nEle;j++){
   vCorr1(j) = mean1[j] + vCorrUnit1(j)*sigma1[j];
   vCorr2(j) = mean2[j] + vCorrUnit2(j)*sigma2[j];
  }
  
  //Place individual vectors into <vectors>
  allCorr1.push_back(vCorr1);
  allCorr2.push_back(vCorr2);
 }//end loop over samples

 //Make covariance matrices
 TMatrixD cov1(nEle,nEle);
 TMatrixD cov2(nEle,nEle);
 cov1(0,0) = corr1(0,0)*(sigX1*sigX1);
 cov1(0,1) = corr1(0,1)*(sigX1*sigY1);
 cov1(1,0) = corr1(1,0)*(sigY1*sigX1);
 cov1(1,1) = corr1(1,1)*(sigY1*sigY1);

 cov2(0,0) = corr2(0,0)*(sigX2*sigX2);
 cov2(0,1) = corr2(0,1)*(sigX2*sigY2);
 cov2(1,0) = corr2(1,0)*(sigY2*sigX2);
 cov2(1,1) = corr2(1,1)*(sigY2*sigY2);

 cout << endl;
 cout << "Correlation 1" << endl;
 corr1.Print();
 cout << "Correlation 2" << endl;
 corr2.Print();

 cout << "Covariance 1" << endl;
 cov1.Print();
 cout << "Covariance 2" << endl;
 cov2.Print();

 TMatrixD covSum1 = cov1+cov2;

 TMatrixD covXY(nEle,nEle);
 TMatrixD corrXY(nEle,nEle);
 double x,y,xy,xAvg,yAvg,xyAvg,xAvg2,yAvg2,xyAvg2,xStd,yStd,xyStd;
 double xSum,ySum,xySum,xSum2,ySum2,xySum2;

 for(int j=0;j<nEle;j++){
  for(int i=0;i<nEle;i++){
   xSum = 0.0;
   ySum = 0.0;
   xySum = 0.0;
   xSum2 = 0.0;
   ySum2 = 0.0;
   xySum2 = 0.0;
    for(Long64_t iVec=0;iVec<nVec;iVec++){
     x = allCorr1.at(iVec)(i) + allCorr2.at(iVec)(j); 
     y = allCorr1.at(iVec)(j) + allCorr2.at(iVec)(i);
     xy = x*y;

     xSum += x;
     ySum += y;
     xySum += xy;

     xSum2 += x*x;
     ySum2 += y*y;
     xySum2 += xy*xy;
    }

   xAvg = xSum/(1.0*nVec);
   yAvg = ySum/(1.0*nVec);
   xyAvg = xySum/(1.0*nVec);

   xAvg2 = xSum2/(1.0*nVec);
   yAvg2 = ySum2/(1.0*nVec);
   xyAvg2 = xySum2/(1.0*nVec);

   xStd = sqrt(xAvg2-xAvg*xAvg);
   yStd = sqrt(yAvg2-yAvg*yAvg);
   xyStd = sqrt(xyAvg2-xyAvg*xyAvg);

   covXY(i,j) = (xyAvg-xAvg*yAvg);
   corrXY(i,j) = covXY(i,j)/((xStd*yStd));
  }
 } 

 cout << "Calculated Covariance:" << endl;
 covXY.Print();
 cout << "Cov1 + Cov2" << endl;
 covSum1.Print();
 cout << "Calculated Correlation:" << endl;
 corrXY.Print();
}





