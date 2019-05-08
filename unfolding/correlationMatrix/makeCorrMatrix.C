

///////////////////////////////////////////////////////////////////////////////////
//This follows the procedure for calculating the covariance matrix as outlined in//
//'Linear Algebra and its Applications' by Lay, Lay, and McDonald on page 427    //
///////////////////////////////////////////////////////////////////////////////////

void makeCorrMatrix()
{
 const int nEle = 3;//number of elements in the vectors
 const int nVec = 4;//number of vectors 
 TH1::SetDefaultSumw2();
 gStyle->SetOptStat(0);
 gStyle->SetPalette(1);
 //All vectors stored in a 2D array
 //There are nVec vectors with nEle elements in each vector
 //double vector[nVec][nEle];
 
 //Defined vector for testing
 double vector[nVec][nEle] = {
  {1,2,1},
  {4,2,13},
  {7,8,1},
  {8,4,5}
};
 
 //Initializing matrices
 TMatrixD matrixB(nEle,nVec);
 TMatrixD matrixBT(nVec,nEle);
 TMatrixD covM(nEle,nEle);
 TMatrixD corrM(nEle,nEle);

 //Initialize random number generator
 TRandom3*random = new TRandom3(time(0));

 //Initialize variables and arrays
 double rndEle,vecSum[nEle],vecAvg[nEle],vecM[nVec][nEle],vecSum2[nEle],vecStd[nEle];
 
 //Set sums to zer 
 for(int i=0;i<nEle;i++){
  vecSum[i] = 0;
  vecSum2[i]=0;
 }

 for(int iEle=0;iEle<nEle;iEle++){
  for(int jVec=0;jVec<nVec;jVec++){
   //rndEle = 8*random->Rndm();//random number
   //vector[jVec][iEle] = rndEle;//random number placed in vector
   vecSum[iEle] += vector[jVec][iEle];//sum of element iEle summed over vectors
   vecSum2[iEle] += vector[jVec][iEle]*vector[jVec][iEle];//sum of squares of elements
  }
  vecAvg[iEle] = vecSum[iEle]/nVec;//vector of averages
  //vector of standard deviations
  vecStd[iEle] = sqrt(vecSum2[iEle]/nVec-vecSum[iEle]*vecSum[iEle]/(nVec*nVec*1.0));
 }
 
 for(int iEle=0;iEle<nEle;iEle++){
  for(int jVec=0;jVec<nVec;jVec++){
   //Calculation of matrix made up of each vector minus the average vector
   vecM[jVec][iEle] = vector[jVec][iEle]-vecAvg[iEle];//vectors minus average vector
   matrixB(iEle,jVec) = vecM[jVec][iEle];
  }
 }

 matrixBT.Transpose(matrixB);//the transpose of matrixB calculated above
 covM = matrixB*matrixBT;//covariance matrix without proper weight
 for(int jEle=0;jEle<nEle;jEle++){
  for(int iEle=0;iEle<nEle;iEle++){
   covM(iEle,jEle) = covM(iEle,jEle)/(nVec);//covariance matrix
   corrM(iEle,jEle) = covM(iEle,jEle)/(vecStd[iEle]*vecStd[iEle]);//correlation matrix
  }
 }

 //Output of all vectors and matrices
 TString vecName;
 TString avgMName;
 for(int j=1;j<nVec+1;j++){
  vecName = "Vector";
  vecName += j;
  vecName += ": (";
  cout << vecName;
  for(int i=0;i<nEle;i++){
   cout << vector[j-1][i];
   if(i!=nEle-1) cout << ", ";
  } 
  cout << ")" << endl;
 }
 //Print average vector
 cout << "Avg Vector: (";
 for(int i=0;i<nEle;i++){
  cout << vecAvg[i];
  if(i!=nEle-1) cout << ", ";
 }
 cout <<")" << endl;

 for(int j=1;j<nVec+1;j++){
  avgMName= "Vector-Avg";
  avgMName+= j;
  avgMName+= ": (";
  cout << avgMName;
  for(int i=0;i<nEle;i++){
   cout << vecM[j-1][i];
   if(i!=nEle-1) cout << ", ";
  } 
  cout << ")" << endl;
 }
 matrixB.Print();
 matrixBT.Print(); 
 cout << "The covariance matrix is:" << endl;
 covM.Print();
 cout << "The correlations matrix is:" << endl;
 corrM.Print();


 //Now calculate covariance in a different way for comparison
 cout << "Starting manual covariance calculation ... " << endl;
 TH2D*hCov = new TH2D("hCov","",3,0,3,3,0,3);
 TH2D*hCorr = new TH2D("hCorr","",3,0,3,3,0,3);
 double x,y,xy,xAvg,yAvg,xyAvg,xAvg2,yAvg2,xyAvg2,xStd,yStd,xyStd,covXY,corrXY;
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
     x = vector[iVec][i];
     y = vector[iVec][j];
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
 
   covXY = (xyAvg-xAvg*yAvg);
   corrXY = covXY/((xStd*yStd));

   cout << "Covariance: " << endl;
   cout << "i = " << i << ", j = " << j << endl;
   cout << covXY << endl;
   cout << endl;

   hCov->SetBinContent(i+1,j+1,covXY);
   hCorr->SetBinContent(i+1,j+1,corrXY);
  }
 }
 TCanvas*canvas = new TCanvas("canvas","",10,10,1200,1000);
 hCov->Draw("colz");
 TCanvas*canvas2 = new TCanvas("canvas2","",10,10,1200,1000);
 hCorr->Draw("colz");
 
}
