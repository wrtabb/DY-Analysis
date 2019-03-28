

///////////////////////////////////////////////////////////////////////////////////
//This follows the procedure for calculating the covariance matrix as outlined in//
//'Linear Algebra and its Applications' by Lay, Lay, and McDonald on page 427    //
///////////////////////////////////////////////////////////////////////////////////

void makeCorrMatrix()
{
 const int nEle = 3;//number of elements in the vectors
 const int nVec = 4;//number of vectors 

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
 covM.Print();
 corrM.Print();
}

