//wrapper to readXMLGaussianMixtureModel

#include "mex.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "GaussianMixtureModel.h"
#include "external/xmlParser2/svlStrUtils.h"
#include "external/xmlParser2/xmlParser.h"


using namespace std;

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) 
{
    
    char *input_buf;
    
    
    if (nrhs!=1 || nlhs>1)
        mexErrMsgTxt("The number of input and/or output arguments os wrong."); 
    
    /* input must be a string */
    if ( mxIsChar(prhs[0]) != 1)
      mexErrMsgTxt("Input must be a string.");

    /* copy the string data from prhs[0] into a C string input_ buf.    */
    input_buf = mxArrayToString(prhs[0]);
    
    //check if file exists
    ifstream inCheck(input_buf);
    if(!inCheck.is_open())
        mexErrMsgTxt("Requested file does not exist");
    inCheck.close();
   
    
    mexPrintf("Reading initial solution from xml file %s\n",input_buf);
    XMLNode xMainNode=XMLNode::openFileHelper(input_buf,"document");
    vector<GaussianMixtureModel*> stack;
    

    //read stack
	int n=xMainNode.nChildNode("GaussianMixtureModel");

    
    for(int ii=0;ii<n;ii++) stack.push_back(new GaussianMixtureModel(xMainNode,ii));
    
    
    mexPrintf("Read %d mixtures in XML file\n",(int)(stack.size()));
    
    
    //allocate memory and prepare struct to return
    const int nfields=18;//it has to match fieldnames
    const char *fieldnames[]={"m","id","dims","scale","nu","beta","alpha","W","nuPrior","betaPrior","alphaPrior","mPrior","WPrior","sigmaDistPrior","splitScore","lineage","parent","svIdx"};       /* pointers to field names */
        
    
    /* create a Nx1 struct matrix for output  */
    plhs[0] = mxCreateStructMatrix(stack.size(), 1, nfields, fieldnames);
    
    
    const size_t vecSize=sizeof(double)*dimsImage;
    const size_t matrixSize=sizeof(double)*dimsImage*dimsImage;
    double *ptr;
    for(unsigned int ii=0;ii<stack.size();ii++)
    {
        
        GaussianMixtureModel *GM=stack[ii];
        
        //set m
        mxArray *field_value_0=mxCreateDoubleMatrix(1,dimsImage,mxREAL);
        ptr=mxGetPr(field_value_0);
        //for(int jj=0;jj<dimsImage;jj++) ptr[jj]=GM->m_k(jj);
        memcpy(ptr,GM->m_k.data(),vecSize);
        mxSetFieldByNumber(plhs[0],ii,0,field_value_0);//field number is 0-indexed
        
        //set id
        mxArray *field_value_1=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(field_value_1) = GM->id;
        mxSetFieldByNumber(plhs[0],ii,1,field_value_1);//field number is 0-indexed
        
        //set dims
        mxArray *field_value_2=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(field_value_2) = dimsImage;
        mxSetFieldByNumber(plhs[0],ii,2,field_value_2);//field number is 0-indexed
        
        //set scale
        mxArray *field_value_3=mxCreateDoubleMatrix(1,dimsImage,mxREAL);
        ptr=mxGetPr(field_value_3);
        for(int jj=0;jj<dimsImage;jj++) ptr[jj]=GaussianMixtureModel::scale[jj];
        mxSetFieldByNumber(plhs[0],ii,3,field_value_3);//field number is 0-indexed
        
        //set nu
        mxArray *field_value_4=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(field_value_4) = GM->nu_k;
        mxSetFieldByNumber(plhs[0],ii,4,field_value_4);//field number is 0-indexed
        
        //set beta
        mxArray *field_value_5=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(field_value_5) = GM->beta_k;
        mxSetFieldByNumber(plhs[0],ii,5,field_value_5);//field number is 0-indexed
        
        //set alpha
        mxArray *field_value_6=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(field_value_6) = GM->alpha_k;
        mxSetFieldByNumber(plhs[0],ii,6,field_value_6);//field number is 0-indexed
        
        //set W
        mxArray *field_value_7=mxCreateDoubleMatrix(dimsImage,dimsImage,mxREAL);
        ptr=mxGetPr(field_value_7);
        memcpy(ptr,GM->W_k.data(),matrixSize);
        mxSetFieldByNumber(plhs[0],ii,7,field_value_7);//field number is 0-indexed
        
        
        //set nuPrior
        mxArray *field_value_8=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(field_value_8) = GM->nu_o;
        mxSetFieldByNumber(plhs[0],ii,8,field_value_8);//field number is 0-indexed
        
        //set betaPrio
        mxArray *field_value_9=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(field_value_9) = GM->beta_o;
        mxSetFieldByNumber(plhs[0],ii,9,field_value_9);//field number is 0-indexed
        
        //set alphaPrior
        mxArray *field_value_10=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(field_value_10) = GM->alpha_o;
        mxSetFieldByNumber(plhs[0],ii,10,field_value_10);//field number is 0-indexed
        
        //set mPrior
        mxArray *field_value_11=mxCreateDoubleMatrix(1,dimsImage,mxREAL);
        ptr=mxGetPr(field_value_11);
        //for(int jj=0;jj<dimsImage;jj++) ptr[jj]=GM->m_o(jj);
        memcpy(ptr,GM->m_o.data(),vecSize);
        mxSetFieldByNumber(plhs[0],ii,11,field_value_11);//field number is 0-indexed
        
        
        //set WPrior
        mxArray *field_value_12=mxCreateDoubleMatrix(dimsImage,dimsImage,mxREAL);
        ptr=mxGetPr(field_value_12);
        memcpy(ptr,GM->W_o.data(),matrixSize);
        mxSetFieldByNumber(plhs[0],ii,12,field_value_12);//field number is 0-indexed
        
        //set sigma dist split score
        mxArray *field_value_13=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(field_value_13) = GM->sigmaDist_o;
        mxSetFieldByNumber(plhs[0],ii,13,field_value_13);//field number is 0-indexed
        
        //set split score
        mxArray *field_value_14=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(field_value_14) = GM->splitScore;
        mxSetFieldByNumber(plhs[0],ii,14,field_value_14);//field number is 0-indexed
        
        //set lineage
        mxArray *field_value_15=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(field_value_15) = GM->lineageId;
        mxSetFieldByNumber(plhs[0],ii,15,field_value_15);//field number is 0-indexed
        
        //set parent
        mxArray *field_value_16=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(field_value_16) = GM->parentId;
        mxSetFieldByNumber(plhs[0],ii,16,field_value_16);//field number is 0-indexed
        
        //set supervoxel idx
        mxArray *field_value_17=mxCreateDoubleMatrix(1,GM->supervoxelIdx.size(), mxREAL);
        ptr=mxGetPr(field_value_17);
        for(size_t aa = 0; aa < GM->supervoxelIdx.size(); aa++)
            ptr[aa] = GM->supervoxelIdx[aa];      
        mxSetFieldByNumber(plhs[0],ii,17,field_value_17);//field number is 0-indexed
    }
    

    //release memory
    mxFree(input_buf);       
    for(int ii=0;ii<n;ii++) delete stack[ii];
    stack.clear();
    
    return;
} 