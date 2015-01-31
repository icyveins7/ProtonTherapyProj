#include "mex.h"
#include<math.h>

/*
 * gamma_C.c
 *
 * Gamma Function (for real numbers)
 *
 * Input with z
 *
 * The calling syntax is:
 *
 *		out = gamma_C(z)
 *
 * This is a MEX-file for MATLAB, compiled with version 8.4.0.150421(R2014b),
 * using Microsoft Visual Studio 2013.
 *
 * If unable to load, download Visual C++ Redistributable Packages for Visual Studio 2013
 * http://www.microsoft.com/en-us/download/details.aspx?id=40784
 *
*/
void gamma_C(double z, double *out){
	*out=tgamma(z);
}

/* The gateway function */
void mexFunction( mwSize nlhs, mxArray *plhs[],
                  mwSize nrhs, const mxArray *prhs[])
{
    double z;                                           
    double *out;            
    
    /* check for proper number of arguments */
    if(nrhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","1 input required.");
    }
//     if(nlhs!=1) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
//     }
    /* make sure the first input argument is scalar */
    if(  !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input z must be type double.");
    }
    
//     /* make sure the second input argument is type double */
//     if( !mxIsDouble(prhs[1]) || 
//          mxIsComplex(prhs[1])) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input b must be type double.");
//     }
//     
//     /* make sure the 3rd input argument is type double */
//     if( !mxIsDouble(prhs[2]) || 
//          mxIsComplex(prhs[2])) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input z must be type double.");
//     }
    
    /* get the value of the scalar inputs  */
    z = mxGetScalar(prhs[0]);
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleScalar(0);
    
    /* get a pointer to the real data in the output matrix */
    out = mxGetPr(plhs[0]);

    /* call the computational routine */
    gamma_C(z, out);
}