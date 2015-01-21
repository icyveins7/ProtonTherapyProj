#include "mex.h"
#include<math.h>

/*
 * Bortfeld_1997dll.c
 *
 * 1F1 Hypergeometric Function (for real numbers)
 *
 * Input with a, b, z, precision(accuracy in decimal places)
 *
 * The calling syntax is:
 *
 *		out = F11(a, b, z, precision)
 *
 * This is a MEX-file for MATLAB, compiled with version 8.4.0.150421(R2014b).
*/
void F11(double a, double b, double z, int precision, double *out){
	//current term
	double term = 1;
	double total = 1;
	for (mwSize k = 1; k >= 0; k++){
		term = ((term*z) / k)*(a + k - 1) / (b + k - 1);
		total = total + term;
		//breaks if term is less than precision or total has overflowed
		if (term < pow(10, -(double)precision) || isinf(total)){
			break;
		}
	}
	*out=total;
}

/* The gateway function */
void mexFunction( mwSize nlhs, mxArray *plhs[],
                  mwSize nrhs, const mxArray *prhs[])
{
    double a, b, z;             
    mwSize precision;                                
    double *out;            
    
    /* check for proper number of arguments */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","4 inputs required.");
    }
//     if(nlhs!=1) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
//     }
    /* make sure the first input argument is scalar */
    if(  !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input a must be type double.");
    }
    
    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input b must be type double.");
    }
    
    /* make sure the 3rd input argument is type double */
    if( !mxIsDouble(prhs[2]) || 
         mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input z must be type double.");
    }
    
    /* get the value of the scalar inputs  */
    a = mxGetScalar(prhs[0]);
    b = mxGetScalar(prhs[1]);
    z = mxGetScalar(prhs[2]);
    precision = mxGetScalar(prhs[3]);

    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleScalar(0);
    
    /* get a pointer to the real data in the output matrix */
    out = mxGetPr(plhs[0]);

    /* call the computational routine */
    F11(a, b, z, precision, out);
}