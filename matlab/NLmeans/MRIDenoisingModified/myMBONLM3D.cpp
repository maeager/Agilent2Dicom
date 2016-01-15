/*

  Modified NL means filter MRONLM

  Michael Eager
  Zhaolin Chen

  Monash Biomedical Imaging
  Monash University, 2015

 */



/* Pierrick Coupe - pierrick.coupe@gmail.com                               */
/* Jose V. Manjon - jmanjon@fis.upv.es                                     */
/* Brain Imaging Center, Montreal Neurological Institute.                  */
/* Mc Gill University                                                      */
/*                                                                         */
/* Copyright (C) 2008 Pierrick Coupe and Jose V. Manjon                    */

/***************************************************************************
 *              3D Adaptive Multiresolution Non-Local Means Filter          *
 * Pierrick Coupe a, Jose V. Manjon, Montserrat Robles , D. Louis Collins   *
 ***************************************************************************/


/*                          Details on ONLM filter                        */
/***************************************************************************
 *  The ONLM filter is described in:                                       *
 *                                                                         *
 *  P. Coup�, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     *
 *  An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic*
 *  Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441, *
 *  Avril 2008                                                             *
 ***************************************************************************/


#include "math.h"
#include "mex.h"
#include <stdlib.h>
#include "matrix.h"

/* Multithreading stuff*/
#ifdef _WIN32
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#endif

#define max(x,y) x < y ? y : x
#define min(x,y) x > y ? y : x
enum gamma_kernel_method {GAMMAMINMAX = 0, GAMMAMULT, GAMMADR, GAMMAEXPDR, GAMMADIFF};

int GAMMAfunction = 0;


typedef struct {
    int rows;
    int cols;
    int slices;
    float * in_image;
    float * means_image;
    float * var_image;
    float * coilsens_image;
    float * estimate;
    float * weight;
    float * label;
    int ini;
    int fin;
    int radiusB;
    int radiusS;
    float sigma;
} myargument;

bool rician = false;

/* Function which compute the weighted average for one block */
void Average_block(float *ima, int x, int y, int z, int neighborhoodsize, float *average, float weight, int sx, int sy, int sz)
{
    int x_pos, y_pos, z_pos;
    bool is_outside;
    int a, b, c, ns, sxy, count, p1, p2;

    ns = 2 * neighborhoodsize + 1;
    sxy = sx * sy;

    count = 0;

    for (c = 0; c < ns; c++)
    {
        z_pos = z + c - neighborhoodsize;
        if ((z_pos < 0) || (z_pos > sz - 1)) z_pos = z;
        for (b = 0; b < ns; b++)
        {
            y_pos = y + b - neighborhoodsize;
            if ((y_pos < 0) || (y_pos > sy - 1)) y_pos = y;
            for (a = 0; a < ns; a++)
            {
                x_pos = x + a - neighborhoodsize;
                if ((x_pos < 0) || (x_pos > sx - 1)) x_pos = x;
#ifdef FP_FAST_FMA
		p2 = fma(z_pos,sxy,fma(y_pos,sx,x_pos));
#else


                p2 = z_pos * (sxy) + (y_pos * sx) + x_pos;
#endif
                if (rician){
#ifdef FP_FAST_FMA
		  average[count]=fma (ima[p2] * ima[p2],weight,average[count]);
#else
                    average[count] = average[count] + ima[p2] * ima[p2] * weight;
#endif
                }else{
#ifdef FP_FAST_FMA
		  average[count]=fma (ima[p2],weight,average[count]);
#else
                    average[count] = average[count] + ima[p2] * weight;
#endif
		}
                count++;
            }
        }
    }

}

/* Function which computes the value assigned to each voxel */
void Value_block(float *Estimate, float *Label, int x, int y, int z, int neighborhoodsize, float *average, float global_sum, int sx, int sy, int sz, float bias)
{
    int x_pos, y_pos, z_pos, p1;
    int ret;
    bool is_outside;
    float value = 0.0;
    float label = 0.0;
    float denoised_value = 0.0;
    int count = 0 ;
    int a, b, c, ns, sxy;

    extern bool rician;

    ns = 2 * neighborhoodsize + 1;
    sxy = sx * sy;


    for (c = 0; c < ns; c++)
    {
        z_pos = z + c - neighborhoodsize;
        if (!((z_pos < 0) || (z_pos > sz - 1)))
            for (b = 0; b < ns; b++)
            {
                y_pos = y + b - neighborhoodsize;
                if (!((y_pos < 0) || (y_pos > sy - 1)))
                    for (a = 0; a < ns; a++)
                    {
                        x_pos = x + a - neighborhoodsize;
                        if (!((x_pos < 0) || (x_pos > sx - 1)))
                        {
#ifdef FP_FAST_FMA
			  p1 = fma(z_pos,sxy,fma(y_pos,sx,x_pos));
#else
                            p1 = z_pos * (sxy) + (y_pos * sx) + x_pos;
#endif
                            value = Estimate[p1];

                            if (rician)
                            {
#ifdef FP_FAST_FMA
			      denoised_value = fma(average[count],1.0f/global_sum,bias);
#else
                                denoised_value  = (average[count] / global_sum) - bias;
#endif
                                if (denoised_value > 0)
                                    denoised_value = sqrt(denoised_value);
                                else {
                                    denoised_value = 0.0;
                                }
                                value = value + denoised_value;
                            }
                            else {
#ifdef FP_FAST_FMA
			      value = fma(average[count],1.0f/global_sum,value);
#else
                                value = value + (average[count] / global_sum);
#endif
                            }

                            //label = Label[p1];
                            Estimate[p1] = value;
                            Label[p1] += 1;
                        }
                        count++;
                    }
            }
    }
}


float distance(float* ima, int x, int y, int z, int nx, int ny, int nz, int f, int sx, int sy, int sz)
{

    float d, acu, distancetotal, inc;
    int i, j, k, ni1, nj1, ni2, nj2, nk1, nk2, kk, sxy, p1, p2;

    sxy = sx * sy;

    acu = 0;
    distancetotal = 0;

    for (k = -f; k <= f; k++)
    {
        nk1 = z + k;
        nk2 = nz + k;
        if (nk1 < 0) nk1 = -nk1;
        if (nk2 < 0) nk2 = -nk2;
        if (nk1 >= sz) nk1 = 2 * sz - nk1 - 1;
        if (nk2 >= sz) nk2 = 2 * sz - nk2 - 1;
        for (j = -f; j <= f; j++)
        {
            nj1 = y + j;
            nj2 = ny + j;
            if (nj1 < 0) nj1 = -nj1;
            if (nj2 < 0) nj2 = -nj2;
            if (nj1 >= sy) nj1 = 2 * sy - nj1 - 1;
            if (nj2 >= sy) nj2 = 2 * sy - nj2 - 1;
            for (i = -f; i <= f; i++)
            {
                ni1 = x + i;
                ni2 = nx + i;
                if (ni1 < 0) ni1 = -ni1;
                if (ni2 < 0) ni2 = -ni2;
                if (ni1 >= sx) ni1 = 2 * sx - ni1 - 1;
                if (ni2 >= sx) ni2 = 2 * sx - ni2 - 1;

                p1 = nk1 * (sxy) + (nj1 * sx) + ni1;
                p2 = nk2 * (sxy) + (nj2 * sx) + ni2;
#ifdef FP_FAST_FMA
		distancetotal = fma((ima[p1] - ima[p2]) , (ima[p1] - ima[p2]),distancetotal);
#else
                distancetotal = distancetotal + ((ima[p1] - ima[p2]) * (ima[p1] - ima[p2]));
#endif     
           acu = acu + 1;
            }
        }
    }
    d = distancetotal / acu;

    return d;

}

float distanceB1(float* ima, float* coilsens, int x, int y, int z, int nx, int ny, int nz, int f, int sx, int sy, int sz)
{

    float d, acu, distancetotal, inc;
    int i, j, k, ni1, nj1, ni2, nj2, nk1, nk2, kk, sxy, p1, p2;

    sxy = sx * sy;

    acu = 0;
    distancetotal = 0;

    for (k = -f; k <= f; k++)
    {
        nk1 = z + k;
        nk2 = nz + k;
        if (nk1 < 0) nk1 = -nk1;
        if (nk2 < 0) nk2 = -nk2;
        if (nk1 >= sz) nk1 = 2 * sz - nk1 - 1;
        if (nk2 >= sz) nk2 = 2 * sz - nk2 - 1;
        for (j = -f; j <= f; j++)
        {
            nj1 = y + j;
            nj2 = ny + j;
            if (nj1 < 0) nj1 = -nj1;
            if (nj2 < 0) nj2 = -nj2;
            if (nj1 >= sy) nj1 = 2 * sy - nj1 - 1;
            if (nj2 >= sy) nj2 = 2 * sy - nj2 - 1;
            for (i = -f; i <= f; i++)
            {
                ni1 = x + i;
                ni2 = nx + i;
                if (ni1 < 0) ni1 = -ni1;
                if (ni2 < 0) ni2 = -ni2;
                if (ni1 >= sx) ni1 = 2 * sx - ni1 - 1;
                if (ni2 >= sx) ni2 = 2 * sx - ni2 - 1;

                p1 = nk1 * (sxy) + (nj1 * sx) + ni1;
                p2 = nk2 * (sxy) + (nj2 * sx) + ni2;
#ifdef FP_FAST_FMA
		distancetotal = fma((ima[p1]/ coilsens[p1] - ima[p2]/ coilsens[p2]) , (ima[p1]/ coilsens[p1] - ima[p2]/ coilsens[p2]),distancetotal);
#else

                distancetotal += (ima[p1] / coilsens[p1] - ima[p2] / coilsens[p2]) * (ima[p1] / coilsens[p1] - ima[p2] / coilsens[p2]);
#endif
                acu = acu + 1;
            }
        }
    }
    d = distancetotal / acu;

    return d;

}


/*calculate gamma */
float gammac(float* coilsens, int ni1, int nj1, int nk1, int ni2, int nj2, int nk2, int sx, int sy, int sz, int method)
{

    float gamma;
    int sxy, p1, p2;

    sxy = sx * sy;

    if (ni1 < 0) ni1 = -ni1;
    if (nj1 < 0) nj1 = -nj1;
    if (ni2 < 0) ni2 = -ni2;
    if (nj2 < 0) nj2 = -nj2;
    if (nk1 < 0) nk1 = -nk1;
    if (nk2 < 0) nk2 = -nk2;

    if (ni1 >= sx) ni1 = 2 * sx - ni1 - 1;
    if (nj1 >= sy) nj1 = 2 * sy - nj1 - 1;
    if (nk1 >= sz) nk1 = 2 * sz - nk1 - 1;
    if (ni2 >= sx) ni2 = 2 * sx - ni2 - 1;
    if (nj2 >= sy) nj2 = 2 * sy - nj2 - 1;
    if (nk2 >= sz) nk2 = 2 * sz - nk2 - 1;

    p1 = nk1 * (sxy) + (nj1 * sx) + ni1;
    p2 = nk2 * (sxy) + (nj2 * sx) + ni2;

    if (method == GAMMAMINMAX) {
        // Max-min method
        gamma = max(coilsens[p1], coilsens[p2]) / min(coilsens[p1] , coilsens[p2]);
    }
    else if (method == GAMMAEXPDR) {
        // Exponential DR method
        gamma = exp(fabs(coilsens[p1] - coilsens[p2]));
    }
    else if (method == GAMMADR) {
        // DR method
        gamma = (1 + (fabs(coilsens[p1] - coilsens[p2])));
    }
    else if (method == GAMMADIFF) {
        // DR method
        gamma = ((fabs(coilsens[p1] - coilsens[p2])));
    }
    return gamma;
}


void* ThreadFunc(void* pArguments)
{
    float  *Estimate, *Label, *Weight, *ima, *means, *coilsens, *variances, *average;
    float bias, sigma, gamma, epsilon, mu1, var1, totalweight, wmax, t1, t2, d, w, h, hh;
    int rows, cols, slices, ini, fin, radiusB, radiusS, init, Ndims; //constants
    int i, j, k, rc, ii, jj, kk, ni, nj, nk; //loop variables
    int p1, p2; //pixel index one and two

    myargument arg;
    arg = *(myargument *) pArguments;

    rows = arg.rows;
    cols = arg.cols;
    slices = arg.slices;
    ini = arg.ini;
    fin = arg.fin;
    ima = arg.in_image;
    means = arg.means_image;
    coilsens = arg.coilsens_image;
    variances = arg.var_image;
    Estimate = arg.estimate;
    Weight = arg.weight;
    Label = arg.label;
    radiusB = arg.radiusB;
    radiusS = arg.radiusS;
    sigma = arg.sigma;

    h = sigma;
    hh = h * h;
    bias = 2 * sigma * sigma;

    /*filter*/
    epsilon = 0.00001;
    mu1 = 0.95;
    var1 = 0.5;
    init = 0;
    rc = rows * cols;

    Ndims = (2 * radiusS + 1) * (2 * radiusS + 1) * (2 * radiusS + 1);

    average = (float *) malloc(Ndims * sizeof(float));

    for (k = ini; k < fin; k += 2)
    {
        for (j = 0; j < rows; j += 2)
        {
            for (i = 0; i < cols; i += 2)
            {
                for (init = 0 ; init < Ndims; init++)
                    average[init] = 0.0;

                /*average=0;*/
                totalweight = 0.0;
                p1 = k * rc + (j * cols) + i;

                if ((means[p1]) > epsilon && (variances[p1] > epsilon))
                {
                    wmax = 0.0;
                    for (kk = -radiusB; kk <= radiusB; kk++)
                    {
                        nk = k + kk;
                        if (nk >= 0 && nk < slices)
                            for (jj = -radiusB; jj <= radiusB; jj++)
                            {
                                nj = j + jj;
                                if (nj >= 0 && nj < rows)
                                    for (ii = -radiusB; ii <= radiusB; ii++)
                                    {
                                        ni = i + ii;
                                        if (ii == 0 && jj == 0 && kk == 0) continue;

                                        if (ni >= 0 && ni < cols)
                                        {
                                            p2 = nk * (rc) + (nj * cols) + ni;

                                            if ((means[p2]) > epsilon && (variances[p2] > epsilon))
                                            {
                                                t1 = (means[p1]) / (means[p2]);
                                                t2 = (variances[p1]) / (variances[p2]);

                                                if (t1 > mu1 && t1 < (1 / mu1) && t2 > var1 && t2 < (1 / var1))
                                                {
                                                    if (GAMMAfunction == GAMMAMULT) {
                                                        d = distanceB1(ima, coilsens, i, j, k, ni, nj, nk, radiusS, cols, rows, slices);
                                                        gamma = 1.0;
                                                    }
                                                    else {
                                                        d = distance(ima, i, j, k, ni, nj, nk, radiusS, cols, rows, slices);
                                                        // Eager Ammendment
                                                        // add coil sensitivty B1 correction factor to weight
                                                        // gamma = max(coilsens[p1], coilsens[p2]) / min(coilsens[p1], coilsens[p2]);
                                                        gamma = gammac(coilsens, i, j, k, ni, nj, nk, cols, rows, slices, GAMMAfunction);
                                                    }

                                                    if (gamma > epsilon) {
                                                        w = exp(-(1 / (gamma * gamma)) * d / hh);
                                                    }
                                                    else {
                                                        w = exp(-d / (hh));
                                                    }
                                                    if (w > wmax) wmax = w;

                                                    Average_block(ima, ni, nj, nk, radiusS, average, w, cols, rows, slices);

                                                    totalweight = totalweight + w;
                                                }
                                            }
                                        }
                                    }
                            }
                    }
                    if (wmax == 0.0) wmax = 1.0;

                    Average_block(ima, i, j, k, radiusS, average, wmax, cols, rows, slices);

                    totalweight = totalweight + wmax;

                    if (totalweight != 0.0)
                        Value_block(Estimate, Label, i, j, k, radiusS, average, totalweight, cols, rows, slices, bias);

                }
                else
                {
                    wmax = 1.0;
                    Average_block(ima, i, j, k, radiusS, average, wmax, cols, rows, slices);
                    totalweight = totalweight + wmax;
                    Value_block(Estimate, Label, i, j, k, radiusS, average, totalweight, cols, rows, slices, bias);

                }
                Weight[ p1 ] = totalweight;
                /*Average_block(ima,i,j,k,radiusS,average,wmax,cols,rows,slices);     */
            }
        }
    }
#ifdef _WIN32
    _endthreadex(0);
#else
    pthread_exit(0);
#endif

    return 0;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*Declarations*/
    mxArray *xData;
    float *ima, *fima, *average, *coilsens;
    mxArray *Mxmeans, *Mxvariances, *Mxweight, *MxEstimate, *MxLabel, *MxCoil;
    float *means, *variances, *Estimate, *Label, *weight;
    mxArray *pv;
    float hSigma, w, totalweight, wmax, d, mean, var, t1, t2, hh, epsilon, mu1, var1, label, estimate;
    int Ndims, i, j, k, ii, jj, kk, ni, nj, nk, p1, radiusB, radiusS, ndim, indice, init, Nthreads, ini, fin, r;
    const int  *dims, *coildims;

    myargument *ThreadArgs;

#ifdef _WIN32
    HANDLE *ThreadList; /* Handles to the worker threads*/
#else
    pthread_t * ThreadList;
#endif

    if (nrhs < 5 || nrhs > 7)
    {
        mexErrMsgIdAndTxt("myMBONLM3D:mexfunction",
                          "Not enough arguments\nUsage: OutputImage=myMBONLM3D(InputImage,searcharea,patcharea,sigma,rician,[Coilsens],[GAMMAfunction])");
    }

    /*Copy input pointer x*/

    /*Get matrix x*/
    if (mxIsSparse(prhs[0]) ||
            mxIsComplex(prhs[0]) ||
            mxIsDouble(prhs[0])  ||
            mxGetNumberOfElements(prhs[0]) == 1 ||
            mxGetNumberOfDimensions(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("myMBONLM:mexfunction", "input1 must be full matrix of real float values.");

    }
    ima = (float*) mxGetPr(prhs[0]);

    ndim = mxGetNumberOfDimensions(prhs[0]);
    dims = mxGetDimensions(prhs[0]);
    mexPrintf("myMBONLM3D: Starting\n");  mexEvalString("drawnow");


    /*Get the patch area, search area and sigma*/
    if (mxIsComplex(prhs[1]) ||
            mxGetNumberOfElements(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("myMBONLM:mexfunction", " patch area must be an integer.");

    }
    radiusB = (int)(mxGetScalar(prhs[1]));

    if (mxIsComplex(prhs[2]) ||
            mxGetNumberOfElements(prhs[2]) != 1) {
        mexErrMsgIdAndTxt("myMBONLM:mexfunction", " search area must be an integer.");
    }
    radiusS = (int)(mxGetScalar(prhs[2]));

    if (mxIsSparse(prhs[3]) ||
            mxIsComplex(prhs[3]) ||
            mxGetNumberOfElements(prhs[3]) != 1) {
        mexErrMsgIdAndTxt("myMBONLM:mexfunction", " sigma must be real value.");
    }
    hSigma = (float)(mxGetScalar(prhs[3]));

    if (mxIsSparse(prhs[4]) ||
            mxIsComplex(prhs[4]) ||
            mxGetNumberOfElements(prhs[4]) != 1) {
        mexErrMsgIdAndTxt("myMBONLM:mexfunction", " rician must be real value or bool.");
    }
    r = (int)(mxGetScalar(prhs[4]));
    if (r > 0) rician = true;

    if (nrhs >= 6) {
        if (mxIsSparse(prhs[5]) ||
                mxIsComplex(prhs[5]) ||
                mxIsDouble(prhs[5])  ||
                mxGetNumberOfElements(prhs[5]) == 1 ||
                mxGetNumberOfDimensions(prhs[5]) != 3)
        {
            mexErrMsgIdAndTxt("myMBONLM:mexfunction", " coil sens must be full matrix of real float values.");
        }
        coildims = mxGetDimensions(prhs[5]);
        if (coildims[0] != dims[0] || coildims[1] != dims[1] || coildims[2] != dims[2])
        {
            mexErrMsgIdAndTxt("myMBONLM:mexfunction", " coil dims must equal input image dims.");
        }
        coilsens = (float*) mxGetPr(prhs[5]);
    }
    else {
        MxCoil = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);
        coilsens = (float*) mxGetPr(MxCoil);
        for (i = 0; i < dims[2] *dims[1] * dims[0]; i++)
            coilsens[i] = 1.0;
    }
    if (nrhs >= 7) {
        if (mxIsSparse(prhs[6]) ||
                mxIsComplex(prhs[6]) ||
                mxGetNumberOfElements(prhs[6]) != 1) {
            mexErrMsgIdAndTxt("myMBONLM:mexfunction", " GAMMAfunction must be between 0 and 4.");
        }

        GAMMAfunction = (int)(mxGetScalar(prhs[6]));

        switch (GAMMAfunction) {
        case GAMMAMULT:
            GAMMAfunction = GAMMAMULT;
            mexPrintf("myMBONLM: GAMMAfunction set to GAMMAMULT.\n");break;
        case GAMMADR:
            GAMMAfunction = GAMMADR;
            mexPrintf("myMBONLM: GAMMAfunction set to GAMMADR.\n");break;
        case GAMMAEXPDR:
            GAMMAfunction = GAMMAEXPDR;
            mexPrintf("myMBONLM: GAMMAfunction set to GAMMAEXPDR.\n");break;
        case GAMMADIFF:
            GAMMAfunction = GAMMADIFF;
            mexPrintf("myMBONLM: GAMMAfunction set to GAMMADIFF.\n");break;
        otherwise:
            GAMMAfunction = GAMMAMINMAX;
            mexPrintf("myMBONLM: GAMMAfunction set to default GAMMAMINMAX.\n");
        }
    }
    else {
        GAMMAfunction = GAMMAMINMAX;
        mexPrintf("myMBONLM: GAMMAfunction set to default (GAMMAMINMAX).\n");
    }




    Ndims = (int)pow((2 * radiusS + 1), ndim);

    /* Allocate memory and assign output pointer  and */
    /* Get a pointer to the data space in our newly allocated memory */
    plhs[0] = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);
    fima = (float*)mxGetPr(plhs[0]);
    if (nlhs > 1) {
        plhs[1] = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);
        weight = (float*) mxGetPr(plhs[1]);
    }
    else {
        Mxweight = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);
        weight = (float*) mxGetPr(Mxweight);
    }

    if (nlhs > 2) {
        plhs[2] = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);
        means = (float*) mxGetPr(plhs[2]);
    }
    else {
        Mxmeans = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);
        means = (float*) mxGetPr(Mxmeans);
    }
    if (nlhs > 3) {
        plhs[3] = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);
        variances = (float*) mxGetPr(plhs[3]);
    }
    else {
        Mxvariances = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);
        variances = (float*) mxGetPr(Mxvariances);
    }
    if (nlhs > 4) {
        plhs[4] = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);
        Estimate = (float*) mxGetPr(plhs[4]);
    }
    else {
        MxEstimate = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);
        Estimate = (float*) mxGetPr(MxEstimate);
    }
    if (nlhs > 5) {
        plhs[5] = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);
        Label = (float*) mxGetPr(plhs[5]);
    }
    else {
        MxLabel = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);
        Label = (float*) mxGetPr(MxLabel);
    }

    average = (float*) malloc(Ndims * sizeof(float));

    for (i = 0; i < dims[2] *dims[1] * dims[0]; i++)
    {
        Estimate[i] = 0.0;
        Label[i] = 0.0;
        fima[i] = 0.0;
    }
    mexPrintf("  myMBONLM3D: Arguments parsed \n");  mexEvalString("drawnow");


    for (k = 0; k < dims[2]; k++)
    {
        for (i = 0; i < dims[1]; i++)
        {
            for (j = 0; j < dims[0]; j++)
            {
                mean = 0;
                indice = 0;
                for (ii = -1; ii <= 1; ii++)
                {
                    ni = i + ii;
                    if (ni < 0) ni = -ni;
                    if (ni >= dims[1]) ni = 2 * dims[1] - ni - 1;
                    for (jj = -1; jj <= 1; jj++)
                    {
                        nj = j + jj;
                        if (nj < 0) nj = -nj;
                        if (nj >= dims[0]) nj = 2 * dims[0] - nj - 1;
                        for (kk = -1; kk <= 1; kk++)
                        {
                            nk = k + kk;
                            if (nk < 0) nk = -nk;
                            if (nk >= dims[2]) nk = 2 * dims[2] - nk - 1;
                            mean = mean + ima[nk * (dims[0] * dims[1]) + (ni * dims[0]) + nj];
                            indice = indice + 1;

                        }
                    }
                }
                mean = mean / indice;
                means[k * (dims[0]*dims[1]) + (i * dims[0]) + j] = mean;
            }
        }
    }

    for (k = 0; k < dims[2]; k++)
    {
        for (i = 0; i < dims[1]; i++)
        {
            for (j = 0; j < dims[0]; j++)
            {
                var = 0;
                indice = 0;
                for (ii = -1; ii <= 1; ii++)
                {
                    ni = i + ii;
                    if (ni >= 0 &&  ni < dims[1])
                        for (jj = -1; jj <= 1; jj++)
                        {
                            nj = j + jj;
                            if (nj >= 0 &&  nj < dims[0])
                                for (kk = -1; kk <= 1; kk++)
                                {
                                    nk = k + kk;
                                    if (nk > 0 &&  nk < dims[2])
                                    {
                                        var = var + (ima[nk * (dims[0] * dims[1]) + (ni * dims[0]) + nj] - means[k * (dims[0] * dims[1]) + (i * dims[0]) + j]) * (ima[nk * (dims[0] * dims[1]) + (ni * dims[0]) + nj] - means[k * (dims[0] * dims[1]) + (i * dims[0]) + j]);
                                        indice = indice + 1;
                                    }
                                }
                        }
                }
                var = var / (indice - 1);
                variances[k * (dims[0]*dims[1]) + (i * dims[0]) + j] = var;
            }
        }
    }

    Nthreads = 24;
    mexPrintf("  myMBONLM3D: Starting threads\n");  mexEvalString("drawnow");

#ifdef _WIN32

    /* Reserve room for handles of threads in ThreadList*/
    ThreadList = (HANDLE*)malloc(Nthreads * sizeof(HANDLE));
    ThreadArgs = (myargument*) malloc(Nthreads * sizeof(myargument));

    for (i = 0; i < Nthreads; i++)
    {
        /* Make Thread Structure   */
        ini = (i * dims[2]) / Nthreads;
        fin = ((i + 1) * dims[2]) / Nthreads;
        ThreadArgs[i].cols = dims[0];
        ThreadArgs[i].rows = dims[1];
        ThreadArgs[i].slices = dims[2];
        ThreadArgs[i].in_image = ima;
        ThreadArgs[i].var_image = variances;
        ThreadArgs[i].means_image = means;
        ThreadArgs[i].coilsens_image = coilsens;
        ThreadArgs[i].estimate = Estimate;
        ThreadArgs[i].weight = weight;
        ThreadArgs[i].label = Label;
        ThreadArgs[i].ini = ini;
        ThreadArgs[i].fin = fin;
        ThreadArgs[i].radiusB = radiusB;
        ThreadArgs[i].radiusS = radiusS;
        ThreadArgs[i].sigma = hSigma;

        ThreadList[i] = (HANDLE)_beginthreadex(NULL, 0, &ThreadFunc, &ThreadArgs[i] , 0, NULL);

    }

    for (i = 0; i < Nthreads; i++) {
        WaitForSingleObject(ThreadList[i], INFINITE);
    }
    for (i = 0; i < Nthreads; i++) {
        CloseHandle(ThreadList[i]);
    }

#else

    /* Reserve room for handles of threads in ThreadList*/
    ThreadList = (pthread_t *) calloc(Nthreads, sizeof(pthread_t));
    ThreadArgs = (myargument*) calloc(Nthreads, sizeof(myargument));

    for (i = 0; i < Nthreads; i++)
    {
        /* Make Thread Structure   */
        ini = (i * dims[2]) / Nthreads;
        fin = ((i + 1) * dims[2]) / Nthreads;
        ThreadArgs[i].cols = dims[0];
        ThreadArgs[i].rows = dims[1];
        ThreadArgs[i].slices = dims[2];
        ThreadArgs[i].in_image = ima;
        ThreadArgs[i].var_image = variances;
        ThreadArgs[i].means_image = means;
        ThreadArgs[i].coilsens_image = coilsens;
        ThreadArgs[i].weight = weight;
        ThreadArgs[i].estimate = Estimate;
        ThreadArgs[i].label = Label;
        ThreadArgs[i].ini = ini;
        ThreadArgs[i].fin = fin;
        ThreadArgs[i].radiusB = radiusB;
        ThreadArgs[i].radiusS = radiusS;
        ThreadArgs[i].sigma = hSigma;
    }

    for (i = 0; i < Nthreads; i++)
    {
        if (pthread_create(&ThreadList[i], NULL, ThreadFunc, &ThreadArgs[i]))
        {
            mexErrMsgIdAndTxt("myMBONLM3D:", " Threads cannot be created\n");

        }
    }

    for (i = 0; i < Nthreads; i++)
    {
        pthread_join(ThreadList[i], NULL);
    }

#endif
    mexPrintf("  myMBONLM3D: Threads finished\n");  mexEvalString("drawnow");

#ifdef _WIN32
    mxFree(ThreadArgs);
    mxFree(ThreadList);
#else
    free(ThreadArgs);
    free(ThreadList);
#endif

    label = 0.0;
    estimate = 0.0;
    /* Aggregation of the estimators (i.e. means computation) */
    for (k = 0; k < dims[2]; k++)
    {
        for (i = 0; i < dims[1]; i++)
        {
            for (j = 0; j < dims[0]; j++)
            {
                p1 = k * (dims[0] * dims[1]) + (i * dims[0]) + j;
                label = Label[p1];
                if (label == 0.0)
                    fima[p1] = ima[p1];
                else
                    fima[p1] = Estimate[p1] / label;
            }
        }
    }
    mexPrintf("myMBONLM3D: Completed.\n");
    return;

}

