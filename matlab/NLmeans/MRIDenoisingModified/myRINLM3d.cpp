/********************
Modified NL means filter for B1 correction

Michael Eager
Zhaolin Chen
Monash Biomedical Imaging

Monash University, 2015
 *******************/


/***************************************************************************                  
/* Jose V. Manjon - jmanjon@fis.upv.es                                     */
/* Universidad Politecnica de Valencia, Spain                              */
/* Pierrick Coupe - pierrick.coupe@gmail.com                               */
/* Brain Imaging Center, Montreal Neurological Institute.                  */
/* Mc Gill University                                                      */
/*                                                                         */
/* Copyright (C) 2010 Jose V. Manjon and Pierrick Coupe                    */

/***************************************************************************
 * New Methods for MRI Denoising based on Sparseness and Self-Similarity    *
 * Jos� V. Manj�n, Pierrick Coup�, Antonio Buades,D. Louis Collins 
 * and Montserrat Robles
 ***************************************************************************/

#include "math.h"
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

/* Multithreading stuff*/
#ifdef _WIN32
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#endif

typedef struct{
  int rows;
  int cols;
  int slices;
  float * in_image;
  float * out_image;
  float * means_image;
  float * ref_image;
  float * pesos;
  float * coilsens_image;
  int ini;
  int fin;
  int radioB;
  int radioS;
  float sigma; 
  int order;
} myargument;

bool rician;

void* ThreadFunc( void* pArguments )
{
  float *ima,*fima,*medias,*pesos,*ref,*coilsens, sigma,w,d,hhh,hh,t1,alfa,x;
  float betaCoil=0.0;
  int ii,jj,kk,ni,nj,nk,i,j,k,ini,fin,rows,cols,slices,p,p1,radioS,radioB,order,rc;    
  extern bool rician;
  /*extern float *table;*/
    
  myargument arg;
  arg=*(myargument *) pArguments;

  rows = arg.rows;    
  cols = arg.cols;
  slices = arg.slices;
  ini = arg.ini;    
  fin = arg.fin;
  ima = arg.in_image;
  fima = arg.out_image;
  medias = arg.means_image;
  ref = arg.ref_image;
  coilsens = arg.coilsens_image;
  pesos = arg.pesos;
  radioB = arg.radioB;
  radioS = arg.radioS;
  sigma = arg.sigma;
  order = arg.order;
        
  hh=2*sigma*sigma;
  alfa=0.5;
  hhh=0.5*sigma*sigma;
  rc=rows*cols;
    
  /* filter*/
        
  if(order==1)
    {        
      for(k=ini;k<fin;k++)
	{
	  for(j=0;j<rows;j++)
	    {
	      for(i=0;i<cols;i++)
		{		
		  p=k*rc+j*cols+i;                      
            
		  if(ima[p]==0) continue;
                
		  for(kk=0;kk<=radioB;kk++)
		    {
		      nk=k+kk;    
		      for(jj=-radioB;jj<=radioB;jj++)
			{
			  nj=j+jj;
			  for(ii=-radioB;ii<=radioB;ii++)
			    {
			      ni=i+ii;								                
			      if(kk==0 && jj<0) continue;             
			      if(kk==0 && jj==0 && ii<=0) continue;  
				
			      if(ni>=0 && nj>=0 && nk>=0 && ni<cols && nj<rows && nk<slices)
				{
				  p1=nk*rc+nj*cols+ni;                        

				  t1 = abs(medias[p]-medias[p1]);  
              
				  if(t1>sigma) continue;

				  w = (ref[p]-ref[p1]);
                     
				  w=(w*w)/hh-1;                     
				  d = (t1*t1)/hhh-1;
              
				  w=w<0?0:w;
				  d=d<0?0:d;
              
				  d=w*0.5 + d*0.5;

				  //Coil sensitivity (B1) correction
				  /*				  
				     betaCoil = (coilsens[p] - coilsens[p1]);
				   */
				  if (betaCoil>0) {
				    d *= (betaCoil*betaCoil);
				  } 
				  
				  if(d<=0) {
				    w=1.0;
				  } else if(d>10) {
				    w=0;
				  } else {
				    w = exp(-d);
				  }
                                                              
				  if(rician>0)
				    {
				      fima[p] += w*ima[p1]*ima[p1];
				      pesos[p]+= w;
                                         
				      fima[p1] += w*ima[p]*ima[p];
				      pesos[p1]+= w;                         
				    }
				  else
				    {
				      fima[p] += w*ima[p1];
				      pesos[p]+= w;
                                         
				      fima[p1] += w*ima[p];
				      pesos[p1]+= w;                         
				    }

				}		 
			    }
			}							
		    }                        			
		}
	    }
	}
    }
  else
    {
      for(k=ini;k<fin;k++)
	{
	  for(j=rows-1;j>=0;j--)
	    {
	      for(i=cols-1;i>=0;i--)
		{		
		  p=k*rc+j*cols+i;                      
		  if(ima[p]==0) continue;
            
		  for(kk=0;kk<=radioB;kk++)
		    {
		      nk=k+kk;
		      for(jj=-radioB;jj<=radioB;jj++)
			{
			  nj=j+jj;			
			  for(ii=-radioB;ii<=radioB;ii++)
			    {
			      ni=i+ii;				                
			      if(kk==0 && jj<0) continue;             
			      if(kk==0 && jj==0 && ii<=0) continue;  
				
			      if(ni>=0 && nj>=0 && nk>=0 && ni<cols && nj<rows && nk<slices)
				{
				  p1=nk*rc+nj*cols+ni;                        

				  t1 = abs(medias[p]-medias[p1]);  
              
				  if(t1>sigma) continue;

				  w = (ref[p] - ref[p1]);
                     
				  w=(w*w)/hh-1;                     
				  d = (t1*t1)/hhh-1;
              
				  w=w<0?0:w;
				  d=d<0?0:d;
              
				  d=w*0.5 + d*0.5;                                                                 
				  //Coil sensitivity (B1) correction
				  //mexPrintf("Beta %f",coilsens[p]);
				  
				  /*
				    betaCoil = (coilsens[p]-coilsens[p1]);
				   */
				      if (betaCoil!=0) {
				      d*= (betaCoil*betaCoil);
				      }
				  
				  if(d<=0) {
				    w=1.0;
				  }else if(d>10) {
				    w=0;
				  }else {
				    w = exp(-d);
				  }
				  if(rician>0)
				    {
				      fima[p] += w*ima[p1]*ima[p1];
				      pesos[p]+= w;
                                         
				      fima[p1] += w*ima[p]*ima[p];
				      pesos[p1]+= w;                         
				    }
				  else
				    {
				      fima[p] += w*ima[p1];
				      pesos[p]+= w;
                                         
				      fima[p1] += w*ima[p];
				      pesos[p1]+= w;                         
				    }
				}		 
			    }
			}							
		    }                        			
		}
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
  mxArray *xData,*pv,*Mxmedias,*Mxpesos,*xref,*MxCoil;
  float *ima, *fima,*lista,*pesos,*ref,*medias,*coilsens;
  float sigma,vr,R,h,w,average,totalweight,wmax,d,media,var,th,hh,hhh,wm,t1,t2,alfa;
  int Ndims,inc,i,j,k,ii,jj,kk,ni,nj,nk,radioB,radioS,ndim,indice,Nthreads,ini,fin,order;
  const int  *dims,*coildims;

  myargument *ThreadArgs;  

#ifdef _WIN32
  HANDLE *ThreadList; /* Handles to the worker threads*/
#else
  pthread_t * ThreadList;
#endif
  mexPrintf("myRINLM: Start mex function \n"); 
  mexEvalString("drawnow");
  if(nrhs != 7)
    {
      mexPrintf("myRINLM3D: Not enough arguments \n");
      mexPrintf( "Usage: OutputImage=myRINLM3D(InputImage,searcharea,patcharea,sigma,ODCTImage,rician,B1CoilSensitivityImage) "); 
      mexEvalString("drawnow");
      return;
    }

  /*Copy input pointer x*/
   mexPrintf("myRINLM: Get input image \n"); 
   mexEvalString("drawnow");
  /*Get matrix x*/
  if ( mxIsSparse(prhs[0]) || 
       mxIsComplex(prhs[0]) || 
       mxIsDouble(prhs[0])  ||
       mxGetNumberOfElements(prhs[0]) == 1 ||
       mxGetNumberOfDimensions(prhs[0]) != 3) 
    { 
      mexErrMsgTxt("myRINLM3d input1 must be full matrix of real float values."); 
      mexEvalString("drawnow");
      return;
    }
  ima = (float*)mxGetPr(prhs[0]);
  ndim = mxGetNumberOfDimensions(prhs[0]);
  dims= mxGetDimensions(prhs[0]);
  mexPrintf("myRINLM: Get params \n"); 
  mexEvalString("drawnow");
  /*Copy input parameters*/
  /*Get the search area*/
  if ( mxIsComplex(prhs[1]) || 
       mxGetNumberOfElements(prhs[1]) != 1) {
      mexErrMsgTxt("myRINLM3d search area must be an integer."); 
      mexEvalString("drawnow");
      return ;
  }
  radioS = (int)(mxGetScalar(prhs[1]));

  if ( mxIsComplex(prhs[2]) || 
       mxGetNumberOfElements(prhs[2]) != 1) {
      mexErrMsgTxt("myRINLM3d patch area must be an integer."); 
      mexEvalString("drawnow");
  }
  radioB = (int)(mxGetScalar(prhs[2]));

  if ( mxIsSparse(prhs[3]) || 
       mxIsComplex(prhs[3]) ||
       mxGetNumberOfElements(prhs[3]) != 1) {
      mexErrMsgTxt("myRINLM3d sigma must be real value."); 
      mexEvalString("drawnow");
  }
  sigma = (float)(mxGetScalar(prhs[3]));

  if ( mxIsSparse(prhs[4]) || 
       mxIsComplex(prhs[4]) || 
       mxIsDouble(prhs[4])  ||
       mxGetNumberOfElements(prhs[4]) == 1 ||
       mxGetNumberOfDimensions(prhs[4]) != 3) {
      mexErrMsgTxt("myRINLM3d ref image must be full matrix of real float values.");
      return;
  }
  ref = (float*)mxGetPr(prhs[4]);

  if ( mxIsSparse(prhs[5]) || 
       mxIsComplex(prhs[5])||
       mxGetNumberOfElements(prhs[5]) != 1 ) {
      mexErrMsgTxt("myRINLM3d rician must be real value or bool."); 
      mexEvalString("drawnow");
  }
  rician = (int)(mxGetScalar(prhs[5]));
    
  h = sigma/3.0f; /* this is due to the removal of part of the noise*/
  mexPrintf("myRINLM: H %f\n",h); 
  mexEvalString("drawnow");
  if (nrhs < 6)
    {
      mexPrintf("myRINLM: Creating neutral B1 correction \n"); 
      mexEvalString("drawnow");

      pv = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL);
      coilsens = (float*)mxGetPr(pv);
      for(k=0;k<dims[2]*dims[1]*dims[0];k++)
	coilsens[k]=1.0;
    }
  else 
    {   
      if ( mxIsSparse(prhs[6]) || 
	   mxIsComplex(prhs[6]) || 
	   mxIsDouble(prhs[6])  ||
	   mxGetNumberOfElements(prhs[6]) == 1 ||
	   mxGetNumberOfDimensions(prhs[6]) != 3) {
	mexErrMsgTxt("myRINLM3d coil sens must be full matrix of real float values.");
	return;
      }     
      coildims = mxGetDimensions(prhs[6]);  
      if (coildims[0]!=dims[0] || coildims[1]!=dims[1] || coildims[2]!=dims[2] )
	{
	  mexErrMsgTxt("myMBONLM coil dims must equal input image dims.");
	  return;
	}
      else
	{ 
	  mexPrintf("myRINLM3d B1 coil sensitivity image correct dimensions.");
	}
      coilsens = (float*)mxGetPr(prhs[6]);
      for(k=0;k<dims[2]*dims[1]*dims[0];k=k+dims[1]*dims[0])
	mexPrintf("%f\t",coilsens[k]); 
    }
  mexPrintf("myRINLM: Creating output matricies \n"); 
  mexEvalString("drawnow");

  /* Allocate memory and assign output pointer */
  /* Get a pointer to the data space in our newly allocated memory */
  plhs[0] = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL);
  fima = (float*) mxGetPr(plhs[0]);
  if (nlhs > 1)
    {
      plhs[1] = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL); 
      medias = (float*) mxGetPr(plhs[1]);
    }
  else
    {
      Mxmedias = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL); 
      medias = (float*) mxGetPr(Mxmedias);
    }
  if (nlhs > 2)
    {
      plhs[2]  = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL); 
      pesos = (float*) mxGetPr(plhs[2]);
    }
  else
    {
      Mxpesos = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL); 
      pesos = (float*) mxGetPr(Mxpesos);
    }
  mexPrintf("myRINLM: Calculating means \n"); 
  mexEvalString("drawnow");

  /* calculate means*/

  for(k=0;k<dims[2];k++)
    {
      for(j=0;j<dims[1];j++)
	{
	  for(i=0;i<dims[0];i++)
	    {
	      media=ref[k*(dims[0]*dims[1])+(j*dims[0])+i];
	      for(ii=-1;ii<=1;ii++)
		{
		  ni=i+ii;
		  if(ni<0) ni=-ni;
		  if(ni>=dims[0]) ni=2*dims[0]-ni-1;              
		  for(jj=-1;jj<=1;jj++)
		    {
		      nj=j+jj;
		      if(nj<0) nj=-nj;
		      if(nj>=dims[1]) nj=2*dims[1]-nj-1;
		      for(kk=-1;kk<=1;kk++)
			{
			  nk=k+kk;               
			  if(nk<0) nk=-nk;
			  if(nk>=dims[2]) nk=2*dims[2]-nk-1;
              
			  if(sqrt(ii*ii+jj*jj+kk*kk)>1)continue;
                          
			  media = media + ref[nk*(dims[0]*dims[1])+(nj*dims[0])+ni];                                       
			}
		    }
		}
	      media=media/8;
	      medias[k*(dims[0]*dims[1])+(j*dims[0])+i]=media;
	    }
	}
    }

  for(k=0;k<dims[2]*dims[1]*dims[0];k++)
    {
      if(rician>0) fima[k]=ima[k]*ima[k];      
      else fima[k]=ima[k];          
      pesos[k]=1.0;
    }

  Nthreads = dims[2]<8?dims[2]:8;
  if ( Nthreads < 1 ) Nthreads = 1;
#ifdef _WIN32

  /* Reserve room for handles of threads in ThreadList*/
  ThreadList = (HANDLE*)malloc(Nthreads* sizeof( HANDLE ));
  ThreadArgs = (myargument*) malloc( Nthreads*sizeof(myargument));

  order=-1;
  for (i=0; i<Nthreads; i++)
    {         
      /* Make Thread Structure*/
      order=-order;
      ini = (i*dims[2])/Nthreads;
      fin = ((i+1)*dims[2])/Nthreads;            
      ThreadArgs[i].cols = dims[0];
      ThreadArgs[i].rows = dims[1];
      ThreadArgs[i].slices = dims[2];
      ThreadArgs[i].in_image = ima;
      ThreadArgs[i].out_image = fima;
      ThreadArgs[i].ref_image = ref;
      ThreadArgs[i].means_image = medias;
      ThreadArgs[i].pesos = pesos;
      ThreadArgs[i].coilsens_image = coilsens;
      ThreadArgs[i].ini=ini;
      ThreadArgs[i].fin=fin;
      ThreadArgs[i].radioB=radioB;
      ThreadArgs[i].radioS=radioS;
      ThreadArgs[i].sigma=h; 
      ThreadArgs[i].order=order; 
    	
      ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &ThreadFunc, &ThreadArgs[i] , 0, NULL );
        
    }
    
  for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
  for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }

#else

  mexPrintf("myRINLM3D: Pthreading %d\n",Nthreads); 
  mexEvalString("drawnow;");
  /* Reserve room for handles of threads in ThreadList*/
  ThreadList = (pthread_t *) calloc(Nthreads,sizeof(pthread_t));
  ThreadArgs = (myargument*) calloc( Nthreads,sizeof(myargument));

  order=-1;
  for (i=0; i<Nthreads; i++)
    {         
      /* Make Thread Structure*/
      order=-order;
      ini=(i*dims[2])/Nthreads;
      fin=((i+1)*dims[2])/Nthreads;            
      ThreadArgs[i].cols=dims[0];
      ThreadArgs[i].rows=dims[1];
      ThreadArgs[i].slices=dims[2];
      ThreadArgs[i].in_image=ima;
      ThreadArgs[i].out_image=fima;
      ThreadArgs[i].ref_image=ref;
      ThreadArgs[i].means_image=medias;
      ThreadArgs[i].pesos=pesos;
      ThreadArgs[i].ini=ini;
      ThreadArgs[i].fin=fin;
      ThreadArgs[i].radioB=radioB;
      ThreadArgs[i].radioS=radioS;
      ThreadArgs[i].sigma=h; 
      ThreadArgs[i].order=order; 
    }
  mexPrintf("myRINLM: Starting threads \n"); 
  mexEvalString("drawnow");
  for (i=0; i<Nthreads; i++)
    {
      if(pthread_create(&ThreadList[i], NULL, ThreadFunc,&ThreadArgs[i]))
	{
	  printf("Threads cannot be created\n"); 
	  mexEvalString("drawnow");
	  return;
	}        
    }
   
  for (i=0; i<Nthreads; i++)
    {
      pthread_join(ThreadList[i],NULL);
    }
#endif
  mexPrintf("myRINLM: Threads completed \n"); 
  mexEvalString("drawnow");

  free(ThreadArgs); 
  free(ThreadList);

  sigma=2*sigma*sigma;
  for(k=0;k<dims[2]*dims[1]*dims[0];k++) 
    {
      if(rician>0)
	{
	  vr=(fima[k]/pesos[k])-sigma;
	  if(vr<0) fima[k]=0;
	  else fima[k]=sqrt(vr);
	  if(ima[k]<=0) fima[k]=0;
	}
      else fima[k]/=pesos[k];
    }
  mexPrintf("myRINLM3D: done \n"); 
  mexEvalString("drawnow");
  return;

}

