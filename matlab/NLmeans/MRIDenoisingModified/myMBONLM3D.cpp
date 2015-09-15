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

typedef struct{
  int rows;
  int cols;
  int slices;
  float * in_image;   
  float * means_image;
  float * var_image;    
  float * estimate;    
  float * label;    
  int ini;
  int fin;
  int radiusB;
  int radiusS;
  float sigma; 
}myargument;

bool rician=false;

/* Function which compute the weighted average for one block */
void Average_block(float *ima,int x,int y,int z,int neighborhoodsize,float *average, float weight, int sx,int sy,int sz)
{
  int x_pos,y_pos,z_pos;
  bool is_outside; 
  int a,b,c,ns,sxy,count;

  ns=2*neighborhoodsize+1;
  sxy=sx*sy;

  count = 0;

  for (c = 0; c<ns;c++)
    {
      for (b = 0; b<ns;b++)
	{
	  for (a = 0; a<ns;a++)
	    {
	
	      is_outside = false;
	      x_pos = x+a-neighborhoodsize;
	      y_pos = y+b-neighborhoodsize;
	      z_pos = z+c-neighborhoodsize;
	
	      if ((z_pos < 0) || (z_pos > sz-1)) is_outside = true;
	      if ((y_pos < 0) || (y_pos > sy-1)) is_outside = true;
	      if ((x_pos < 0) || (x_pos > sx-1)) is_outside = true;
		
	      if(rician)
                {
		  if (is_outside)
		    average[count] = average[count] + ima[z*(sxy)+(y*sx)+x]*ima[z*(sxy)+(y*sx)+x]*weight;
		  else	
		    average[count] = average[count] + ima[z_pos*(sxy)+(y_pos*sx)+x_pos]*ima[z_pos*(sxy)+(y_pos*sx)+x_pos]*weight;
                }
	      else
                {
		  if (is_outside)
		    average[count] = average[count] + ima[z*(sxy)+(y*sx)+x]*weight;
		  else	
		    average[count] = average[count] + ima[z_pos*(sxy)+(y_pos*sx)+x_pos]*weight;
                }
                
	      count++;
	    }
	}
    }
}

/* Function which computes the value assigned to each voxel */
void Value_block(float *Estimate, float *Label,int x,int y,int z,int neighborhoodsize,float *average, float global_sum, int sx,int sy,int sz,float bias)
{
  int x_pos,y_pos,z_pos;
  int ret;
  bool is_outside;
  float value = 0.0;
  float label = 0.0;
  float denoised_value =0.0;
  int count=0 ;
  int a,b,c,ns,sxy;

  extern bool rician;

  ns=2*neighborhoodsize+1;
  sxy=sx*sy;


  for (c = 0; c<ns;c++)
    {
      for (b = 0; b<ns;b++)
	{
	  for (a = 0; a<ns;a++)
	    {		
	      is_outside = false;
	      x_pos = x+a-neighborhoodsize;
	      y_pos = y+b-neighborhoodsize;
	      z_pos = z+c-neighborhoodsize;
	
	      if ((z_pos < 0) || (z_pos > sz-1)) is_outside = true;
	      if ((y_pos < 0) || (y_pos > sy-1)) is_outside = true;
	      if ((x_pos < 0) || (x_pos > sx-1)) is_outside = true;
	      if (!is_outside)
		{
		
		  value = Estimate[z_pos*(sxy)+(y_pos*sx)+x_pos];      
                    
		  if(rician)
                    {
                      denoised_value  = (average[count]/global_sum) - bias;                    
                      if (denoised_value > 0)
			denoised_value = sqrt(denoised_value);
		      else denoised_value = 0.0;
                      value = value + denoised_value;
                    }
		  else value = value + (average[count]/global_sum);                    
                    
		  label = Label[(x_pos + y_pos*sx + z_pos*sxy)];
		  Estimate[z_pos*(sxy)+(y_pos*sx)+x_pos] = value;
		  Label[(x_pos + y_pos*sx + z_pos *sxy)] = label +1;		
		}
	      count++;
	    }
	}
    }
}


float distance(float* ima,int x,int y,int z,int nx,int ny,int nz,int f,int sx,int sy,int sz)
{

  float d,acu,distancetotal,inc;
  int i,j,k,ni1,nj1,ni2,nj2,nk1,nk2,kk,sxy;

  sxy=sx*sy;

  acu=0;
  distancetotal=0;
	
  for(k=-f;k<=f;k++)
    {
      for(j=-f;j<=f;j++)
	{
	  for(i=-f;i<=f;i++)
	    {
	      ni1=x+i;
	      nj1=y+j;
	      nk1=z+k;
	      ni2=nx+i;
	      nj2=ny+j;
	      nk2=nz+k;
			
	      if(ni1<0) ni1=-ni1;
	      if(nj1<0) nj1=-nj1;
	      if(ni2<0) ni2=-ni2;
	      if(nj2<0) nj2=-nj2;
	      if(nk1<0) nk1=-nk1;
	      if(nk2<0) nk2=-nk2;
			
	      if(ni1>=sx) ni1=2*sx-ni1-1;
	      if(nj1>=sy) nj1=2*sy-nj1-1;
	      if(nk1>=sz) nk1=2*sz-nk1-1;
	      if(ni2>=sx) ni2=2*sx-ni2-1;
	      if(nj2>=sy) nj2=2*sy-nj2-1;
	      if(nk2>=sz) nk2=2*sz-nk2-1;
			
	      distancetotal = distancetotal + ((ima[nk1*(sxy)+(nj1*sx)+ni1]-ima[nk2*(sxy)+(nj2*sx)+ni2])*(ima[nk1*(sxy)+(nj1*sx)+ni1]-ima[nk2*(sxy)+(nj2*sx)+ni2]));
	      acu=acu + 1;
	    }
	}
    }

  d=distancetotal/acu;

  return d;

}

void* ThreadFunc( void* pArguments )
{
  float bias,*Estimate,*Label,*ima,*means,*variances,*average,sigma,epsilon,mu1,var1,totalweight,wmax,t1,t2,d,w,h,hh;
  int rows,cols,slices,ini,fin,radiusB,radiusS,init,i,j,k,rc,ii,jj,kk,ni,nj,nk,Ndims;
    
  myargument arg;
  arg=*(myargument *) pArguments;

  rows=arg.rows;    
  cols=arg.cols;
  slices=arg.slices;
  ini=arg.ini;    
  fin=arg.fin;
  ima=arg.in_image;    
  means=arg.means_image;  
  variances=arg.var_image;     
  Estimate=arg.estimate;
  Label=arg.label;    
  radiusB=arg.radiusB;
  radiusS=arg.radiusS;
  sigma=arg.sigma;
               
  h=sigma;
  hh=h*h;
  bias=2*sigma*sigma;
    
  /*filter*/
  epsilon = 0.00001;
  mu1 = 0.95;
  var1 = 0.5;
  init = 0;
  rc=rows*cols;

  Ndims = (2*radiusS+1)*(2*radiusS+1)*(2*radiusS+1);

  average=(float*)malloc(Ndims*sizeof(float));

  for(k=ini;k<fin;k+=2)
    for(j=0;j<rows;j+=2)
      for(i=0;i<cols;i+=2)
	{
	  for (init=0 ; init < Ndims; init++) average[init]=0.0;
			
	  /*average=0;*/

	  totalweight=0.0;					
		
	  if ((means[k*rc+(j*cols)+i])>epsilon && (variances[k*rc+(j*cols)+i]>epsilon))
	    {
	      wmax=0.0;
				
	      for(kk=-radiusB;kk<=radiusB;kk++)
		{
		  for(jj=-radiusB;jj<=radiusB;jj++)
		    {
		      for(ii=-radiusB;ii<=radiusB;ii++)
			{
			  ni=i+ii;
			  nj=j+jj;
			  nk=k+kk;

			  if(ii==0 && jj==0 && kk==0) continue; 
				
			  if(ni>=0 && nj>=0 && nk>=0 && ni<cols && nj<rows && nk<slices)
			    {					
				
			      if ((means[nk*(rc)+(nj*cols)+ni])> epsilon && (variances[nk*rc+(nj*cols)+ni]>epsilon))
				{
				
				  t1 = (means[k*rc+(j*cols)+i])/(means[nk*rc+(nj*cols)+ni]);  
				  t2 = (variances[k*rc+(j*cols)+i])/(variances[nk*rc+(nj*cols)+ni]);
	
				  if(t1>mu1 && t1<(1/mu1) && t2>var1 && t2<(1/var1))
				    {                 
										
				      d=distance(ima,i,j,k,ni,nj,nk,radiusS,cols,rows,slices);
	
				      w = exp(-d/(hh));
	
				      if(w>wmax) wmax = w;
										
				      Average_block(ima,ni,nj,nk,radiusS,average,w,cols,rows,slices);
										
									
				      totalweight = totalweight + w;
				    }
				}
			
			
			
			    }
			}
		    }
							
		}
				
	      if(wmax==0.0) wmax=1.0;
						
	      Average_block(ima,i,j,k,radiusS,average,wmax,cols,rows,slices);
					
	      totalweight = totalweight + wmax;
					
					 
	      if(totalweight != 0.0)
		Value_block(Estimate,Label,i,j,k,radiusS,average,totalweight,cols,rows,slices,bias);
				
	    }
	  else 
            {           
              wmax=1.0;  
	      Average_block(ima,i,j,k,radiusS,average,wmax,cols,rows,slices);	
              totalweight = totalweight + wmax;
	      Value_block(Estimate,Label,i,j,k,radiusS,average,totalweight,cols,rows,slices,bias);
            }
            
	  /*Average_block(ima,i,j,k,radiusS,average,wmax,cols,rows,slices);		*/
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
  float *ima, *fima,*average,*coil;
  mxArray *Mxmeans, *Mxvariances, *MxEstimate, *MxLabel,*MxCoil;
  float *means, *variances, *Estimate, *Label;
  mxArray *pv;
  float hSigma,w,totalweight,wmax,d,mean,var,t1,t2,hh,epsilon,mu1,var1,label,estimate;
  int Ndims,i,j,k,ii,jj,kk,ni,nj,nk,radiusB,radiusS,ndim,indice,init,Nthreads,ini,fin,r;
  const int  *dims,*coildims;

  myargument *ThreadArgs;  

#ifdef _WIN32
  HANDLE *ThreadList; /* Handles to the worker threads*/
#else
  pthread_t * ThreadList;
#endif

  if(nrhs < 5 || nrhs > 6)
    {
      mexPrintf("myMBONLM3D: Not enough arguments \n");
      mexPrintf( "Usage: OutputImage=myMBONLM3D(InputImage,searcharea,patcharea,sigma,rician,[Coilsens]) ");
      return;
    }

  /*Copy input pointer x*/

  /*Get matrix x*/
  if ( mxIsSparse(prhs[0]) || 
       mxIsComplex(prhs[0]) || 
       mxIsDouble(prhs[0])  ||
       mxGetNumberOfElements(prhs[0]) == 1 ||
       mxGetNumberOfDimensions(prhs[0]) != 3) {
      mexErrMsgTxt("myMBONLM input1 must be full matrix of real float values.");
      return;
  }
  ima = (float*) mxGetPr(prhs[0]);

  ndim = mxGetNumberOfDimensions(prhs[0]);
  dims = mxGetDimensions(prhs[0]);

  /*Get the patch area, search area and sigma*/
  if ( mxIsComplex(prhs[1]) || 
       mxGetNumberOfElements(prhs[1]) != 1) {
      mexErrMsgTxt("myMBONLM patch area must be an integer.");
      return ;
  }
  radiusB = (int)(mxGetScalar(prhs[1]));

  if ( mxIsComplex(prhs[2]) || 
       mxGetNumberOfElements(prhs[2]) != 1) {
      mexErrMsgTxt("myMBONLM search area must be an integer.");
  }
  radiusS = (int)(mxGetScalar(prhs[2]));
  
  if ( mxIsSparse(prhs[3]) || 
       mxIsComplex(prhs[3]) ||
       mxGetNumberOfElements(prhs[3]) != 1) {
      mexErrMsgTxt("myMBONLM sigma must be real value.");
  }
  hSigma = (float)(mxGetScalar(prhs[3]));
  
  if ( mxIsSparse(prhs[4]) || 
       mxIsComplex(prhs[4])||
       mxGetNumberOfElements(prhs[4]) != 1 ) {
      mexErrMsgTxt("myMBONLM rician must be real value or bool.");
  }
  r = (int)(mxGetScalar(prhs[4]));
  if(r>0) rician=true;

  if (nrhs == 6){
     if ( mxIsSparse(prhs[5]) || 
       mxIsComplex(prhs[5]) || 
       mxIsDouble(prhs[5])  ||
       mxGetNumberOfElements(prhs[5]) == 1 ||
       mxGetNumberOfDimensions(prhs[5]) != 3) {
      mexErrMsgTxt("myMBONLM coil sens must be full matrix of real float values.");
      return;
     }     
     coildims = mxGetDimensions(prhs[5]);
     if (coildims[0]!=dims[0] || coildims[1]!=dims[1] || coildims[2]!=dims[2] ){
       mexErrMsgTxt("myMBONLM coil dims must equal input image dims.");
      return;
     }     
     coil = (float*)mxGetPr(prhs[5]);
  }else{
    MxCoil = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL);
    coil = (float*) mxGetPr(MxCoil);
  }


  Ndims = (int)pow((2*radiusS+1),ndim);

  /* Allocate memory and assign output pointer  and */
  /* Get a pointer to the data space in our newly allocated memory */
  plhs[0] = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL);
  fima = (float*)mxGetPr(plhs[0]);
  if (nlhs>1){
    plhs[1] = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL); 
    means = (float*)mxGetPr(plhs[1]);
  }else{
    Mxmeans = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL); 
    means = (float*)mxGetPr(Mxmeans);
  } 
  if (nlhs>2){
    plhs[2] == mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL); 
    variances = (float*)mxGetPr(plhs[2]);
  }else{
    Mxvariances = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL); 
    variances = (float*)mxGetPr(Mxvariances);
  }
  if (nlhs>3){
    plhs[3] = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL); 
    Estimate = (float*)mxGetPr(plhs[3]);
  }else{
    MxEstimate = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL); 
    Estimate = (float*)mxGetPr(MxEstimate);
  }
  if (nlhs>4){
    plhs[4] = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL);
    Label = (float*)mxGetPr(plhs[4]);
  }else{
    MxLabel = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS, mxREAL);
    Label = (float*)mxGetPr(MxLabel);
  }

  average=(float*)malloc(Ndims*sizeof(float));

  for (i = 0; i < dims[2] *dims[1] * dims[0];i++)
    {
      Estimate[i] = 0.0;
      Label[i] = 0.0;
      fima[i] = 0.0;
    }


  for(k=0;k<dims[2];k++)
    {
      for(i=0;i<dims[1];i++)
	{
	  for(j=0;j<dims[0];j++)
	    {
	      mean=0;
	      indice=0;
	      for(ii=-1;ii<=1;ii++)
		{
		  for(jj=-1;jj<=1;jj++)
		    {
		      for(kk=-1;kk<=1;kk++)
			{
			  ni=i+ii;
			  nj=j+jj;		   		  
			  nk=k+kk;
						
			  if(ni<0) ni=-ni;
			  if(nj<0) nj=-nj;
			  if(nk<0) nk=-nk;
			  if(ni>=dims[1]) ni=2*dims[1]-ni-1;
			  if(nj>=dims[0]) nj=2*dims[0]-nj-1;
			  if(nk>=dims[2]) nk=2*dims[2]-nk-1;
						
								
			  mean = mean + ima[nk*(dims[0]*dims[1])+(ni*dims[0])+nj];
			  indice=indice+1;                
					
			}
		    }
		}
	      mean=mean/indice;
	      means[k*(dims[0]*dims[1])+(i*dims[0])+j]=mean;
	    }
	}
    }

  for(k=0;k<dims[2];k++)
    {
      for(i=0;i<dims[1];i++)
	{
	  for(j=0;j<dims[0];j++)
	    {
	      var=0;
	      indice=0;
	      for(ii=-1;ii<=1;ii++)
		{
		  for(jj=-1;jj<=1;jj++)
		    {
		      for(kk=-1;kk<=1;kk++)
			{
			  ni=i+ii;
			  nj=j+jj;		   		  
			  nk=k+kk;		   		  
			  if(ni>=0 && nj>=0 && nk>0 && ni<dims[1] && nj<dims[0] && nk<dims[2])
			    {
			      var = var + (ima[nk*(dims[0]*dims[1]) + (ni*dims[0])+nj] - means[k*(dims[0]*dims[1])+(i*dims[0])+j]) * (ima[nk*(dims[0]*dims[1])+(ni*dims[0])+nj] - means[k*(dims[0]*dims[1])+(i*dims[0])+j]);
			      indice = indice+1;
			    }
			}
		    }
		}
	      var=var/(indice-1);
	      variances[k*(dims[0]*dims[1])+(i*dims[0])+j]=var;
	    }
	}
    }

  Nthreads = 8;


#ifdef _WIN32

  /* Reserve room for handles of threads in ThreadList*/
  ThreadList = (HANDLE*)malloc(Nthreads* sizeof( HANDLE ));
  ThreadArgs = (myargument*) malloc( Nthreads*sizeof(myargument));

  for (i=0; i<Nthreads; i++)
    {         
      /* Make Thread Structure   */
      ini=(i*dims[2])/Nthreads;
      fin=((i+1)*dims[2])/Nthreads;            
      ThreadArgs[i].cols=dims[0];
      ThreadArgs[i].rows=dims[1];
      ThreadArgs[i].slices=dims[2];
      ThreadArgs[i].in_image=ima;	
      ThreadArgs[i].var_image=variances;
      ThreadArgs[i].means_image=means;  
      ThreadArgs[i].estimate=Estimate;
      ThreadArgs[i].label=Label;    
      ThreadArgs[i].ini=ini;
      ThreadArgs[i].fin=fin;
      ThreadArgs[i].radiusB=radiusB;
      ThreadArgs[i].radiusS=radiusS;
      ThreadArgs[i].sigma=hSigma;    
    	
      ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &ThreadFunc, &ThreadArgs[i] , 0, NULL );
        
    }
    
  for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
  for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }
    
#else

  /* Reserve room for handles of threads in ThreadList*/
  ThreadList = (pthread_t *) calloc(Nthreads,sizeof(pthread_t));
  ThreadArgs = (myargument*) calloc( Nthreads,sizeof(myargument));

  for (i=0; i<Nthreads; i++)
    {         
      /* Make Thread Structure   */
      ini=(i*dims[2])/Nthreads;
      fin=((i+1)*dims[2])/Nthreads;            
      ThreadArgs[i].cols=dims[0];
      ThreadArgs[i].rows=dims[1];
      ThreadArgs[i].slices=dims[2];
      ThreadArgs[i].in_image=ima;	
      ThreadArgs[i].var_image=variances;
      ThreadArgs[i].means_image=means;  
      ThreadArgs[i].estimate=Estimate;
      ThreadArgs[i].label=Label;    
      ThreadArgs[i].ini=ini;
      ThreadArgs[i].fin=fin;
      ThreadArgs[i].radiusB=radiusB;
      ThreadArgs[i].radiusS=radiusS;
      ThreadArgs[i].sigma=hSigma;  
    }
    	
  for (i=0; i<Nthreads; i++)
    {
      if(pthread_create(&ThreadList[i], NULL, ThreadFunc,&ThreadArgs[i]))
	{
	  printf("Threads cannot be created\n");
	  exit(1);
	}        
    }
   
  for (i=0; i<Nthreads; i++)
    {
      pthread_join(ThreadList[i],NULL);
    }

#endif

  free(ThreadArgs); 
  free(ThreadList);


  label = 0.0;
  estimate = 0.0;
  /* Aggregation of the estimators (i.e. means computation) */
  for (k = 0; k < dims[2]; k++ )
    {
      for (i = 0; i < dims[1]; i++ )
	{
	  for (j = 0; j < dims[0]; j++ )
	    {
	      label = Label[k*(dims[0]*dims[1])+(i*dims[0])+j];
	      if (label == 0.0)
		{
		  fima[k*(dims[0]*dims[1])+(i*dims[1])+j] = ima[k*(dims[0]*dims[1])+(i*dims[0])+j];
	
		}
	      else
		{
		  estimate = Estimate[k*(dims[0]*dims[1])+(i*dims[0])+j];
		  estimate = (estimate/label);
		  fima[k*(dims[0]*dims[1])+(i*dims[0])+j]=estimate;
	
		}
	    }
	}
    }
 
  return;

}

