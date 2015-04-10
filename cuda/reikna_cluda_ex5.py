import numpy as np
from numpy.linalg import norm

from reikna import cluda
from reikna.cluda import functions, dtypes

N = 128
dtype = np.complex64
ftype = np.float32
api = cluda.ocl_api()
thr = api.Thread.create()

sigma=np.ones(3) 
FACTOR = (((np.sqrt(2*np.pi)**3)*np.prod(sigma)))

program = thr.compile("""
KERNEL void gauss_kernel(
    GLOBAL_MEM ${ctype} *dest,
    GLOBAL_MEM ${ctype} *src)
{
  const ulong x = get_global_id(0); //const SIZE_T size1 = get_local_size(0);
  const SIZE_T y = get_local_id(0); //const SIZE_T size1 = get_local_size(1);
  const SIZE_T z = get_group_id(0); //const SIZE_T size2 = get_num_groups(0);
  //const SIZE_T k = get_global_id(0);const SIZE_T dim3 = get_global_size(0);
  const SIZE_T dim1= %d;
  const SIZE_T dim2= %d;
  const SIZE_T dim3= %d;                    
  ${ftype} sigma[3];
  sigma[0]=%f;sigma[1]=%f;sigma[2]=%f;
  ${ftype} factor = %f;            
  const double TWOPISQ = 19.739208802178716; //6.283185307179586;  //2*3.141592;
  const ${ftype} SQRT2PI = 2.5066282746;
  const double CUBEDSQRT2PI = 15.749609945722419;
  const ulong idx = x; //z * dim2 * dim1  + y * dim1  + x;
                      //(SIZE_T)floor((${ftype})(z * dim2 * dim1  + y * dim1  + x )/(${ftype})dim1*dim2*dim3);
  //${ftype} i = ((${ftype})(x) - floor((${ftype})(dim1)/2.0))/(${ftype})(dim1);
  //${ftype} j = ((${ftype})(y) - floor((${ftype})(dim2)/2.0))/(${ftype})(dim2);
  //${ftype} k = ((${ftype})(z) - floor((${ftype})(dim3)/2.0))/(${ftype})(dim3);
  ${ftype} i = (${ftype})((x / dim3) / dim2);
      i = (i - (${ftype})floor((${ftype})(dim1)/2.0))/(${ftype})(dim1);
  ${ftype} j = (${ftype})(x / dim3);
      if((SIZE_T)j > dim2) {j=(${ftype})fmod(j, (${ftype})dim2);};
      j = (j - (${ftype})floor((${ftype})(dim2)/2.0f))/(${ftype})(dim2);
  //Account for large global index (stored as ulong) before performing modulus
  double pre_k=fmod((double)(x) , (double) dim3);
  ${ftype} k = (${ftype}) pre_k;
      k = (k - (${ftype})floor((${ftype})(dim3)/2.0f))/(${ftype})(dim3);

  ${ftype} weight = exp(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  //${ftype} weight = expm1(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]))+1;
  //${ftype} weight= ${exp}(-TWOPISQ*((i*i)*sigma[0]*sigma[0] + (j*j)*sigma[1]*sigma[1] + (k*k)*sigma[2]*sigma[2]));
  
  dest[idx].x = src[idx].x * weight;
  dest[idx].y = src[idx].y * weight;  //(${ftype})k; //
  
}
""" % (N,N,N,sigma[0],sigma[1],sigma[2],FACTOR), render_kwds=dict(ctype=dtypes.ctype(dtype),
                                ftype=dtypes.ctype(ftype),
    exp=functions.exp(ftype)),fast_math=True)

gauss_kernel = program.gauss_kernel

r1 = np.ones((N,N,N)).astype(ftype)   #/N
r2 = np.ones((N,N,N)).astype(ftype)   #/N
a = r1 + 1j * r2
b = r1 - 1j * r2
a_dev = thr.to_device(a)
#b_dev = thr.to_device(b)
#c_dev= thr.to_device(b.ravel())
#sigma_dev = thr.to_device(sigma)
dest_dev = thr.empty_like(a_dev)


# (np.pi).astype(np.float32),
#gauss_kernel(dest_dev, a_dev, sigma_dev, local_size=N,global_size=N*N*N)
gauss_kernel(dest_dev, a_dev, local_size=N,global_size=N*N*N)
dest_out = dest_dev.get()
print dest_out
#max_dest = 


#dest_numpy = (1.0/np.sqrt(2*np.pi))*np.exp(-2*np.pi*(a*a + b*b))
#print(norm(np.reshape(dest_dev.get(),(N,N-1,N+2)) - dest ) / norm(dest) <= 1e-6)
