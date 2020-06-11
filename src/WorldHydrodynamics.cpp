#include "World.h"

#ifdef GPU

void mycudaMalloc(void ** ptr, size_t size){
  std::cout<<cudaGetErrorString(cudaMalloc(ptr, size))<<std::endl;
}

__global__ void kernelFillMobilityMatrix(double * d_R, double * d_A, double * d_MM, int n)
{
  // n = 3 * noOfParticles
  
  // \todo GPU version only works with Hydrodynamic Interactions
  
  //Calculate the upper triangle of Hydrodynamic Interaction Mobility Matrix (symmetric)
  // total n(n+1) entries
  /* convert 1D thread # to 2D upper triangle matrix indices
   * 
   * initially all labels are flipped:
   * k begins in lower right of matrix and fills it from right to left
   * 0 <= k < n(n+1)/2
   * k is a triangular number and i is its root (with some offsetting by 1 required):
   * https://en.wikipedia.org/wiki/Triangular_number#Triangular_roots_and_tests_for_triangular_numbers
   *
   * 0 < i <= n (to avoid % 0) (bottom to top)
   * 0 <= j < n (right to left)
   * i and j are then flipped
   *
   */
  int k = ((n*n+n)/2)-blockIdx.x *blockDim.x -threadIdx.x-1;
  if ( k >= 0){
    int i = (sqrt((double) 8*k+1) + 1) / 2;
    int j = (k+1) % ((i*i-i)/2+1);
    i = n - i;
    j = n - j - 1;
    int pi = i/3;
    int pj = j/3;
    int x = i % 3;
    int y = j % 3;
    if (pi==pj)
      {
	if (x == y)
	  d_MM[i2d(i, j, n)] = 1.0 / d_A[pi];
	else
	  d_MM[i2d(i, j, n)] = 0.0;
      }
    else
      {
	double rijxy;
	double magsq;
	double mag;
	
	rijxy = (d_R[pi*3+x]-d_R[pj*3+x]) * (d_R[pi*3+y]-d_R[pj*3+y]);
	magsq = pow(d_R[pi*3]- d_R[pj*3], 2)
	  + pow(d_R[pi*3+1] - d_R[pj*3+1], 2)
	    + pow(d_R[pi*3+2]- d_R[pj*3+2], 2);
	mag = sqrt(magsq);
	
	//RPY
	if (mag >2*d_A[pi])
	  {
	      if (x==y)
		d_MM[i2d(i, j, n)] = (0.75/mag) * (1.0 + rijxy / magsq +(2.0*d_A[pi]*d_A[pi]/magsq)*(1.0/3.0 - rijxy / magsq));
	      else
		d_MM[i2d(i, j, n)] = 0.75/(mag*magsq) * (1.0-2.0*d_A[pi] * d_A[pi]/magsq) * rijxy;
	  }
	else
	  {
	    if (x==y)
	      d_MM[i2d(i, j, n)] = (1.0/d_A[pi]) * (1.0-9.0*mag/(32.0*d_A[pi]) + 3.0/(32.0*d_A[pi]) * rijxy / mag);
	    else
	      d_MM[i2d(i, j, n)] = 3.0 *rijxy / (32.0 * d_A[pi]*d_A[pi]*mag) ;
	  }
      }
  }
}
#endif


void World::fillMobilityMatrix()
{
  Vector3D rij;
  TYPE_FLOAT mag_rij;
  Matrix3 m, mtemp;
  //fill the upper triangle of the mobility matrix (col indexed)
  for (int i = 0 ; i < noOfParticles ; i++)
    for (int j = i ; j < noOfParticles; j++)
      {
	if (i==j)
	  {
	    for (int x = 0 ; x < 3 ; x ++){
	      for (int y = x ; y < 3 ; y++){
		if (x == y)
		  blas_mm[i2d(i*3+x, i*3+y, blas_n)] = 1.0 / particles[i]->getHRadius();
		else
		  blas_mm[i2d(i*3+x, i*3+y, blas_n)] = 0.0;
	      }
	    }
	  }
	else
	  {
	    rij = particles[i]->r - particles[j]->r;
	    mag_rij = rij.magnitude();
	    TYPE_FLOAT ai = particles[i]->getHRadius();
	    TYPE_FLOAT aj = particles[j]->getHRadius();
	    TYPE_FLOAT amax = std::max(ai, aj);
	    TYPE_FLOAT amin = std::min(ai, aj);
	    
	    /*
	    //Oseen
	    m = Matrix3(rij/mag_rij, rij/mag_rij);
	    m = m + Matrix3(1);
	    MobilityMatrix[i][j] = m * (0.75/mag_rij);
	    */
	    //RPY
	    if (mag_rij >ai+aj)
	      {
		m = Matrix3(rij/mag_rij, rij/mag_rij) * (1-(ai*ai+aj*aj)/(mag_rij*mag_rij));
	    	mtemp = Matrix3(1.0)*(1+(ai*ai+aj*aj)/(3*mag_rij*mag_rij));
		m = m + mtemp;
		for ( int x = 0 ; x < 3 ; x++ )
		  for ( int y = 0 ; y < 3 ; y++ )
		    blas_mm[i2d(i*3+x, j*3+y, blas_n)] = m.get(x, y) * (0.75/mag_rij);
	      }
	    else if (mag_rij < amax - amin)
	      {
		m = Matrix3(1.0/amax);
		for ( int x = 0 ; x < 3 ; x++ )
		  for ( int y = 0 ; y < 3 ; y++ )
		    blas_mm[i2d(i*3+x, j*3+y, blas_n)] = m.get(x, y);

	      }
	    else
	      {
		m = Matrix3(rij/mag_rij, rij/mag_rij)
		  * ((3*pow(pow(ai-aj, 2)-pow(mag_rij, 2), 2))/(32*pow(mag_rij,3)));
		mtemp = Matrix3(1.0)
		  *((16*pow(mag_rij,3)*(ai+aj)-pow(pow(ai-aj,2)+3*pow(mag_rij,2),2))/(32*pow(mag_rij, 3)));

		//		m= Matrix3((1.0-9.0*mag_rij/(32.0*a)));
		//mtemp = Matrix3(rij/mag_rij, rij/mag_rij);
		//mtemp = mtemp * (3.0*mag_rij/(32.0*a));
		m = m + mtemp;
		for ( int x = 0 ; x < 3 ; x++ )
		  for ( int y = 0 ; y < 3 ; y++ )
		    blas_mm[i2d(i*3+x, j*3+y, blas_n)] = m.get(x, y) * (1.0/(ai*aj));
	      }
	  }
      }  
}

void World::fillCholeskyMatrix()
{
  //convention seems to be different a[i,j] =  a[j*a_dim1 + i]
  dpotrf(blas_upper,
	 &blas_n,
	 blas_mm,
	 &blas_n,
	 &blas_cholesky_info);
}
