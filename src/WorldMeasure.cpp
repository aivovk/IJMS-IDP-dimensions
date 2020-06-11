#include "World.h"


TYPE_FLOAT World::averageEndToEndDistanceSquared()
{
  TYPE_FLOAT l_sq = 0;
  for(unsigned int i = 0 ; i < polymers.size() ; i++)
    {
      l_sq += polymers[i]->endToEndDistanceSquared();
    }
  return l_sq / polymers.size();
}

/// \todo does not work for multiple polymers, yet...
/// \todo merge with Polymer::averageBondLengthSquared()
TYPE_FLOAT World::averageBondLength()
{
  TYPE_FLOAT l = 0;
  Vector3D diff;
  
  for (int i = 1 ; i < noOfParticles ; i++)
    {
      diff = particles[i]->r - particles[i-1]->r;
      l += diff.magnitude();
    }
  
  return l/(noOfParticles-1);
}

/// \todo merge with Polymer::averageBondLengthSquared()
TYPE_FLOAT World::averageBondLengthSquared()
{
  TYPE_FLOAT l = 0;
  for (int i = 1 ; i < noOfParticles ; i++)
    l += (particles[i]->r - particles[i-1]->r).magnitudeSquared();
  return l/(noOfParticles-1);
}

TYPE_FLOAT World::averageSeparation(int sep)
{
  TYPE_FLOAT l = 0;
  Vector3D diff;
  for (int i = sep ; i < noOfParticles ; i++)
    {
      diff = particles[i]->r - particles[i-sep]->r;
      l += diff.magnitude();
    }

  return l/(noOfParticles - sep);
}

TYPE_FLOAT World::averageSeparationSquared(int sep)
{
  TYPE_FLOAT l = 0;
  Vector3D diff;
  for (int i = sep ; i < noOfParticles ; i++)
    {
      diff = particles[i]->r - particles[i-sep]->r;
      l += diff.magnitudeSquared();
    }

  return l/(noOfParticles - sep);
}

Vector3D World::averageEndToEndVector()
{
  Vector3D l = Vector3D();
  for(unsigned int i = 0 ; i < polymers.size() ; i++)
    {
      l += polymers[i]->endToEndVector();
    }
  return l / polymers.size();
}
Vector3D World::centreOfMass() const
{
  Vector3D com = Vector3D();
  for(int i = 0 ; i < noOfParticles ; i++)
    {
      com += particles[i]-> r;
    }
  return com / noOfParticles;
}

Vector3D World::endToEndVector()
{
  Vector3D re = particles[noOfParticles-1]->r-particles[0]->r;
  return re;
}


TYPE_FLOAT World::rgSquared()
{
  TYPE_FLOAT rg_sq = 0;
  for(unsigned int i = 0 ; i < polymers.size() ; i++)
    {
      rg_sq += polymers[i]->radiusOfGyrationSquared();
    }
  return rg_sq / polymers.size();
}

//returns the inverse of the hydrodynamic radius (Kirkwood approx.) of the
//current conformation normalized by the size of the first bead
TYPE_FLOAT World::rkInverse()
{
  TYPE_FLOAT rki = 0;
  TYPE_FLOAT eps = 0;
  for(int i = 0 ; i < noOfParticles ; i++){
    eps += 1.0/particles[i]->getHRadius();
    for(int j = i + 1; j< noOfParticles; j++)
      {
	//times 2 since only doing 1/2 the sum
	rki += 2.0/sqrt(separationSquared(i, j));
      }
  }
  return (rki +eps) / (noOfParticles * noOfParticles);
}
/**
 * Calculates the kirkwood approximation (short time diffusion coeffcient) using
 * the RPY tensor 
 *
 *\return inverse of the hydrodynamic radius of the current conformation
 * normalized by the size of the first bead
 * 
 * this gives almost the same result as above but mobility matrix is overwritten
 * by BLAS methods, so DO NOT USE
 */
TYPE_FLOAT World::kirkwoodApproximation()
{
  //matrix is overwritten by cholesky so would either need to call this before (add a flag?)
  //rerun the kernel?
  //add a calculation kernel?
  //or do it on CPU since only every third diagonal is needed and this is not called very often
  //only differs from rkInverse when r<2*a and not significantly
  //trace of diagonal elements is same as Oseen which is what rhInverse does
  TYPE_FLOAT rhi = 0;
  //sum diagonals
  for (int i = 0 ; i < blas_n ; i++)
    rhi += blas_mm[i2d(i, i, blas_n)];
  //rest of the elements
  for (int i = 0 ; i < noOfParticles ; i++)
    for (int j = i + 1 ; j < noOfParticles ; j++)
      for (int k = 0 ; k < 3 ; k++)
	rhi += 2 * blas_mm[i2d(3*i+k, 3*j+k, blas_n)];
  rhi /= (3 * noOfParticles * noOfParticles);
  return rhi;
}

/*
Vector3D World::fixmanCorrection()
{
  Vector3D fc = Vector3D();
  if (WorldSettings::typeHydrodynamics == WorldSettings::HYDRO_NONE)
    fillMobilityMatrix(WorldSettings::RPY);
  for (int i = 0; i < noOfParticles; i++)
    {
      for (int j = 0 ; j < noOfParticles ; j++)
	{
	  fc += MobilityMatrix[i][j] * particles[j]->dr ;
	}
    }
  return fc;
}*/

/*
 * Eigenvalues of Gyration Tensor
 */
Vector3D World::asphericity()
{
  Matrix3 GyrationTensor = Matrix3(0);
  Vector3D com = centreOfMass();

  // symmetrical
  for (int i = 0 ; i < noOfParticles ; i++){
    GyrationTensor = GyrationTensor + Matrix3(particles[i]->r - com, particles[i]->r - com);
  }
  GyrationTensor = GyrationTensor * ( 1.0 / noOfParticles );

  double rgt[9];
  for (int i = 0 ; i < 9 ; i++){
    rgt[i] = GyrationTensor.get(i/3, i%3);
  }  
  
  double eigen[3];
  LAPACKE_dsyevd (LAPACK_ROW_MAJOR, //int matrix_layout, (doesn't matter if symmetrical and everything is filled)
		  'N', //char jobz,
		  'U', //char uplo,
		  3, //lapack_int n,
		  rgt, //double* a,
		  3, //lapack_int lda,
		  eigen); //double* w);
  //eigen is in ascending order
  //TYPE_FLOAT asphericity = eigen[2] - 0.5 * (eigen[1]+eigen[0]);
    
  return Vector3D(eigen[0],eigen[1],eigen[2]);
}

TYPE_FLOAT World::separationSquared(int i, int j)
{
  return (particles[i]->r - particles[j]->r).magnitudeSquared();
}
