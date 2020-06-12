#ifndef WORLD_H
#define WORLD_H

#include <mkl.h>
#include <sys/time.h>
#include <time.h>

#ifdef GPU
#include <cuda.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cuda_runtime_api.h>
#endif

#include "Polymer.h"
#include "Vector3D.h"
#include "Matrix3.h"
#include "CubeSpace.h"
#include "NormalDistribution.h"
#include "PeriodicBoundary.h"
#include "WorldSettings.h"

// indexing 2D BLAS arrays
#define i2d(i,j,ld) (((j)*(ld))+(i))

int timeval_subtract (double *result, struct timeval *x, struct timeval *y);

class World
{
 public:
  World();
  ~World();

  void draw() const;
  TYPE_FLOAT simulate();

  //measure
  TYPE_FLOAT averageBondLength(); ///< \todo average over polymers
  TYPE_FLOAT averageBondLengthSquared(); ///< \todo average over polymers
  
  TYPE_FLOAT averageEndToEndDistanceSquared(); ///< averaged over polymers
  Vector3D endToEndVector(); ///< all particles
  Vector3D averageEndToEndVector(); ///< averaged over polymers

  TYPE_FLOAT rgSquared(); ///< averaged over polymers \todo separate versions
  Vector3D asphericity(); ///< eigenvalues of gyration tensor of all particles

  TYPE_FLOAT rkInverse(); ///< all particles
  TYPE_FLOAT kirkwoodApproximation(); ///< DO NOT USE
  //Vector3D fixmanCorrection();
   
  Vector3D centreOfMass() const; ///< all particles
  
  void resetDisplacementSquared();
  TYPE_FLOAT displacementSquared(int i); ///< of one monomer in last step
  TYPE_FLOAT displacementSquaredOfCentreOfMass(); ///< since reset called
  
  TYPE_FLOAT averageSeparation(int sep);
  TYPE_FLOAT averageSeparationSquared(int sep);
  TYPE_FLOAT separationSquared(int i, int j);
  void updateContactMap(TYPE_FLOAT realDT);
  void updateSeparationIJ(TYPE_FLOAT realDT);
  void updateMonomerDensity(TYPE_FLOAT realDT);
  void updateRadialMonomerDensity(TYPE_FLOAT realDT);
  std::vector< std::vector<TYPE_FLOAT> > * getContactMap(){return &contactMap;};
  std::vector<TYPE_FLOAT> * getSeparationIJ(){return &separationIJ;};
  std::vector<TYPE_FLOAT> * getSeparationSquaredIJ(){return &separationIJ;};
  std::vector<TYPE_FLOAT> * getMonomerDensity(){return &histogramMonomerDensity;};
  std::vector<TYPE_FLOAT> * getNanoparticleDensity(){return &histogramNanoparticleDensity;};
  std::vector<TYPE_FLOAT> * getRadialMonomerDensity(){return &radialDensityHistogram;};
  
  void savePositions();
  void writeContactMap(std::string filename);
  void writeSeparationIJ(std::string filename);
  void writeState(std::string filename);
  void writeStateBinary(std::string filename);
  void writeStateBinaryFast(std::string filename);
  void loadState(std::string filename);
  void loadStateBinary(std::string filename);
  bool isStateSaved(){return stateSaved;};

  int getStep() {return step;};
  TYPE_FLOAT getTime() {return time;};
  int getNoOfParticles() {return noOfParticles;};
 protected:
 private:
  void initializePolymer();
  void fillMobilityMatrix();
  void fillCholeskyMatrix();
  
  // Cell Lists for repulsive and cohesive forces
  CubeSpace csRepulsive;
  CubeSpace csCohesive;
  //CubeSpace csCharged;

  std::vector<Polymer *> polymers;
  std::vector< std::set< int > > * neighboursRepulsive;
  std::vector< std::set< int > > * neighboursCohesive;
  //std::vector< std::set< int > > * neighboursCharged;

  //Hydrodynamic Interaction Matrices
  std::vector< std::vector< Matrix3 > > MobilityMatrix;
  std::vector< std::vector< Matrix3 > > CholeskyMatrix;

  //pointers for BLAS
  
  // do not need to initialize to 0 since only using upper triangle not reading
  // whole thing
  double * blas_mm; // will store both mobility matrix and its Choleksy
		    // decomposition
  double * blas_deterministic_force;
  double * blas_noise_force;
  double * blas_displacement;
  double * blas_r;
  double * blas_a;

  #ifdef GPU
  double * d_mm;
  double * d_deterministic_force;
  double * d_noise_force;
  double * d_displacement;
  double * d_r;
  double * d_a;
  cusolverDnHandle_t d_solver_handle;
  int d_work_size;
  int * d_cholesky_info;  
  double * d_work;
  cublasHandle_t d_cublas_handle;
  #endif
  
  char * blas_upper;
  int blas_n;
  int blas_cholesky_info;
  int blas_one;
  double blas_fone;
  int blas_zero;
  double blas_fzero;
  
  int noOfParticles;
  std::vector<Particle *> particles; ///< pointers to all particles
  std::vector<Vector3D> previous; ///< positions of particles in the previous step
  std::vector<Vector3D> noise; ///< random forces for the current step

  //pointers to only hydrophobic particles
  std::vector<Particle *> hParticles;
  //std::vector<Particle *> cParticles;
  
  //separate lists of charged and hydrophobic particle indices
  std::vector<int> charged;
  std::vector<int> hydrophobic;

  //contact map
  std::vector< std::vector<TYPE_FLOAT> > contactMap;

  //inter-particle distances
  int noOfSepEntries;
  std::vector<TYPE_FLOAT> separationIJ;
  std::vector<TYPE_FLOAT> separationSquaredIJ;

  //density
  std::vector<TYPE_FLOAT> histogramMonomerDensity;
  std::vector<TYPE_FLOAT> histogramNanoparticleDensity;
  std::vector<TYPE_FLOAT> radialDensityHistogram;
  TYPE_FLOAT binSize;

  //useful \todo explain
  std::vector<Vector3D> startPosition;
  Vector3D startPositionOfCentreOfMass;
  int step;
  int savedStep;
  int previousStep;
  TYPE_FLOAT time;
  TYPE_FLOAT savedTime;
  int checkpointNumber;
       
  bool stateSaved;
  std::vector<Vector3D> savedPositions;

};

#endif // WORLD_H
