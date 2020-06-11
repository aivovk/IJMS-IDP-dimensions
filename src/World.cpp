#include "World.h"

// stdout a 2D BLAS array
void print_matrix(TYPE_FLOAT * A, int n)
{
  for (int i = 0 ; i < n ; i++){
    for (int j = 0 ; j < n ; j++){
      std::cout<<A[i2d(i,j,n)]<<", ";
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
}

// stdout a 1D array
void print_vector(TYPE_FLOAT * x, int n)
{
  for (int i = 0 ; i < n ; i++)
    std::cout<<x[i]<<", ";
  std::cout<<std::endl<<std::endl;
}


/**
   Calculate elapsed time (see TIMING macros)
   Code from: https://www.gnu.org/software/libc/manual/html_node/Elapsed-Time.html
*/ 
int timeval_subtract (double *result, struct timeval *x, struct timeval *y) {
    struct timeval result0;

    /* Perform the carry for the later subtraction by updating y. */
    if (x->tv_usec < y->tv_usec) {
        int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
        y->tv_usec -= 1000000 * nsec;
        y->tv_sec += nsec;
    }
    if (x->tv_usec - y->tv_usec > 1000000) {
        int nsec = (y->tv_usec - x->tv_usec) / 1000000;
        y->tv_usec += 1000000 * nsec;
        y->tv_sec -= nsec;
    }

    /* Compute the time remaining to wait.
     tv_usec is certainly positive. */
    result0.tv_sec = x->tv_sec - y->tv_sec;
    result0.tv_usec = x->tv_usec - y->tv_usec;
    *result = ((double)result0.tv_usec)/1e6 + (double)result0.tv_sec;

    /* Return 1 if result is negative. */
    return x->tv_sec < y->tv_sec;
}

void World::initializePolymer()
{
  WorldSettings::typeForceExternal = WorldSettings::EXT_NONE;

  WorldSettings::pbc = PeriodicBoundary(false,
                                        false,
                                        false,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0);
    
  std::ifstream file;
  file.open(WorldSettings::nupFile.c_str());

  if(!file.is_open())
    {
      std::cerr<<"Unable to read NUP sequence file: "<<WorldSettings::nupFile<<std::endl;
      return;
    }

  std::string line;
  getline(file, line);


  bool random = true;
  bool fixedStart = false;
  Vector3D start = WorldSettings::offset * Vector3D(1,1,1);
  Vector3D orientation = Vector3D(1,0,0);
  polymers.push_back(new Polymer(line.size(),
				 line,
				 random,
				 fixedStart,
				 start,
				 orientation,
				 &particles));
}


World::World()
{  
  stateSaved = false;
  time = 0;
  step = 0;
  checkpointNumber = 0;
  
  std::cout<<"adding polymer "<<std::endl;
  initializePolymer();
  
  noOfParticles = particles.size();
  std::cout<<"Total number of particles: "<<noOfParticles<<std::endl;
  WorldSettings::noOfMonomers = noOfParticles;

  for (int i = 0 ; i < noOfParticles ; i++)
    {
      MobilityMatrix.push_back(std::vector< Matrix3 >());
      CholeskyMatrix.push_back(std::vector< Matrix3 >());
      for (int j = 0 ; j < noOfParticles ; j++)
        {
	  MobilityMatrix[i].push_back(Matrix3(0));;
  	  CholeskyMatrix[i].push_back(Matrix3(0));
	}
    }
 
  //initialize vars for BLAS
  blas_upper = new char[1];
  blas_upper[0] = 'U';
  
  blas_one = 1;
  blas_fone = 1.0;
  blas_zero = 0;
  blas_fzero = 0.0;
  
  blas_n = 3 * noOfParticles;
  blas_mm = (double*) malloc(blas_n * blas_n *sizeof(double));
  blas_deterministic_force = (double *) malloc(blas_n *sizeof(double));
  blas_noise_force = (double *) malloc(blas_n * sizeof(double));
  blas_displacement = (double *) malloc(blas_n * sizeof(double));
  // monomer positions
  blas_r = (double *) malloc(blas_n * sizeof(double));
  // hydrodynamic radii
  blas_a = (double *) malloc(noOfParticles * sizeof(double));
  for (int i = 0 ; i < noOfParticles ; i++)
    blas_a[i] = particles[i]->getHRadius();
  
  #ifdef GPU
  //allocate memory on GPU
  mycudaMalloc((void **) &d_mm, blas_n * blas_n * sizeof(double));
  mycudaMalloc((void **) &d_deterministic_force, blas_n * sizeof(double));
  mycudaMalloc((void **) &d_noise_force,  blas_n * sizeof(double));
  mycudaMalloc((void **) &d_displacement, blas_n * sizeof(double));
  mycudaMalloc((void **) &d_r, blas_n * sizeof(double));
  mycudaMalloc((void **) &d_a, noOfParticles * sizeof(double));
  cudaMemcpy(d_a, blas_a, noOfParticles * sizeof(double), cudaMemcpyHostToDevice);
  
  //for cholesky decomposition on GPU
  mycudaMalloc((void **) &d_cholesky_info, sizeof(int));
  cusolverDnCreate(&d_solver_handle);	
  cusolverDnDpotrf_bufferSize(d_solver_handle, CUBLAS_FILL_MODE_UPPER, blas_n, d_mm, blas_n, &d_work_size);
  mycudaMalloc((void **) &d_work, d_work_size*sizeof(double));

  //for cuBLAS on gpu
  cublasCreate(&d_cublas_handle);
  #endif

  // initialize helper lists
  for (int i = 0 ; i < noOfParticles ; i++)
    {
      TYPE_FLOAT charge = particles[i]->getCharge();
      TYPE_FLOAT hydrophobicity = particles[i]->getHydrophobicity();
      if (charge != 0)
	charged.push_back(i);
      if (hydrophobicity != 0)
	hydrophobic.push_back(i);

      noise.push_back(Vector3D());

      // \todo these will be incorrect once a checkpoint is loaded
      savedPositions.push_back(particles[i]->r - WorldSettings::offsetVector);
      // previous is set in the loadStateBinary function
      previous.push_back(particles[i]->r);
    }
  
  //initialize contact map
  /*contactMap.resize(noOfParticles);
  for (int i = 0 ; i < noOfParticles ; i++)
    {
      contactMap[i].resize(noOfParticles);
      for (int j = 0 ; j < noOfParticles ; j++)
	{
	  contactMap[i][j] = 0;
	}
	}*/
  //initialize separationIJ
  noOfSepEntries = 0;
  separationIJ.resize(noOfParticles);
  separationSquaredIJ.resize(noOfParticles);
  for (int i = 0 ; i < noOfParticles ; i++){
    separationIJ[i] = 0;
    separationSquaredIJ[i] = 0;
  }

  loadStateBinary(WorldSettings::checkpointFile);
  
  if(WorldSettings::coulombStrength == 0)
    {
      charged.clear();
    }
  if(WorldSettings::hydrophobicStrength == 0)
    {
      hydrophobic.clear();
    }

  startPositionOfCentreOfMass = Vector3D();
  for ( int i = 0 ; i < noOfParticles ; i++)
    {
      startPosition.push_back(particles[i]->r);
      startPositionOfCentreOfMass += startPosition[i];
    }
  startPositionOfCentreOfMass = startPositionOfCentreOfMass/noOfParticles;
  
  
  for (unsigned int i = 0 ; i < hydrophobic.size() ; i++)
    hParticles.push_back(particles[hydrophobic[i]]);

  //for (unsigned int i = 0 ; i < charged.size() ; i++)
  //  cParticles.push_back(particles[charged[i]]);
  
  csRepulsive = CubeSpace(2 * (int) sqrt(noOfParticles),
			  &particles,
			  WorldSettings::aaProperties->maxLJRadius()*2 
			  + WorldSettings::neighbourListBuffer * WorldSettings::bondLength);
  csCohesive = CubeSpace((int) sqrt(noOfParticles),
			 &hParticles,
			 WorldSettings::aaProperties->maxLJRadius()*2 * 
			 WorldSettings::hydrophobicCutoff +
			 WorldSettings::neighbourListBuffer * WorldSettings::bondLength);
  //csCharged = CubeSpace(40, &rC, WorldSettings::bondLength * WorldSettings::debyeLength*3);

  

  csRepulsive.update();
  csCohesive.update();
  //csCharged.update();
  neighboursRepulsive = csRepulsive.getNeighbours();
  neighboursCohesive = csCohesive.getNeighbours();
  //neighboursCharged = csCharged.getNeighbours();
  csRepulsive.clear();
  csCohesive.clear();
  //csCharged.clear();
}

World::~World()
{
  free(blas_mm);
  delete [] blas_upper;
  free(blas_deterministic_force);
  free(blas_noise_force);
  free(blas_displacement);
  free(blas_r);
  free(blas_a);

  #ifdef GPU
  cudaFree(d_mm);
  cudaFree(d_deterministic_force);
  cudaFree(d_noise_force);
  cudaFree(d_displacement);
  cudaFree(d_r);
  cudaFree(d_a);
  cusolverDnDestroy(d_solver_handle);
  cublasDestroy(d_cublas_handle);
  #endif

  /// \todo memcheck
}

TYPE_FLOAT World::simulate()
{
  #ifdef TIMING
  // timing vars
  double restime;
  struct timeval  tdr0, tdr1, tdr2, tdr3;
  gettimeofday (&tdr0, NULL);
  #endif
  
  /* compute mobility matrix */
  //HI and hydrophobic forces (good solvent) turned off during eqSteps
  if (WorldSettings::typeHydrodynamics != WorldSettings::HYDRO_NONE
      && step >= WorldSettings::eqSteps){
    for (int i = 0 ; i < noOfParticles ; i++){
      for (int x = 0 ; x < 3 ; x++)
	blas_r[i*3+x] = particles[i]->r.get(x);
    }

    #ifdef GPU
    cudaMemcpy(d_r, blas_r, blas_n * sizeof(double), cudaMemcpyHostToDevice);
    int DIM = (blas_n*blas_n+blas_n)/2; // sum always even
    int blockSize = 512;
    int gridSize = DIM/blockSize + (DIM % blockSize > 0);
    cudaDeviceSynchronize();
    kernelFillMobilityMatrix<<<gridSize,blockSize>>>(d_r, d_a, d_mm, blas_n);
    #else
    // CPU version
    fillMobilityMatrix();
    #endif
    
    #ifdef TIMING
    #ifdef GPU
    cudaDeviceSynchronize();
    #endif
    gettimeofday (&tdr1, NULL);
    timeval_subtract (&restime, &tdr1, &tdr0);
    std::cout<<"mobility matrix took "<<1000*restime<<" ms"<<std::endl;
    #endif
  }

  /* generate noise terms */
  #ifdef TIMING
  gettimeofday (&tdr0, NULL);
  #endif

  for ( int i = 0 ; i < noOfParticles ; i++ )
    {
      noise[i] = noiseTerm(WorldSettings::sqrt_dt);
    }  
  #ifdef TIMING
  gettimeofday (&tdr1, NULL);
  timeval_subtract (&restime, &tdr1, &tdr0);
  std::cout<<"rng took "<<1000*restime<<" ms"<<std::endl;
  #endif

  /* compute neighbour lists */
  #ifdef TIMING
  gettimeofday (&tdr0, NULL);
  #endif
  if (step%WorldSettings::skipNeighbourUpdate == 0)
    {
      csRepulsive.update();
      neighboursRepulsive = csRepulsive.getNeighbours();
      csRepulsive.clear();

      if (step >= WorldSettings::eqSteps - WorldSettings::skipNeighbourUpdate){
	csCohesive.update();
	neighboursCohesive = csCohesive.getNeighbours();
	csCohesive.clear();
      }

      //csCharged.update();
      //neighboursCharged = csCharged.getNeighbours();
      //csCharged.clear();
    }
  #ifdef TIMING
  gettimeofday (&tdr1, NULL);
  timeval_subtract (&restime, &tdr1, &tdr0);
  std::cout<<"neighbour list build took "<<1000*restime<<" ms"<<std::endl;
  #endif

  /* external/charged forces */
  #ifdef TIMING
  gettimeofday (&tdr0, NULL);
  #endif
  Vector3D force;
  TYPE_FLOAT hStrength;
  TYPE_FLOAT avgSize;
  std::vector<int>::iterator iter;
  std::vector<int>::iterator iterH;

  //external force (not used)
  for (int i = 0; i < noOfParticles; i++)
    {
      particles[i]->dr = Vector3D();// + forceExternal(particles[i]->r, (particles[i]->getSize()+1)/2);
    }
  
  //electrostatic interactions
  for(unsigned int i = 0 ; i < charged.size() ; i++)
    {
      //the charged list is in ascending order
      for(unsigned int j = i + 1 ; j < charged.size() ; j++)
	{
	  if(charged[j] > charged[i] + 1)
	    {
	      force = forceCharge(WorldSettings::pbc.convertR(particles[charged[i]]->r
							      -particles[charged[j]]->r,
							      WorldSettings::bondLength),
				 particles[charged[i]]->getCharge() * particles[charged[j]]->getCharge());
	      particles[charged[i]]->dr += force;
	      particles[charged[j]]->dr -= force;
	    }
	}
      
    }
  #ifdef TIMING
  gettimeofday (&tdr1, NULL);
  timeval_subtract (&restime, &tdr1, &tdr0);
  std::cout<<"external/charged/misc forces took "<<1000*restime<<" ms"<<std::endl;
  #endif

  /* hydrophobic (cohesive) forces */
  #ifdef TIMING
  gettimeofday (&tdr0, NULL);
  #endif
  //equilibration happens in a good solvent with no HI
  if (step >= WorldSettings::eqSteps){
    for(unsigned int i = 0 ; i < neighboursCohesive -> size() ; i++)
      {
	for (iterH = (*neighboursCohesive)[i].begin();
	     iterH != (*neighboursCohesive)[i].end();
	     iterH++)
	  {
	    // only non-bonded monomers
	    if (hydrophobic[*iterH] > hydrophobic[i] + 1
		|| hydrophobic[*iterH] < hydrophobic[i] - 1)
	      {
		hStrength = sqrt(particles[hydrophobic[i]]->getHydrophobicity() 
				 * particles[hydrophobic[*iterH]]->getHydrophobicity());
		
		avgSize = (particles[hydrophobic[i]]->getLJRadius() 
			   + particles[hydrophobic[*iterH]]->getLJRadius())
		  / WorldSettings::bondLength;
		
		
		force = hStrength *
		  forceCohesive(WorldSettings::pbc.convertR(particles[hydrophobic[i]]->r
								 -particles[hydrophobic[*iterH]]->r,
								 WorldSettings::bondLength),
				     avgSize);
		
		particles[hydrophobic[i]]->dr += force;
		particles[hydrophobic[*iterH]]->dr -= force;
	      }	  
	  }
      }
  }
  #ifdef TIMING
  gettimeofday (&tdr1, NULL);
  timeval_subtract (&restime, &tdr1, &tdr0);
  std::cout<<"hydrophobic forces took "<<1000*restime<<" ms"<<std::endl;
  #endif
  
  /* repulsive forces */
  #ifdef TIMING
  gettimeofday (&tdr0, NULL);
  #endif

  for (int i = 0; i < noOfParticles; i++)
    {
      for (iter = (*neighboursRepulsive)[i].begin(); iter != (*neighboursRepulsive)[i].end(); iter++)
	{
	  if (*iter > i + 1 || *iter < i - 1)
	    {
	      avgSize = (particles[i]->getLJRadius() +
			 particles[*iter]->getLJRadius())/WorldSettings::bondLength;
	    }
	  else
	    {
	      avgSize = WorldSettings::LJSize/WorldSettings::bondLength;
	    }
	  //if(*iter != i) //this check is already done when building the neighbour list
	  force = forceRepulsive(WorldSettings::pbc.convertR(particles[i]->r
							       -particles[*iter]->r,
							       WorldSettings::bondLength),
				   avgSize);
	  
	  particles[i]->dr += force;
	  particles[*iter]->dr -= force;
	}
    }
  #ifdef TIMING
  gettimeofday (&tdr1, NULL);
  timeval_subtract (&restime, &tdr1, &tdr0);
  std::cout<<"lj repulsion took "<<1000*restime<<" ms"<<std::endl;
  #endif

  /* bond forces */
  #ifdef TIMING
  gettimeofday (&tdr0, NULL);
  #endif
  for (unsigned int i = 0; i < polymers.size(); i++)
    {
      polymers[i]->simulate(step*WorldSettings::dt);
    }
  #ifdef TIMING
  gettimeofday (&tdr1, NULL);
  timeval_subtract (&restime, &tdr1, &tdr0);
  std::cout<<"bond forces took "<<1000*restime<<" ms"<<std::endl;
  #endif
   
  /* update BLAS arrays for matrix vector multiplication */
  #ifdef TIMING
  gettimeofday (&tdr2, NULL);
  #endif
  for ( int i = 0 ; i < noOfParticles ; i++){
    for ( int j = 0 ; j < 3 ; j++){
      blas_deterministic_force[i*3 + j] = particles[i]->dr.get(j);
      blas_noise_force[i*3 + j] = noise[i].get(j);
    }
  }
  
  //Brownian dynamics
  double realDT = WorldSettings::dt; /// \todo unnecessary
  
  if (WorldSettings::typeHydrodynamics == WorldSettings::HYDRO_NONE
      || step < WorldSettings::eqSteps){
    for (int i = 0 ; i < blas_n ; i++)
      //TODO cublas daxpy
      blas_displacement[i] = blas_deterministic_force[i]*
	realDT/blas_a[i/3] +
	blas_noise_force[i]/sqrt(blas_a[i/3]);
    
  }else{ //Hydrodynamic interactions

    // GPU code starts here
    #ifdef GPU

    //copy memory to GPU
    cudaMemcpy(d_deterministic_force, blas_deterministic_force, blas_n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_noise_force, blas_noise_force, blas_n * sizeof(double), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();

    //  d_displacement <- dt * d_mm  * d_deterministic_force + 0 * d_displacement
    cublasDsymv(d_cublas_handle, 
		CUBLAS_FILL_MODE_UPPER,
                blas_n, 
		&realDT,
                d_mm,
		blas_n,
                d_deterministic_force, 
		1,
		&blas_fzero,
                d_displacement, 
		1);    
    
    #ifdef TIMING
    cudaDeviceSynchronize();
    gettimeofday (&tdr0, NULL);
    #endif

    // Cholesky decomposition of mobility matrix
    cusolverDnDpotrf(d_solver_handle,
		     CUBLAS_FILL_MODE_UPPER,
		     blas_n,
		     d_mm,
		     blas_n,
		     d_work,
		     d_work_size,
		     d_cholesky_info);
        
    #ifdef TIMING
    cudaDeviceSynchronize();
    gettimeofday (&tdr1, NULL);
    timeval_subtract (&restime, &tdr1, &tdr0);
    std::cout<<"cholesky decomposition took "<<1000*restime<<" ms"<<std::endl;
    #endif

    // d_noise_force <- d_mm * d_noise_force
    cublasDtrmv(d_cublas_handle,
		CUBLAS_FILL_MODE_UPPER,
		CUBLAS_OP_T,
		CUBLAS_DIAG_NON_UNIT,
		blas_n, 
		d_mm,
		blas_n,
		d_noise_force,
		1);

    // d_displacement <- 1 * d_noise_force + d_displacement
    cublasDaxpy(d_cublas_handle,
		blas_n,
		&blas_fone,
		d_noise_force,
		1,
		d_displacement,
		1);

    cudaDeviceSynchronize();

    cudaMemcpy(blas_displacement, d_displacement, blas_n * sizeof(double), cudaMemcpyDeviceToHost);

    //Don't need this? copying to CPU
    cudaDeviceSynchronize();

    // GPU CODE ENDS HERE
    #else

    /// CPU Hydrodynamic Interactions
    
    //print_vector(blas_deterministic_force, blas_n);
    //print_matrix(blas_mm, blas_n);
    
    // deterministic displacement
    // y <- alpha*A*x + beta*y
    cblas_dsymv(CblasColMajor,
                CblasUpper,
                blas_n,
                realDT, // alpha
                blas_mm,                      // A
                blas_n,
                blas_deterministic_force,     // x
                1,
                0,                       // beta
                blas_displacement, //y
                1);
    //print_vector(blas_displacement, blas_n);
    
    #ifdef TIMING
    gettimeofday (&tdr0, NULL);
    #endif
    
    fillCholeskyMatrix();
    
    #ifdef TIMING
    gettimeofday (&tdr1, NULL);
    timeval_subtract (&restime, &tdr1, &tdr0);
    std::cout<<"cholesky decomposition took "<<1000*restime<<" ms"<<std::endl;
    #endif
    // n <- B*n
    cblas_dtrmv(CblasColMajor,
                CblasUpper,
                CblasTrans,
                CblasNonUnit,
                blas_n,
                blas_mm,             // B
                blas_n,
                blas_noise_force,     // n
                1);
    //print_vector(blas_noise_force, blas_n);
    cblas_daxpy(blas_n,
                1.0,
                blas_noise_force,
                1,
                blas_displacement,
                1);
    //print_vector(blas_displacement, blas_n);
    #endif
  }
  // finish updating positions:
  // - copy to particle classes
  // - check fixed positions and max displacements
  previousStep = step;
  for (int i = 0; i < noOfParticles; i++)
    {
      if(!(particles[i]->isFixed()))
	{
	  previous[i] = particles[i]->r;
	  
	  particles[i]->r += Vector3D(blas_displacement[i*3],
				      blas_displacement[i*3+1],
				      blas_displacement[i*3+2]);
	  
	  Vector3D change = particles[i]->r - previous[i];
	  TYPE_FLOAT change_mag = change.magnitude();
	  if (change_mag > WorldSettings::MAX_FORCE || std::isnan(change_mag)||std::isinf(change_mag)){
	    std::cout<<"MAX"<<std::endl;
	    particles[i]->r = previous[i]+WorldSettings::MAX_FORCE* change/change_mag;
	    WorldSettings::errorCountLJ++;
	  }
	}
    }
  #ifdef TIMING
  gettimeofday (&tdr3, NULL);
  timeval_subtract (&restime, &tdr3, &tdr2);
  std::cout<<"updating particle positions (with cholesky) took "<<1000*restime<<" ms"<<std::endl;
  #endif

  /// \todo unnecessary?
  for (unsigned int i = 0 ; i < hydrophobic.size() ; i++)
    hParticles[i] = particles[hydrophobic[i]];
  //for (unsigned int i = 0 ; i < charged.size() ; i++)
  //  cParticles[i] = particles[charged[i]];
  
  step++;
  time += realDT;
  return realDT;
}

void World::draw() const
{
#ifdef GL
  GLUquadricObj* pQuadric = gluNewQuadric();
  // center the polymer
  Vector3D c = centreOfMass();
  glTranslatef(-c.x,
	       -c.y,
	       -c.z);

  //draw polymer bonds
  glColor3f(1,1,1);
  for(unsigned int i = 0 ; i < polymers.size(); i++)
    polymers[i]->draw();

  // draw particles
  int numSlices = 10;
  int numStacks = 10;

  for (int i = 0; i < noOfParticles; i++)
    {
      GLfloat radius = 0.5 * particles[i]->getHRadius();
      if (particles[i]->getType() == NANO_PARTICLE)
	glColor3f(1,0,1);
      else if (particles[i]->getHydrophobicity() > 0)
	glColor3f(0,1,0);
      else if (particles[i]->getCharge() > 0)
	glColor3f(1,0,0);
      else if (particles[i]->getCharge() < 0)
	glColor3f(0,0,1);
      else
	  glColor3f(1,1,1);

      glPushMatrix();
      glTranslatef(particles[i]->r.x, 
		   particles[i]->r.y, 
		   particles[i]->r.z);
      gluSphere(pQuadric,radius,numSlices,numStacks);
      glPopMatrix();
    }
  free(pQuadric);
#endif
}

/* measure methods */

void World::updateContactMap(TYPE_FLOAT realDT)
{
  /* brute force
  for(int i = 0 ; i < noOfParticles ; i++)
    for(int j = i + 1 ; j < noOfParticles ; j++)
      if(WorldSettings::pbc.convertR(particles[i]->r-particles[j]->r).magnitudeSquared() 
	 < WorldSettings::bondLength * WorldSettings::bondLength * 4 * 4)
	contactMap[i][j] = contactMap[i][j] + realDT;
  */
  std::vector<int>::iterator iter;
  //using LJ neighbour list (shorter distance)
  for (int i = 0; i < noOfParticles; i++)
    {
      for (iter = (*neighboursRepulsive)[i].begin(); iter != (*neighboursRepulsive)[i].end(); iter++)
	{
	  if(WorldSettings::pbc.convertR(particles[i]->r-particles[*iter]->r, 
					 WorldSettings::bondLength).magnitudeSquared() 
	     < WorldSettings::bondLength * WorldSettings::bondLength * 1.54 * 1.54)
	    contactMap[i][*iter] = contactMap[i][*iter] + realDT;
	}
    }
}

void World::updateSeparationIJ(TYPE_FLOAT realDT)
{
  noOfSepEntries += 1;
  for (int sep = 1 ; sep < noOfParticles ; sep++){
    separationIJ[sep] += averageSeparation(sep);
    separationSquaredIJ[sep] += averageSeparationSquared(sep);
  }
}

void World::updateMonomerDensity(TYPE_FLOAT realDT)
{
  unsigned int z; // profile will be in the z direction;
  binSize = WorldSettings::bondLength;
  for(int i = 0 ; i < noOfParticles ; i++)
    {
      z = (particles[i]->r.z - WorldSettings::offset) / binSize;
      if(z >= 0)
	{
	  if (particles[i]->getType() == AMINO_ACID)
	    {
	      if (z >= histogramMonomerDensity.size())
		{
		  histogramMonomerDensity.resize(z + 1, 0);
		}
	      histogramMonomerDensity[z] += realDT;
	    }
	  else if ( particles[i]->getType() == NANO_PARTICLE)
	    {
	      if (z >= histogramNanoparticleDensity.size())
		{
		  histogramNanoparticleDensity.resize(z + 1, 0);
		}
	      histogramNanoparticleDensity[z] += realDT;
	    }
	}
      
    }
}

void World::updateRadialMonomerDensity(TYPE_FLOAT realDT)
{
  unsigned int rad;
  binSize = WorldSettings::bondLength;
  for(int i = 0 ; i < noOfParticles ; i++)
    {
      rad = Vector3D(WorldSettings::offset-particles[i]->r.x,
		     WorldSettings::offset-particles[i]->r.y, 
		     0).magnitude()/binSize;
      if(rad >= radialDensityHistogram.size())
        {
	  radialDensityHistogram.resize(rad+1, 0);
        }
      radialDensityHistogram[rad] += realDT;
    }
}

void World::resetDisplacementSquared()
{
  startPositionOfCentreOfMass = Vector3D();
  for ( int i = 0 ; i < noOfParticles ; i++)
    {
      startPosition[i] = particles[i]->r;
      startPositionOfCentreOfMass += startPosition[i];
    }
  startPositionOfCentreOfMass = startPositionOfCentreOfMass/noOfParticles;
}
TYPE_FLOAT World::displacementSquared(int i)
{
  //TYPE_FLOAT d_sq = (particles[i]->r - startPosition[i]).magnitudeSquared();
  //
  TYPE_FLOAT d_sq = (particles[i]->r - previous[i]).magnitudeSquared();
  return d_sq;
}
TYPE_FLOAT World::displacementSquaredOfCentreOfMass()
{
  Vector3D r_cm = Vector3D();
  for ( int i = 0 ; i < noOfParticles ; i++)
    r_cm += particles[i]->r;
  r_cm = r_cm/noOfParticles;
  return (r_cm - startPositionOfCentreOfMass).magnitudeSquared();
}

void World::savePositions()
{
  for ( int i = 0 ; i < noOfParticles ; i++)
    savedPositions[i] = particles[i]->r - WorldSettings::offsetVector;
  savedTime = time;
  savedStep = step;
}
