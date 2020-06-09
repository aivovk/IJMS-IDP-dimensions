#include "NormalDistribution.h"

NormalDistribution::NormalDistribution()
{
  // seed with epoch time
  // causes problems if two runs are started at the same second
  init(time(NULL));
}

NormalDistribution::NormalDistribution(unsigned int seed)
{
  init(seed);
}

void NormalDistribution::init(unsigned int seed)
{
  const gsl_rng_type * typeRNG;
  
  gsl_rng_env_setup();
  
  typeRNG = gsl_rng_mt19937;
  generator = gsl_rng_alloc (typeRNG);
  gsl_rng_set(generator, seed);     
}

NormalDistribution::~NormalDistribution()
{
  gsl_rng_free(generator);
}

TYPE_FLOAT NormalDistribution::generateUnitNormal()
{
  return myGSL::gsl_ran_gaussian(generator, 1);
}
