#include "NormalDistribution.h"

//const int NormalDistribution::TR1 = 1;
const int NormalDistribution::GSL = 2;


NormalDistribution::NormalDistribution(int type_)
{
  // seed with epoch time
  // causes problems if two runs are started at the same second
  init(type_, time(NULL));
}

NormalDistribution::NormalDistribution(int type_, unsigned int seed)
{
  init(type_, seed);
}

void NormalDistribution::init(int type_, unsigned int seed)
{
  type = type_;
  //if ( type == TR1 )
    //{
      //eng.seed(seed);
      //normal = Normal (0, 1); //unit normal
      //normal.reset();
      
    //}
  if ( type == GSL )
    {
      const gsl_rng_type * T;
      
      gsl_rng_env_setup();
      
      //T = gsl_rng_default;
      T = gsl_rng_mt19937; //faster than TR1
      r = gsl_rng_alloc (T);
      gsl_rng_set(r, seed);     
    }
}

NormalDistribution::~NormalDistribution()
{
  if ( type == GSL )
    gsl_rng_free(r);
}

TYPE_FLOAT NormalDistribution::unitNormal()
{
  //if ( type == TR1 )
  //  return normal(eng); //pick a random number
  if ( type == GSL )
    return myGSL::gsl_ran_gaussian(r, 1);
  return 0;
}
