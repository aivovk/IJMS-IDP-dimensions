#ifndef NORMALDISTRIBUTION_H
#define NORMALDISTRIBUTION_H

#include <ctime>
#include <gsl/gsl_rng.h>

// avoids conflicts with MKL
namespace myGSL
{
#include <gsl/gsl_randist.h>
}

#include "Vector3D.h"

/*
  Wrapper for random number generation
 */
class NormalDistribution
{
 public:
  NormalDistribution();
  NormalDistribution(unsigned int seed);
  void init(unsigned int seed);
  ~NormalDistribution();

  // generate a random number from the unit normal distribution
  TYPE_FLOAT generateUnitNormal();

 protected:
 private:
  gsl_rng * generator;
};

#endif // NORMALDISTRIBUTION_H
