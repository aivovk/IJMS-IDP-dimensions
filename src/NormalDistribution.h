#ifndef NORMALDISTRIBUTION_H
#define NORMALDISTRIBUTION_H

#include <ctime>
//#include <tr1/random>
#include <gsl/gsl_rng.h>

// avoids conflicts with MKL
namespace myGSL
{
#include <gsl/gsl_randist.h>
}

#include "Vector3D.h"

// generator and distribution from TR1/random
//typedef std::tr1::ranlux64_base_01 Engine;
//typedef std::tr1::normal_distribution<TYPE_FLOAT> Normal;

/*
  Wrapper for random number generation
 */
class NormalDistribution
{
 public:
  // constants to define which library to use
  //static const int TR1;
  static const int GSL;

  NormalDistribution(int type_);
  NormalDistribution(int type_, unsigned int seed);
  void init(int type_, unsigned int seed);
  ~NormalDistribution();

  // generate a random number from the unit normal distribution
  TYPE_FLOAT unitNormal();

 protected:
 private:
  // GSL or TR1/random
  int type;

  // TR1/random random number generator
  //Engine eng;
  //Normal normal;

  // GSL random number generator
  gsl_rng * r;
};

#endif // NORMALDISTRIBUTION_H
