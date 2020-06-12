#include "Polymer.h"

Polymer::Polymer()
{
}

Polymer::Polymer(int noOfMonomers,
		 std::string sequence,
		 bool random,
                 bool fixedStart,
                 Vector3D start,
                 Vector3D orientation,
                 std::vector<Particle *> * particles)
{
  this->noOfMonomers = noOfMonomers;
  
  monomers.resize(noOfMonomers);
   
  // calculate initial positions and create Particle array
  for (int i = noOfMonomers - 1; i > -1; i--)
    {
      Vector3D position;
      
      bool fixed = false;
      if (i==0 && fixedStart)
	fixed = true;

      Particle * next = NULL;
      if (i < noOfMonomers - 1)
	next = &monomers[i+1];

      // Particles arranged in a straight line
      if (WorldSettings::initialCondition == WorldSettings::IC_LINE) {
	position = start + i * WorldSettings::bondLength * orientation;
      } else {
	if ( i == noOfMonomers - 1 ) {
	  position = start; ///< \todo start is position of last monomer
	} else {
	  // Self Avoiding Walk
	  /// \todo check for correctness after rewrite
	  /// \todo reimplement counter for impossible configurations
	  if (WorldSettings::initialCondition == WorldSettings::IC_SAW) {
	    bool valid;
	    do {
	      valid = true;
	      Vector3D step = noiseTerm(1);
	      position = monomers[i+1].r +
		WorldSettings::bondLength * step / step.magnitude();

	      for (int j = 2 ; valid && j + i < noOfMonomers ; j++)
		valid = isSelfAvoiding(monomers[i + j].r - position);
	    } while (!valid);
	  }
	  
	  // Random Walk
	  if (WorldSettings::initialCondition == WorldSettings::IC_RW) {
	    Vector3D step = noiseTerm(1);
	    position = monomers[i+1].r +
	      WorldSettings::bondLength * step / step.magnitude();
	  }
	}
      }
      monomers[i] = Particle(position,
			  AMINO_ACID,
			  sequence[i],
			  next,
			  fixed);
    }
  for(int i = 0; i < noOfMonomers; i++)
    {  
      particles -> push_back(&(monomers[i])); 
    }
}

Polymer::~Polymer()
{
}

TYPE_FLOAT Polymer::endToEndDistanceSquared() const
{
  return (monomers[0].r - monomers[noOfMonomers-1].r).magnitudeSquared();
}

Vector3D Polymer::endToEndVector() const
{
  if (noOfMonomers <= 1)
    return Vector3D();
  return monomers[noOfMonomers-1].r - monomers[0].r;
}

Vector3D Polymer::centreOfMass() const
{
  Vector3D rCOM;
  for (int i = 0; i < noOfMonomers; i++)
    {
      rCOM += monomers[i].r;
    }
  return rCOM/noOfMonomers;
}

TYPE_FLOAT Polymer::radiusOfGyrationSquared() const
{
  TYPE_FLOAT rGsquared = 0;
  Vector3D rCOM = centreOfMass();
  
  for (int i = 0; i < noOfMonomers; i++)
    {
      rGsquared += (monomers[i].r - rCOM).magnitudeSquared();
    }
  return rGsquared/noOfMonomers;
}

/// \todo duplicate function in World
TYPE_FLOAT Polymer::averageBondLengthSquared() const
{
  if (noOfMonomers <= 1)
    return 0;
  TYPE_FLOAT sumOfSquaredBondLength = 0;
  for (int i = 1; i < noOfMonomers; i++)
    {
      Vector3D bondVector = monomers[i].r - monomers[i-1].r;
      sumOfSquaredBondLength += bondVector.magnitudeSquared();
    }
  return sumOfSquaredBondLength / (noOfMonomers - 1);
}

// check for fixed particle occurs in World, when positions are updated
void Polymer::simulate(TYPE_FLOAT t)
{
  if (noOfMonomers > 1)
    {
      monomers[0].dr += forceBond(monomers[0].r - monomers[1].r);
      for (int i = 1 ; i < noOfMonomers - 1 ; i++)
	{
	  monomers[i].dr += forceBond(monomers[i].r - monomers[i+1].r) 
	    + forceBond(monomers[i].r - monomers[i-1].r);
        }
      monomers[noOfMonomers - 1].dr += forceBond(monomers[noOfMonomers - 1].r
					      - monomers[noOfMonomers - 2].r);
    }
}

/// \todo should use Particle LJ radii
bool Polymer::isSelfAvoiding(Vector3D r) const
{
  if (r.magnitude() > WorldSettings::bondLength)
    return true;
  return false;
}

void Polymer::draw() const {
#ifdef GL
  GLUquadricObj* pQuadric = gluNewQuadric();
  int numSlices = 10;
  int numStacks = 10;
  GLfloat bondRadius = 0.2;

  for (int i = 0; i < noOfMonomers - 1; i++)
    {
      glPushMatrix();
      glTranslatef(monomers[i].r.x, 
		   monomers[i].r.y, 
		   monomers[i].r.z);
      
      Vector3D bondVector = monomers[i+1].r-monomers[i].r;
      
      /* Draw a cylinder along the direction of the bond vector
       *
       * Before the call to gluCylinder, we need to rotate from (0,0,1) to the
       * bondVector
       * 
       * Naive implementation breaks in several cases, so used:
       * https://github.com/curran/renderCyliner/blob/master/renderCylinder.c
       */
      
      GLfloat rotAngle;// = acos(bondVector.z/bondVector.magnitude()) * 180.0 / M_PI;
      Vector3D rotVector(- bondVector.y / bondVector.z,
			 bondVector.x / bondVector.z,
			 0);
      
      if (fabs(bondVector.z) < 1.0e-3) {
	rotAngle = acos( bondVector.x / bondVector.magnitude() ) * 180.0 / M_PI; 
	if ( bondVector.y <= 0.0 )
	  rotAngle = - rotAngle;
	glRotatef(90.0, 0.0, 1, 0.0);
	glRotatef(rotAngle, -1.0, 0.0, 0.0); 
      }
      else {
	rotAngle = acos( bondVector.z / bondVector.magnitude() ) * 180.0 / M_PI;
	if ( bondVector.z <= 0.0 )
	  rotAngle = -rotAngle;
	glRotatef(rotAngle, rotVector.x, rotVector.y, rotVector.z);
      }
	
      gluCylinder(pQuadric,
		  bondRadius,
		  bondRadius,
		  bondVector.magnitude(),
		  numSlices,
		  numStacks);

      glPopMatrix();
    }
  free(pQuadric);
#endif
}
