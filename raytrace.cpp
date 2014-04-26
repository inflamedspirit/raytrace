#include <iostream>
#include <fstream>
#include <stdio.h> //needed?
#include <stdlib.h>//needed?
#include <cmath>
#include <algorithm> // for max
#include <cstdlib>   // for atoi

#define E 2.71
#define PI 3.14
#define URAND() (double)rand()/(double)RAND_MAX

#define dt 0.1
#define hw 0.01 //really should be hw/c, since this is used as a momentum
#define gw 0.05

#define omega0_global 1
#define omega_global .9
#define C 1

using namespace std;

class Particle
{
public:
  double x[3];          ///< Position
  double v[3];          ///< Velocity
  double internal_time; ///< The internal time is the time since last absorbtion/emission
  bool state;           ///< Ground or Excited state
  double tau;           ///< Emission time constant (needed?)

  Particle(double _x, double _y, double _z, double _vx, double _vy, double _vz); //constructor
  Particle& operator-=(const Particle& rhs);                                     //subtraction op

  void update(void);      ///< Runs every time step, updates everything
  void info(ofstream*);   ///< Prints particle information (for getting data)
  void printinfo(void);   ///< Prints particle information to cout

  void checkEmit(void);   ///< Randomly determines if should emit, then accounts for momentum change
  void checkAbsorb(void); ///< Randomly determines if should absorb, then accounts for momentum change
  double Lorentzian(double omega0, double omega); ///<
  double ZeemanShift(int jtransition); ///< Describes the shift in energy of the given state
  double DopplerShift(int lasernumber); ///< Describes the shift in observed frequency of the given beam
  double P_emit(double);   ///< The probability of emitting, depends on time
  double P_absorb(double); ///< The probability of absorbing, depends on time
  void interact(Particle); ///< Allows a particle to interact with another particle
};

Particle::Particle(double _x, double _y, double _z, double _vx, double _vy, double _vz)
{
  x[0]=_x;
  x[1]=_y;
  x[2]=_z;
  v[0]=_vx;
  v[1]=_vy;
  v[2]=_vz;
  internal_time=0;
  tau = 1;
}

void Particle::update(void)
{
  x[0]+=v[0]*dt;
  x[1]+=v[1]*dt;
  x[2]+=v[2]*dt;
  internal_time+=dt;
  
  if(state){
    checkEmit();
  } else {
    checkAbsorb();
  }

}

void Particle::printinfo(void)
{
  cout << x[0] << " " << x[1] << " " << x[2] << " "
       << v[0] << " " << v[1] << " " << v[2] << " "
       << sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) << " "
       << sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) << " "
       << (int)state << "\n";
}

void Particle::info(ofstream * file)
{
  (*file) << x[0] << " " << x[1] << " " << x[2] << " "
       << v[0] << " " << v[1] << " " << v[2] << " "
       << sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) << " "
       << sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) << " "
       << (int)state << "\n";
}


void Particle::checkEmit(void)
{
    // Based on some decay rate, emit photons in a random direction
    // and get a hw kick in the opposite direction
    if( URAND() < P_emit(internal_time) )
      {
	state=false;
	double theta = URAND()*PI;
	double phi = URAND()*2*PI;
	v[0] += hw*cos(phi)*sin(theta);
	v[1] += hw*sin(phi)*sin(theta);
	v[2] += hw*cos(theta);
	internal_time = 0.0f;
      }
}

void Particle::checkAbsorb(void)
{

  if( URAND() < P_absorb(internal_time) )
      {
	state=true;
	internal_time = 0.0f;

	//We want a probability distribution with 6 outcomes
	//(or 12 outcomes when we consider magnetic splitting)
	//The probability will be proportional to the linewidth
	//of the particular transition

	double P_beam[6];
	double P_total=0;

	//For each beam we have a different detuning
	//labeling is in the order x, y, z, -x, -y, -z,
	for(int i=0; i<6; i++){
	  P_beam[i] = Lorentzian( omega0_global + ZeemanShift(0), omega_global + DopplerShift(i));
	  P_total += P_beam[i];
	}
	//normalize probabilities
	for(int i=0; i<6; i++){
	  P_beam[i] = P_beam[i]/P_total;
	}

/*
	//debug sum
	P_total=0;
	for(int i=0; i<6; i++){
	  P_total += P_beam[i];
	}
	//debugprint
	cout << "DEBUG:" << P_total << " ";
	for(int i=0; i<6; i++){
	  cout << P_beam[i] << " ";
	}
	cout << "\n";
*/

	//turn into running sum
	for(int i=1; i<6; i++){
	  P_beam[i] += P_beam[i-1];
	}


	
	double r = URAND();
	//pick one direction randomly
	if(r<P_beam[0]){
	  v[0]-=hw;
	}else if(r<P_beam[1]){
	  v[1]-=hw;
	}else if(r<P_beam[2]){
	  v[2]-=hw;
	}else if(r<P_beam[3]){
	  v[0]+=hw;
	}else if(r<P_beam[4]){
	  v[1]+=hw;
	}else if(r<P_beam[5]){
	  v[2]+=hw;
	}else{
	  cout << "ERROR! Probability!\n";
	}

	//Simple and temporary method that absorbs photons from
	//the most tuned laser (according to "zeeman" shifts)
	//and the most tuned laser (according to "doppler" shifts)       
	//This is not physically accurate in any way.
	/*
	double poses[6];
	for(int i=0; i<3; i++){
	  poses[i] = x[i] + v[i];    // we might want a scaling factor?
	  poses[i+3] = -x[i] - v[i];
	}
	double cur_max=poses[0];
	int cur_i=0;
	for(int i=1; i<6; i++){
	  if( poses[i] > cur_max ){
	    cur_i = i;
	    cur_max = poses[i];
	  }
	}
	if( cur_i < 3 ){
	  v[cur_i] -= hw;
	} else {
	  v[cur_i-3] += hw;
	}
	*/
      }
}

double Particle::Lorentzian(double omega0, double omega)
{
  return  omega0/(pow(omega0-omega,2.0)+pow(omega0/2.0,2.0));
}

double Particle::DopplerShift(int lasernumber) //v,lasernumber
{
  switch(lasernumber){
  case 0:
    return v[0]*omega0_global/C;
    break;
  case 1:
    return v[1]*omega0_global/C;
    break;
  case 2:
    return v[2]*omega0_global/C;
    break;
  case 3:
    return -v[0]*omega0_global/C;
    break;
  case 4:
    return -v[1]*omega0_global/C;
    break;
  case 5:
    return -v[2]*omega0_global/C;
    break;
  }

}

double Particle::ZeemanShift(int jtransition) // x,t,transition
{
  switch(jtransition){
  case 0:
    return 0;
    break;
  case 1:
    return 0;
    break;
      }
}

double Particle::P_emit(double internal_time)
{
  return 1-pow(E,-internal_time/tau);
}

double Particle::P_absorb(double internal_time)
{
  return 1.0;
}

void Particle::interact(Particle otherp)
{
  //lets do... orbits?

  /*
    double r32 = pow(
    abs((double)
    (pow((x[0]-otherp.x[0]),2.0)+
    pow((x[1]-otherp.x[1]),2.0)+
    pow((x[2]-otherp.x[2]),2.0))
    )
    ,3.0/2.0);
    double deltax[3];
    
    for(int i=0; i<3; i++){
    deltax[i] = (otherp.x[i]-x[i]);
    v[i] += gw*deltax[i]/r32;
    }
  */
  
  
  //try spring instead
  double deltax[3];

  for(int i=0; i<3; i++){
    deltax[i] = (otherp.x[i]-x[i]);
    v[i] += gw*deltax[i];
  }
  
}


//subtraction functionality
Particle& Particle::operator-=(const Particle& rhs)
{
  // actual addition of rhs to *this
  for(int i=0;i<3;i++){
    this->x[i] = this->x[i] - rhs.x[i];
    this->v[i] = this->v[i] - rhs.v[i];
  }
  this->state = false;
  this->internal_time = 0;
  return *this;
}
inline Particle operator-(Particle lhs, const Particle& rhs)
{
  lhs -= rhs;
  return lhs;
}



int main( int argc, char** argv)
{
  if(argc < 2){
    cout << "Number of steps?\n";
    return 0;
  }

  //for data output
  ofstream ofstream_p1;
  ofstream ofstream_p2;
  ofstream ofstream_pd;
  ofstream_p1.open("data1");
  ofstream_p2.open("data2");
  ofstream_pd.open("datad");

  #define NUM_OF_PARTICLES 1000
  Particle * Ps[NUM_OF_PARTICLES];
  for(int i=0; i<NUM_OF_PARTICLES; i++){
    Ps[i] = new Particle((URAND()-0.5)*2*5,(URAND()-0.5)*2*5,(URAND()-0.5)*2*5,0.0,0.0,0.0);
  }

  Particle p1(0,0,0,1,0,0);
  Particle p2(0,0,0,0,0,0);

  for(int i=0; i<atoi(argv[1]); i++){

    //update particles
    p1.update();
    p2.update();


    //add an interaction between particles
    //p1.interact(p2);
    //    p2.interact(p1);

    //print info to files
    p1.info(&ofstream_p1);
    p2.info(&ofstream_p2);

    /*
    (p1-p2).info(&ofstream_pd);
    for(int i=0; i<NUM_OF_PARTICLES; i++){
      Ps[i]->update();
      cout << Ps[i]->x[0] << " " << Ps[i]->x[1] << " " << Ps[i]->x[2] << " ";
    }
    cout << "\n";
    */

  }
   
  ofstream_p1.close();
  ofstream_p2.close();

}
