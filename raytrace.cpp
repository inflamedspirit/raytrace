#include <iostream>
#include <fstream>
#include <stdio.h> //needed?
#include <stdlib.h>//needed?
#include <cmath>
#include <algorithm> // for max
#include <cstdlib>   // for atoi

#define E 2.71
#define PI 3.14
#define URAND() ((double)rand()/(double)RAND_MAX)

#define dt 0.1
#define hw 0.01 //really should be hw/c, since this is used as a momentum
#define gw 0.05

#define omega0_global 1
#define omega_global .9
#define C 1

#define ATOM_X 0
#define ATOM_Y 0
#define ATOM_Z 0
#define STEP_SIZE 0.001
#define TOLERANCE 0.01
#define TIME_MAX 100000
using namespace std;

class Surface
{
public:
  double (*func)(double, double); //surface function def

  Surface(double (*func_)(double,double));
};

Surface::Surface(double (*func_)(double,double)){
  func = func_;
}

//surfaces
double interface1(double x, double y)
{
  double z = 1;
  return z;
}


class Ray
{
public:
  double x[3];          ///< Position
  double phi,theta;
  double phase; 
  double intensity;

  Ray(double _x, double _y, double _z, double _phase, double _intensity); //constructor
  Ray(double _x, double _y, double _z, double _phi, double _theta, double _phase, double _intensity);


  void info(ofstream*);   ///< Prints particle information (for getting data)
  void printinfo(void);   ///< Prints particle information to cout

  /*
  Ray& operator-=(const Ray& rhs);                                     //subtraction op
  void update(void);      ///< Runs every time step, updates everything
  */

};

Ray::Ray(double _x, double _y, double _z, double _phi, double _theta, double _phase, double _intensity)
{
  x[0]=_x;
  x[1]=_y;
  x[2]=_z;
  phi = _phi;
  theta = _theta;
  phase = _phase;
  intensity = _intensity;
}

Ray::Ray(double _x, double _y, double _z, double _phase, double _intensity)
{
  x[0]=_x;
  x[1]=_y;
  x[2]=_z;
  phi = URAND()*2*PI;
  theta = URAND()*PI;
  phase = _phase;
  intensity = _intensity;
}

void Ray::printinfo(void)
{
  cout << x[0] << " " << x[1] << " " << x[2] << " "
       << phi  << " " << theta << " "
       << intensity << " " << phase << " "
       << sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) << " "
       << "\n";
}

void Ray::info(ofstream * file)
{
  cout << x[0] << " " << x[1] << " " << x[2] << " "
       << phi  << " " << theta << " "
       << intensity << " " << phase << " "
       << sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) << " "
       << "\n";
}

class Vect
{
public:
  double x1[3];
  double x2[3];
  Vect(Ray ray, Surface surf);
  void info(ofstream*);   ///< Prints particle information (for getting data)
  void printinfo(void);   ///< Prints particle information to cout


};

Vect::Vect(Ray ray, Surface surf)
{
  double dx = cos(ray.phi)*sin(ray.theta)*STEP_SIZE;
  double dy = cos(ray.phi)*sin(ray.theta)*STEP_SIZE;
  double dz = cos(ray.theta)*STEP_SIZE;

  x1[0] = ray.x[0];
  x1[1] = ray.x[1];
  x1[2] = ray.x[2];

  x2[0] = ray.x[0];
  x2[1] = ray.x[1];
  x2[2] = ray.x[2];

  double diff = abs(x2[2]-surf.func(x2[0],x2[1]));
  int time = 0;

  while(diff > TOLERANCE && time < TIME_MAX){
    time++;
    x2[0]+=dx;
    x2[1]+=dy;
    x2[2]+=dz;

    diff = abs(x2[2]-surf.func(x2[0],x2[1]));
    if(time == TIME_MAX){
      x2[0] = ray.x[0];
      x2[1] = ray.x[1];
      x2[2] = ray.x[2];
      return;
    }
  }

}


void Vect::printinfo(void)
{
  cout << x1[0] << " " << x1[1] << " " << x1[2] << " "
       << x2[0] << " " << x2[1] << " " << x2[2] << " "
       << "\n";
}

void Vect::info(ofstream * file)
{
  cout << x1[0] << " " << x1[1] << " " << x1[2] << " "
       << x2[0] << " " << x2[1] << " " << x2[2] << " "
       << "\n";
}



/*
//subtraction functionality
Ray& Ray::operator-=(const Ray& rhs)
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
inline Ray operator-(Ray lhs, const Ray& rhs)
{
  lhs -= rhs;
  return lhs;
}
*/




int main( int argc, char** argv)
{
  if(argc < 2){
    cout << "Number of steps?\n";
    return 0;
  }

  //for data output
  //  ofstream ofstream1;
  //  ofstream1.open("data1");

  //make the surface
  Surface s1(interface1);

  //make the ray  Ray p1(x,y,z,phi,theta,phase,intensity);
  //  Ray r1(0.0,0.0,0.0,0.0,(2*PI/360.0)*4.0,0.0,1.0);
  //  Ray r2(0.0,0.0,0.0,0.0,(2*PI/360.0)*8.0,0.0,1.0);

  #define NUM_OF_RAYS 1000
  Ray * rays[NUM_OF_RAYS];
  Vect * vects[NUM_OF_RAYS];
  for(int i=0; i<NUM_OF_RAYS; i++){
    rays[i] = new Ray(0.0,0.0,0.0,(URAND()*2.0*PI),(URAND()*PI),0.0,1.0);
    (vects[i]) = new Vect(*(rays[i]),interface1);
    vects[i]->printinfo();
  }

   
  //  ofstream1.close();

}
