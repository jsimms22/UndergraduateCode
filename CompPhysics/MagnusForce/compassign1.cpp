#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

using namespace std;

const double g = 9.80665;	//Acceleration due to gravity near 
			  	//the surface of the Earth m/s^2
const double m = .04593;  	//Mass of a golf ball, in kg
const double tstep = .01;	//Size of the time step, in seconds
const double r = .0427;		//Radius of a golf ball
const int n = 500;		//Size of the arrays for data storage
				//Values do not need to saved in an array
				//It would be more efficient to write out
				//After the value is calculated to a file
const int q = n - 2;		//Iterating max for Verlet loop
const double pi = 4 * atan(1);	//Approximate value of PI
const double k = 1.225 * pi * r * r; //Drag coef * area * density of fluid
				     //Density of air is 1

/*Function to calculate total acceleration due to
  gravity, drag, and magnus effect
  INPUTS:
  integer i: Current iterated value
  double m: Mass of the object, golf ball
  double g: Acceleration due to gravity
  double x[n], y[n], z[n]: Position arrays of the ball
  double Vx[n], Vy[n], Vz[n]: Velocity arrays of the ball
  double Vwx, Vwy, Vwz: Velocity of the wind
  double Vxnet[n], Vynet[n], Vznet[n]: Net velocity arrays of the ball/wind
  double Wx, Wy, Wz: Angulare velocities of the ball
  double Ax[n], Ay[n], Az[n]: Acceleration arrays of the ball

  OUTPUTS:
  double Ax[n], Ay[n], Az[n]: Acceleration arrays of the ball 

  MATH:
  Magnus Force = Cl * rho * A * (wxv)^2
	Cl: coefficient of lift
	rho: density of the fluid the object is passing through
	A: area of the object
	wxv: angular velocity crossed with velocity of the object
  Drag Force = k * (v.v) * u
	k: coeficient of drag * A * rho; 
		k/m = .25 s^-1
	v.v: velocity of the object dotted with itself
	u: unit vector of velocity
*/
double Accel(int i, double m, double g,
             double x[n], double y[n], double z[n],
	     double Vx[n], double Vy[n], double Vz[n], 
             double Vwx, double Vwy,double Vwz, 
             double Vxnet[n], double Vynet[n], double Vznet[n],
             double Wx, double Wy, double Wz, 
             double Ax[n], double Ay[n], double Az[n], double Vmag) 
{
  if (Vynet[i+1] >= 0) 
  {
    Ax[i+1] = (((.47 * k * (pow(Vxnet[i+1],2) + pow(Vynet[i+1],2) + pow(Vznet[i+1],2)) / m) * (Vxnet[i+1] / Vmag))) + ((.25 * m * k * (Wy * Vznet[i+1] - Vynet[i+1] * Wz)) / m) + Ax[i+1] + 0;

    Ay[i+1] = (((.47 * k * (pow(Vxnet[i+1],2) + pow(Vynet[i+1],2) + pow(Vznet[i+1],2)) / m) * (Vynet[i+1] / Vmag))) + g + ((-.25 * m * k * (Wx * Vznet[i+1] - Wz * Vxnet[i+1])) / m) + Ay[i+1];

    Az[i+1] = (((.47 * k * (pow(Vxnet[i+1],2) + pow(Vynet[i+1],2) + pow(Vznet[i+1],2)) / m) * (Vznet[i+1] / Vmag))) + 0 + ((.25 * m * k * (Wx * Vynet[i+1] - Wy * Vxnet[i+1])) / m) + Az[i+1];
  }
  if (Vynet[i+1] < 0) 
  {
    Ax[i+1] = (((.47 * k * (pow(Vxnet[i+1],2) + pow(Vynet[i+1],2) + pow(Vznet[i+1],2)) / m) * (Vxnet[i+1] / Vmag))) + 0 - ((.25 * m * k * (Wy * Vznet[i+1] - Vynet[i+1] * Wz)) / m) + Ax[i+1];

    Ay[i+1] = (((-.47 * k * (pow(Vxnet[i+1],2) + pow(Vynet[i+1],2) + pow(Vznet[i+1],2)) / m) * (Vynet[i+1] / Vmag))) + g - ((.25 * m * k * (Wx * Vznet[i+1] - Wz * Vxnet[i+1])) / m) + Ay[i+1];

    Az[i+1] = (((.47 * k * (pow(Vxnet[i+1],2) + pow(Vynet[i+1],2) + pow(Vznet[i+1],2)) / m) * (Vznet[i+1] / Vmag))) + 0 + ((.25 * m * k * (Wx * Vynet[i+1] - Wy * Vxnet[i+1])) / m) + Az[i+1];
  }
  return Ax[n], Ay[n], Az[n];
}

/*Function to calculate force
  INPUTS:
  integer i: Current iterated value
  double m: mass of the object, golf ball
  double Ax[n], Ay[n], Az[n]: Acceleration arrays of the ball
  double Fx[n], Fy[n], Fz[n]: Force arrays of the ball
  
  OUTPUTS:
  double Fx[n], Fy[n], Fz[n]: Force arrays of the ball    
*/
double Force(int i, double m, 
	     double Ax[n], double Ay[n], double Az[n], 
	     double Fx[n], double Fy[n], double Fz[n]) 
{
  Fx[i+1] = m * Ax[i+1];
  Fy[i+1] = m * Ay[i+1];
  Fz[i+1] = m * Az[i+1];
  return Fx[n], Fy[n], Fz[n];
}

/*Function for writing out, and plotting saved data
  INPUTS:
  integer i: Current iterated value
  double x[n], y[n], z[n]: Position arrays of the ball
  double Ax[n], Ay[n], Az[n]: Acceleration arrays of the ball
  double Fx[n], Fy[n], Fz[n]: Force arrays of the ball

  OUTPUTS:
  datafile data: data.dat file where desired arrays are printed to
  gnuplot graph: 2D graph for position coordinates to be printed
*/
void writeData(int i, double x[n], double y[n], double z[n], 
	       double Vxnet[n], double Vynet[n], double Vznet[n], 
	       double Ax[n], double Ay[n], double Az[n])
{
  int ii = i + 1;
  string filename = "golfBall.dat";

  ofstream dataFile;
  dataFile.open(filename.c_str());

  for (i = 0; i < ii; i++)
  {
    dataFile << x[i] << " " << y[i] << /*" " << z[i] <<*/ "\n";
    //dataFile << Vxnet[i] << " " << Vynet[i] << " " << Vznet[i] << endl;
    //dataFile << Ax[i] << " " << Ay[i] << " " << Az[i] << endl;
  }
  dataFile.close();
    
  string plot = "graph.gnuplot";

  ofstream gnuIchooseYou;
  gnuIchooseYou.open(plot.c_str());

  gnuIchooseYou << "set terminal png medium size 760,360 background '#FFFFFF'" << "\n";
  gnuIchooseYou << "set xlabel \"x (meters)\"" << "\n";
  gnuIchooseYou << "set ylabel \"y (meters)\"" << "\n";
  //gnuIchooseYou << "set zlabel \"z (meters)\"" << "\n";
  gnuIchooseYou << "set title \"Distance Graph - 1600 rpm\"" << "\n";
  gnuIchooseYou << "set key off" << "\n";
  gnuIchooseYou << "set output \"tmp.png\"" << "\n";
  gnuIchooseYou << "plot \"golfBall.dat\"" << " with lines" << "\n";
  gnuIchooseYou.close();
  string pngFilename = "distancegraph3";
  system("gnuplot graph.gnuplot");
  rename("tmp.png", pngFilename.c_str());
  string write = "display " + pngFilename;
  system(write.c_str());
}

int main() 
{
  double x[n], y[n], z[n]; 	//Cartesian position coordinates, m
  double Vx[n], Vy[n], Vz[n];   //Velocity of the object, m/s
  double Vwx, Vwy, Vwz; 	//Velocity of the wind, m/s
  double Vxnet[n], Vynet[n], Vznet[n]; //Net velocity of wind and object, m/s
  double Wx, Wy, Wz; 		//Angular velocity for Magnus Effect, rpm
  double Ax[n], Ay[n], Az[n];   //Acceleration of the object, m/s^2
  double Fx[n], Fy[n], Fz[n];   //Forces, N
  double Vmag;			//Magnitude of the velocity vector, ||m/s||

  /*.611 radians is approx 35 degrees*/

  /*Initial position coordinates*/
  x[0] = 0; y[0] = 0; z[0] = 0;
  /*----------------------------*/
  
  /*Initial velocity of the object*/ 
  Vx[0] = (500 * cos(.611)); 
  Vy[0] = (500 * sin(.611)); 
  Vz[0] = 0;
  /*------------------------------*/
  
  /*Initial velocity of the wind*/
  Vwx = 0; 
  Vwy = 0; 
  Vwz = 0;
  /*----------------------------*/

  /*Net velocities of the object and ball*/
  Vxnet[0] = Vx[0] + Vwx; 
  Vynet[0] = Vy[0] + Vwy; 
  Vznet[0] = Vz[0] + Vwz;
  /*-------------------------------------*/

  /*Magnitude of the net velocity vectors*/
  Vmag = abs(sqrt(pow(Vxnet[0],2) 
         + pow(Vynet[0],2) 
         + pow(Vznet[0],2)));
  /*-------------------------------------*/

  /*Initial acceleration due to gravity*/
  Ax[0] = 0; Ay[0] = g; Az[0] = 0;
  /*-----------------------------------*/

  /*Initial angular velocity of the ball*/
  /*Currently set to 1600, units are rpm*/
  Wx = 0*(376.99 * cos(.611)) * 2; 
  Wy = 0*(376.99 * sin(.611)) * 2; 
  Wz = 0;
  /*------------------------------------*/

  /*Euler's method to approximate second time step*/
  Accel(-1, m, g, x, y, z, Vx, Vy, Vz, Vwx, Vwy, Vwz, Vxnet, Vynet, Vznet, 
        Wx, Wy, Wz, Ax, Ay, Az, Vmag);
  Force(-1, m, Az, Ay, Az, Fx, Fy, Fz);

  Vx[1] = Vx[0] - Ax[0] * tstep;
  Vy[1] = Vy[0] - Ay[0] * tstep;
  Vz[1] = Vz[0] - Az[0] * tstep;

  x[1] = (Vx[0] * tstep) - (.5 * Ax[0] * pow(tstep, 2)) + x[0];
  y[1] = (Vy[0] * tstep) - (.5 * Ay[0] * pow(tstep, 2)) + y[0];
  z[1] = (Vz[0] * tstep) - (.5 * Az[0] * pow(tstep, 2)) + z[0];

  Vxnet[1] = Vx[1] + Vwx;
  Vynet[1] = Vy[1] + Vwy;
  Vznet[1] = Vz[1] + Vwz;
  /*---------------------------------------------*/

  /*Verlet's method to approximate all subsequent time steps*/
  int i;
  for (i = 0; i <= q; i++) 
  {
    /*Update the magnitude of the velocity vectors*/
    /*Using previous time step's net velocities   */
    Vmag = abs(sqrt(pow(Vxnet[i+1],2) 
	   + pow(Vynet[i+1],2) 
	   + pow(Vznet[i+1],2)));
    /*--------------------------------------------*/
    
    /*Calculate accelerations and forces*/
    Accel(i, m, g, x, y, z, Vx, Vy, Vz, Vwx, Vwy, Vwz, Vxnet, Vynet, Vznet, 
          Wx, Wy, Wz, Ax, Ay, Az, Vmag);
    Force(i, m, Ax, Ay, Az, Fx, Fy, Fz);
    /*----------------------------------*/

    /*Calculate next time step's values using Verlet's method*/
    x[i+2] = 2 * x[i+1] - x[i] - Ax[i+1] * pow(tstep,2);
    y[i+2] = 2 * y[i+1] - y[i] - Ay[i+1] * pow(tstep,2);
    z[i+2] = 2 * z[i+1] - z[i] - Az[i+1] * pow(tstep,2);
        
    Vx[i+2] = (x[i+2] - x[i+1]) / tstep;
    Vy[i+2] = (y[i+2] - y[i+1]) / tstep;
    Vz[i+2] = (z[i+2] - z[i+1]) / tstep;

    Vxnet[i+2] = Vx[i+2] + Vwx;
    Vynet[i+2] = Vy[i+2] + Vwy;
    Vznet[i+2] = Vz[i+2] + Vwz;
    /*-------------------------------------------------------*/

    /*Output to terminal for debugging purposes*/
    /*cout << "position step " << i << " in cartesian: " 
         << x[i] << " " << y[i] << " " << z[i] << endl;
    */
    /*-----------------------------------------*/
    
    /*Conditional statement for ending once the object hits the ground*/ 
    if (y[i] < 0) 
    {
      break;
    }
    /*----------------------------------------------------------------*/
  }

  writeData(i, x, y, z, Vxnet, Vynet, Vznet, Ax, Ay, Az);

  return 0;
}
