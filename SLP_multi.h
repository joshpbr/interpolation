#include "interpt.h"
#include <TRandom3.h>

double SLP_multi(int pmt, double x_target, double y_target)
{


  //get slp response
double coordinates [3][16];

  //get x
for( int i = 0; i<16; i++ )
{
    if (i%4==0)
        coordinates [0][i]=-3.0-x_target;
    if (i%4==1)
        coordinates [0][i]=-1.0-x_target;
    if (i%4==2)
        coordinates [0][i]=1.0-x_target;
    if (i%4==3)
        coordinates [0][i]=3.0-x_target;
}

//get y
for(int i=0; i<16; i++)
{
    if(i/4==0)
    coordinates [1][i]=3.0-y_target;

    if(i/4==1)
    coordinates [1][i]=1.0-y_target;

    if(i/4==2)
    coordinates [1][i]=-1.0-y_target;

    if(i/4==3)
    coordinates [1][i]=-3.0-y_target;
}


//get Z
float z=0;
float slp_response [16];  //Things I changed
if (x_target>4||x_target<-4||y_target>4||y_target<-4)
{
    for(int i=0; i<16; i++)
    {
        coordinates [2][i]=0;
        slp_response [i]= coordinates[2][i];
    
    }
}
else
{
for(int i=0; i<16; i++)  //Things I changed
{
 z=interpt(coordinates[1][i],coordinates[0][i]);
 coordinates[2][i]=z;
 slp_response [i]= coordinates[2][i];

}
}

static Double_t erpd[16];

for(int n = 0; n < 16; n++) 
{
  erpd[n] = slp_response [n]; 
}

return erpd[pmt];

}

