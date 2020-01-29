#include "interpt.h"
#include <TRandom3.h>

double SLP_multi(int pmt, double x_target, double y_target)
{
//#include <iostream>
//#include <cmath>

//double x_target=-3.75;
//double y_target=-2.765;


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
/*
int a=0;
float coordinates [3][81]; //Things i changed
        for (float x=-4.0; x<=4.0; x=x+1.0)
        {

            for (float y=-4.0; y<=4.0; y=y+1.0)
            {
            coordinates [0][a] = x-x_target;
            coordinates [1][a] = y-y_target;
            a++;
            }
        }
*/


//get Z
float z=0;
float slp_response [16];  //Things I changed
if (x_target>4||x_target<-4||y_target>4||y_target<-4)
{
    for(int i=0; i<16; i++)
    {
        coordinates [2][i]=0;
        slp_response [i]= coordinates[2][i];
       // std::cout<<coordinates[0][i]+x_target<<" "<<coordinates[1][i]+y_target<<" "<<coordinates[2][i]; //things I changed
       //std::cout<<"\n";
    }
}
else
{
for(int i=0; i<16; i++)  //Things I changed
{
 z=interpt(coordinates[1][i],coordinates[0][i]);
 coordinates[2][i]=z;
 slp_response [i]= coordinates[2][i];
// std::cout<<coordinates[0][i]+x_target<<" "<<coordinates[1][i]+y_target<<" "<<coordinates[2][i]; //things I changed
// std::cout<<"\n";
}
}

//Fill Histo

/*
ChannelEnergy = new TH2D("RPD_TotalChannelEnergy","Total Energy Deposited, MeV",4,0,4,4,0,4);
ChannelEnergy->GetXaxis()->SetBinLabel(4,"Column1");
ChannelEnergy->GetXaxis()->SetBinLabel(3,"Column2");
ChannelEnergy->GetXaxis()->SetBinLabel(2,"Column3");
ChannelEnergy->GetXaxis()->SetBinLabel(1,"Column4");

ChannelEnergy->GetYaxis()->SetBinLabel(4,"Row1");
ChannelEnergy->GetYaxis()->SetBinLabel(3,"Row2");
ChannelEnergy->GetYaxis()->SetBinLabel(2,"Row3");
ChannelEnergy->GetYaxis()->SetBinLabel(1,"Row4");

ChannelEnergy->SetOption("textcolz");
*/
TRandom3 r;
static Double_t erpd[16];

for(int n = 0; n < 16; n++) 
{
  erpd[n] = slp_response [n]; //Z values replace PMT_1.....
}

/*
//First Row (Top)
ChannelEnergy->Fill(0.5,3.5,erpd[0]);//Total Event Energy For channel 1
ChannelEnergy->Fill(1.5,3.5,erpd[1]);//Total Event Energy For channel 2
ChannelEnergy->Fill(2.5,3.5,erpd[2]);//Total Event Energy For channel 3
ChannelEnergy->Fill(3.5,3.5,erpd[3]);//Total Event Energy For channel 4

//Second Row (Second From Top)
ChannelEnergy->Fill(0.5,2.5,erpd[4]);//Total Event Energy For channel 5
ChannelEnergy->Fill(1.5,2.5,erpd[5]);//Total Event Energy For channel 6
ChannelEnergy->Fill(2.5,2.5,erpd[6]);//Total Event Energy For channel 7
ChannelEnergy->Fill(3.5,2.5,erpd[7]);//Total Event Energy For channel 8

//Third Row (Third From Top)
ChannelEnergy->Fill(0.5,1.5,erpd[8]);//Total Event Energy For channel 9
ChannelEnergy->Fill(1.5,1.5,erpd[9]);//Total Event Energy For channel 10
ChannelEnergy->Fill(2.5,1.5,erpd[10]);//Total Event Energy For channel 11
ChannelEnergy->Fill(3.5,1.5,erpd[11]);//Total Event Energy For channel 12

//Fourth Row (Bottom Row)
ChannelEnergy->Fill(0.5,0.5,erpd[12]);//Total Event Energy For channel 13
ChannelEnergy->Fill(1.5,0.5,erpd[13]);//Total Event Energy For channel 14
ChannelEnergy->Fill(2.5,0.5,erpd[14]);//Total Event Energy For channel 15
ChannelEnergy->Fill(3.5,0.5,erpd[15]);//Total Event Energy For channel 16


ChannelEnergy->Draw();
*/
return erpd[pmt];







/*
    slp_response [0]= PMT_1_EdistPerHit.GetEntries();
    slp_response [1]= PMT_2_EdistPerHit.GetEntries();
    slp_response [2]= PMT_3_EdistPerHit.GetEntries();
    slp_response [3]= PMT_4_EdistPerHit.GetEntries();
    slp_response [4]= PMT_5_EdistPerHit.GetEntries();
    slp_response [5]= PMT_6_EdistPerHit.GetEntries();
    slp_response [6]= PMT_7_EdistPerHit.GetEntries();
    slp_response [7]= PMT_8_EdistPerHit.GetEntries();
    slp_response [8]= PMT_9_EdistPerHit.GetEntries();
    slp_response [9]= PMT_10_EdistPerHit.GetEntries();
    slp_response [10]= PMT_11_EdistPerHit.GetEntries();
    slp_response [11]= PMT_12_EdistPerHit.GetEntries();
    slp_response [12]= PMT_13_EdistPerHit.GetEntries();
    slp_response [13]= PMT_14_EdistPerHit.GetEntries();
    slp_response [14]= PMT_15_EdistPerHit.GetEntries();
    slp_response [15]= PMT_16_EdistPerHit.GetEntries();
*/

}

