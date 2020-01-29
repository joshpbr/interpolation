#include <TRandom.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <TString.h>
#include <sstream>
#include "SLP_multi.h"
#include <TH2.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include "alglibinternal.h"
#include "alglibmisc.h"
#include "ap.h"
#include "dataanalysis.h"
#include "diffequations.h"
#include "fasttransforms.h"
#include "integration.h"
#include "interpolation.h"
#include "linalg.h"
#include "optimization.h"
#include "solvers.h"
#include "specialfunctions.h"
#include "statistics.h"
#include "stdafx.h"

int Josh_Optimized_Code_UnweightedVersion() {
	Int_t experiments=5000;
std::ofstream resultfile("SL_altitude_results.txt");  //the results are here, the first column is the actual rp, the second is the calculated one, the last column is the variable being measured, in this case the y position.  THis structure is preserved across output files from different experiments 
Double_t output [5000][3];
for (Int_t f=0; f<experiments; f++)
{

cout << "experiment # " << f << " of " << experiments <<  endl;

/*SL_GAUSSIAN_GENERATOR
//Created by Jeffrey Zboray, Summer 2015
//This program generates random sets of coordinates in the cartesian plane
//The points are generated in a two dimensional gaussian distribution
//This program includes functionalities to transform the coordinate pairs linearly in the cartesian plane, as well as to rotate the set of coordinates
//The fully transformed coordinate set is outputted as "pgun_coordinate_set.csv"*/
/*
//#include <iostream>
#include <TRandom.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <TString.h>
#include <sstream>
*/

    Double_t mu=0;

//Creates StdDev for  x coordinates
    TRandom3 xsigma_mt (0);
    Double_t xsigma=xsigma_mt.Uniform(0.0, 2.0);

//Particles is the number of coordinate pairs created

    Int_t particles=120; //<---------VALUE YOU MAY WISH TO CHANGE
    Double_t xy [particles][2];

//variables used later to verify coordinate set generated is in a gaussian 
//distribution
    Int_t xstddev1=0;
    Int_t xstddev2=0;
    Int_t xstddev3=0;

    Int_t ystddev1=0;
    Int_t ystddev2=0;
    Int_t ystddev3=0;

	// xsigma=xsigma*2;
	//ysigma=ysigma*2;

    for (Int_t i=0; i<particles; i++)
    {
//generates x coordinates along a normal distribution
		TRandom3 xcoor(0);
		Double_t x=xcoor.Gaus(mu, xsigma);
	//	x=x*2;//<<<<<<<<<<<<<<<<<<<<<Important number
		xy [i][0]= x;

/*Verifies x-values are generated along a normal distribution, not vital to program operation*/
        Double_t number=xy [i][0];
    if (number<=(mu+xsigma)&&number>=(mu-xsigma))
        xstddev1++;
    if ((number<(mu-xsigma)&&number>=(mu-2*xsigma))||(number>(mu+xsigma)&&number<=mu+2*xsigma))
        xstddev2++;
    if ((number<(mu-2*xsigma)||number>mu+2*xsigma))
        xstddev3++;
    }

//Creates StdDev for Y coordinates
    TRandom3 ysigma_mt (0);
    Double_t ysigma=ysigma_mt.Uniform(0.0, 2.0);


    for (Int_t i=0; i<particles; i++)
    {
//generates Y values on normal distribution
		TRandom3 ycoor (0);
		Double_t y=ycoor.Gaus(mu, ysigma);
	//	y=y*2;//<<<<<<<<<<<<important number
       		 xy [i][1]=y;
        Double_t number=xy [i][1];
/*Verifies y values are correctly generated along normal distribution*/
    if (number<=(mu+ysigma)&&number>=(mu-ysigma))
        ystddev1++;
    if ((number<(mu-ysigma)&&number>=(mu-2*ysigma))||(number>(mu+ysigma)&&number<=mu+2*ysigma))
        ystddev2++;
    if ((number<(mu-2*ysigma)||number>mu+2*ysigma))
        ystddev3++;
    }



//TRANSFORMATIONS:

//Rotation Matrix
//generates angle of rotation
    const Double_t pi=3.141592653589793;
    const Double_t twopi=pi*2;

    TRandom3 tmt (0);
    Double_t thetashift=tmt.Uniform(0.0, twopi);


    Double_t coordinate_set [particles][2];


    for (Int_t i=0; i<particles; i++)
    {
    coordinate_set [i][0]=xy [i][0]*cos(thetashift)-xy [i][1]*sin(thetashift);
    coordinate_set [i][1]=xy [i][1]* cos(thetashift)+xy [i][0]*sin(thetashift);
    }

//Linear Transformation


	TRandom3 xmt (0);
    Double_t xshift=xmt.Uniform(-4.0, 4.0);


	TRandom3 ymt(0);
	Double_t yshift=ymt.Uniform(-4.0,4.0);

Double_t xavg=0;
Double_t yavg=0;
/*for(Int_t i=0; i<particles; i++)
{
    xavg=xavg+coordinate_set [i][0];
    yavg=yavg+coordinate_set [i][1];
}
    xavg=xavg/particles;
    yavg=yavg/particles;*/
//center fixing
//    yshift=0.68682;
//    xshift=0.561124;


    for (Int_t i=0; i<particles; i++)
    {
    coordinate_set [i][0]=(coordinate_set [i][0]+ xshift);
    coordinate_set [i][1]=(coordinate_set [i][1]+ yshift);
    }

//Reaction plane calculation
xavg=0;
yavg=0;
for(Int_t i=0; i<particles; i++)
{
    xavg=xavg+coordinate_set [i][0];
    yavg=yavg+coordinate_set [i][1];
}
    xavg=xavg/particles;
    yavg=yavg/particles;
    
    Double_t reaction_plane=atan(yavg/xavg);
    //std::cout<<"Center="<<xavg<<" , "<<yavg<<"\n";
    //std::cout<<"Uncorrected Reaction Plane="<<reaction_plane<<endl;

//corrects for atan's preference for answers <pi/2
if (xavg<0 && yavg>0)
      reaction_plane=reaction_plane+pi;
if (xavg<0 && yavg<0)
      reaction_plane=reaction_plane+pi;
if (xavg>0 && yavg<0)
      reaction_plane=reaction_plane+twopi;

//std::cout<<"Reaction Plane="<<reaction_plane;



//csv creation
    std::ofstream newfile("pgun_coordinate_set.txt");
   // newfile.open("pgun_coordinate_set.csv");


    for(Int_t i=0; i<particles; i++)
        {
//        std::cout<<coordinate_set [i][0]<<" , "<<coordinate_set [i][1]<<"\n";
        newfile<<coordinate_set [i][0]<<" , "<<coordinate_set [i][1]<<"\n";
        }

    newfile.close();
/*
//macro writer
//macmaker
 for (int i=0; i<particles; i++)
    {
    std::stringstream converter;
    converter<<i;
    std::string n = converter.str();
    std::string macfilename_string = "condor"+n+"_gauss_pgun.mac";
    std::string root_file_string="condor"+n+"_gauss_pgun.root";
    char* condor_root_file = root_file_string.c_str();
    char* macfilename = macfilename_string.c_str();
//    std::cout<<"\nmacfilename= "<<macfilename;
    std::ofstream macfile(macfilename);

   // macfile.open(macfilename);
    macfile<<"/test/histo/setRootName "<<condor_root_file<<"\n";//<--CHANGE
    macfile<<"/run/verbose 2 \n";
    macfile<<"/event/verbose 0 \n";
    macfile<<"/gps/source/multiplevertex true\n\n";

    macfile<<"/gps/number 1\n";
    macfile<<"/gps/pos/type Point \n";
    macfile<<"/gps/particle neutron \n";
    macfile<<"/gps/direction 0. 0. 1.\n";
    macfile<<"/gps/energy 2.55 TeV\n";
    macfile<<"/tracking/verbose 0\n";
    macfile<<"/gps/position "<<coordinate_set [i][0]<<" "<<coordinate_set [i][1]<<" -6 cm\n";
    macfile<<"\n\n/run/beamOn 1";
    macfile.close();
    
    string jdlfilename_string = "gauss_condor"+n+".jdl";
    char* jdlfilename = jdlfilename_string.c_str();
  //  cout<<"\njdlfilename= "<<jdlfilename;
    ofstream jdlfile(jdlfilename);
    jdlfile<<"universe = vanilla\n";
    jdlfile<<"Executable = condor-executable.sh\n";
    jdlfile<<"getenv = True\n";
    jdlfile<<"requirements = (OpSys == \"LINUX\")\n";
	jdlfile<<"should_transfer_files = NO\n";
	jdlfile<<"Requirements = TARGET.FileSystemDomain == \"privnet\"\n";
	jdlfile<<"Output = /data/users/jzboray/ZDCTDR_07_15/logs/simple_$(cluster)_$(process).stdout\n";
	jdlfile<<"Error  = /data/users/jzboray/ZDCTDR_07_15/logs/simple_$(cluster)_$(process).stderr\n";
	jdlfile<<"Log    = /data/users/jzboray/ZDCTDR_07_15/logs/simple_$(cluster)_$(process).condor\n";
	jdlfile<<"Arguments = simple_$(cluster) $(process) HitDevel /data/users/jzboray/ZDCTDR_07_15/logs "<< "/data/users/jzboray/ZDCTDR_07_15/holder/"<<macfilename<<" 1234567 1\n";
	jdlfile<<"Queue 1\n";
	jdlfile.close();
    }
*/
//SL feeding
	Double_t Z_hadd[16]={0};
	Double_t Z_temp[16]={0};
	/*std::cout<<endl;*/
	for(Int_t i=0; i<particles; i++)
        {
	//std::cout<<"COOORINATES: "<<coordinate_set [i][0]<<","<<coordinate_set [i][1]<<endl;
		for(Int_t a=0; a<16; a++)
		{
		Z_temp[a]=SLP_multi(a, coordinate_set [i][0] , coordinate_set [i][1]);
		Z_hadd[a]=Z_hadd[a]+Z_temp[a];
		//std::cout<<"\t"<<Z_hadd[16]<<endl;
		}
        }

    

TRandom3 r;
Double_t erpd[16]={0};


for(Int_t q = 0; q < 16; q++)
{
  erpd[q] = Z_hadd[q];
}

//Fill histo
//Courtesy of Christopher Ferraioli
TH2D* ChannelEnergy;
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
ChannelEnergy->SaveAs("ShowerLibraryPlot_temp.jpg");


//CALCULATE REACTION PLANE
//Courtesy of Christopher Ferraioli


//////sets up variables we're going to use
Double_t phi[4][4]={0};
Double_t posx[4]={-3.0,-1.0,1.0,3.0};
Double_t posy[4]={3.0,1.0,-1.0,-3.0};
Double_t xrp[4][4]={0};
Double_t yrp[4][4]={0};


/////sets up x y coordinates of center of each quartz segment
for(Int_t i = 0; i < 4; i++)
{
  for(Int_t j = 0; j < 4; j++)
  {
  xrp[i][j]=posx[j];
  yrp[i][j]=posy[i];
  }
}

///////changes from coordinates to angle
for(Int_t i = 0; i < 4; i++)
{
  for(Int_t j = 0; j < 4; j++)
  {
  phi[i][j]=atan2(yrp[i][j],xrp[i][j]);
  //phi[i][j]=atan2(1.0,-1.0);
  }
}


/////get pmt count
Int_t pmthits[4][4]={0};
pmthits[0][0]=erpd[0];
pmthits[0][1]=erpd[1];
pmthits[0][2]=erpd[2];
pmthits[0][3]=erpd[3];
pmthits[1][0]=erpd[4];
pmthits[1][1]=erpd[5];
pmthits[1][2]=erpd[6];
pmthits[1][3]=erpd[7];
pmthits[2][0]=erpd[8];
pmthits[2][1]=erpd[9];
pmthits[2][2]=erpd[10];
pmthits[2][3]=erpd[11];
pmthits[3][0]=erpd[12];
pmthits[3][1]=erpd[13];
pmthits[3][2]=erpd[14];
pmthits[3][3]=erpd[15];




////for event plane 1
Int_t p = 1;

////calculate x and y vectors
Double_t xvecpmt = 0;
Double_t yvecpmt = 0;
Double_t xvecrpd = 0;
Double_t yvecrpd = 0;
for(Int_t i = 0; i < 4; i++)
{
  for(Int_t j = 0; j < 4; j++)
  {
  xvecpmt+=pmthits[i][j]*cos(p*phi[i][j]);
  yvecpmt+=pmthits[i][j]*sin(p*phi[i][j]);
  //xvecrpd+=rpdenergy[i][j]*cos(p*phi[i][j]);
  //yvecrpd+=rpdenergy[i][j]*sin(p*phi[i][j]);
  }
}



////calculate psi
Double_t psipmt = 0;
Double_t psirpd = 0;
psipmt=(atan(yvecpmt/xvecpmt))/p;
psirpd=(atan(yvecrpd/xvecrpd))/p;



//corrects for atan's preference for answers <pi/2
if (xavg<0 && yavg>0)
      psipmt=psipmt+pi;
if (xavg<0 && yavg<0)
      psipmt=psipmt+pi;
if (xavg>0 && yavg<0)
      psipmt=psipmt+twopi;
if (reaction_plane>1.2 && psipmt<-1.0)
      {
	psipmt=psipmt+pi;
      }
if (reaction_plane>4.6 && psipmt>7.0)
      {
        psipmt=psipmt-pi;
      }
if (reaction_plane>4.2 && psipmt<2.0)
        psipmt=psipmt+pi; //This was slice 2
if (reaction_plane<2.0 && psipmt>4.0)
        psipmt=psipmt-pi;

//added to fix the mystery slice
if (psipmt<reaction_plane-1.5 && reaction_plane<2 && reaction_plane>0.5)
        psipmt=psipmt+pi;



//cout<<"Psi_1 for PMT distro is "<<psipmt<<"\n"<<endl;
//cout<<"SL_altitude2.0_plotter.C"<<endl;
//cout<<"Experiment number "<<f<<"\n"<<endl;

cout << 0 ;
//make results array
output [f][0] = reaction_plane;
output [f][1] = psipmt;
output [f][2] = yshift;
//ends master loop
}

cout << 0 ;
TCanvas *c1 = new TCanvas("c1","c", 800, 1000);
c1->SetGrid();

Double_t input_rp[5000]={0.};
Double_t output_rp[5000]={0.};

for (Int_t j=0; j<experiments; j++)
{
 input_rp[j] = output [j][0];
 output_rp[j] = output [j][1];
} 
TGraph *gr = new TGraph(experiments, input_rp, output_rp);


  std::ofstream newfile("rp_coordinate_setW.txt");

    for(Int_t i=0; i<experiments; i++)
        {
        if(std::isnan(output_rp[i])==true)
        continue;
        
        if(std::isnan(input_rp[i])==true)
        continue;
        
        else
        newfile<<input_rp[i]<<" , "<<output_rp[i]<<"\n";
        }

    newfile.close();

cout << 0 ;
TLine *line = new TLine(0.0,0.0,6.3,6.3);
gr->SetTitle("Calculated Reaction Plane vs Input Reaction Plane");
gr->GetXaxis()->SetTitle("Input Reaction Plane (Radians)");
gr->GetYaxis()->SetTitle("Output Reaction Plane (Radians)");
gr->Draw("A*");
line->SetLineColor(kBlue);
line->Draw();
c1->SaveAs("UWrpd.pdf");







//results output txt creation
  //  std::ofstream resultfile("SL_altitude_results.txt");
    //resultfile.open("SL_altitude_results.txt");
    for(Int_t k=0; k<experiments; k++)
        {
        resultfile<<output [k][0]<<","<<output [k][1]<<","<<output [k][2]<<"\n";
//"SL_altitude_results.txt")
         }
//    resultfile.close();
         return 0;
}
