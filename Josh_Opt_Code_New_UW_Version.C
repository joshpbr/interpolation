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
#include <iostream>

int Josh_Opt_Code_New_UW_Version() {
	Int_t experiments = 1;
	Double_t Widths = 1;
	std::ofstream resultfile("SL_altitude_results.txt");
	Double_t output [1][10];

	for (Int_t f = 0; f<experiments; f++)
	{
		cout << "experiment # " << f << " of " << experiments << endl;

		Double_t mu = 0;

		//Creating StdDev for X coords
		TRandom3 xsigma_mt (0);
		Double_t xsigma = xsigma_mt.Uniform(0.0, 2.0);

		//Particles are number of coord pairs created
		Int_t particles = 100;
		Double_t xy [particles][2];

		//Variables used to verify if coords are in gaussian dist
		Int_t xstddev1 = 0;
		Int_t xstddev2 = 0;
		Int_t xstddev3 = 0;

		Int_t ystddev1 = 0;
		Int_t ystddev2 = 0;
		Int_t ystddev3 = 0;

		for (Int_t i = 0; i<particles; i++)
		{
			TRandom3 xcoor(0);
			Double_t x = xcoor.Gaus(mu, Widths);
			xy [i][0] = x;
		}


		TRandom3 ysigma_mt (0);
		Double_t ysigma = ysigma_mt.Uniform(0.0, 2.0);

		for (Int_t i = 0; i<particles; i++) 
		{
			TRandom3 ycoor(0);
			Double_t y = ycoor.Gaus(mu, Widths);
			xy [i][1] = y;
		}	

		Double_t Simga_p_x = 0;
		Double_t Simga_p_y = 0;

		for (Int_t i = 0; i<particles; i++)
		{
			Simga_p_x = Simga_p_x + xy[i][0];
		}

		for (Int_t i = 0; i<particles; i++)
		{
			Simga_p_y = Simga_p_y + xy[i][1];
		}

		Double_t Q_p_x = 0;
		Double_t Q_p_y = 0;

		Q_p_x = Simga_p_x / particles;
		Q_p_y = Simga_p_y / particles;

		//Transformations

		//Rotation

			const Double_t pi = 3.141592653589792;
			const Double_t twopi = pi * 2;

			TRandom3 tmt (0);
			Double_t theta_shift = tmt.Uniform(0.0, twopi);

			Double_t coordinate_set [particles][2];

			for (Int_t i = 0; i<particles; i++)
			{
				coordinate_set [i][0] = xy [i][0] * cos(theta_shift) - xy [i][1] * sin(theta_shift);
				coordinate_set [i][1] = xy [i][1] * cos(theta_shift) + xy [i][0] * sin(theta_shift);
			}

		//Linear Transform
		
		TRandom3 xmt (0);
		Double_t x_shift =  2;

		TRandom3 ymt (0);
		Double_t y_shift = 0;

		Double_t x_avg = 0;
		Double_t y_avg = 0;

		for (Int_t i = 0; i<particles; i++)
		{
			coordinate_set [i][0] = (coordinate_set [i][0] + x_shift);
			coordinate_set [i][1] = (coordinate_set [i][1] + y_shift);
		}

		//Reaction Plane Calc

		x_avg = 0;
		y_avg = 0;

		for (Int_t i = 0; i<particles; i++)
		{
			x_avg = x_avg + coordinate_set [i][0];
			y_avg = y_avg + coordinate_set [i][1];
		}

		x_avg = x_avg / particles;
		y_avg = y_avg / particles;

		Double_t reaction_plane = atan(y_avg / x_avg);

		//Corrects for atan preference

		
		if (x_avg<0 && y_avg>0)
			reaction_plane = reaction_plane + pi;
		if (x_avg<0 && y_avg<0)
			reaction_plane = reaction_plane + pi;
		if (x_avg>0 && y_avg<0)
			reaction_plane = reaction_plane + twopi;
		

		//SL Feeding
		Double_t Z_hadd[16] = {0};
		Double_t Z_temp[16] = {0};

		for (Int_t i = 0; i<particles; i ++)
		{
			for (Int_t a = 0; a<16; a++)
			{
				Z_temp[a] = SLP_multi(a, coordinate_set [i][0], coordinate_set [i][1]);
				Z_hadd[a] = Z_hadd[a] + Z_temp[a];
			}
		}

		TRandom3 r;
		Double_t erpd[16] = {0};

		for (Int_t q = 0; q<16; q++)
		{
			erpd[q] = Z_hadd[q];
		}

		//Calc Reaction plane

		Double_t phi[4][4] = {0};
		Double_t posx[4] = {-3.0, -1.0, 1.0, 3.0};
		Double_t posy[4] = {3.0, 1.0, -1.0, -3.0};
		Double_t xrp[4][4] = {0};
		Double_t yrp[4][4] = {0};

		//sets xy coords to center of each quartz seg
		for (Int_t i = 0; i<4; i++)
		{
			for (Int_t j = 0; j<4; j++)
			{
				xrp[i][j] = posx[j];
				yrp[i][j] = posy[i];
			}
		}

		//Changes from coords to angle
		for (Int_t i = 0; i<4; i++)
		{
			for (Int_t j = 0; j<4; j++)
			{
				phi[i][j] = atan2(yrp[i][j], xrp[i][j]);
			}
		}

		//Get Pmt count
		Int_t pmthits[4][4] = {0};
		pmthits[0][0] = erpd[0]*1;
		pmthits[0][1] = erpd[1]*1;
		pmthits[0][2] = erpd[2]*1;
		pmthits[0][3] = erpd[3]*1;
		pmthits[1][0] = erpd[4]*1;
		pmthits[1][1] = erpd[5]*1;
		pmthits[1][2] = erpd[6]*1;
		pmthits[1][3] = erpd[7]*1;
		pmthits[2][0] = erpd[8]*1;
		pmthits[2][1] = erpd[9]*1;
		pmthits[2][2] = erpd[10]*1;
		pmthits[2][3] = erpd[11]*1;
		pmthits[3][0] = erpd[12]*1;
		pmthits[3][1] = erpd[13]*1;
		pmthits[3][2] = erpd[14]*1;
		pmthits[3][3] = erpd[15]*1;

		/*Double_t sigma_e = 0;

		for (Int_t i = 0; i<16; i++)
		{
			sigma_e = sigma_e + erpd[i];
		}

		Double_t n_x_sigma_e = 0;
		Double_t n_y_sigma_e = 0;
		Double_t n_x = 3.5;
		
		for (Int_t i = 0; i<4; i++)
		{
			Double_t n_y = 3.5;

			for (Int_t j = 0; j<4; j++)
			{
				n_x_sigma_e = n_x_sigma_e + (pmthits[i][j] * n_x);
				n_y_sigma_e = n_y_sigma_e + (pmthits[i][j] * n_y);
				n_y = n_y - 2.0;
			}
			n_x = n_x - 2.0;
		}

		Double_t Q_x_vec = 0;012;
		Double_t Q_y_vec = 0;
		Q_x_vec = n_x_sigma_e / sigma_e;
		Q_y_vec = n_y_sigma_e / sigma_e;
		
		*/


		// For event plane 1
		Int_t p = 1;

		Double_t x_vec_pmt = 0;
		Double_t y_vec_pmt = 0;
		Double_t x_vec_rpd = 0;
		Double_t y_vec_rpd = 0;

		for (Int_t i = 0; i<4; i++)
		{
			for (Int_t j = 0; j<4; j++)
			{
				x_vec_pmt += pmthits[i][j] * cos(p*phi[i][j]);
				y_vec_pmt += pmthits[i][j] * sin(p*phi[i][j]);
			}
		}

		//Calc psi
		Double_t psi_pmt = 0;
		Double_t psi_rpd = 0;
		psi_pmt = (atan(y_vec_pmt / x_vec_pmt)) / p;
		psi_rpd = (atan(y_vec_rpd / x_vec_rpd)) / p;


		//corrects for atans preference
		
		if (x_avg<0 && y_avg>0)
			psi_pmt = psi_pmt + pi;
		if (x_avg<0 && y_avg<0)
			psi_pmt = psi_pmt + pi;
		if (x_avg>0 && y_avg<0)
			psi_pmt = psi_pmt + twopi;
		if (reaction_plane>1.2 && psi_pmt <- 1.0)
		{
			psi_pmt = psi_pmt + pi;
		}
		if (reaction_plane>4.6 && psi_pmt>7.0)
		{
			psi_pmt = psi_pmt - pi;
		}
		if (reaction_plane>4.2 && psi_pmt<2.0)
			psi_pmt = psi_pmt + pi;
		if (reaction_plane<2.0 && psi_pmt>4.0)
			psi_pmt = psi_pmt - pi;
		if (psi_pmt<reaction_plane-1.5 && reaction_plane<2 && reaction_plane>0.5)
			psi_pmt = psi_pmt + pi;
			

		//make results array
		output [f][0] = reaction_plane;
		output [f][1] = psi_pmt;
		


TH2* NRPDBlockDiagram = new TH2D("NRPD Block Diagram", "Total Energy Deposited", 4, 0, 4, 4, 0, 4);

////////////////////////////////////////////START THE BLOCK DIAGRAMS///////////////////////////////////////////////////////

NRPDBlockDiagram->GetXaxis()->SetBinLabel(4, "Column4");
NRPDBlockDiagram->GetXaxis()->SetBinLabel(3, "Column3");
NRPDBlockDiagram->GetXaxis()->SetBinLabel(2, "Column2");
NRPDBlockDiagram->GetXaxis()->SetBinLabel(1, "Column1");

NRPDBlockDiagram->GetYaxis()->SetBinLabel(4, "Row4");
NRPDBlockDiagram->GetYaxis()->SetBinLabel(3, "Row3");
NRPDBlockDiagram->GetYaxis()->SetBinLabel(2, "Row2");
NRPDBlockDiagram->GetYaxis()->SetBinLabel(1, "Row1");

NRPDBlockDiagram->SetOption("textcolz");


//First Row (Top)
NRPDBlockDiagram->Fill(3.5, 3.5, Z_hadd[3]/100); //Total Event Energy For channel 7
NRPDBlockDiagram->Fill(3.5, 2.5, Z_hadd[7]/100); //Total Event Energy For channel 3
NRPDBlockDiagram->Fill(3.5, 1.5, Z_hadd[11]/100);//Total Event Energy For channel 16
NRPDBlockDiagram->Fill(3.5, 0.5, Z_hadd[15]/100);//Total Event Energy For channel 12

//Second Row (Second From Top)
NRPDBlockDiagram->Fill(2.5, 3.5, Z_hadd[2]/100); //Total Event Energy For channel 6
NRPDBlockDiagram->Fill(2.5, 2.5, Z_hadd[6]/100); //Total Event Energy For channel 2
NRPDBlockDiagram->Fill(2.5, 1.5, Z_hadd[10]/100);//Total Event Energy For channel 13
NRPDBlockDiagram->Fill(2.5, 0.5, Z_hadd[14]/100); //Total Event Energy For channel 9

//Third Row (Third From Top)
NRPDBlockDiagram->Fill(1.5, 3.5, Z_hadd[1]/100); //Total Event Energy For channel 5
NRPDBlockDiagram->Fill(1.5, 2.5, Z_hadd[5]/100); //Total Event Energy For channel 1
NRPDBlockDiagram->Fill(1.5, 1.5, Z_hadd[9]/100);//Total Event Energy For channel 14
NRPDBlockDiagram->Fill(1.5, 0.5, Z_hadd[13]/100);//Total Event Energy For channel 10

//Fourth Row (Bottom Row)
NRPDBlockDiagram->Fill(0.5, 3.5, (Z_hadd[0]/100)); //Total Event Energy For channel 8
NRPDBlockDiagram->Fill(0.5, 2.5, (Z_hadd[4]/100)); //Total Event Energy For channel 4
NRPDBlockDiagram->Fill(0.5, 1.5, (Z_hadd[8]/100));//Total Event Energy For channel 15
NRPDBlockDiagram->Fill(0.5, 0.5, (Z_hadd[12]/100));//Total Event Energy For channel 11

//////////////////////////////////////////////END THE BLOCK DIAGRAMS///////////////////////////////////////////////////////

TCanvas* b1 = new TCanvas();
NRPDBlockDiagram->SetTitle(Form("Run %d NRPD Block Diagram, All Cent", particles));
gStyle->SetOptFit(); //Allows you to optimize your histogram upon opening in the TBrowser.
gStyle->SetOptStat(0); //Gets rid of stupid white box from upper right corner
NRPDBlockDiagram->SetMarkerSize(1.8); //Changes the font size within the block diagram.
NRPDBlockDiagram->Draw();
b1->SaveAs(Form("NRPD_BlockDiagram_AllCent_Run%d.pdf", experiments));

std::ofstream newfile("energies.txt");

	
		newfile<< Z_hadd[0]<<" , "<<Z_hadd[1]<<" , "<<Z_hadd[2]<<" , "<<Z_hadd[3]<<" , "<<Z_hadd[4]<<" , "<<Z_hadd[5]<<" , "<<Z_hadd[6]<<" , "<<Z_hadd[7]<<" , "<<Z_hadd[8]<<" , "<<Z_hadd[9]<<" , "<<Z_hadd[10]<<" , "<<Z_hadd[11]<<" , "<<Z_hadd[12]<<" , "<<Z_hadd[13]<<" , "<<Z_hadd[14]<<" , "<<Z_hadd[15]<<"\n";
		newfile.close();

		}

TCanvas *c1 = new TCanvas("c1","c", 800, 1000);
c1->SetGrid();

Double_t input_rp[20000]={0.};
Double_t output_rp[20000]={0.};

std::ofstream newfile("output.txt");



    for(int i=0; i<experiments; i++)
        {
        newfile<< output [i][0]<<" , "<< output [i][1]<<"\n";
        }
        newfile.close();



for (Int_t j=0; j<experiments; j++)
{
 input_rp[j] = output [j][0];
 output_rp[j] = output [j][1];
} 
TGraph *gr = new TGraph(experiments, input_rp, output_rp);

TLine *line = new TLine(0.0,0.0,6.3,6.3);
gr->SetTitle("Calculated Reaction Plane vs Input Reaction Plane");
gr->GetXaxis()->SetTitle("Input Reaction Plane (Radians)");
gr->GetYaxis()->SetTitle("Output Reaction Plane (Radians)");
gr->Draw("A*");
line->SetLineColor(kBlue);
line->Draw();
c1->SaveAs("UWrpd.pdf");


}
