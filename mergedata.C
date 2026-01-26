//Macro that Analyse the root files
// 
// Before calling this macro the files 
//    DigOut.root   
//
// files should be created by running
//    raw_to_root.C

#include <stdio.h>
#include <fcntl.h>
#include <TTree.h>
#include <TFile.h>
#include <TObject.h>
#include <TNtuple.h>
#include "Riostream.h"

#include <math.h>
#include "TMath.h"
#include <TRandom.h>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"

#include <string.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TGraph.h"
#include "TLine.h"

#include "ROOT/RDataFrame.hxx" //to add tree RDataFrame.hxx
#include <iostream>

#include "WaveFormFunctions.C"

// class Pulse : public TObject {
// public:
//   float fadqTime;
//   float fMax;
//   float fInt;
//   float fWidth;
//   float fT0_30;
//   float fBase;  s
//   float fMaxBin;
//   float fRCharge;
 
//   Pulse() { };
//   ClassDef(Pulse,1);
// };


void mergedata(){
  TFile* file = new TFile("Treed3.root");
  TChain *tree = new TChain("tree");
  // Open the ROOT file
  tree->Add("Treed3.root");

  //tree variables
  int run, event; 
  Pulse *B0 = new Pulse();   
  Pulse *B1 = new Pulse();
  Pulse *B2 = new Pulse();
  Pulse *B3 = new Pulse();  
  
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("event",&event);  

  tree->SetBranchAddress("B0",&B0);  //ch1
  tree->SetBranchAddress("B1",&B1);  //ch2
  tree->SetBranchAddress("B2",&B2);  //ch3 
  tree->SetBranchAddress("B3",&B3);  //ch4-Trig 


  Long64_t nEvents = tree->GetEntries();
  int iEvent;   
  // int const noCh = 4;
  // bool adqComplete[noCh];
  // for(int j=0; j<noCh; j++){
  //     adqComplete[j]=false;
  // }

  Int_t i, iChannels = bins_per_record+1;

  //histogramas para los valores màximos (amplitudes)
  auto hB0Max = new TH1F("hB0Max","Amplitude Distribution ch1; Amplitude, V; Entries",150,-0.05,0.14);
  auto hB1Max = new TH1F("hB1Max","Amplitude Distribution ch2; Amplitude, V; Entries",150,-0.05,0.14);
  auto hB2Max = new TH1F("hB2Max","Amplitude Distribution ch3; Amplitude, V; Entries",150,-0.05,0.14);
  auto hB3Max = new TH1F("hB3Max","Amplitude Distribution ch4; Amplitude, V; Entries",150,-1.,1.1);
  //para la dist de ancho de pulso
  auto hB0Width = new TH1F("hB0Width", "PulseWidth Distribution ch1; Time (ns); Entries", 150, 0, 600);
  auto hB1Width = new TH1F("hB0Width", "PulseWidth Distribution ch2; Time (ns); Entries", 150, 0, 600);
  auto hB2Width = new TH1F("hB0Width", "PulseWidth Distribution ch3; Time (ns); Entries", 150, 0, 600);
  auto hB3Width = new TH1F("hB0Width", "PulseWidth Distribution ch4; Time (ns); Entries", 150, 63, 64);
  //para persistencias
  TH2F *h_Per0= (TH2F*)file->Get("h_Pers"); 
  TH2F *h_Per1= (TH2F*)file->Get("h_PersTwo"); 
  TH2F *h_Per2= (TH2F*)file->Get("h_PersThree"); 
  TH2F *h_Per3= (TH2F*)file->Get("h_PersFour"); 
  //para waveforms
  TH1F *h_wf1 = (TH1F*)file->Get("h_WF1");
  TH1F *h_wf2 = (TH1F*)file->Get("h_WF2");
  TH1F *h_wf3 = (TH1F*)file->Get("h_WF3");
  TH1F *h_wf4 = (TH1F*)file->Get("h_WF4");
  //dist de diferencia de tiempos
  auto hDiffTime = new TH1F("hDiffTime ", "(A) T0_30 - T0_CFD; Time (ns); Entries", 60, -1.2,0.2);
  auto hDiffTime2 = new TH1F("hDiffTime2 ", "(B) T0_30 - T0_CFD; Time (ns); Entries", 60, -1.2,0.2);
  auto hDiffTime3 = new TH1F("hDiffTime3 ", "(C) T0_30 - T0_CFD; Time (ns); Entries", 60, -1.2,0.2);

  auto hDiffTime4 = new TH1F("T0_CFD", "t(C3)-t(B3), uncorrected; Time (ns); Entries", 25, -6, 6);
  auto hDiffTime5 = new TH1F("T0_CFD", "t(B3)-t(A3), uncorrected; Time (ns); Entries", 25, -6, 6);

  //for correction
  auto h_u0 = new TProfile("h_u[0]","t(C3) - t(B3) vs Width(C3); Width, ns;dT,ns",150.,0.,250,-20.,20.);
  auto h_u1 = new TProfile("h_u[1]","t(B3) - t(C3) vs Width(B3); Width, ns;dT,ns",150.,0.,250,-20.,20.);
  auto h_t = new TH1F("t'=T0_CFD-m*width", "t'(C3)-t'(B3), corrected; Time (ns); Entries", 40, -6, 6);
  

  //auto c4 = new TCanvas("c4", "Waveforms");x
  //c4->Divide(2,2);

  //llenado de los histogramas de distribuciòn de amplitudes bajo la condiciòn
  for (iEvent=0; iEvent<nEvents; iEvent++) {
    tree->GetEntry(iEvent);

    hB0Max->Fill(B0->fMax);
    hB1Max->Fill(B1->fMax);
    hB2Max->Fill(B2->fMax);
    hB3Max->Fill(B3->fMax);

    hB0Width->Fill(B0->fWidth);
    hB1Width->Fill(B1->fWidth);
    hB2Width->Fill(B2->fWidth);
    hB3Width->Fill(B3->fWidth);

    if (B0->fMax > 0.1 && B1->fMax > 0.1 && B2->fMax > 0.1 
    && B0->fT0_30*bin_width>90 && B0->fT0_30*bin_width<100 && B1->fT0_30*bin_width>90 && B1->fT0_30*bin_width<100 && B2->fT0_30*bin_width>90 && B2->fT0_30*bin_width<100
    && B0->fT0CFD*bin_width>90 && B0->fT0CFD*bin_width<100 && B1->fT0CFD*bin_width>90 && B1->fT0CFD*bin_width<100 && B2->fT0CFD*bin_width>90 && B2->fT0CFD*bin_width<100
    && B0->fWidth > 0 && B1->fWidth > 0 && B2->fWidth > 0 ){

    hDiffTime->Fill(B0->fT0_30*bin_width - B0->fT0CFD*bin_width);  
    hDiffTime2->Fill(B1->fT0_30*bin_width - B1->fT0CFD*bin_width);
    hDiffTime3->Fill(B2->fT0_30*bin_width - B2->fT0CFD*bin_width);
    hDiffTime4->Fill(B2->fT0CFD*bin_width - B1->fT0CFD*bin_width);
    hDiffTime5->Fill(B1->fT0CFD*bin_width - B0->fT0CFD*bin_width);

    h_u0->Fill(B2->fWidth*bin_width, B2->fT0CFD*bin_width - B1->fT0CFD*bin_width);
    h_u1->Fill(B1->fWidth*bin_width, B1->fT0CFD*bin_width - B0->fT0CFD*bin_width);
    

    }

  }  

  //display amplitude histograms


  auto c3 = new TCanvas("c3", "Amplitude Distribution", 10,10,1400,700);
  c3->Divide(2,2);

  c3->cd(1);
  gPad->SetLogy();  
  hB0Max->SetFillStyle(3001);
  hB0Max->SetFillColor(kRed);
  hB0Max->Draw();
  c3->cd(2);
  gPad->SetLogy();
  hB1Max->SetFillStyle(3001);
  hB1Max->SetFillColor(kRed);
  hB1Max->Draw();
  c3->cd(3);
  gPad->SetLogy();
  hB2Max->SetFillStyle(3001);
  hB2Max->SetFillColor(kRed);
  hB2Max->Draw();
  c3->cd(4);
  gPad->SetLogy();
  hB3Max->SetFillStyle(3001);
  hB3Max->SetFillColor(kRed);
  hB3Max->Draw();

  c3->Draw();

  auto c1 = new TCanvas("c1", "Pulse Width Distribution");
  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetLogy();
  hB0Width->Draw();
  c1->cd(2);
  gPad->SetLogy();
  hB1Width->Draw();
  c1->cd(3);
  gPad->SetLogy();
  hB2Width->Draw();
  c1->cd(4);
  gPad->SetLogy();
  hB3Width->Draw();


  auto c2 = new TCanvas("c2", "Persistencias");
  c2->Divide(2,2);
  c2->cd(1);
  gPad -> SetLogz();
  h_Per0->Draw();
  c2->cd(2);
  gPad -> SetLogz();
  h_Per1->Draw();
  c2->cd(3);
  gPad -> SetLogz();
  h_Per2->Draw();
  c2->cd(4);
  gPad -> SetLogz();
  h_Per3->Draw();

  auto c5 = new TCanvas("c5", "Time Difference Distribution");
  c5->Divide(2,3);
  c5->cd(1);
  hDiffTime->Draw();
  c5->cd(2);
  hDiffTime2->Draw();
  c5->cd(3);
  hDiffTime3->Draw();
  c5->cd(4);
  hDiffTime4->Fit("gaus");
  hDiffTime4->Draw();
  c5->cd(5);
  hDiffTime5->Fit("gaus");
  hDiffTime5->Draw();


 

  auto c6 = new TCanvas("c6", "DeltaT vs Width");
  c6->Divide(2,2);
  c6->cd(1);
  h_u0->Fit("pol1", "R", "", 15, 50);
  h_u0->Draw();
  c6->cd(2);
  h_u1->Fit("pol1", "R", "", 15, 58);
  h_u1->Draw();
  c6->cd(3);
  hDiffTime4->Fit("gaus");
  hDiffTime4->Draw();

//now correction
    TF1 *fitFunc = h_u0->GetFunction("pol1");
    Double_t slope = fitFunc->GetParameter(1); //pendiente 
    std::cout << "Pendiente m1=" << slope << std::endl;
  
    TF1 *fitFunc2 = h_u1->GetFunction("pol1");
    Double_t slope2 = fitFunc2->GetParameter(1); //pendiente
    std::cout << "Pendiente m2=" << slope2 << std::endl;

    for (iEvent=0; iEvent<nEvents; iEvent++) {
      tree->GetEntry(iEvent);

      if (B0->fMax > 0.1 && B1->fMax > 0.1 && B2->fMax > 0.1 
        && B0->fT0_30*bin_width>90 && B0->fT0_30*bin_width<100 && B1->fT0_30*bin_width>90 && B1->fT0_30*bin_width<100 && B2->fT0_30*bin_width>90 && B2->fT0_30*bin_width<100
        && B0->fT0CFD*bin_width>90 && B0->fT0CFD*bin_width<100 && B1->fT0CFD*bin_width>90 && B1->fT0CFD*bin_width<100 && B2->fT0CFD*bin_width>90 && B2->fT0CFD*bin_width<100
        && B0->fWidth > 0 && B1->fWidth > 0 && B2->fWidth > 0 ){

          double_t T0B2 = B2->fT0CFD*bin_width - B2->fWidth*bin_width*(-0.06);
          double_t T0B1 = B1->fT0CFD*bin_width - B1->fWidth*bin_width*(-0.06);
          h_t->Fill(T0B2 - T0B1);
        }
    }
  c6->cd(4);
  h_t->Fit("gaus");
  h_t->Draw();


}