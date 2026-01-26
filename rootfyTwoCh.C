// Before to start this macro you should check:
// 1. The input file should be without headers
// 2. the output from the Tektronix sometimes is -inf, we chosen the value and remplaced it for 0.0  in order
// 3. The input should be without headers and the "," should be remplaced by " "

// This program were made to explore quickly the
// output from a digitizer
// To run this program you should write:
// root wfDisplay\(60\)
// where 60 is the number of events to explore


#include <stdio.h>
#include <fcntl.h>
#include <TTree.h>
#include <TFile.h>
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

#include <iostream>
#include <fstream>
#include <string>

#include <stdlib.h> //Atof funtion

#include "WaveFormFunctions.C"

void rootfyTwoCh(){
// void rootfy(int noWFfrom, int noWFto, int iSet, TString output){ 
    int noWFfrom; int noWFto; int iSet; TString output;
    //rootfy.C\(1,1001,1,\"data1_1\"\)
    //noWFfrom 1,1001,2001,2500
    //iSet the adq number 1-11
    //output file format dataiSet_1,2,3

    TString mainfile ="d3";
   gROOT->Reset();
   Bool_t eof = false;
   
//    TString mainfile=Form("waveforms%01d",iSet);
   Int_t iChannels = bins_per_record+1;
   Int_t hour=0, minute=0,second=0;
   Float_t time=1.,adqTime0=1.,adqTime1=1.,adqTime2=1.,adqTime3=1., wf0=1., wf1=1.,wf2=1., wf3=1.;
   
  int iBLfrom=0.05/bin_width;     // time window for baseline subrtaction
  int iBLto=50./bin_width;
  int iPULSEfrom=51./bin_width;  // time window for pulse
  int iPULSEto=500./bin_width;

  //for trigger pulse
  int iBLfromT =5./bin_width;
  int iBLtoT = 159./bin_width;
  int iPULSEfromT = 160./bin_width;
  int iPULSEtoT = 290./bin_width;

  float fBase=1.;     // average baseline determined beforehand   
  float fBaseCFD=1.;  // average baseline determined beforehand

//Construction the histograms for the waveforms
  TH1F *h_WF1 = new TH1F("h_WF1","h_WF1;Time [ns]; Voltage [V]",iChannels,0,bin_width*(bins_per_record-1));
  TH1F *h_WF2 = new TH1F("h_WF2","h_WF2;Time [ns]; Voltage [V]",iChannels,0,bin_width*(bins_per_record-1));
  TH1F *h_WF3 = new TH1F("h_WF3","h_WF3;Time [ns]; Voltage [V]",iChannels,0,bin_width*(bins_per_record-1));
  TH1F *h_WF4 = new TH1F("h_WF4","h_WF4;Time [ns]; Voltage [V]",iChannels,0,bin_width*(bins_per_record-1));

  h_WF1->SetStats(0);    h_WF2->SetStats(0);    h_WF3->SetStats(0);    h_WF4->SetStats(0);


  TH2F *h_Pers = new TH2F("h_Pers","h_Pers;Time [ns]; Voltage [V]",iChannels-1,0,bin_width*(bins_per_record-1),
                                      300,-0.05,0.15);
  TH2F *h_PersTwo = new TH2F("h_PersTwo","h_PersTwo;Time [ns]; Voltage [V]",iChannels-1,0,bin_width*(bins_per_record-1),
                                      300,-0.05,0.15);
  TH2F *h_PersThree = new TH2F("h_PersThree","h_PersThree;Time [ns]; Voltage [V]",iChannels-1,0,bin_width*(bins_per_record-1),
                                      300,-0.05,0.15);
  TH2F *h_PersFour = new TH2F("h_PersFour","h_PersFour;Time [ns]; Voltage [V]",iChannels-1,0,bin_width*(bins_per_record-1),
                                      75,-1,0.2);
    
  TH1F *h_WFatt1 = new TH1F("h_WFatt1","h_WFatt1;Time [ns]; Voltage [V]",iChannels,0,bin_width*(bins_per_record-1));
  TH1F *h_WFatt2 = new TH1F("h_WFatt2","h_WFatt2;Time [ns]; Voltage [V]",iChannels,0,bin_width*(bins_per_record-1));
  TH1F *h_WFatt3 = new TH1F("h_WFatt3","h_WFatt3;Time [ns]; Voltage [V]",iChannels,0,bin_width*(bins_per_record-1));
  TH1F *h_WFatt4 = new TH1F("h_WFatt4","h_WFatt4;Time [ns]; Voltage [V]",iChannels,0,bin_width*(bins_per_record-1));
  h_WFatt1->SetStats(0); h_WFatt1->SetLineColor(kRed); h_WFatt2->SetStats(0); h_WFatt3->SetStats(0); h_WFatt4->SetStats(0);
  h_WFatt2->SetLineColor(kRed); h_WFatt3->SetLineColor(kRed); h_WFatt4->SetLineColor(kRed);


  TH1F *h_WFCFD1 = new TH1F("h_WFCFD1","h_WFCFD1;Time [ns]; Voltage [V]",iChannels,0,bin_width*(bins_per_record-1));
  TH1F *h_WFCFD2 = new TH1F("h_WFCFD2","h_WFCFD2;Time [ns]; Voltage [V]",iChannels,0,bin_width*(bins_per_record-1));
  TH1F *h_WFCFD3 = new TH1F("h_WFCFD3","h_WFCFD3;Time [ns]; Voltage [V]",iChannels,0,bin_width*(bins_per_record-1));
  TH1F *h_WFCFD4 = new TH1F("h_WFCFD4","h_WFCFD4;Time [ns]; Voltage [V]",iChannels,0,bin_width*(bins_per_record-1));
  h_WFCFD1->SetStats(0); h_WFCFD1->SetLineColor(kGreen); h_WFCFD2->SetStats(0); h_WFCFD3->SetStats(0); h_WFCFD4->SetStats(0);
  h_WFCFD2->SetLineColor(kGreen); h_WFCFD3->SetLineColor(kGreen); h_WFCFD4->SetLineColor(kGreen);


  OSC_Record fWaveForm0,fWaveForm1,fWaveForm2,fWaveForm3; ;
  OSC_Record fWaveFormBL0,fWaveFormBL1,fWaveFormBL2,fWaveFormBL3;
  OSC_Record fCFDWaveForm0,fCFDWaveForm1,fCFDWaveForm2,fCFDWaveForm3;
  OSC_Record fCFDWaveFormBL0,fCFDWaveFormBL1,fCFDWaveFormBL2,fCFDWaveFormBL3;

  OSC_Record *pfCFDWaveForm[4];
  OSC_Record *pfCFDWaveFormBL[4];
  OSC_Record *pfWaveForm[4];
  OSC_Record *pfWaveFormBL[4];

  pfWaveForm[0]    = &fWaveForm0;
  pfWaveFormBL[0]  = &fWaveFormBL0;
  pfWaveForm[1]    = &fWaveForm1;
  pfWaveFormBL[1]  = &fWaveFormBL1;
  pfWaveForm[2]    = &fWaveForm2;
  pfWaveFormBL[2]  = &fWaveFormBL2;
  pfWaveForm[3]    = &fWaveForm3;
  pfWaveFormBL[3]  = &fWaveFormBL3;

  pfCFDWaveForm[0]    = &fCFDWaveForm0;
  pfCFDWaveFormBL[0]  = &fCFDWaveFormBL0;
  pfCFDWaveForm[1]    = &fCFDWaveForm1;
  pfCFDWaveFormBL[1]  = &fCFDWaveFormBL1;
  pfCFDWaveForm[2]    = &fCFDWaveForm2;
  pfCFDWaveFormBL[2]  = &fCFDWaveFormBL2;
  pfCFDWaveForm[3]    = &fCFDWaveForm3;
  pfCFDWaveFormBL[3]  = &fCFDWaveFormBL3;
  
    Pulse *B0 = new Pulse();
    Pulse *B1 = new Pulse();
    Pulse *B2 = new Pulse();
    Pulse *B3 = new Pulse();

  Pulse *pPulse[4];

  pPulse[0] = B0;
  pPulse[1] = B1;
  pPulse[2] = B2;
  pPulse[3] = B3;

  int iEvent;
  TTree tTree("tree","Parameters of pulses in Oscilloscope input channels");
  int run = 1;
  tTree.Branch("run",&run,"i1/I");//''
  tTree.Branch("event",&iEvent,"i1/I");
  tTree.Branch("B0",&B0,10000,1);  
  tTree.Branch("B1",&B1,10000,1);  
  tTree.Branch("B2",&B2,10000,1);
  tTree.Branch("B3",&B3,10000,1);
  
  TString waveNo;

  //cout << "here" <<endl;
  int adqEvents;
//   TString mainfile="BrTh1on84LY_1";
   TString recordName1="C1Trace";
//    TString recordName2=recordName1; 
   TString recordName2="C2Trace";
   TString recordName3="C3Trace";
   TString recordName4="C4Trace";

  adqEvents=1900; //max
  
  TCanvas *c0 = new TCanvas("c0","WaveForm",10,10,1600,900);
  gPad->SetGrid();
  c0->Divide(2,2);
  
//    cout << "here2" <<endl;

  ifstream in[adqEvents];  
  ifstream inTwo[adqEvents];  
  ifstream inThree[adqEvents];
  ifstream inFour[adqEvents];

   int auxCounter=-1;
   for(int iFile=1; iFile<adqEvents; iFile++){
    auxCounter++;   
       if(iFile<10) {waveNo = Form("0000%01d", iFile);}
       else if(iFile<100){waveNo = Form("000%01d", iFile);}
       else if(iFile<1000){waveNo = Form("00%01d", iFile);}
       else if(iFile<10000){waveNo = Form("0%01d", iFile);}
       else if(iFile<100000){waveNo = Form("%01d", iFile);}
    TString input1 = recordName1 + waveNo+ ".txt";
    TString input2 = recordName2 + waveNo+ ".txt";
    TString input3 = recordName3 + waveNo+ ".txt";
    TString input4 = recordName4 + waveNo+ ".txt";

    TString inputfile1=mainfile+"/"+input1; 
    TString inputfile2=mainfile+"/"+input2; 
    TString inputfile3=mainfile+"/"+input3;
    TString inputfile4=mainfile+"/"+input4;
//     cout << "rootty:" << inputfile1<<endl;
//     cout << "rootty:" << inputfile2<<endl;
    in[iFile].open(inputfile1);
    inTwo[iFile].open(inputfile2);
    inThree[iFile].open(inputfile3);
    inFour[iFile].open(inputfile4);
    //adqTime START
    //first
    jumpToLine(in[iFile],4);
         in[iFile].ignore(14);
         in[iFile] >> hour;
         in[iFile].ignore();
         in[iFile] >> minute;
         in[iFile].ignore();
         in[iFile] >> second;
    adqTime0= (hour*60.*60.)+(minute*60.) + second;
    //second
    jumpToLine(inTwo[iFile],4);
         inTwo[iFile].ignore(14);
         inTwo[iFile] >> hour;
         inTwo[iFile].ignore();
         inTwo[iFile] >> minute;
         inTwo[iFile].ignore();
         inTwo[iFile] >> second;
    adqTime1= (hour*60.*60.)+(minute*60.) + second;

    //third
    jumpToLine(inThree[iFile],4);
         inThree[iFile].ignore(14);
         inThree[iFile] >> hour;
         inThree[iFile].ignore();
         inThree[iFile] >> minute;
         inThree[iFile].ignore();
         inThree[iFile] >> second;
    adqTime2= (hour*60.*60.)+(minute*60.) + second;

    //fourth
    jumpToLine(inFour[iFile],4);
         inFour[iFile].ignore(14);
         inFour[iFile] >> hour;
         inFour[iFile].ignore();
         inFour[iFile] >> minute;
         inFour[iFile].ignore();
         inFour[iFile] >> second;
    adqTime3= (hour*60.*60.)+(minute*60.) + second;

    //adqTime Finish
    jumpToLine(in[iFile],6);
    jumpToLine(inTwo[iFile],6);
    jumpToLine(inThree[iFile],6);
    jumpToLine(inFour[iFile],6);
    
//     if(in[iFile].eof())  break;        
    if(!in[iFile].is_open())  break;       


      float att=0.3;
      int delay=40; //2*0.5ns (cada bin son 0.5ns)

    //DELAY - ATT - SIGMA - MEAN D3

//  2, 0.3 shows bad response

// ch1   8 - 0.3 - 0.8716 - 1.262
// ch2   8 - 0.3 - 0.8897 - 1.445
// ch3   8 - 0.3 - 0.9631 - 1.863
//
// ch1   6 - 0.3 - 1.03 - 1.822
// ch2   6 - 0.3 - 1.081 - 2.105  (worse with 0.2)
// ch3   6 - 0.3 - 1.114 - 2.584

//ch1    10 - 0.3 - 0.7064 - 0.8426
//ch2    10 - 0.3 - 0.6933 - 0.9455
//ch3    10 - 0.3 - 0.7982 - 1.269

//ch1    20 - 0.3 - 0.2383 - -0.05741
//ch2    20 - 0.3 - 0.2427 - -0.0391
//ch3    20 - 0.3 - 0.2624 - 0.07253

//ch1   40 - 0.3 - 0.1667 - -0.2613
//ch2   40 - 0.3 - 0.1743 - -0.2757
//ch3   40 - 0.3 - 0.1571 - -0.2486

//DELAY-ATT-SIGMA-MEAN 

      for(int i=0; i<iChannels; i++){
//         in[auxCounter] >> time >> wf0;
//         inTwo[auxCounter] >> time >> wf1;

         in[auxCounter] >> time;
         in[auxCounter].ignore();  in[auxCounter] >> wf0;

         inTwo[auxCounter] >> time;
         inTwo[auxCounter].ignore();  inTwo[auxCounter] >> wf1;

         inThree[auxCounter] >> time;
         inThree[auxCounter].ignore();  inThree[auxCounter] >> wf2;

         inFour[auxCounter] >> time;
         inFour[auxCounter].ignore();  inFour[auxCounter] >> wf3;

        h_WF1->SetBinContent(i,wf0);
        h_WF2->SetBinContent(i,wf1);  
        h_WF3->SetBinContent(i,wf2);
        h_WF4->SetBinContent(i,wf3);
        h_Pers->Fill(i*bin_width,wf0);
        h_PersTwo->Fill(i*bin_width,wf1);
        h_PersThree->Fill(i*bin_width,wf2);
        h_PersFour->Fill(i*bin_width,wf3);

        fWaveForm0.data[i]=wf0;
        fWaveForm1.data[i]=wf1;
        fWaveForm2.data[i]=wf2;
        fWaveForm3.data[i]=wf3; 
        
        double_t attWf0=att*wf0;
        double_t attWf1=att*wf1;
        double_t attWf2=att*wf2;
        double_t attWf3=att*wf3;

        h_WFatt1->SetBinContent(i-delay,attWf0);
        h_WFatt2->SetBinContent(i-delay,attWf1);
        h_WFatt3->SetBinContent(i-delay,attWf2);
        h_WFatt4->SetBinContent(i-delay,attWf3);
      }

      for(int i=0; i<iChannels; i++){
        //cout << wf0 << " ***" << attWf0 << " ***" << h_WFatt1->GetBinContent(i)  <<endl;

        double_t CFDwf1=h_WF1->GetBinContent(i) - h_WFatt1->GetBinContent(i);
        double_t CFDwf2=h_WF2->GetBinContent(i) - h_WFatt2->GetBinContent(i);
        double_t CFDwf3=h_WF3->GetBinContent(i) - h_WFatt3->GetBinContent(i );
        double_t CFDwf4=h_WF4->GetBinContent(i) - h_WFatt4->GetBinContent(i );

        h_WFCFD1->SetBinContent(i, CFDwf1);
        h_WFCFD2->SetBinContent(i, CFDwf2);
        h_WFCFD3->SetBinContent(i, CFDwf3);
        h_WFCFD4->SetBinContent(i, CFDwf4);

        fCFDWaveForm0.data[i]=h_WFCFD1->GetBinContent(i);
        fCFDWaveForm1.data[i]=h_WFCFD2->GetBinContent(i);
        fCFDWaveForm2.data[i]=h_WFCFD3->GetBinContent(i);
        fCFDWaveForm3.data[i]=h_WFCFD4->GetBinContent(i);        
        
      }

        //c8->cd(1);
        //h_WFatt1->Draw();
      //START DRAW the wfs
            TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);

      if (iFile<50){
            c0->cd(1);
            h_WF1->Draw();
            h_WFatt1->Draw("same");
            h_WFCFD1->Draw("same");
            //gPad->Modified(); 
            //gPad->Update();                 
            c0->cd(2); 
            h_WF2->Draw();
            h_WFatt2->Draw("same");
            h_WFCFD2->Draw("same");
            //gPad->Modified(); 
            //gPad->Update();     
            c0->cd(3);
            h_WF3->Draw();
            h_WFatt3->Draw("same");
            h_WFCFD3->Draw("same");
            //gPad->Modified(); 
            //gPad->Update();     
            c0->cd(4);
            h_WF4->Draw(); 
            h_WFatt4->Draw("same");
            h_WFCFD4->Draw("same");
            gPad->Modified(); 
            gPad->Update();        
            leg->AddEntry(h_WF4,"Original","l");
            leg->AddEntry(h_WFatt4,"Attenuated","l");
            leg->AddEntry(h_WFCFD4,"CFD","l");
            leg->Draw();
      }
      //END DRAW

        //first
        pPulse[0]->fadqTime=adqTime0;
        fBase = GetBaseLine(pfWaveForm[0], iBLfrom, iBLto);
        pPulse[0]->fBase;
        SubtractBaseLine(pfWaveForm[0], pfWaveFormBL[0], fBase );
        SubtractBaseLine(pfCFDWaveForm[0], pfCFDWaveFormBL[0], fBase );
//         InvertWaveForm(pfWaveFormBL[0], pfWaveFormBL[0]);
        pPulse[0]->fMaxBin = GetPeakPosition(pfWaveFormBL[0], iPULSEfrom , iPULSEto );
        pPulse[0]->fMax = GetPeak(pfWaveFormBL[0], iPULSEfrom , iPULSEto );
        pPulse[0]->fInt = GetIntegral(pfWaveFormBL[0], iPULSEfrom , iPULSEto , kFALSE);
        pPulse[0]->fWidth = GetTdcWidth(pfWaveFormBL[0], iPULSEfrom , iPULSEto ,0.10);
        pPulse[0]->fRCharge = GetRCharge(pfWaveFormBL[0], iPULSEfrom , iPULSEto );
        pPulse[0]->fT0_30 = GetFrontThresholdPosition(pfWaveFormBL[0], iPULSEfrom , iPULSEto , 0.3);
        pPulse[0]->fT0CFD = GetZeroCross(pfCFDWaveFormBL[0], iPULSEfrom , iPULSEto);
        
        //second
        pPulse[1]->fadqTime=adqTime1;
        fBase = GetBaseLine(pfWaveForm[1], iBLfrom, iBLto);
        pPulse[1]->fBase;
        SubtractBaseLine(pfWaveForm[1], pfWaveFormBL[1], fBase );
        SubtractBaseLine(pfCFDWaveForm[1], pfCFDWaveFormBL[1], fBase);
//         InvertWaveForm(pfWaveFormBL[1], pfWaveFormBL[1]);
        pPulse[1]->fMaxBin = GetPeakPosition(pfWaveFormBL[1], iPULSEfrom , iPULSEto );
        pPulse[1]->fMax = GetPeak(pfWaveFormBL[1], iPULSEfrom , iPULSEto );
        pPulse[1]->fInt = GetIntegral(pfWaveFormBL[1], iPULSEfrom , iPULSEto , kFALSE);
        pPulse[1]->fWidth = GetTdcWidth(pfWaveFormBL[1], iPULSEfrom , iPULSEto ,0.10);
        pPulse[1]->fRCharge = GetRCharge(pfWaveFormBL[1], iPULSEfrom , iPULSEto );
        pPulse[1]->fT0_30 = GetFrontThresholdPosition(pfWaveFormBL[1], iPULSEfrom , iPULSEto , 0.3);  
        pPulse[1]->fT0CFD = GetZeroCross(pfCFDWaveFormBL[1], iPULSEfrom , iPULSEto);
        
        //third
        pPulse[2]->fadqTime=adqTime2;
        fBase = GetBaseLine(pfWaveForm[2], iBLfrom, iBLto);
        pPulse[2]->fBase;
        SubtractBaseLine(pfWaveForm[2], pfWaveFormBL[2], fBase );
        SubtractBaseLine(pfCFDWaveForm[2], pfCFDWaveFormBL[2], fBase );
//         InvertWaveForm(pfWaveFormBL[2], pfWaveFormBL[2]);
        pPulse[2]->fMaxBin = GetPeakPosition(pfWaveFormBL[2], iPULSEfrom , iPULSEto );
        pPulse[2]->fMax = GetPeak(pfWaveFormBL[2], iPULSEfrom , iPULSEto );
        pPulse[2]->fInt = GetIntegral(pfWaveFormBL[2], iPULSEfrom , iPULSEto , kFALSE);
        pPulse[2]->fWidth = GetTdcWidth(pfWaveFormBL[2], iPULSEfrom , iPULSEto ,0.10);
        pPulse[2]->fRCharge = GetRCharge(pfWaveFormBL[2], iPULSEfrom , iPULSEto );
        pPulse[2]->fT0_30 = GetFrontThresholdPosition(pfWaveFormBL[2], iPULSEfrom , iPULSEto , 0.3);
        pPulse[2]->fT0CFD = GetZeroCross(pfCFDWaveFormBL[2], iPULSEfrom , iPULSEto);
        
        //fourth
        pPulse[3]->fadqTime=adqTime3;
        fBase = GetBaseLineT(pfWaveForm[3], iBLfromT, iBLtoT);
        pPulse[3]->fBase;
        SubtractBaseLine(pfWaveForm[3], pfWaveFormBL[3], fBase );
        InvertWaveForm(pfWaveFormBL[3], pfWaveFormBL[3]);
        pPulse[3]->fMaxBin = GetPeakPosition(pfWaveFormBL[3], iPULSEfromT , iPULSEtoT );
        pPulse[3]->fMax = GetPeak(pfWaveFormBL[3], iPULSEfromT , iPULSEtoT );
        pPulse[3]->fInt = GetIntegral(pfWaveFormBL[3], iPULSEfromT , iPULSEtoT , kFALSE);
        pPulse[3]->fWidth = GetTdcWidth(pfWaveFormBL[3], iPULSEfromT , iPULSEtoT ,0.50);
        pPulse[3]->fRCharge = GetRCharge(pfWaveFormBL[3], iPULSEfromT , iPULSEtoT );
        pPulse[3]->fT0_30 = GetFrontThresholdPosition(pfWaveForm[3], iPULSEfromT , iPULSEtoT , 0.3);

        tTree.Fill();
      if(iFile<5){
          cout<< "Max=" << GetPeak(pfWaveFormBL[3], iPULSEfrom , iPULSEto ) 
          << "arrivalTime=" << 0.5 * GetFrontThresholdPosition(pfWaveFormBL[1], iPULSEfrom , iPULSEto , 0.3) <<endl;
      }        
//         cout << "---------wf no" << iFile << endl;
        in[auxCounter].close();     
        inTwo[auxCounter].close();     
        inThree[auxCounter].close();
        inFour[auxCounter].close();
   }
   
     //TCanvas *c1 = new TCanvas("c1","WaveForm",10,10,1400,700);
     //gPad->SetGrid();
     //h_WF1->Draw("same"); 
     //h_WF2->Draw(); 
//     TCanvas *c2 = new TCanvas("c2","WaveForm",10,10,1400,700);
//     h_Pers->Draw("colz");

  TString outputFile = "Tree"+mainfile+".root";
//   TString outputFile = output+".root";
  TFile*  file = new TFile(outputFile,"recreate");
  tTree.Write();

  //h_WF1->Write();
  //h_WF2->Write();
  //h_WF3->Write();
  //h_WF4->Write();
  h_Pers->Write();
  h_PersTwo->Write();
  h_PersThree->Write();
  h_PersFour->Write();
  h_WFCFD1->Write();

  //h_WFatt1->Write();
  file->Close(); 
}

