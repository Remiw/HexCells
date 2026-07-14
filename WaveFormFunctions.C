// Notes for the correct implementation of this funtions:
// 1. Always check if the variables bins_per_record and bin_width are in accordance with your data
// 2. For the funtion GetBaseLine look at the region and choose the correct values for bl_from and bl_to

#include "TH1.h"
#include "TF1.h"

void jumpToLine(std::istream& os, int n){
	// Clear error flags, just in case.
	os.clear();	
	// Start reading from the beginning of the file.
	os.seekg(0, std::ios::beg);
	// Skip to line n.
	for (int i = 1; i < n; ++i)
	{
		os.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
}


class Pulse : public TObject {
public:
  float fadqTime;
  float fMax;
  float fInt;
  float fWidth;
  float fT0_30;
  float fBase;  
  float fMaxBin;
  float fRCharge;
  float fT0CFD;
 
  Pulse() { };
  ClassDef(Pulse,1);
};
 
const Int_t nChannels = 4;
//two channels
int const bins_per_record = 1003-1; // number of points per waveform, the value of the Record Length:3125
float const bin_width = 0.5; // here is the Sample Interval, better if you use nanoseconds (ns) 3.20000000e-10

//one channel
 //int const bins_per_record = 2002-1; // number of points per waveform, the value of the Record Length:3125
//float const bin_width = 0.5/2.; // here is the Sample Interval, better if you use nanoseconds (ns) 3.20000000e-10


float const bl_from = -0.01; // choose the MIN value for the baseline region
float const bl_to = 0.01; // choose the MAX value for the baseline region

typedef struct {
  std::array<float,bins_per_record> data;
} OSC_Record;

int CopyWaveForm(OSC_Record* pWFin, OSC_Record* pWFout) {
  for (int i=0; i<bins_per_record; i++) pWFout->data[i] = pWFin->data[i];
  return 0;
}

void InvertWaveForm(OSC_Record* pWFin, OSC_Record* pWFout) {
  for (int i=0; i<bins_per_record; i++) pWFout->data[i] = -pWFin->data[i];
}

int SumUpWaveForms(OSC_Record* pWFin_1, OSC_Record* pWFin_2, OSC_Record* pWFout) {
  for (int i=0; i<bins_per_record; i++) pWFout->data[i] = pWFin_1->data[i] + pWFin_2->data[i];
  return 0;
}


int CheckChannelRange(int base_from, int base_to) {
  if(base_from < 0)    { cout << "unexpected first bin for base " << base_from << endl; return -1; }
  if(base_from > bins_per_record) { cout << "unexpected first bin for base " << base_from << endl; return -1; }
  if(base_to > bins_per_record)   { cout << "unexpected last bin for base " << base_to << endl; return -1; }
  if(base_to < 0)      { cout << "unexpected last bin for base " << base_to << endl; return -1; }
  if(base_from > base_to) { 
    cout << "unexpected interval for base: " << base_from << " .. " << base_to << endl; 
    return -1;
  } 
  return 0;
}

float GetBaseLine(OSC_Record* pWFin, int base_from, int base_to) {

  if(CheckChannelRange(base_from, base_to) != 0) return -100000.;

  int iBaseFrom = base_from;
  int iBaseTo = base_to;

  TH1F  *hTmp = new TH1F("hTmp","Pedestal",100,bl_from,bl_to);
  int i;
  float base = 0.;
  for (i=iBaseFrom; i<=iBaseTo; i++) { hTmp->Fill(pWFin->data[i]); }

  Double_t iBin = hTmp->GetMean();
  hTmp->Delete();

  return iBin;

}

void SubtractBaseLine(OSC_Record* pWFin, OSC_Record* pWFout, float base) {
  for (int i=0; i<bins_per_record; i++) pWFout->data[i] = pWFin->data[i] - base;
  return;
}

float SubtractBaseLine(OSC_Record* pWFin, OSC_Record* pWFout, int base_from, int base_to) {

  float base = GetBaseLine(pWFin, base_from, base_to);
  for (int i=0; i<bins_per_record; i++) pWFout->data[i] = pWFin->data[i] - base;
  return base;
}


float GetIntegral(OSC_Record* rec, int from, int to, bool only_positive = kFALSE) {

  // if only_positive is true, negative part of the pulse is not added to the sum
  // in this case one has to expect slight additional baseline shift

  if(CheckChannelRange(from, to) != 0) return -1000000.;

  float sum = 0.;
  for (int i=from; i<=to; i++) { 
    sum = sum + rec->data[i];
    if(only_positive && (rec->data[i] < 0)) sum = sum - rec->data[i];
  }
  return sum;

}

float GetRCharge(OSC_Record* rec, int from, int to) {
//If the integral is correctly computed this value will be zero

  if(CheckChannelRange(from, to) != 0) return -1000000.;

  float sum = 0.;
  for (int i=from; i<=to; i++) { 
    sum = sum + rec->data[i];
    sum = sum - rec->data[i];
  }
  return sum;

}


int GetPeakPosition(OSC_Record* rec, int from, int to, int sign = 1) {

  int iFrom = from;
  int iTo = to;

  int iMinPosition = 0;
  int iMaxPosition = 0;

  float max = -5000.;
  float min =  5000.;

  for (int i=iFrom; i<=iTo; i++) { 
    if (max<rec->data[i]) { max = rec->data[i]; iMaxPosition = i; }
    if (min>rec->data[i]) { min = rec->data[i]; iMinPosition = i; }
  }

  if(sign < 0) return iMinPosition;
  return iMaxPosition; 

}


float GetPeak(OSC_Record* rec, int from, int to, int sign = 1) {

  int iPeak = GetPeakPosition(rec, from, to, sign);

  return rec->data[iPeak];

}


float GetFrontThresholdPosition(OSC_Record *rec, int from, int to, float threshold) {

  // WaveForm is expected to have a positive pulse and a substracted baseline
  // before calling use SubstractBaseline() and InvertWaveForm() if needed

  // threshold is a fraction from the peak value, should be from 0.05 to 0.95

  if((threshold < 0.05) || (threshold > 0.95)) {
    cout << "GetFrontThresholdPosition() error. Unexpected threshold value " << threshold << endl;
    return -1.;
  }

  int   iPeak = GetPeakPosition(rec, from, to);
  float fThreshold = rec->data[iPeak]*threshold;

  float front;
  int   i, iAbove, iBelow;

  iBelow = from;
  for (i=iPeak; i>=from; i--) {
    if ((rec->data[i] < fThreshold) && (rec->data[i+1] >= fThreshold)) {
      iBelow = i;
      break;
    }
  }
  iAbove = iPeak;
  for (i=iPeak; i>=from; i--) {
    if ((rec->data[i-1] <= fThreshold) && (rec->data[i] > fThreshold)) {
      iAbove = i;
      break;
    }
  }
  front = iBelow + (iAbove-iBelow)*(fThreshold - rec->data[iBelow])/(rec->data[iAbove]-rec->data[iBelow]);

  return front;
}

float GetTailThresholdPosition(OSC_Record *rec, int from, int to, float threshold) {

  // WaveForm is expected to have a positive pulse and a substracted baseline
  // before calling use SubstractBaseline() and InvertWaveform() if needed 
 
  // threshold is a fraction from the peak value, should be from 0.05 to 0.95


  if((threshold < 0.05) || (threshold > 0.95)) {
    cout << "GetTailThresholdPosition() error. Unexpected threshold value " << threshold << endl;
    return -1000000.;
  }

  int   iPeak = GetPeakPosition(rec, from, to);
  float fThreshold = rec->data[iPeak]*threshold;

  float tail; 
  int   i, iAbove, iBelow;

  iBelow = to;
  for (i=iPeak; i<=to; i++) { 
    if ((rec->data[i] < fThreshold) && (rec->data[i-1] >= fThreshold)) {
      iBelow = i;
      break;
    }
  }  
  iAbove = iPeak;
  for (i=iPeak; i<=to; i++) { 
    if ((rec->data[i+1] <= fThreshold) && (rec->data[i] > fThreshold)) {
      iAbove = i;
      break;
    }
  }  
  tail = iAbove + (iBelow-iAbove)*(rec->data[iAbove]-fThreshold)/(rec->data[iAbove]-rec->data[iBelow]);
  
  return tail;
}

float GetTdcWidth(OSC_Record *rec, int from, int to, float threshold) {

  float fMax = GetPeak(rec, from, to);
  float fFraction = threshold/fMax;
  if(fFraction > 0.95) return -1000.;
  if(fFraction < 0.05) fFraction = 0.05;

  float fStart = GetFrontThresholdPosition(rec, from, to, fFraction);
  float fStop = GetTailThresholdPosition(rec, from, to, fFraction);

  return fStop - fStart;

}

//now for TRIGGER PULSE (BL RANGES CHANGE)
        //fourth
        //pPulse[3]->fadqTime=adqTime3;
        //fBase = GetBaseLine(pfWaveForm[3], iBLfromT, iBLtoT);
        //pPulse[3]->fBase;
        //SubtractBaseLine(pfWaveForm[3], pfWaveFormBL[3], fBase );
        //InvertWaveForm(pfWaveFormBL[3], pfWaveFormBL[3]);
        //pPulse[3]->fMaxBin = GetPeakPosition(pfWaveFormBL[3], iPULSEfromT , iPULSEtoT );
        //pPulse[3]->fMax = GetPeak(pfWaveFormBL[3], iPULSEfromT , iPULSEtoT );
        //pPulse[3]->fInt = GetIntegral(pfWaveFormBL[3], iPULSEfromT , iPULSEtoT , kFALSE);
        //pPulse[3]->fWidth = GetTdcWidth(pfWaveFormBL[3], iPULSEfromT , iPULSEtoT ,0.50);
        //pPulse[3]->fRCharge = GetRCharge(pfWaveFormBL[3], iPULSEfromT , iPULSEtoT );
        //pPulse[3]->fT0_30 = GetFrontThresholdPosition(pfWaveFormBL[3], iPULSEfromT , iPULSEtoT , 0.3);

float const bl_fromT = 0.05; // choose the MIN value for the baseline region
float const bl_toT = 159; // choose the MAX value for the baseline region

float GetBaseLineT(OSC_Record* pWFin, int base_from, int base_to) {

  if(CheckChannelRange(base_from, base_to) != 0) return -100000.;

  int iBaseFrom = base_from;
  int iBaseTo = base_to;

  TH1F  *hTmp = new TH1F("hTmp","Pedestal",500,bl_fromT,bl_toT);
  int i;
  float base = 0.;
  for (i=iBaseFrom; i<=iBaseTo; i++) { hTmp->Fill(pWFin->data[i]); }

  Double_t iBin = hTmp->GetMean();
  hTmp->Delete();

  return iBin;

}

int GetZeroCross(OSC_Record* rec, int from, int to) {
  int iFrom = from;
  int iTo = to - 10 ;//Deberíamos asegurarnos de que no se salga del array
  int ZeroCrossPosition = -1000; 

  for (int i = iFrom; i <= iTo; i++) {
    if (rec->data[i] < 0 && rec->data[i + 1] > 0
          && rec->data[i + 2] > rec->data[i+1] 
          && rec->data[i + 3] > rec->data[i + 2]
          && rec->data[i + 4] > rec->data[i + 3] 
          && rec->data[i + 5] > rec->data[i + 4]
          && rec->data[i + 6] >  rec->data[i + 5] 
          && rec->data[i + 7] > rec->data[i + 6]
          && rec->data[i + 8] > rec->data[i + 7] 
          && rec->data[i + 9] > rec->data[i + 8] 
          && rec->data[i + 10] > rec->data[i + 9]) {
          ZeroCrossPosition = i + 1;  // Devuelve el número de bin donde ocurre el cruce
          return ZeroCrossPosition;
          //cout << "El cruce de 0 ocurre en el bin: " << ZeroCrossPosition << endl;
    }
  }

}


// nt GetPeakPosition(OSC_Record* rec, int from, int to, int sign = 1) {

//   int iFrom = from;
//   int iTo = to;

//   int iMinPosition = 0;
//   int iMaxPosition = 0;

//   float max = -5000.;
//   float min =  5000.;

//   for (int i=iFrom; i<=iTo; i++) { 
//     if (max<rec->data[i]) { max = rec->data[i]; iMaxPosition = i; }
//     if (min>rec->data[i]) { min = rec->data[i]; iMinPosition = i; }
//   }

//   if(sign < 0) return iMinPosition;
//   return iMaxPosition; 

// }