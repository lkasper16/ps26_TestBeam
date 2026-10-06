#define trdclass_ps26_cxx
#include "trdclass_ps26.h"
#include "PlotLib.C"
#include "GNN/gnn_model.h"
#include "GNN/gnn_model.cpp"
#include "GNN/toGraph.cpp"

#define NPRT 1000
#define USE_TRK
#define MAX_PRINT 10
//
#define USE_GNN  1
#define USE_FIT  1
////#define USE_CLUST 1
#define USE_PULSE 0
//
////#define USE_125_RAW
//#define USE_250_PULSE
#define MAX_CLUST 500
#define MAX_NODES 100
#define USE_MAXPOS 1
//
#define SAVE_TRACK_HITS
#define SAVE_PDF
//#define WRITE_CSV
#define DEBUG 0
#define USE_TRD_EXT_TRACK 0
//
//#define GAIN_CALIB
//
int timeSwitchRun1=6155; //RunNum where timing window was changed to 250
int timeSwitchRun2=6210; //RunNum where timing window was changed back to 200
double firstTimeWin=200.;
double secondTimeWin=250.;
int noRadList[] = {8265,8266,8267,8268,8269,8270};
int listSize = sizeof(noRadList) / sizeof(noRadList[0]);
int argonRunStop = 8223;
int argonRunStart = 8271;

//-- For single evt clustering display, uncomment BOTH:
//#define SHOW_EVTbyEVT
//#define SHOW_EVT_DISPLAY

void WriteToCSV(std::ofstream &csvFile, float v1, float v2, float v3, float v4, float v5, float v6, float v7, float v8, float v9, float v10) {
  csvFile<<v1<<","<<v2<<","<<v3<<","<<v4<<","<<v5<<","<<v6<<","<<v7<<","<<v8<<","<<v9<<","<<v10<<std::endl;
}

//===================================
//      TRD DAQ Channel Mapping
//===================================

// -- Triple-GEMTRD mapping --
int Get3GEMChan(int ch, int slot, int runNum) {
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  float dchan = cardChannel+cardNumber*24+(slot-3)*72.;
  if ((slot==6 && ch>23) || slot==7 || slot==8 || (slot==9 && ch<48)) {
    #ifdef GAIN_CALIB
      if (dchan-240.==0 || dchan-240.==18 || dchan-240.==47 || dchan-240.==78 || dchan-240.==97 || dchan-240.>222 || dchan-240.<24 || dchan-240.==135 || dchan-240.==119) { return -1; } 
    #else
      if (dchan-240.==227. || dchan-240.==220. || dchan-240.==221. || dchan-240.==226. || dchan-240.==224. ||dchan-240.==16. ||  dchan-240.==223. || dchan-240.==117. || dchan-240.==119. || dchan-240.==120. || dchan-240.==131. || dchan-240.==116. || dchan-240.==129. || dchan-240.==135. || dchan-240.==138. || dchan-240.==130. || dchan-240.==180. || dchan-240.==184. || dchan-240.==181. || dchan-240.==183. || dchan-240.==217. || dchan-240.==219. || dchan-240.==229. || dchan-240.==231. || dchan-240.==191. || dchan-240.==78. || dchan-240.<24. || dchan-240.==216. || dchan-240.==218.) { return -1; }
    #endif
    else { return dchan-240.; }
  }
  return -1;
}

// -- MMG1TRD mapping --
int GetMMG1Chan(int ch, int slot, int runNum) {
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  float dchan = cardChannel+cardNumber*24+(slot-3)*72.;
  if (slot==3 || slot==4 || slot==5 || (slot==6 && ch<24)) {
    #ifdef GAIN_CALIB
      if (dchan==16. || dchan==1. || dchan==23. || dchan==45. || dchan==98. || dchan==99. || dchan==101. || (dchan>=119. && dchan<=130.) || dchan==149. || dchan==151. || (dchan>=224. && dchan<=236.)) { return -1; }
    #else
      if (dchan==16. || dchan==1. || dchan==23. || dchan==45. || dchan==98. || dchan==99. || dchan==101. || (dchan>=119. && dchan<=130.) || dchan==15. || dchan==149. || dchan==151. || (dchan>=224. && dchan<=236.)) { return -1; }
    #endif
    else { return dchan; }
  }
  return -1;
}

// -- uRWELL+GEM mapping --
int GetRWELLXChan(int ch, int slot, int runNum) {
  if (slot>10) slot=(slot-2);
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  cardChannel = 23-cardChannel; //new style electronics have inverted polarity
  float dchan = cardChannel+cardNumber*24+(slot-3)*72.;
    if (slot==10 || (slot==11 && ch<48)) {
        if ( (dchan-504.)==97. || (dchan-504.)==98. || (dchan-504.)==99. || (dchan-504.)==100. ) { return -1; }
        else { return dchan-504.; }
    }
  return -1;
}

int GetRWELLYChan(int ch, int slot, int runNum) {
  if (slot>10) slot=(slot-2);
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  cardChannel = 23-cardChannel; //new style electronics have inverted polarity
  float dchan = cardChannel+cardNumber*24+(slot-3)*72.;
    if ((slot==11 && ch>47) || (slot==12 && ch<48)) {
      /*if ((dchan-600.)==42.) { return -1; } else {*/ return dchan-600.;//} //-624.
    }
  return -1;
}

//============ END DAQ Channel Mapping ============

//============ ADC Gain Calibrations ============
float Get3GEMCalib(float amp, int ch, int runNum) {
  float calibrated_amp = amp;
  if ((runNum>=6295 && runNum<=6297) || (runNum>=6315 && runNum<=6317)) { //LG Cu
    double gemCoefficients[] ={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.664062,0.664062,0.712891,0.820312,0.849609,0.800781,0.78125,0.74707,0.756836,0.805664,0.732422,0.712891,0.717773,0.800781,0.766602,0.712891,0.722656,0.74707,0.703125,0.708008,0.761719,1,1,1,1,1,0.688477,1,1,0.874023,1,0.727539,0.498047,1,0.81543,1,0.581055,1,1,1,1,1,1,1,1,0.888672,0.878906,0.810547,0.854492,0.864258,0.874023,0.927734,0.854492,0.800781,0.800781,0.81543,0.820312,0.800781,0.825195,0.834961,0.869141,0.839844,0.874023,0.844727,0.893555,0.883789,0.874023,0.878906,0.795898,0.844727,0.844727,0.859375,0.869141,0.849609,0.81543,0.820312,0.888672,0.864258,0.795898,0.839844,0.908203,0.864258,0.825195,0.844727,0.805664,0.786133,0.849609,0.786133,0.869141,0.854492,0.90332,1,0.849609,0.834961,0.874023,0.898438,0.864258,0.859375,0.893555,0.820312,0.898438,1,0.830078,0.917969,0.893555,0.874023,0.869141,0.830078,0.932617,0.805664,0.893555,0.834961,0.878906,0.859375,0.864258,0.869141,0.888672,0.864258,0.844727,0.869141,0.859375,0.839844,0.893555,0.864258,0.834961,0.878906,0.859375,0.820312,0.805664,0.800781,0.839844,0.874023,0.854492,0.81543,0.869141,0.820312,0.864258,0.908203,0.864258,0.810547,0.849609,0.81543,0.805664,0.810547,0.800781,0.805664,0.834961,1,0.800781,0.791016,0.810547,0.825195,1,0.825195,0.805664,1,0.727539,0.78125,0.849609,0.742188,0.742188,0.830078,0.830078,0.908203,0.825195,0.791016,1,0.737305,0.776367,0.766602,0.683594,0.756836,0.776367,0.732422,0.839844,0.874023,0.810547,0.766602,0.756836,0.737305,0.810547,0.771484,0.722656,0.742188,0.737305,1,1,1,0.839844,1,0.717773,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    calibrated_amp = amp/gemCoefficients[ch];
  }
  else if ((runNum>=6301 && runNum<=6304) || (runNum>=6318 && runNum<=6320)) { //HG Cu
    double gemCoefficients[240] = {0.942383,1,1,0.893555,0.810547,0.581055,0.810547,0.834961,0.874023,0.888672,0.844727,0.859375,0.830078,0.834961,0.756836,0.620117,1,0.898438,0.90332,0.874023,0.854492,0.869141,0.859375,0.869141,1,0.859375,0.830078,0.883789,0.888672,0.917969,0.878906,0.913086,0.9375,0.864258,0.844727,0.90332,0.888672,0.908203,0.961914,0.942383,0.874023,0.874023,0.874023,0.893555,0.90332,0.839844,0.830078,0.878906,1,1,1,1,1,0.839844,1,1,0.97168,1,0.795898,0.541992,1,0.922852,1,0.634766,1,1,1,1,1,1,1,1,0.913086,0.952148,0.849609,0.913086,0.957031,0.90332,0.952148,0.922852,0.893555,0.869141,0.893555,0.869141,0.854492,0.893555,0.888672,0.913086,0.893555,0.9375,0.898438,0.947266,0.927734,0.913086,0.932617,0.839844,0.888672,0.883789,0.893555,0.908203,0.898438,0.878906,0.883789,0.942383,0.917969,0.839844,0.888672,0.957031,0.908203,0.883789,0.922852,0.893555,0.859375,0.913086,0.849609,0.908203,0.869141,0.917969,1,0.883789,0.898438,0.927734,0.922852,0.893555,0.888672,0.913086,0.849609,0.927734,1,0.878906,0.957031,0.942383,0.917969,0.922852,0.859375,0.966797,0.898438,0.961914,0.883789,0.947266,0.913086,0.913086,0.922852,0.947266,0.9375,0.917969,0.90332,0.917969,0.878906,0.9375,0.893555,0.869141,0.908203,0.90332,0.869141,0.878906,0.869141,0.893555,0.927734,0.893555,0.859375,0.942383,0.883789,0.917969,0.976562,0.927734,0.893555,0.90332,0.874023,0.844727,0.864258,0.849609,0.864258,0.878906,1,0.908203,0.922852,0.932617,0.9375,1,0.878906,0.849609,1,0.922852,0.854492,0.922852,0.908203,0.898438,0.917969,0.874023,0.942383,0.859375,0.898438,1,0.839844,0.913086,0.849609,0.839844,0.927734,0.952148,0.917969,0.893555,0.908203,0.898438,0.908203,0.878906,0.859375,0.898438,0.849609,0.913086,0.839844,0.839844,1,1,1,0.898438,1,0.97168,0.927734,0.893555,0.883789,0.908203,0.952148,0.932617,1,0.927734,0.869141,0.893555,0.864258,0.913086,0.874023,0.878906,0.888672,0.864258,0.874023,0.864258,0.874023,0.888672,0.917969,0.97168};
    calibrated_amp = amp/gemCoefficients[ch];
  }
    else if ((runNum>=6381 && runNum<=6385) || (runNum>=6392 && runNum<=6394)) { //LG Al
    double gemCoefficients[240] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.771484,0.722656,0.791016,0.78125,0.595703,0.644531,0.810547,0.844727,0.81543,0.786133,0.751953,0.703125,0.74707,0.708008,0.776367,0.849609,1,1,1,1,1,0.805664,1,1,1,1,0.683594,0.498047,1,0.859375,1,0.600586,1,1,1,1,1,1,1,1,0.874023,0.874023,0.795898,0.854492,1,0.859375,0.883789,0.849609,0.854492,0.810547,0.839844,0.820312,0.791016,0.834961,0.839844,0.854492,0.81543,0.859375,0.839844,0.878906,0.844727,0.844727,0.874023,0.830078,0.81543,0.839844,0.844727,0.859375,0.844727,0.820312,0.844727,0.874023,0.878906,0.805664,0.830078,0.908203,0.864258,0.820312,0.810547,0.805664,0.776367,0.844727,0.776367,0.834961,0.800781,0.864258,1,0.854492,0.830078,0.854492,0.869141,0.839844,0.834961,0.869141,1,1,1,0.820312,0.893555,0.878906,0.854492,0.844727,0.795898,0.908203,0.849609,0.90332,0.830078,0.888672,0.849609,0.869141,0.888672,0.913086,0.878906,0.854492,0.844727,0.854492,0.849609,0.908203,0.849609,0.839844,0.864258,0.849609,0.830078,0.825195,0.81543,0.854492,0.878906,0.869141,0.834961,0.874023,0.810547,0.859375,0.917969,0.90332,0.859375,0.849609,0.839844,0.805664,0.849609,0.81543,0.800781,0.810547,1,0.825195,0.859375,0.869141,0.825195,1,0.81543,0.81543,1,0.893555,0.795898,0.81543,0.825195,0.830078,0.830078,0.737305,0.854492,0.786133,0.820312,1,0.742188,0.825195,0.737305,0.78125,0.834961,0.791016,0.810547,0.78125,0.776367,0.795898,0.830078,0.78125,0.81543,0.834961,0.776367,0.810547,0.776367,0.737305,1,1,1,0.756836,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    calibrated_amp = amp/gemCoefficients[ch];
  }
  else if ((runNum>=6389 && runNum<=6391) || (runNum>=6385 && runNum<=6388)) { //HG Al
    double gemCoefficients[240] = {1,1,1,1,1,1,1,1,1,1,1,1,1,0.820312,0.742188,0.625,1,0.888672,0.90332,0.859375,0.834961,0.864258,0.844727,0.839844,1,0.874023,0.854492,0.878906,0.859375,0.90332,0.864258,0.888672,0.927734,0.869141,0.830078,0.878906,0.859375,0.883789,0.932617,0.917969,0.874023,0.898438,0.874023,0.878906,0.913086,0.849609,0.839844,0.883789,1,1,1,1,1,0.830078,1,1,1,1,0.776367,0.537109,1,0.913086,1,0.634766,1,1,1,1,1,1,1,1,0.917969,0.952148,0.844727,0.917969,1,0.908203,0.957031,0.9375,0.898438,0.864258,0.898438,0.874023,0.854492,0.893555,0.893555,0.927734,0.898438,0.932617,0.893555,0.9375,0.932617,0.913086,0.932617,0.893555,0.888672,0.878906,0.893555,0.913086,0.893555,0.859375,0.874023,0.9375,0.917969,0.839844,0.878906,0.947266,0.893555,0.864258,0.913086,0.883789,0.839844,0.893555,0.834961,0.90332,0.869141,0.917969,1,0.888672,0.878906,0.908203,0.90332,0.883789,0.888672,0.922852,1,1,1,0.888672,0.961914,0.932617,0.908203,0.917969,0.869141,0.966797,0.893555,0.957031,0.878906,0.932617,0.90332,0.908203,0.917969,0.957031,0.9375,0.913086,0.908203,0.913086,0.888672,0.947266,0.898438,0.864258,0.898438,0.90332,0.874023,0.874023,0.864258,0.888672,0.9375,0.913086,0.864258,0.9375,0.888672,0.922852,0.97168,0.922852,0.898438,0.908203,0.874023,0.844727,0.874023,0.854492,0.864258,0.878906,1,0.898438,0.917969,0.922852,0.9375,1,0.898438,0.864258,1,0.947266,0.874023,0.927734,0.917969,0.908203,0.942383,0.893555,0.957031,0.869141,0.90332,1,0.854492,0.917969,0.859375,0.849609,0.922852,0.957031,0.927734,0.888672,0.913086,0.913086,0.913086,0.874023,0.869141,0.913086,0.849609,0.893555,0.825195,0.834961,1,1,1,0.917969,1,0.961914,0.917969,0.849609,0.883789,0.908203,0.922852,0.893555,1,0.922852,0.839844,0.859375,0.820312,0.908203,0.869141,0.834961,0.820312,0.839844,0.830078,0.854492,0.874023,1,0.893555,1};
    calibrated_amp = amp/gemCoefficients[ch];
  }
  if (ch<0 || ch>=240) return amp;
  return calibrated_amp;
}

float GetMMGCalib(float amp, int ch, int runNum) {
  float calibrated_amp = amp;
  if ((runNum>=6295 && runNum<=6297) || (runNum>=6315 && runNum<=6317)) { //LG Cu
    double mmgCoefficients[240] = {1,0.81543,1,0.849609,0.917969,0.805664,1,0.917969,1,0.849609,1,0.820312,0.771484,0.771484,0.893555,0.849609,1,0.791016,0.839844,0.849609,0.791016,0.849609,0.795898,0.673828,0.771484,0.722656,0.849609,0.859375,0.708008,0.761719,0.97168,0.854492,0.917969,0.864258,0.869141,0.844727,0.874023,0.869141,0.854492,0.942383,0.927734,0.90332,0.922852,0.834961,0.883789,1,0.90332,1,0.869141,0.859375,0.844727,0.834961,0.898438,0.898438,0.869141,0.9375,0.908203,0.90332,0.932617,0.908203,0.883789,0.874023,0.898438,0.947266,0.878906,0.869141,0.859375,0.834961,0.820312,0.898438,0.917969,1,0.859375,0.869141,0.883789,0.883789,0.878906,0.878906,0.874023,0.859375,0.9375,0.859375,0.859375,0.90332,0.839844,0.849609,0.932617,0.844727,0.859375,0.874023,0.869141,1,1,0.834961,0.825195,0.859375,0.864258,0.888672,1,0.864258,0.854492,0.97168,0.81543,0.825195,0.849609,0.844727,0.795898,0.81543,0.795898,0.839844,0.825195,0.830078,0.81543,0.78125,0.825195,0.927734,0.810547,0.922852,0.786133,0.825195,1,0.859375,1,0.869141,1,0.839844,1,0.869141,1,0.874023,0.820312,0.859375,0.874023,0.810547,0.776367,0.800781,1,0.825195,0.849609,0.805664,0.771484,0.830078,0.78125,0.90332,0.751953,0.795898,0.810547,0.81543,0.810547,0.786133,0.81543,0.805664,0.771484,0.878906,0.908203,0.810547,0.756836,0.756836,0.751953,0.805664,0.78125,0.791016,0.756836,0.849609,0.834961,0.795898,0.805664,0.761719,0.708008,0.742188,0.761719,0.791016,0.883789,0.922852,0.898438,0.771484,0.834961,0.839844,0.693359,0.742188,1,1,0.639648,0.605469,1,1,0.800781,1,1,1,1,1,1,1,0.708008,1,0.732422,0.839844,0.761719,0.737305,0.708008,0.688477,0.751953,0.722656,0.708008,0.74707,0.717773,0.727539,0.805664,0.81543,0.791016,0.708008,0.766602,0.791016,0.761719,1,0.771484,0.742188,0.654297,0.727539,0.74707,0.761719,0.761719,0.727539,1,0.742188,1,0.708008,0.668945,0.712891,0.678711,1,1,0.800781,0.795898,1,0.791016,0.703125,0.708008,0.742188};
    calibrated_amp = amp/mmgCoefficients[ch];
  }
  else if ((runNum>=6301 && runNum<=6304) || (runNum>=6318 && runNum<=6320)) { //HG Cu
    double mmgCoefficients[240] = {1,0.839844,1,0.878906,0.927734,0.854492,1,0.932617,1,0.878906,1,0.913086,0.888672,0.825195,0.957031,0.913086,1,0.854492,0.883789,0.878906,0.820312,0.878906,0.849609,0.732422,0.839844,0.78125,0.883789,0.893555,0.742188,0.776367,0.981445,0.996094,0.976562,0.927734,0.922852,0.898438,0.927734,0.917969,0.908203,0.976562,0.976562,0.97168,0.97168,0.908203,0.9375,1,0.9375,0.991211,0.932617,0.927734,0.908203,0.883789,0.9375,0.947266,0.898438,0.952148,0.922852,0.961914,0.961914,0.9375,0.917969,0.917969,0.9375,0.976562,0.947266,0.932617,0.917969,0.908203,0.888672,0.947266,0.966797,1,0.908203,0.922852,0.922852,0.966797,0.947266,0.913086,0.913086,0.90332,0.976562,0.97168,0.9375,0.952148,0.9375,0.976562,0.976562,0.957031,0.9375,0.9375,0.908203,1,1,0.932617,0.908203,0.961914,0.947266,0.947266,1,0.961914,0.922852,0.981445,0.878906,0.922852,0.947266,0.913086,0.859375,0.917969,0.869141,0.922852,0.908203,0.97168,0.942383,0.839844,0.869141,0.942383,0.888672,0.961914,0.927734,0.869141,1,0.90332,1,0.942383,1,0.898438,1,0.9375,1,0.952148,0.888672,0.932617,0.922852,0.908203,0.898438,0.913086,1,0.942383,0.913086,0.883789,0.888672,0.908203,0.917969,0.932617,0.888672,0.893555,0.917969,0.883789,0.864258,0.90332,0.913086,0.908203,0.908203,0.917969,0.90332,0.888672,0.869141,0.90332,0.927734,0.854492,0.9375,0.898438,0.888672,0.9375,0.913086,0.878906,0.883789,0.864258,0.908203,0.864258,0.830078,0.869141,0.834961,0.869141,0.874023,0.927734,0.888672,0.859375,0.844727,0.932617,1,1,0.81543,0.913086,1,1,0.957031,1,1,1,1,1,1,1,0.90332,1,0.854492,0.932617,0.874023,0.9375,0.883789,0.810547,0.830078,0.942383,0.844727,0.854492,0.859375,0.917969,0.966797,0.97168,0.917969,0.893555,0.908203,0.952148,0.898438,1,0.898438,0.844727,0.874023,0.888672,0.859375,0.874023,0.917969,0.888672,1,0.844727,1,0.854492,0.849609,0.898438,0.893555,0.864258,0.898438,0.888672,0.888672,1,0.883789,0.932617,0.864258,0.874023};
    calibrated_amp = amp/mmgCoefficients[ch];
  }
  else if ((runNum>=6381 && runNum<=6385) || (runNum>=6392 && runNum<=6394)) { //LG Al
    double mmgCoefficients[240] = {1,0.742188,1,0.649414,0.693359,0.74707,1,0.825195,1,0.708008,1,0.820312,0.771484,0.688477,0.693359,0.673828,1,0.737305,0.668945,0.722656,0.727539,0.751953,0.771484,0.74707,0.805664,0.698242,0.810547,0.830078,0.678711,0.712891,0.922852,0.966797,0.947266,0.864258,0.883789,0.839844,0.854492,0.874023,0.849609,0.947266,0.942383,0.942383,0.97168,0.839844,0.913086,1,0.859375,1,0.849609,0.844727,0.849609,0.820312,0.869141,0.888672,0.830078,0.957031,0.927734,0.957031,0.957031,0.913086,0.888672,0.878906,0.908203,0.957031,0.898438,0.864258,0.854492,0.869141,0.820312,0.893555,0.917969,1,0.830078,0.854492,0.869141,0.908203,0.869141,0.849609,0.859375,0.849609,0.90332,0.81543,0.849609,0.888672,0.854492,0.9375,0.966797,0.883789,0.913086,0.913086,0.908203,1,1,0.927734,0.898438,0.883789,0.942383,0.927734,1,0.957031,0.922852,0.976562,0.888672,0.942383,0.961914,0.932617,0.810547,0.854492,0.786133,0.849609,0.844727,0.859375,0.810547,0.771484,0.771484,0.869141,0.791016,0.927734,0.834961,0.805664,1,0.810547,1,0.864258,1,0.849609,1,0.869141,1,0.859375,1,0.712891,0.805664,0.864258,0.830078,0.869141,1,0.776367,0.830078,0.766602,0.786133,0.825195,0.805664,0.727539,0.795898,0.722656,0.786133,0.708008,0.722656,0.800781,0.761719,0.830078,0.830078,0.810547,0.854492,0.78125,0.698242,0.712891,0.756836,0.78125,0.78125,0.883789,0.766602,0.732422,0.74707,0.595703,0.776367,0.59082,0.595703,0.776367,0.693359,0.722656,0.634766,0.703125,0.703125,0.776367,0.820312,0.664062,0.668945,0.649414,1,1,0.717773,0.649414,1,1,0.625,1,1,1,1,1,1,1,0.703125,1,0.65918,0.786133,0.712891,0.732422,0.664062,0.761719,0.664062,0.678711,0.644531,0.59082,0.727539,0.776367,0.678711,0.737305,0.629883,0.795898,0.795898,0.581055,0.668945,1,0.820312,0.78125,0.795898,0.639648,0.722656,0.678711,0.698242,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    calibrated_amp = amp/mmgCoefficients[ch];
  }
  else if ((runNum>=6389 && runNum<=6391) || (runNum>=6385 && runNum<=6388)) { //HG Al
    double mmgCoefficients[240] = {1,0.825195,1,0.874023,0.898438,0.839844,1,0.942383,1,0.854492,1,0.90332,0.854492,0.834961,0.898438,0.883789,1,0.844727,0.878906,0.883789,0.800781,0.893555,0.844727,0.830078,0.854492,0.786133,0.859375,0.859375,0.732422,0.78125,0.976562,0.976562,0.976562,0.922852,0.922852,0.893555,0.922852,0.922852,0.90332,0.976562,0.976562,0.97168,0.97168,0.917969,0.942383,1,0.9375,0.90332,0.932617,0.932617,0.908203,0.869141,0.9375,0.952148,0.893555,0.957031,0.913086,0.957031,0.961914,0.932617,0.913086,0.908203,0.9375,0.976562,0.942383,0.927734,0.917969,0.898438,0.869141,0.942383,0.966797,1,0.917969,0.927734,0.90332,0.957031,0.957031,0.917969,0.917969,0.908203,0.976562,0.97168,0.932617,0.947266,0.9375,0.966797,0.976562,0.913086,0.9375,0.932617,0.917969,1,1,0.942383,0.888672,0.888672,0.947266,0.952148,1,0.947266,0.917969,0.976562,0.869141,0.908203,0.942383,0.893555,0.869141,0.908203,0.878906,0.913086,0.893555,0.966797,0.9375,0.830078,0.864258,0.932617,0.869141,0.957031,0.917969,0.869141,1,0.893555,1,0.947266,1,0.90332,1,0.9375,1,0.952148,1,0.922852,0.913086,0.913086,0.888672,0.893555,1,0.922852,0.893555,0.869141,0.883789,0.869141,0.913086,0.917969,0.869141,0.849609,0.90332,0.854492,0.834961,0.898438,0.878906,0.874023,0.869141,0.888672,0.922852,0.878906,0.864258,0.883789,0.874023,0.849609,0.898438,0.913086,0.839844,0.908203,0.854492,0.864258,0.874023,0.864258,0.888672,0.898438,0.844727,0.869141,0.825195,0.859375,0.820312,0.869141,0.854492,0.883789,0.849609,0.917969,1,1,0.776367,0.874023,1,1,0.961914,1,1,1,1,1,1,1,0.898438,1,0.869141,0.869141,0.830078,0.869141,0.859375,0.800781,0.800781,0.917969,0.883789,0.839844,0.771484,0.883789,0.922852,0.957031,0.883789,0.854492,0.883789,0.932617,0.878906,1,0.898438,0.81543,0.854492,0.800781,0.820312,0.74707,0.883789,0.810547,1,0.805664,1,0.849609,0.820312,0.830078,0.78125,0.81543,0.825195,0.834961,0.81543,1,0.825195,0.869141,0.756836,0.786133};
    calibrated_amp = amp/mmgCoefficients[ch];
  }
  if (ch<0 || ch>=240) return amp;
  return calibrated_amp;
}

//============ END ADC Gain Calibrations ============

//--UNDO SRS strip-to-mm conversion from JANA2
/*int CalcAPVChannel(int peakPosition) {
  double pitch = 0.4;
  return (peakPosition - (-0.5*(102.4 - pitch))) / pitch;
}
*/
void trdclass_ps26::Loop() {
  
  if (fChain == 0) return;

  //==================================================================================================
  //            B o o k    H i s t o g r a m s
  //==================================================================================================
  
  gErrorIgnoreLevel = kBreak; // Suppress warning messages from empty chi^2 fit data
  TList *HistList = new TList();
  #ifdef WRITE_CSV
    std::ofstream csvFile("EventByEvent.csv");
  #endif
  //-----------------  (canvas 0) Event Display ----------
  #ifdef SHOW_EVT_DISPLAY
    //-----------------  canvas 2 FPGA Display ----------
    char c2Title[256];
    sprintf(c2Title,"FPGA_Event_Display_Run=%d",RunNum);
    TCanvas *c2 = new TCanvas("FPGA",c2Title,1000,100,1500,1000);
    c2->Divide(5,3); c2->cd(1);
  #endif
  
  //iRunNum-based timing window change
  double histTime = -1;
  if (RunNum<=timeSwitchRun1 || RunNum>=timeSwitchRun2) { histTime = firstTimeWin; } else { histTime = secondTimeWin; }
  // Track fit in time
  TF1 fx1("fx1","pol1",40,150);
  TF1 fx2("fx2","pol1",40,150);
  f125_fit = new TH2F("f125_fit","Triple GEM-TRD Track Fit; Time Response (*8ns) ; X Channel",(int)histTime,0.5,histTime+0.5,240,-0.5,239.5);
  mmg1_f125_fit = new TH2F("mmg1_f125_fit","MMG1-TRD Track Fit; Time Response (*8ns) ; X Channel",(int)histTime,0.5,histTime+0.5,240,-0.5,239.5);
  urw_f125_fit = new TH2F("urw_f125_fit","uRWELL-TRD Track Fit; Time Response (*8ns) ; X Channel",(int)histTime,0.5,histTime+0.5,120,-0.5,119.5);
  //-- TRD - GEMTRKR alignment --------
  double xx1=-37., yy1=-55.,  xx2=53., yy2=44.;
  double aa=(yy2-yy1)/(xx2-xx1);
  double bb=yy1-aa*xx1;
  TF1 ftrk("ftrk","[0]*x+[1]",-55.,55.);
  ftrk.SetParameter(0,aa);
  ftrk.SetParameter(1,bb);
  TF1 ftrkr("ftrkr","(x-[1])/[0]",0.,255.);
  ftrkr.SetParameter(0,aa);
  ftrkr.SetParameter(1,bb);
  //---- Define Z Positions [mm] -----
  float z1 = 0.;
  float z2 = 1031.08;
  float zmmg1 = 278.12;
  float zgem = 554.5;
  float zurw = 830.93;
  //---- Define Y Positions [mm] -----
  float y1 = 0.;
  float y2 = -3.90;
  float ymmg1 = -2.84;
  float ygem = -4.68;
  float yurw = -3.81;
  //---- Define X Positions [mm] -----
  /*
  float x1 = 0.;
  float x2 = 121.90;
  float xmmg1 = 33.22;
  float xgem = 65.57;
  float xurw = 97.99;
  */
  //---- Set ADC Thresholds ----
  float TGEM_THRESH=155.; //200
  float MMG1_THRESH=130.; //110
  float URW_THRESH=140.; //150
  float TRKR_THRESH=600.; //700.
  float FE55_THRESH=500.;
  
  hcount = new TH1D("hcount","Count",3,0,3);  HistList->Add(hcount);
  hcount->SetStats(0);   hcount->SetFillColor(38);   hcount->SetMinimum(1.);
  #if ROOT_VERSION_CODE > ROOT_VERSION(6,0,0)
    hcount->SetCanExtend(TH1::kXaxis);
  #else
    hcount->SetBit(TH1::kCanRebin);
  #endif
  
  htgem_nhits = new TH1F("hgem_nhits","Triple GEM-TRD X Pulses (fADC)",90,0,90);  HistList->Add(htgem_nhits);
  hmmg1_nhits = new TH1F("hmmg1_nhits","MMG1-TRD X Pulses (fADC)",90,0,90);  HistList->Add(hmmg1_nhits);
  hurw_nxhits = new TH1F("hurw_nxhits","uRWELL-TRD X Pulses (fADC)",75,0,75);  HistList->Add(hurw_nxhits);
  hurw_nyhits = new TH1F("hurw_nyhits","uRWELL-TRD Y Pulses (fADC)",75,0,75);  HistList->Add(hurw_nyhits);
  hgt1_nhits = new TH1F("hgt1_nhits","GEM-TRKR1 Pulses (SRS)",12,0,12);  HistList->Add(hgt1_nhits);
  hgt2_nhits = new TH1F("hgt2_nhits","GEM-TRKR2 Pulses (SRS)",12,0,12);  HistList->Add(hgt2_nhits);
  
  htgem_tmp_nhits = new TH1F("hgem_tmp_nhits","Triple GEM-TRD X Pulses (fADC) Pre-Cuts",90,0,90);  HistList->Add(htgem_tmp_nhits);
  hmmg1_tmp_nhits = new TH1F("hmmg1_tmp_nhits","MMG1-TRD X Pulses (fADC) Pre-Cuts",90,0,90);  HistList->Add(hmmg1_tmp_nhits);
  hurw_tmp_nxhits = new TH1F("hurw_tmp_nxhits","uRWELL-TRD X Pulses (fADC) Pre-Cuts",75,0,75);  HistList->Add(hurw_tmp_nxhits);
  
  cout<<"**************************RunNum="<<RunNum<<endl;
  int nx0=100;
  int mfac=30; //25; //60
  int ufac=15;
  int ny0=256; //240;
  int uny0=128;
  double Ymin=0.;    double Ymax=ny0*0.4;     double uYmax=uny0*0.8;
  double Xmin=0.;    double Xmax=30.; //double Xmax=26.;
  mhevt  = new TH2F("mhevt","MMG1-TRD Event Display; z pos [mm]; y pos [mm]",nx0+mfac,Xmin,Xmax,ny0,Ymin,Ymax); mhevt->SetStats(0); mhevt->SetMaximum(10.); mhevt->SetMinimum(0.);
  mhevtc = new TH2F("mhevtc","Clustering; FADC bins; MMG1 strips",nx0+mfac,-0.5,(nx0+mfac)-0.5,ny0,-0.5,ny0-0.5);  mhevtc->SetStats(0);   mhevtc->SetMinimum(0.07); mhevtc->SetMaximum(40.);
  mhevtf = new TH2F("mhevtf","MMG1: Clusters for FPGA; z pos [mm]; y pos [mm]",nx0+mfac,Xmin,Xmax,ny0,Ymin,Ymax);  mhevtf->SetStats(0); mhevtf->SetMaximum(10.);
  
  uhevt  = new TH2F("uhevt","URW-TRD Event Display; z pos [mm]; y pos [mm]",nx0+ufac,Xmin,Xmax,uny0,Ymin,uYmax); uhevt->SetStats(0); uhevt->SetMaximum(10.); uhevt->SetMinimum(0.);
  uhevtc = new TH2F("uhevtc","Clustering; FADC bins; URW strips",nx0+ufac,-0.5,(nx0+ufac)-0.5,uny0,-0.5,uny0-0.5);  uhevtc->SetStats(0);   uhevtc->SetMinimum(0.07); uhevtc->SetMaximum(40.);
  uhevtf = new TH2F("uhevtf","URW: Clusters for FPGA; z pos [mm]; y pos [mm]",nx0+ufac,Xmin,Xmax,uny0,Ymin,uYmax);  uhevtf->SetStats(0); uhevtf->SetMaximum(10.);
  
  hevt = new TH2F("hevt","GEM-TRD Event display; z pos [time *8ns]; y pos [chan #]",nx0,Xmin,Xmax,ny0,Ymin,Ymax); hevt->SetStats(0); hevt->SetMaximum(10.); hevt->SetMinimum(0.);
  hevtc = new TH2F("hevtc","Clustering; FADC bins; GEM strips",nx0,-0.5,nx0-0.5,ny0,-0.5,ny0-0.5); hevtc->SetStats(0);   hevtc->SetMinimum(0.07); hevtc->SetMaximum(40.);
  hevtf = new TH2F("hevtf","GEM: Clusters for FPGA; z pos [mm]; y pos [mm]",nx0,Xmin,Xmax,ny0,Ymin,Ymax);  hevtf->SetStats(0); hevtf->SetMaximum(10.);
  #if (USE_PULSE>0)
    hevtk = new TH2F("hevtk","Event display; z pos [mm]; y pos [mm]",nx0,Xmin,Xmax,ny0,Ymin,Ymax); /*hevtk->SetStats(0);*/ hevtk->SetMaximum(10.);
    hevtck = new TH2F("hevtck","Clustering; FADC bins; GEM strips",nx0,-0.5,nx0-0.5,ny0,-0.5,ny0-0.5);
  #endif
  
  //////
  f125_el = new TH1F("f125_el","Triple GEM-TRD f125 Peak Amp; ADC Amplitude ; Counts ",450,0.,4500.);          HistList->Add(f125_el);
  f125_el_max = new TH1F("f125_el_max","GEM-TRD f125 Max Amp; ADC Amplitude ; Counts ",450,0.,4500.);          HistList->Add(f125_el_max);
  f125_el_max_late = new TH1F("f125_el_max_late","GEM-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",450,0.,4500.);           HistList->Add(f125_el_max_late);
  mmg1_f125_el = new TH1F("mmg1_f125_el","MMG1-TRD f125 Peak Amp; ADC Amplitude ; Counts ",450,0.,4500.);       HistList->Add(mmg1_f125_el);
  urw_f125_el_x = new TH1F("urw_f125_el_x","uRWELL-TRD X f125 Peak Amp; ADC Amplitude ; Counts ",410,0.,4100.);       HistList->Add(urw_f125_el_x);
  urw_f125_el_y = new TH1F("urw_f125_el_y","uRWELL-TRD Y f125 Peak Amp; ADC Amplitude ; Counts ",410,0.,4100.);       HistList->Add(urw_f125_el_y);
  mmg1_f125_el_max = new TH1F("mmg1_f125_el_max","MMG1-TRD f125 Max Amp; ADC Amplitude ; Counts ",450,0.,4500.);           HistList->Add(mmg1_f125_el_max);
  mmg1_f125_el_max_late = new TH1F("mmg1_f125_el_max_late","MMG1-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",450,0.,4500.);           HistList->Add(mmg1_f125_el_max_late);
  urw_f125_el_xmax = new TH1F("urw_f125_el_xmax","uRWELL-TRD X f125 Max Amp; ADC Amplitude ; Counts ",410,0.,4100.);           HistList->Add(urw_f125_el_xmax);
  urw_f125_el_xmax_late = new TH1F("urw_f125_el_xmax_late","uRWELL-TRD X f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",410,0.,4100.);           HistList->Add(urw_f125_el_xmax_late);
  urw_f125_el_ymax = new TH1F("urw_f125_el_ymax","uRWELL-TRD Y f125 Max Amp; ADC Amplitude ; Counts ",410,0.,4100.);           HistList->Add(urw_f125_el_ymax);

  // --- SRS ---
  hgemtrkr_1_max_xch = new TH1F("hgemtrkr_1_max_xch"," GEM-TRKR1 Max X Position ; X Chan [mm]",256,-0.2,102.2);  HistList->Add(hgemtrkr_1_max_xch);
  hgemtrkr_1_max_xamp = new TH1F("hgemtrkr_1_max_xamp"," GEM-TRKR1 Max X Amp ; ADC Amp ",410,0.,4100.);  HistList->Add(hgemtrkr_1_max_xamp);
  hgemtrkr_2_max_xch = new TH1F("hgemtrkr_2_max_xch"," GEM-TRKR2 Max X Position ; X Chan [mm]",256,-0.2,102.2);  HistList->Add(hgemtrkr_2_max_xch);
  hgemtrkr_2_max_xamp = new TH1F("hgemtrkr_2_max_xamp"," GEM-TRKR2 Max X Amp ; ADC Amp ",410,0.,4100.);  HistList->Add(hgemtrkr_2_max_xamp);
  hgemtrkr_max_xcorr = new TH2F("hgemtrkr_max_xcorr","Max X Correlation for GEM-TRKRs; GEM-TRKR1 X [mm]; GEM-TRKR2 X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(hgemtrkr_max_xcorr);
  hgemtrkr_max_ycorr = new TH2F("hgemtrkr_max_ycorr","Max Y Correlation for GEM-TRKRs; GEM-TRKR1 Y [mm]; GEM-TRKR2 Y [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(hgemtrkr_max_ycorr);
  
  hgemtrkr_1D_xcorr = new TH1F("hgemtrkr_1D_xcorr"," Corrected XCorr for GEM-TRKRs ; GEM-TRKR2 X - GEM-TRKR1 X [mm]",100,-20.,20.);  HistList->Add(hgemtrkr_1D_xcorr);
hgemtrkr_1D_ycorr = new TH1F("hgemtrkr_1D_ycorr"," Corrected YCorr for GEM-TRKRs ; GEM-TRKR2 Y - GEM-TRKR1 Y [mm]",100,-20.,20.);  HistList->Add(hgemtrkr_1D_ycorr);
  
  //--GEMTracker 1
  hgemtrkr_1_peak_xy = new TH2F("hgemtrkr_1_peak_xy","GEM-TRKR1 Peak X-Y Correlation; Peak X [mm]; Peak Y [mm] ",256,-0.2,102.2,256,-0.2,102.2);    hgemtrkr_1_peak_xy->SetStats(0); HistList->Add(hgemtrkr_1_peak_xy);
  hgemtrkr_1_max_xy = new TH2F("hgemtrkr_1_max_xy","GEM-TRKR1 X-Y Correlation for Max Hits; Peak X [mm]; Peak Y [mm] ",256,-0.2,102.2,256,-0.2,102.2);    hgemtrkr_1_max_xy->SetStats(0); HistList->Add(hgemtrkr_1_max_xy);
  hgemtrkr_1_peak_x = new TH1F("hgemtrkr_1_peak_x"," GEM-TRKR1 Peak X Position; X [mm] ",256,-0.2,102.2);  HistList->Add(hgemtrkr_1_peak_x);
  hgemtrkr_1_peak_y = new TH1F("hgemtrkr_1_peak_y"," GEM-TRKR1 Peak Y Position; Y [mm] ",256,-0.2,102.2);  HistList->Add(hgemtrkr_1_peak_y);
  hgemtrkr_1_peak_x_height = new TH1F("hgemtrkr_1_peak_x_height"," GEM-TRKR1 Peak Amplitudes in X; ADC Value ",410,0.,4100.);  HistList->Add(hgemtrkr_1_peak_x_height);
  hgemtrkr_1_peak_y_height = new TH1F("hgemtrkr_1_peak_y_height"," GEM-TRKR1 Peak Amplitudes in Y ; ADC Value ",410,0.,4100.);  HistList->Add(hgemtrkr_1_peak_y_height);
  //--GEMTracker 2
  hgemtrkr_2_peak_xy = new TH2F("hgemtrkr_2_peak_xy","GEM-TRKR2 Peak X-Y Correlation; Peak X [mm]; Peak Y [mm] ",256,-0.2,102.2,256,-0.2,102.2); hgemtrkr_2_peak_xy->SetStats(0); HistList->Add(hgemtrkr_2_peak_xy);
  hgemtrkr_2_max_xy = new TH2F("hgemtrkr_2_max_xy","GEM-TRKR2 X-Y Correlation for Max Hits; Peak X [mm]; Peak Y [mm] ",256,-0.2,102.2,256,-0.2,102.2); hgemtrkr_2_max_xy->SetStats(0);   HistList->Add(hgemtrkr_2_max_xy);
  hgemtrkr_2_peak_x = new TH1F("hgemtrkr_2_peak_x"," GEM-TRKR2 Peak X Position; X [mm] ",256,-0.2,102.2);  HistList->Add(hgemtrkr_2_peak_x);
  hgemtrkr_2_peak_y = new TH1F("hgemtrkr_2_peak_y"," GEM-TRKR2 Peak Y Position ; Y [mm] ",256,-0.2,102.2);  HistList->Add(hgemtrkr_2_peak_y);
  hgemtrkr_2_peak_x_height = new TH1F("hgemtrkr_2_peak_x_height"," GEM-TRKR2 Peak Amplitudes in X; ADC Value ",410,0.,4100.);  HistList->Add(hgemtrkr_2_peak_x_height);
  hgemtrkr_2_peak_y_height = new TH1F("hgemtrkr_2_peak_y_height"," GEM-TRKR2 Peak Amplitudes in Y ; ADC Value ",410,0.,4100.);  HistList->Add(hgemtrkr_2_peak_y_height);
  
  mmg1_peak_y = new TH1F("mmg1_peak_y"," MMG1-TRD Peak Y Position (SRS) ; Y [mm] ",256,-0.2,102.2);  HistList->Add(mmg1_peak_y);
  hmmg1_peak_y_height = new TH1F("hmmg1_peak_y_height"," MMG1-TRD Peak Amplitudes in Y ; ADC Value ",410,0.,4100.);  HistList->Add(hmmg1_peak_y_height);
  tgem_peak_y = new TH1F("tgem_peak_y","Triple GEM-TRD Peak Y Position (SRS) ; Y [mm] ",256,-0.2,102.2);  HistList->Add(tgem_peak_y);
  htgem_peak_y_height = new TH1F("htgem_peak_y_height","Triple GEM-TRD Peak Amplitudes in Y ; ADC Value ",410,0.,4100.);  HistList->Add(htgem_peak_y_height);
  hgemtrkr_1_tgem = new TH2F("hgemtrkr_1_tgem","GEM-TRKR1 & Triple GEM-TRD Y Correlation; Triple GEM-TRD Y [mm]; GEM-TRKR1 Y [mm] ",256,-0.2,102.2,256,-0.2,102.2); hgemtrkr_1_tgem->SetStats(0); HistList->Add(hgemtrkr_1_tgem);
  hgemtrkr_1_mmg1 = new TH2F("hgemtrkr_1_mmg1","GEM-TRKR1 & MMG1 Y Correlation; MMG1-TRD Y [mm]; GEM-TRKR1 Y [mm] ",256,-0.2,102.2,256,-0.2,102.2); hgemtrkr_1_mmg1->SetStats(0); HistList->Add(hgemtrkr_1_mmg1);
  
  //--External Tracking
  f125_el_tracker_hits = new TH1F("f125_el_tracker_hits","GEM-TRD Track Extr. Hits; X Chan [mm]",128,-0.2,102.2);   HistList->Add(f125_el_tracker_hits);  f125_el_tracker_hits->Sumw2();
  f125_el_tracker_eff = new TH1F("f125_el_tracker_eff","GEM-TRD Track Eff.; X Chan [mm]",128,-0.2,102.2);   HistList->Add(f125_el_tracker_eff); f125_el_tracker_eff->Sumw2();
  tgem_residuals = new TH1F("tgem_residuals","Triple GEM-TRD Residual Hits; X Chan [mm] (actual - expected)",125,-25.,25.);   HistList->Add(tgem_residuals);
  tgem_residualscorr = new TH1F("tgem_residualscorr","Triple GEM-TRD Residual Hits WITH CORR.; X Chan [mm] (actual - expected)",125,-25.,25.);   HistList->Add(tgem_residualscorr);
  tgem_residual_ch = new TH2F("tgem_residual_ch","Triple GEM-TRD Residual Hits vs Chan; X Chan [mm] (actual); X Chan [mm] (actual - expected)",175, 11.5, 81.5,35,-7.5,6.5); tgem_residual_ch->SetStats(0); HistList->Add(tgem_residual_ch);
  tgem_residual_chcorr = new TH2F("tgem_residual_chcorr","Triple GEM-TRD Residual Hits vs Chan WITH CORR.; X Chan [mm] (actual); X Chan [mm] (actual - expected)",175, 11.5, 81.5,35,-7.5,6.5); tgem_residual_chcorr->SetStats(0); HistList->Add(tgem_residual_chcorr);
  mmg1_f125_el_tracker_hits = new TH1F("mmg1_f125_el_tracker_hits","MMG1-TRD Track Extr. Hits; X Chan [mm]",128,-0.2,102.2);   HistList->Add(mmg1_f125_el_tracker_hits);  mmg1_f125_el_tracker_hits->Sumw2();
  mmg1_f125_el_tracker_eff = new TH1F("mmg1_f125_el_tracker_eff","MMG1-TRD Track Eff.; X Chan [mm]",128,-0.2,102.2);   HistList->Add(mmg1_f125_el_tracker_eff);  mmg1_f125_el_tracker_eff->Sumw2();
  mmg1_residuals = new TH1F("mmg1_residuals","MMG1-TRD Residual Hits; X Chan [mm] (actual - expected)",125,-25.,25.);   HistList->Add(mmg1_residuals);
  mmg1_residualscorr = new TH1F("mmg1_residualscorr","MMG1-TRD Residual Hits WITH CORR.; X Chan [mm] (actual - expected)",125,-25.,25.);   HistList->Add(mmg1_residualscorr);
  mmg1_residual_ch = new TH2F("mmg1_residual_ch","MMG1-TRD Residual Hits vs Chan; X Chan [mm] (actual); X Chan [mm] (actual - expected)",125, 10.5, 60.5,35,-7.5,6.5); mmg1_residual_ch->SetStats(0); HistList->Add(mmg1_residual_ch);
  mmg1_residual_chcorr = new TH2F("mmg1_residual_chcorr","MMG1-TRD Residual Hits vs Chan WITH CORR.; X Chan [mm] (actual); X Chan [mm] (actual - expected)",125, 10.5, 60.5,35,-7.5,6.5); mmg1_residual_chcorr->SetStats(0); HistList->Add(mmg1_residual_chcorr);
  urw_f125_x_tracker_hits = new TH1F("urw_f125_x_tracker_hits","uRWELL-TRD X Track Extr. Hits; X Chan [mm]",64,-0.4,102.);   HistList->Add(urw_f125_x_tracker_hits); urw_f125_x_tracker_hits->Sumw2();
  urw_f125_x_tracker_eff = new TH1F("urw_f125_x_tracker_eff","uRWELL-TRD X Track Eff.; X Chan [mm]",64,-0.4,102.);   HistList->Add(urw_f125_x_tracker_eff); urw_f125_x_tracker_eff->Sumw2();
  urw_x_residuals = new TH1F("urw_x_residuals","uRWELL-TRD X Residual Hits; X Chan [mm] (actual - expected)",65,-26.,26.);   HistList->Add(urw_x_residuals);
  urw_x_residualscorr = new TH1F("urw_x_residualscorr","uRWELL-TRD X Residual Hits WITH CORR.; X Chan [mm] (actual - expected)",65,-26.,26.);   HistList->Add(urw_x_residualscorr);
  urw_x_residual_ch = new TH2F("urw_x_residual_ch","uRWELL-TRD Residual Hits vs Chan; X Chan [mm] (actual); X Chan [mm] (actual - expected)",95, 12.5, 88.5,15,-4.5,7.5); urw_x_residual_ch->SetStats(0); HistList->Add(urw_x_residual_ch);
  urw_x_residual_chcorr = new TH2F("urw_x_residual_chcorr","uRWELL-TRD Residual Hits vs Chan WITH CORR.; X Chan [mm] (actual); X Chan [mm] (actual - expected)",95, 12.5, 88.5,15,-6.5,5.5); urw_x_residual_chcorr->SetStats(0); HistList->Add(urw_x_residual_chcorr);
  
  //Efficiency Plots
  p_tgem_eff = new TEfficiency("p_tgem_eff","Triple GEM-TRD Tracking Efficiency",128,-0.2,102.2); p_tgem_eff = 0;  HistList->Add(p_tgem_eff);
  p_mmg1_eff = new TEfficiency("p_mmg1_eff","MMG1-TRD Tracking Efficiency",128,-0.2,102.2); p_mmg1_eff = 0;  HistList->Add(p_mmg1_eff);
  p_urw_eff = new TEfficiency("p_urw_eff","uRWELL-TRD Tracking Efficiency",64,-0.4,102.); p_urw_eff = 0;  HistList->Add(p_urw_eff);
  
  //-- Amplitude Histos --
  f125_el_amp2d = new TH2F("f125_el_amp2d","Triple GEM-TRD ADC Amp in Time (After Selections); Time Response (*8ns) ; X Channel ",(int)histTime,0.5,histTime+.5,240,-0.5,239.5); f125_el_amp2d->SetStats(0); HistList->Add(f125_el_amp2d);
  f125_el_amp2d_max = new TH2F("f125_el_amp2d_max","Triple GEM-TRD Max ADC Amp in Time; Time Response (*8ns) ; X [mm] ",(int)histTime,0.5,histTime+.5,256,-0.2,102.2); f125_el_amp2d_max->SetStats(0); HistList->Add(f125_el_amp2d_max);
  f125_xVSamp = new TH2F("f125_xVSamp","Triple GEM-TRD X Channel vs ADC Amp; X Channel; ADC Value",240,-0.5,239.5,410,0.,4100.); f125_xVSamp->SetStats(0); HistList->Add(f125_xVSamp);
  f125_xVSamp_max = new TH2F("f125_xVSamp_max","Triple GEM-TRD X Channel vs Max ADC Amp; X Channel; ADC Value",240,-0.5,239.5,410,0.,4100.); f125_xVSamp_max->SetStats(0); HistList->Add(f125_xVSamp_max);
  f125_timeVSamp = new TH2F("f125_timeVSamp","Triple GEM-TRD X ADC Amp in Time; Time Response (*8ns); ADC Value",(int)histTime,0.5,histTime+.5,410,0.,4100.); f125_timeVSamp->SetStats(0); HistList->Add(f125_timeVSamp);
  f125_timeVSamp_max = new TH2F("f125_timeVSamp_max","Triple GEM-TRD X Max ADC Amp in Time; Time Response (*8ns); ADC Value",(int)histTime,0.5,histTime+.5,410,0.,4100.); f125_timeVSamp_max->SetStats(0); HistList->Add(f125_timeVSamp_max);

  urw_f125_x_amp2d = new TH2F("urw_f125_x_amp2d","uRWELL-TRD X ADC Amp in Time (After Selections); Time Response (*8ns) ; X Channel ",(int)histTime,0.5,histTime+.5,120,-0.5,119.5); urw_f125_x_amp2d->SetStats(0); HistList->Add(urw_f125_x_amp2d);
  urw_f125_x_amp2d_max = new TH2F("urw_f125_x_amp2d_max","uRWELL-TRD X Max ADC Amp in Time; Time Response (*8ns) ; X [mm] ",(int)histTime,0.5,histTime+.5,128,-0.2,102.2); urw_f125_x_amp2d_max->SetStats(0); HistList->Add(urw_f125_x_amp2d_max);
  urw_f125_xVSamp = new TH2F("urw_f125_xVSamp","uRWELL-TRD X Channel vs ADC Amp; X Channel; ADC Value",120,-0.5,119.5,410,0.,4100.); urw_f125_xVSamp->SetStats(0); HistList->Add(urw_f125_xVSamp);
  urw_f125_xVSamp_max = new TH2F("urw_f125_xVSamp_max","uRWELL-TRD X Channel vs Max ADC Amp; X Channel; ADC Value",120,-0.5,119.5,410,0.,4100.); urw_f125_xVSamp_max->SetStats(0); HistList->Add(urw_f125_xVSamp_max);
  urw_f125_x_timeVSamp = new TH2F("urw_f125_x_timeVSamp","uRWELL-TRD X ADC Amp in Time (After Selections); Time Response (*8ns); ADC Value",(int)histTime,0.5,histTime+.5,410,0.,4100.); urw_f125_x_timeVSamp->SetStats(0); HistList->Add(urw_f125_x_timeVSamp);
  urw_f125_x_timeVSamp_max = new TH2F("urw_f125_x_timeVSamp_max","uRWELL-TRD X Max ADC Amp in Time; Time Response (*8ns); ADC Value",(int)histTime,0.5,histTime+.5,410,0.,4100.); urw_f125_x_timeVSamp_max->SetStats(0); HistList->Add(urw_f125_x_timeVSamp_max);

  mmg1_f125_el_amp2d = new TH2F("mmg1_f125_el_amp2d","MMG1-TRD ADC Amp in Time (After Selections); Time Response (*8ns) ; X Channel ",(int)histTime,0.5,histTime+.5,240,-0.5,239.5); mmg1_f125_el_amp2d->SetStats(0); HistList->Add(mmg1_f125_el_amp2d);
  mmg1_f125_el_amp2d_max = new TH2F("mmg1_f125_el_amp2d_max","MMG1-TRD Max ADC Amp in Time; Time Response (*8ns) ; X [mm] ",(int)histTime,0.5,histTime+.5,256,-0.2,102.2); mmg1_f125_el_amp2d_max->SetStats(0); HistList->Add(mmg1_f125_el_amp2d_max);
  mmg1_f125_xVSamp = new TH2F("mmg1_f125_xVSamp","MMG1-TRD X Channel vs ADC Amp; X Channel; ADC Value",240,-0.5,239.5,410,0.,4100.); mmg1_f125_xVSamp->SetStats(0); HistList->Add(mmg1_f125_xVSamp);
  mmg1_f125_xVSamp_max = new TH2F("mmg1_f125_xVSamp_max","MMG1-TRD X Channel vs Max ADC Amp; X Channel; ADC Value",240,-0.5,239.5,410,0.,4100.); mmg1_f125_xVSamp_max->SetStats(0); HistList->Add(mmg1_f125_xVSamp_max);
  mmg1_f125_timeVSamp = new TH2F("mmg1_f125_timeVSamp","MMG1-TRD X ADC Amp in Time; Time Response (*8ns); ADC Value",(int)histTime,0.5,histTime+.5,410,0.,4100.); mmg1_f125_timeVSamp->SetStats(0); HistList->Add(mmg1_f125_timeVSamp);
  mmg1_f125_timeVSamp_max = new TH2F("mmg1_f125_timeVSamp_max","MMG1-TRD X Max ADC Amp in Time; Time Response (*8ns); ADC Value",(int)histTime,0.5,histTime+.5,410,0.,4100.); mmg1_f125_timeVSamp_max->SetStats(0); HistList->Add(mmg1_f125_timeVSamp_max);

  urw_f125_y_amp2d = new TH2F("urw_f125_y_amp2d","uRWELL-TRD Y ADC Amp in Time; Time Response (*8ns) ; Y Channel ",(int)histTime,0.5,histTime+.5,120,-0.5,119.5); urw_f125_y_amp2d->SetStats(0); HistList->Add(urw_f125_y_amp2d);
  urw_f125_y_amp2d_max = new TH2F("urw_f125_y_amp2d_max","uRWELL-TRD Y Max ADC Amp in Time; Time Response (*8ns) ; Y [mm] ",(int)histTime,0.5,histTime+.5,128,-0.2,102.2); urw_f125_y_amp2d_max->SetStats(0); HistList->Add(urw_f125_y_amp2d_max);
  urw_f125_yVSamp = new TH2F("urw_f125_yVSamp","uRWELL-TRD Y Channel vs ADC Amp; Y Channel; ADC Value",120,-0.5,119.5,410,0.,4100.); urw_f125_yVSamp->SetStats(0); HistList->Add(urw_f125_yVSamp);
   urw_f125_yVSamp_max = new TH2F("urw_f125_yVSamp_max","uRWELL-TRD Y Channel vs Max ADC Amp; Y Channel; ADC Value",120,-0.5,119.5,410,0.,4100.); urw_f125_yVSamp_max->SetStats(0); HistList->Add(urw_f125_yVSamp_max);
  urw_f125_y_timeVSamp = new TH2F("urw_f125_y_timeVSamp","uRWELL-TRD Y ADC Amp in Time (After Selections); Time Response (*8ns); ADC Value",(int)histTime,0.5,histTime+.5,410,0.,4100.); urw_f125_y_timeVSamp->SetStats(0); HistList->Add(urw_f125_y_timeVSamp);
  urw_f125_y_timeVSamp_max = new TH2F("urw_f125_y_timeVSamp_max","uRWELL-TRD Y Max ADC Amp in Time; Time Response (*8ns); ADC Value",(int)histTime,0.5,histTime+.5,410,0.,4100.); urw_f125_y_timeVSamp_max->SetStats(0); HistList->Add(urw_f125_y_timeVSamp_max);
  
  f125_el_amp2ds = new TH2F("f125_el_amp2ds","Triple GEM-TRD ADC Amp in Time (Before Selections); Time Response (*8ns) ; X Channel ",(int)histTime,0.5,histTime+.5,240,-0.5,239.5); f125_el_amp2ds->SetStats(0); HistList->Add(f125_el_amp2ds);
  urw_f125_x_amp2ds = new TH2F("urw_f125_x_amp2ds","uRWELL-TRD X ADC Amp in Time (Before Selections); Time Response (*8ns) ; X Channel ",(int)histTime,0.5,histTime+.5,120,-0.5,119.5); urw_f125_x_amp2ds->SetStats(0); HistList->Add(urw_f125_x_amp2ds);
  mmg1_f125_el_amp2ds = new TH2F("mmg1_f125_el_amp2ds","MMG1-TRD ADC Amp in Time (Before Selections); Time Response (*8ns) ; X Channel ",(int)histTime,0.5,histTime+.5,240,-0.5,239.5); mmg1_f125_el_amp2ds->SetStats(0); HistList->Add(mmg1_f125_el_amp2ds);
  urw_f125_y_amp2ds = new TH2F("urw_f125_y_amp2ds","uRWELL-TRD Y ADC Amp in Time (Before Selections); Time Response (*8ns) ; Y Channel ",(int)histTime,0.5,histTime+.5,120,-0.5,119.5); urw_f125_y_amp2ds->SetStats(0); HistList->Add(urw_f125_y_amp2ds);

  //-- 2D Hit Displays --
  htgem_xy = new TH2F("htgem_xy","Triple GEM-TRD X-Y Hit Display; X (fADC) [mm]; Y (SRS) [mm] ",256,-0.2,102.2,256,-0.2,102.2); htgem_xy->SetStats(0);  HistList->Add(htgem_xy);
  htgem_max_xy = new TH2F("htgem_max_xy","Triple GEM-TRD X-Y Max Hit Display; X (fADC) [mm]; Y (SRS) [mm] ",256,-0.2,102.2,256,-0.2,102.2); htgem_max_xy->SetStats(0);  HistList->Add(htgem_max_xy);
  hmmg1_xy = new TH2F("hmmg1_xy","MMG1-TRD X-Y Hit Display; X (fADC) [mm]; Y (SRS) [mm] ",256,-0.2,102.2,256,-0.2,102.2); hmmg1_xy->SetStats(0);  HistList->Add(hmmg1_xy);
  hmmg1_max_xy = new TH2F("hmmg1_max_xy","MMG1-TRD X-Y Max Hit Display; X (fADC) [mm]; Y (SRS) [mm] ",256,-0.2,102.2,256,-0.2,102.2); hmmg1_max_xy->SetStats(0);  HistList->Add(hmmg1_max_xy);
  hurw_xy = new TH2F("hurw_xy","uRWELL-TRD X-Y Hit Display; X (fADC) [mm]; Y (fADC) [mm] ",128,-0.4,102.,128,-0.4,102.); hurw_xy->SetStats(0);  HistList->Add(hurw_xy);
  hurw_max_xy = new TH2F("hurw_max_xy","uRWELL-TRD X-Y Max Hit Display; X (fADC) [mm]; Y (fADC) [mm] ",128,-0.4,102.,128,-0.4,102.); hurw_max_xy->SetStats(0);  HistList->Add(hurw_max_xy);

  //-- X,Y Correlations --
  hmmg1_tgem_ydiff = new TH1F("hmmg1_tgem_ydiff","Y Position Diff. Between MMG1 & GEM; GEM Y(SRS) - MMG-1 Y(SRS) [mm]",100,-20.,20.);    HistList->Add(hmmg1_tgem_ydiff);
  hmmg1_urw_xdiff = new TH1F("hmmg1_urw_xdiff","X Position Diff. Between MMG1 & uRWell; MMG-1 X - uRWell X [mm]",100,-20.,20.);    HistList->Add(hmmg1_urw_xdiff);
  hmmg1_tgem_xdiff = new TH1F("hmmg1_tgem_xdiff","X Position Diff. Between MMG1 & GEM; MMG-1 X - GEM X [mm]",100,-20.,20.);    HistList->Add(hmmg1_tgem_xdiff);
  htgem_urw_xdiff = new TH1F("htgem_urw_xdiff","X Position Diff. Between GEM & uRWell; GEM X - uRWell X [mm]",100,-20.,20.);    HistList->Add(htgem_urw_xdiff);
  htgem_timeDiff = new TH2F("htgem_timeDiff","Time of Max Hit Diff. vs Channel for GEM; GEM X Channel# ; (Max Pulse Time - All Pulse Times) (*8ns)",240,-0.5,239.5,120,-121,119);    HistList->Add(htgem_timeDiff);
  hmmg1_timeDiff = new TH2F("hmmg1_timeDiff","Time of Max Hit Diff. vs Channel for MMG; MMG X Channel# ; (Max Pulse Time - All Pulse Times) (*8ns)",240,-0.5,239.5,120,-121,119);    HistList->Add(hmmg1_timeDiff);
  hurw_timeDiff = new TH2F("hurw_timeDiff","Time of Max Hit Diff. vs Channel for uRWELL; uRWELL X Channel# ; (Max Pulse Time - All Pulse Times) (*8ns)",120,-0.5,119.5,120,-121,119);    HistList->Add(hurw_timeDiff);
  htgem_2DPulseMultiplicity = new TH2F("htgem_2DPulseMultiplicity","Pulse Multiplicity per Channel for GEM; GEM X Channel# ; N fADC Pulses",240,-0.5,239.5,20,-0.5,19.5);    HistList->Add(htgem_2DPulseMultiplicity);
  hmmg1_2DPulseMultiplicity = new TH2F("hmmg1_2DPulseMultiplicity","Pulse Multiplicity per Channel for MMG; MMG X Channel# ; N fADC Pulses",240,-0.5,239.5,20,-0.5,19.5);    HistList->Add(hmmg1_2DPulseMultiplicity);
  hurw_2DPulseMultiplicity = new TH2F("hurw_2DPulseMultiplicity","Pulse Multiplicity per Channel for uRWELL; uRWELL X Channel# ; N fADC Pulses",120,-0.5,119.5,20,-0.5,19.5);    HistList->Add(hurw_2DPulseMultiplicity);
  htgem_2DPulseVsChan = new TH2F("htgem_2DPulseVsChan","Pulse Multiplicity per Channel Weighted for GEM; GEM X Channel# ; N fADC Pulses, Weighted",240,-0.5,239.5,150,-0.5,149.5);    HistList->Add(htgem_2DPulseVsChan);
  hmmg1_2DPulseVsChan = new TH2F("hmmg1_2DPulseVsChan","Pulse Multiplicity per Channel Weighted for MMG; MMG X Channel# ; N fADC Pulses, Weighted",240,-0.5,239.5,150,-0.5,149.5);    HistList->Add(hmmg1_2DPulseVsChan);
  hurw_2DPulseVsChan = new TH2F("hurw_2DPulseVsChan","Pulse Multiplicity per Channel Weighted for uRWELL; uRWELL X Channel# ; N fADC Pulses, Weighted",120,-0.5,119.5,150,-0.5,149.5);    HistList->Add(hurw_2DPulseVsChan);
  
  tgem_mmg1_xcorr = new TH2F("tgem_mmg1_xcorr","Triple GEM-TRD X Correlation With MMG1-TRD; Triple GEM-TRD X [mm]; MMG1-TRD X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(tgem_mmg1_xcorr);
  urw_tgem_xcorr = new TH2F("urw_tgem_xcorr","Triple GEM-TRD X Correlation With uRWELL-TRD; uRWELL-TRD X [mm]; Triple GEM-TRD X [mm]",128,-0.4,102.,256,-0.2,102.2); urw_tgem_xcorr->SetStats(0); HistList->Add(urw_tgem_xcorr);
  urw_mmg1_xcorr = new TH2F("urw_mmg1_xcorr","MMG1-TRD X Correlation With uRWELL-TRD; uRWELL-TRD X [mm];  MMG1-TRD X [mm]",128,-0.4,102.,256,-0.2,102.2); urw_mmg1_xcorr->SetStats(0); HistList->Add(urw_mmg1_xcorr);
  tgem_mmg1_ycorr = new TH2F("tgem_mmg1_ycorr","Triple GEM-TRD Y Correlation With MMG1-TRD; Triple GEM-TRD Y [mm]; MMG1-TRD Y [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(tgem_mmg1_ycorr);
  tgem_mmg1_max_xcorr = new TH2F("tgem_mmg1_max_xcorr","Triple GEM-TRD Max X Correlation With MMG1-TRD; Triple GEM-TRD X [mm]; MMG1-TRD X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(tgem_mmg1_max_xcorr);
  tgem_urw_max_xcorr = new TH2F("tgem_urw_max_xcorr","Triple GEM-TRD Max X Correlation With uRWELL-TRD; uRWELL-TRD X [mm]; Triple GEM-TRD X [mm]",128,-0.4,102.,256,-0.2,102.2); HistList->Add(tgem_urw_max_xcorr);
  urw_mmg1_max_xcorr = new TH2F("urw_mmg1_max_xcorr","uRWELL-TRD Max X Correlation With MMG1-TRD; uRWELL-TRD X [mm]; MMG1-TRD X [mm] ",128,-0.4,102.,256,-0.2,102.2); HistList->Add(urw_mmg1_max_xcorr);
  
  htgem_trdTrackCorr = new TH2F("htgem_trdTrackCorr","Triple GEM-TRD X Correlation With Matched Track; Triple GEM-TRD X [mm]; GEM-TRD Matched Ext. Track Extrap. X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(htgem_trdTrackCorr);
  hmmg1_trdTrackCorr = new TH2F("hmmg1_trdTrackCorr","MMG-TRD X Correlation With Matched Track; MMG-TRD X [mm]; MMG-TRD Matched Ext. Track Extrap. X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(hmmg1_trdTrackCorr);
  hurw_trdTrackCorr = new TH2F("hurw_trdTrackCorr","uRWELL-TRD X Correlation With Matched Track; uRWELL-TRD X [mm]; uRWELL-TRD Matched Ext. Track Extrap. X [mm] ",128,-0.4,102.,128,-0.4,102.); HistList->Add(hurw_trdTrackCorr);
  
  tgem_gt1_xcorr = new TH2F("tgem_gt1_xcorr","Triple GEM-TRD X Correlation With GEMTRKR-1; Triple GEM-TRD X [mm]; GEMTRKR-1 X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(tgem_gt1_xcorr);
  tgem_gt2_xcorr = new TH2F("tgem_gt2_xcorr","Triple GEM-TRD X Correlation With GEMTRKR-2; Triple GEM-TRD X [mm]; GEMTRKR-2 X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(tgem_gt2_xcorr);
  mmg1_gt1_xcorr = new TH2F("mmg1_gt1_xcorr","MMG1-TRD X Correlation With GEMTRKR-1; MMG1-TRD X [mm]; GEMTRKR-1 X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(mmg1_gt1_xcorr);
  mmg1_gt2_xcorr = new TH2F("mmg1_gt2_xcorr","MMG1-TRD X Correlation With GEMTRKR-2; MMG1-TRD X [mm]; GEMTRKR-2 X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(mmg1_gt2_xcorr);
  urw_gt1_xcorr = new TH2F("urw_gt1_xcorr","uRWell-TRD X Correlation With GEMTRKR-1; uRWell-TRD X [mm]; GEMTRKR-1 X [mm] ",128,-0.4,102.,256,-0.2,102.2); HistList->Add(urw_gt1_xcorr);
  urw_gt2_xcorr = new TH2F("urw_gt2_xcorr","uRWell-TRD X Correlation With GEMTRKR-2; uRWell-TRD X [mm]; GEMTRKR-2 X [mm] ",128,-0.4,102.,256,-0.2,102.2); HistList->Add(urw_gt2_xcorr);
  
  #if (USE_CLUST>0)
  hgemClusterDiff_el = new TH1F("hgemClusterDiff_el","GEM Cluster Distance from MMG1 Track; Distance [mm]",200,-20.,20.); HistList->Add(hgemClusterDiff_el);
  hmmg1ClusterDiff_el = new TH1F("hmmg1ClusterDiff_el","MMG1 Cluster Distance from GEM Track; Distance [mm]",200,-20.,20.); HistList->Add(hmmg1ClusterDiff_el);
  hurwClusterDiff_el = new TH1F("hurwClusterDiff_el","URW Cluster Distance from GEM Track; Distance [mm]",200,-20.,20.); HistList->Add(hurwClusterDiff_el);
  hClusterMaxdEdx_el = new TH1F("hClusterMaxdEdx_el","GEM Max Cluster Energy; Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hClusterMaxdEdx_el);
  hClusterTotaldEdx_el = new TH1F("hClusterTotaldEdx_el","GEM Total Cluster Energy; Total Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hClusterTotaldEdx_el);
  hmmg1ClusterMaxdEdx_el = new TH1F("hmmg1ClusterMaxdEdx_el","MMG1 Max Cluster Energy; Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hmmg1ClusterMaxdEdx_el);
  hmmg1ClusterTotaldEdx_el = new TH1F("hmmg1ClusterTotaldEdx_el","MMG1 Total Cluster Energy; Total Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hmmg1ClusterTotaldEdx_el);
  hurwClusterMaxdEdx_el = new TH1F("hurwClusterMaxdEdx_el","URW Max Cluster Energy; Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hurwClusterMaxdEdx_el);
  hurwClusterTotaldEdx_el = new TH1F("hurwClusterTotaldEdx_el","URW Total Cluster Energy; Total Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hurwClusterTotaldEdx_el);
  #endif
  
  hgemPulseDiff_el = new TH1F("hgemPulseDiff_el","GEM Pulse Distance from MMG1 Track; (Max MMG Pulse Pos) - (All GEM Pulse Positions) [mm]",150,-30.,30.); HistList->Add(hgemPulseDiff_el);
  hmmg1PulseDiff_el = new TH1F("hmmg1PulseDiff_el","MMG1 Pulse Distance from GEM Track; (All MMG Pulse Positions) - (Max GEM Pulse Pos) [mm]",150,-30.,30.); HistList->Add(hmmg1PulseDiff_el);
  hurwPulseDiff_el = new TH1F("hurwPulseDiff_el","uRWELL Pulse Distance from GEM Track; (Max GEM Pulse Pos) - (All uRW Pulse Positions) [mm]",75,-30.,30.); HistList->Add(hurwPulseDiff_el);
  hurwPulseDiff_mmg = new TH1F("hurwPulseDiff_mmg","uRWELL Pulse Distance from MMG Track; (Max MMG Pulse Pos) - (All uRW Pulse Positions) [mm]",75,-30.,30.); HistList->Add(hurwPulseDiff_mmg);
  hTrackDiff = new TH1F("hTrackDiff","Tracker 1 & 2 X Pos Difference; Distance [mm]",125,-25.,25.); HistList->Add(hTrackDiff);
  
  //----- Gain Calib ----
  #ifdef GAIN_CALIB
  hmmg1GainSum = new TH1F("hmmg1GainSum","MMG1 Pulse Cluster Sum; ADC Sum; Counts ",300,0.,15000); HistList->Add(hmmg1GainSum);
  hmmg1GainSpread = new TH1F("hmmg1GainSpread","MMG1 Pulse Cluster Spread; N Strips Fired; Counts ",12,-0.5,11.5); HistList->Add(hmmg1GainSpread);
  hmmg1GainMultiplicity = new TH2F("hmmg1GainMultiplicity","MMG1 2D Pulse Cluster Sum vs Strips; N Strips Fired; ADC Sum",12,-0.5,11.5,300,0.,15000); HistList->Add(hmmg1GainMultiplicity);
  hmmg12DGain = new TH2F("hmmg12DGain","MMG1 2D Pulse Cluster Amp. vs Position; MMG1 Channel# ; ADC Amplitude",240,-0.5,239.5,100,0.,4960); HistList->Add(hmmg12DGain);
  hmmg12DMaxGain = new TH2F("hmmg12DMaxGain","MMG1 2D Pulse Cluster Max Sum vs Max Position; MMG1 Max Hit Channel# ; Cluster ADC Sum",240,-0.5,239.5,300,0.,15000); HistList->Add(hmmg12DMaxGain);
  hmmg12DMaxGainSingle = new TH2F("hmmg12DMaxGainSingle","MMG1 2D Pulse Cluster Max Amp. vs Max Position; MMG1 Max Hit Channel# ; Cluster Max ADC Amp",240,-0.5,239.5,100,0.,4960); HistList->Add(hmmg12DMaxGainSingle);
  hmmg12DGainLateTime = new TH2F("hmmg12DGainLateTime","MMG1 2D Pulse Cluster Max Amp. vs Time; Drift Time (*8ns); ADC Amplitude",200,0.5,200.5,100,0.,4960); HistList->Add(hmmg12DGainLateTime);
  hmmg1GainMaxLate = new TH1F("hmmg1GainMaxLate","MMG1 Pulse Cluster Max Amp. in Latest Time; ADC Amplitude",100,0.,4960); HistList->Add(hmmg1GainMaxLate);
  hmmg12DGainSumLateTime = new TH2F("hmmg12DGainSumLateTime","MMG1 2D Pulse Cluster Sum vs Time; Drift Time (*8ns); ADC Amplitude Sum",200,0.5,200.5,300,0.,15000); HistList->Add(hmmg12DGainSumLateTime);
  hmmg1GainSumLate = new TH1F("hmmg1GainSumLate","MMG1 Pulse Cluster Sum Amp. in Latest Time; ADC Amplitude Sum",300,0.,15000); HistList->Add(hmmg1GainSumLate);
  
  hgemGainSum = new TH1F("hgemGainSum","Triple-GEM Pulse Cluster Sum; ADC Sum; Counts ",440,0.,22000); HistList->Add(hgemGainSum);
  hgemGainSpread = new TH1F("hgemGainSpread","Triple-GEM Pulse Cluster Spread; N Strips Fired; Counts ",12,-0.5,11.5); HistList->Add(hgemGainSpread);
  hgemGainMultiplicity = new TH2F("hgemGainMultiplicity","Triple-GEM 2D Pulse Cluster Sum vs Strips; N Strips Fired; ADC Sum",12,-0.5,11.5,440,0.,22000); HistList->Add(hgemGainMultiplicity);
  hgem2DGain = new TH2F("hgem2DGain","Triple-GEM 2D Pulse Cluster Amp. vs Position; GEM Channel# ; ADC Amplitude",240,-0.5,239.5,100,0.,4960); HistList->Add(hgem2DGain);
  hgem2DMaxGain = new TH2F("hgem2DMaxGain","Triple-GEM 2D Pulse Cluster Max Sum vs Max Position; GEM Max Hit Channel# ; Cluster ADC Sum",240,-0.5,239.5,440,0.,22000); HistList->Add(hgem2DMaxGain);
  hgem2DMaxGainSingle = new TH2F("hgem2DMaxGainSingle","Triple-GEM 2D Pulse Cluster Max Amp. vs Max Position; GEM Max Hit Channel# ; Cluster Max ADC Amp",240,-0.5,239.5,100,0.,4960); HistList->Add(hgem2DMaxGainSingle);
  hgem2DGainLateTime = new TH2F("hgem2DGainLateTime","Triple-GEM 2D Pulse Cluster Max Amp. vs Time; Drift Time (*8ns); ADC Amplitude",200,0.5,200.5,100,0.,4960); HistList->Add(hgem2DGainLateTime);
  hgemGainMaxLate = new TH1F("hgemGainMaxLate","Triple-GEM Pulse Cluster Max Amp. in Latest Time; ADC Amplitude",100,0.,4960); HistList->Add(hgemGainMaxLate);
  hgem2DGainSumLateTime = new TH2F("hgem2DGainSumLateTime","Triple-GEM 2D Pulse Cluster Sum vs Time; Drift Time (*8ns); ADC Amplitude Sum",200,0.5,200.5,440,0.,22000); HistList->Add(hgem2DGainSumLateTime);
  hgemGainSumLate = new TH1F("hgemGainSumLate","Triple-GEM Pulse Cluster Sum Amp. in Latest Time; ADC Amplitude Sum",440,0.,22000); HistList->Add(hgemGainSumLate);
  
  hurwXGainSum = new TH1F("hurwXGainSum","uRWELL X-Plane Pulse Cluster Sum; ADC Sum; Counts ",300,0.,15000); HistList->Add(hurwXGainSum);
  hurwXGainSpread = new TH1F("hurwXGainSpread","uRWELL X-Plane Pulse Cluster Spread; N Strips Fired; Counts ",12,-0.5,11.5); HistList->Add(hurwXGainSpread);
  hurwXGainMultiplicity = new TH2F("hurwXGainMultiplicity","uRWELL X-Plane 2D Pulse Cluster Sum vs Strips; N Strips Fired; ADC Sum",12,-0.5,11.5,300,0.,15000); HistList->Add(hurwXGainMultiplicity);
  hurwX2DGain = new TH2F("hurwX2DGain","uRWELL X-Plane 2D Pulse Cluster Amp. vs Position; uRWELL X Channel# ; ADC Amplitude",120,-0.5,119.5,100,0.,4960); HistList->Add(hurwX2DGain);
  hurwX2DMaxGain = new TH2F("hurwX2DMaxGain","uRWELL X-Plane 2D Pulse Cluster Max Sum vs Max Position; uRWELL Max Hit X Channel# ; Cluster ADC Sum",120,-0.5,119.5,300,0.,15000); HistList->Add(hurwX2DMaxGain);
  hurwX2DMaxGainSingle = new TH2F("hurwX2DMaxGainSingle","uRWELL X-Plane 2D Pulse Cluster Max Amp. vs Max Position; uRWELL Max Hit X Channel# ; Cluster Max ADC Amp",120,-0.5,119.5,100,0.,4960); HistList->Add(hurwX2DMaxGainSingle);
  hurwX2DGainLateTime = new TH2F("hurwX2DGainLateTime","uRWELL X-Plane 2D Pulse Cluster Max Amp. vs Time; Drift Time (*8ns); ADC Amplitude",200,0.5,200.5,100,0.,4960); HistList->Add(hurwX2DGainLateTime);
  hurwXGainMaxLate = new TH1F("hurwXGainMaxLate","uRWELL X-Plane Pulse Cluster Max Amp. in Latest Time; ADC Amplitude",100,0.,4960); HistList->Add(hurwXGainMaxLate);
  hurwX2DGainSumLateTime = new TH2F("hurwX2DGainSumLateTime","uRWELL X-Plane 2D Pulse Cluster Sum vs Time; Drift Time (*8ns); ADC Amplitude Sum",200,0.5,200.5,300,0.,15000); HistList->Add(hurwX2DGainSumLateTime);
  hurwXGainSumLate = new TH1F("hurwXGainSumLate","uRWELL X-Plane Pulse Cluster Sum Amp. in Latest Time; ADC Amplitude Sum",300,0.,15000); HistList->Add(hurwXGainSumLate);
  
  hurwYGainSum = new TH1F("hurwYGainSum","uRWELL Y-Plane Pulse Cluster Sum; ADC Sum; Counts ",300,0.,15000); HistList->Add(hurwYGainSum);
  hurwYGainSpread = new TH1F("hurwYGainSpread","uRWELL Y-Plane Pulse Cluster Spread; N Strips Fired; Counts ",12,-0.5,11.5); HistList->Add(hurwYGainSpread);
  hurwYGainMultiplicity = new TH2F("hurwYGainMultiplicity","uRWELL Y-Plane 2D Pulse Cluster Sum vs Strips; N Strips Fired; ADC Sum",12,-0.5,11.5,300,0.,15000); HistList->Add(hurwYGainMultiplicity);
  hurwY2DGain = new TH2F("hurwY2DGain","uRWELL Y-Plane 2D Pulse Cluster Amp. vs Position; uRWELL Y Channel# ; ADC Amplitude",120,-0.5,119.5,100,0.,4960); HistList->Add(hurwY2DGain);
  hurwY2DMaxGain = new TH2F("hurwY2DMaxGain","uRWELL Y-Plane 2D Pulse Cluster Max Sum vs Max Position; uRWELL Max Hit Y Channel# ; Cluster ADC Sum",120,-0.5,119.5,300,0.,15000); HistList->Add(hurwY2DMaxGain);
  hurwY2DMaxGainSingle = new TH2F("hurwY2DMaxGainSingle","uRWELL Y-Plane 2D Pulse Cluster Max Amp. vs Max Position; uRWELL Max Hit Y Channel# ; Cluster Max ADC Amp",120,-0.5,119.5,100,0.,4960); HistList->Add(hurwY2DMaxGainSingle);
  hurwY2DGainLateTime = new TH2F("hurwY2DGainLateTime","uRWELL Y-Plane 2D Pulse Cluster Max Amp. vs Time; Drift Time (*8ns); ADC Amplitude",200,0.5,200.5,100,0.,4960); HistList->Add(hurwY2DGainLateTime);
  hurwYGainMaxLate = new TH1F("hurwYGainMaxLate","uRWELL Y-Plane Pulse Cluster Max Amp. in Latest Time; ADC Amplitude",100,0.,4960); HistList->Add(hurwYGainMaxLate);
  hurwY2DGainSumLateTime = new TH2F("hurwY2DGainSumLateTime","uRWELL Y-Plane 2D Pulse Cluster Sum vs Time; Drift Time (*8ns); ADC Amplitude Sum",200,0.5,200.5,300,0.,15000); HistList->Add(hurwY2DGainSumLateTime);
  hurwYGainSumLate = new TH1F("hurwYGainSumLate","uRWELL Y-Plane Pulse Cluster Sum Amp. in Latest Time; ADC Amplitude Sum",300,0.,15000); HistList->Add(hurwYGainSumLate);
  #endif
  //=========================================
  
  TFile* fHits;
  #ifdef SAVE_TRACK_HITS
    #if ANALYZE_MERGED
    char hitsFileName[256]; sprintf(hitsFileName, "RootOutput/ps26/merged/trd_singleTrackHits_Run_%06d_%06dEntries.root",RunNum,nEntries);
    #else
    char hitsFileName[256]; sprintf(hitsFileName, "RootOutput/ps26/trd_singleTrackHits_Run_%06d.root", RunNum);
    #endif
    fHits = new TFile(hitsFileName, "RECREATE");
    tgem_zHist = new TH1F("tgem_zHist", "tgem_zHist", 20, 40., 160.);
    mmg1_zHist = new TH1F("mmg1_zHist", "mmg1_zHist", 20, 50., 170.);
    urw_zHist = new TH1F("urw_zHist", "urw_zHist", 20, 40., 160.);
    //-- Triple GEM-TRD
    EVENT_VECT_GEM = new TTree("gem_hits","Triple-GEM TTree with single track hit info");
    EVENT_VECT_GEM->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_GEM->Branch("nhit",&tgem_nhit,"tgem_nhit/I");
    EVENT_VECT_GEM->Branch("nclu",&tgem_nclu,"tgem_nclu/I");
    EVENT_VECT_GEM->Branch("ntracks",&tgem_ntracks,"tgem_ntracks/I");
    EVENT_VECT_GEM->Branch("xpos",&tgem_xpos);
    EVENT_VECT_GEM->Branch("zpos",&tgem_zpos);
    EVENT_VECT_GEM->Branch("dedx",&tgem_dedx);
    EVENT_VECT_GEM->Branch("parID",&tgem_parID);
    EVENT_VECT_GEM->Branch("zHist",&tgem_zHist_vect);
    EVENT_VECT_GEM->Branch("xposc",&clu_xpos);
    EVENT_VECT_GEM->Branch("zposc",&clu_zpos);
    EVENT_VECT_GEM->Branch("dedxc",&clu_dedx);
    EVENT_VECT_GEM->Branch("widthc",&clu_width);
    EVENT_VECT_GEM->Branch("xposc_max",&clu_xpos_max,"clu_xpos_max/f");
    EVENT_VECT_GEM->Branch("zposc_max",&clu_zpos_max,"clu_zpos_max/f");
    EVENT_VECT_GEM->Branch("dedxc_max",&clu_dedx_max,"clu_dedx_max/f");
    EVENT_VECT_GEM->Branch("widthc_max",&clu_width_max,"clu_width_max/f");
    EVENT_VECT_GEM->Branch("dedxc_tot",&clu_dedx_tot,"clu_dedx_tot/f");
    EVENT_VECT_GEM->Branch("xch_max",&tgem_xch_max,"tgem_xch_max/f");
    EVENT_VECT_GEM->Branch("amp_max",&tgem_amp_max,"tgem_amp_max/f");
    EVENT_VECT_GEM->Branch("time_max",&tgem_time_max,"tgem_time_max/f");
    EVENT_VECT_GEM->Branch("chi2",&tgem_chi2cc);
    EVENT_VECT_GEM->Branch("Fint",&tgem_integral);
    EVENT_VECT_GEM->Branch("a0",&a0);
    EVENT_VECT_GEM->Branch("a1",&a1);
    
    //-- MMG1-TRD
    EVENT_VECT_MMG1 = new TTree("mmg1_hits","MMG1 TTree with single track hit info");
    EVENT_VECT_MMG1->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_MMG1->Branch("nhit",&mmg1_nhit,"mmg1_nhit/I");
    EVENT_VECT_MMG1->Branch("nclu",&mmg1_nclu,"mmg1_nclu/I");
    EVENT_VECT_MMG1->Branch("ntracks",&mmg1_ntracks,"mmg1_ntracks/I");
    EVENT_VECT_MMG1->Branch("xpos",&mmg1_xpos);
    EVENT_VECT_MMG1->Branch("zpos",&mmg1_zpos);
    EVENT_VECT_MMG1->Branch("dedx",&mmg1_dedx);
    EVENT_VECT_MMG1->Branch("parID",&mmg1_parID);
    EVENT_VECT_MMG1->Branch("zHist",&mmg1_zHist_vect);
    EVENT_VECT_MMG1->Branch("xposc",&mmg1_clu_xpos);
    EVENT_VECT_MMG1->Branch("zposc",&mmg1_clu_zpos);
    EVENT_VECT_MMG1->Branch("dedxc",&mmg1_clu_dedx);
    EVENT_VECT_MMG1->Branch("widthc",&mmg1_clu_width);
    EVENT_VECT_MMG1->Branch("xposc_max",&mmg1_clu_xpos_max,"mmg1_clu_xpos_max/f");
    EVENT_VECT_MMG1->Branch("zposc_max",&mmg1_clu_zpos_max,"mmg1_clu_zpos_max/f");
    EVENT_VECT_MMG1->Branch("dedxc_max",&mmg1_clu_dedx_max,"mmg1_clu_dedx_max/f");
    EVENT_VECT_MMG1->Branch("widthc_max",&mmg1_clu_width_max,"mmg1_clu_width_max/f");
    EVENT_VECT_MMG1->Branch("dedxc_tot",&mmg1_clu_dedx_tot,"mmg1_clu_dedx_tot/f");
    EVENT_VECT_MMG1->Branch("xch_max",&mmg1_xch_max,"mmg1_xch_max/f");
    EVENT_VECT_MMG1->Branch("amp_max",&mmg1_amp_max,"mmg1_amp_max/f");
    EVENT_VECT_MMG1->Branch("time_max",&mmg1_time_max,"mmg1_time_max/f");
    EVENT_VECT_MMG1->Branch("chi2",&mmg1_chi2cc);
    EVENT_VECT_MMG1->Branch("Fint",&mmg1_integral);
    EVENT_VECT_MMG1->Branch("a0",&mmg1_a0);
    EVENT_VECT_MMG1->Branch("a1",&mmg1_a1);
    
     //-- uRWELL-TRD
    EVENT_VECT_URW = new TTree("urw_hits","uRWELL TTree with single track hit info");
    EVENT_VECT_URW->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_URW->Branch("nhit",&urw_nhit,"urw_nhit/I");
    EVENT_VECT_URW->Branch("nyhit",&urw_nyhit,"urw_nyhit/I");
    EVENT_VECT_URW->Branch("nclu",&urw_nclu,"urw_nclu/I");
    EVENT_VECT_URW->Branch("ntracks",&urw_ntracks,"urw_ntracks/I");
    EVENT_VECT_URW->Branch("xpos",&urw_xpos);
    EVENT_VECT_URW->Branch("zpos",&urw_zpos);
    EVENT_VECT_URW->Branch("dedx",&urw_dedx);
    EVENT_VECT_URW->Branch("parID",&urw_parID);
    EVENT_VECT_URW->Branch("zHist",&urw_zHist_vect);
    EVENT_VECT_URW->Branch("xposc",&urw_clu_xpos);
    EVENT_VECT_URW->Branch("zposc",&urw_clu_zpos);
    EVENT_VECT_URW->Branch("dedxc",&urw_clu_dedx);
    EVENT_VECT_URW->Branch("widthc",&urw_clu_width);
    EVENT_VECT_URW->Branch("xposc_max",&urw_clu_xpos_max,"urw_clu_xpos_max/f");
    EVENT_VECT_URW->Branch("zposc_max",&urw_clu_zpos_max,"urw_clu_zpos_max/f");
    EVENT_VECT_URW->Branch("dedxc_max",&urw_clu_dedx_max,"urw_clu_dedx_max/f");
    EVENT_VECT_URW->Branch("widthc_max",&urw_clu_width_max,"urw_clu_width_max/f");
    EVENT_VECT_URW->Branch("dedxc_tot",&urw_clu_dedx_tot,"urw_clu_dedx_tot/f");
    EVENT_VECT_URW->Branch("xch_max",&urw_xch_max,"urw_xch_max/f");
    EVENT_VECT_URW->Branch("amp_max",&urw_amp_max,"urw_amp_max/f");
    EVENT_VECT_URW->Branch("time_max",&urw_time_max,"urw_time_max/f");
    EVENT_VECT_URW->Branch("chi2",&urw_chi2cc);
    EVENT_VECT_URW->Branch("Fint",&urw_integral);
    EVENT_VECT_URW->Branch("a0",&urw_a0);
    EVENT_VECT_URW->Branch("a1",&urw_a1);
    
  #endif
  
  TStopwatch timer;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes=0, nb=0;
  if (MaxEvt>0) nentries=MaxEvt;  //-- limit number of events for test
  
  //==================================================================================================
  //                      E v e n t    L o o p
  //==================================================================================================
  printf("===============  Begin Event Loop - 1st evt=%lld, Last evt=%lld =============== \n",FirstEvt,MaxEvt);
  timer.Start();
  Long64_t jentry=0;
  
  for (jentry=FirstEvt; jentry<nentries; jentry++) { //-- Event Loop --
    Count("EVT");
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (!(jentry%NPRT))
      printf("------- evt=%llu  f125_raw_count=%llu f125_pulse_count=%llu srs_peak_count=%llu \n",jentry,f125_wraw_count, f125_pulse_count, gem_peak_count);
    event_num=jentry;
    
    bool match=false, match_mmg1=false, match_urw=false;
    bool trackFound_gem=false, trackFound_mmg1=false, trackFound_urw=false;
    bool GoodYCorr=false, GoodXCorr=false;
    
    tgem_nhit=0;
    mmg1_nhit=0;
    urw_nhit=0;
    urw_nyhit=0;
    tgem_nclu=0;
    mmg1_nclu=0;
    urw_nclu=0;
    tgem_tmp_nhit=0;
    mmg1_tmp_nhit=0;
    urw_tmp_nhit=0;
    
    //-- Triple GEM-TRD
    tgem_xpos.clear();
    tgem_zpos.clear();
    tgem_dedx.clear();
    tgem_parID.clear();
    tgem_zHist->Reset();
    tgem_zHist_vect.clear();
    clu_xpos.clear();
    clu_zpos.clear();
    clu_dedx.clear();
    clu_width.clear();
    clu_xpos_max=0;
    clu_zpos_max=0;
    clu_dedx_max=0;
    clu_width_max=0;
    clu_dedx_tot=0;
    tgem_amp_max=0;
    tgem_xch_max=0;
    tgem_time_max=0;
    tgem_chi2cc.clear();
    tgem_integral.clear();
    double a0 = -1., a1=-1.;
    
    //-- MMG1-TRD
    mmg1_xpos.clear();
    mmg1_zpos.clear();
    mmg1_dedx.clear();
    mmg1_parID.clear();
    mmg1_zHist->Reset();
    mmg1_zHist_vect.clear();
    mmg1_clu_xpos.clear();
    mmg1_clu_zpos.clear();
    mmg1_clu_dedx.clear();
    mmg1_clu_width.clear();
    mmg1_clu_xpos_max=0;
    mmg1_clu_zpos_max=0;
    mmg1_clu_dedx_max=0;
    mmg1_clu_width_max=0;
    mmg1_clu_dedx_tot=0;
    mmg1_amp_max=0;
    mmg1_xch_max=0;
    mmg1_time_max=0;
    mmg1_chi2cc.clear();
    mmg1_integral.clear();
    double a0_mmg1 = -1., a1_mmg1=-1.;
    
    //-- uRWELL-TRD
    urw_xpos.clear();
    urw_zpos.clear();
    urw_dedx.clear();
    urw_parID.clear();
    urw_zHist->Reset();
    urw_zHist_vect.clear();
    urw_clu_xpos.clear();
    urw_clu_zpos.clear();
    urw_clu_dedx.clear();
    urw_clu_width.clear();
    urw_clu_xpos_max=0;
    urw_clu_zpos_max=0;
    urw_clu_dedx_max=0;
    urw_clu_width_max=0;
    urw_clu_dedx_tot=0;
    urw_amp_max=0;
    urw_xch_max=0;
    urw_time_max=0;
    urw_chi2cc.clear();
    urw_integral.clear();
    double a0_urw = -1., a1_urw=-1.;
    
    //==================================================================================================
    //                    Process SRS data
    //==================================================================================================
    
    //===========================================================
    //  GEMTracker (SRS) Correlations with TRD Prototypes
    //===========================================================
    
    ULong64_t gt_1_idx_x = 0, gt_1_idx_y=0, gt_2_idx_x = 0, gt_2_idx_y = 0, mmg1_idx_y=0, tgem_idx_y=0;
    int gt1_nhit=0, gt2_nhit=0;
    
    double gemtrkr1_xamp_max=-1., gemtrkr1_xch_max=-1000.;
    double gemtrkr1_yamp_max=-1., gemtrkr1_ych_max=-1000.;
    double gemtrkr2_xamp_max=-1., gemtrkr2_xch_max=-1000.;
    double gemtrkr2_yamp_max=-1., gemtrkr2_ych_max=-1000.;
    double tgem_yamp_max=-1., tgem_ych_max=-1000.;
    double mmg1_yamp_max=-1., mmg1_ych_max=-1000.;
    
    double gemtrkr_1_peak_pos_y[gem_peak_count];
    double gemtrkr_1_peak_pos_x[gem_peak_count];
    double gemtrkr_1_peak_x_height[gem_peak_count];
    double gemtrkr_1_peak_y_height[gem_peak_count];
    double gemtrkr_2_peak_pos_y[gem_peak_count];
    double gemtrkr_2_peak_pos_x[gem_peak_count];
    double gemtrkr_2_peak_x_height[gem_peak_count];
    double gemtrkr_2_peak_y_height[gem_peak_count];
    double mmg1_peak_pos_y[gem_peak_count];
    double mmg1_peak_y_height[gem_peak_count];
    double tgem_peak_pos_y[gem_peak_count];
    double tgem_peak_y_height[gem_peak_count];
    
    #if USE_TRD_EXT_TRACK
    //Skip SRS info
    #else
    for (ULong64_t i=0; i<gem_peak_count; i++) {
      gemtrkr_1_peak_pos_y[i] = -1000;
      gemtrkr_1_peak_pos_x[i] = -1000;
      gemtrkr_1_peak_x_height[i] = -1000;
      gemtrkr_1_peak_y_height[i] = -1000;
      gemtrkr_2_peak_pos_y[i] = -1000;
      gemtrkr_2_peak_pos_x[i] = -1000;
      gemtrkr_2_peak_x_height[i] = -1000;
      gemtrkr_2_peak_y_height[i] = -1000;
      mmg1_peak_pos_y[i] = -1000;
      mmg1_peak_y_height[i] = -1000;
      tgem_peak_pos_y[i] = -1000;
      tgem_peak_y_height[i] = -1000;
    }
   
    for (ULong64_t i=0; i<gem_peak_count; i++) { //-- SRS Peaks Loop
      
      if (gem_peak_plane_name->at(i) == "GEMTR1X") {
        gemtrkr_1_peak_x_height[gt_1_idx_x] = gem_peak_height->at(i);
        if (gemtrkr_1_peak_x_height[gt_1_idx_x]>TRKR_THRESH) {
          gemtrkr_1_peak_pos_x[gt_1_idx_x] = gem_peak_real_pos->at(i);
          if (gemtrkr_1_peak_pos_x[gt_1_idx_x]<=0) gemtrkr_1_peak_pos_x[gt_1_idx_x]+=51.2; else if (gemtrkr_1_peak_pos_x[gt_1_idx_x]>0) gemtrkr_1_peak_pos_x[gt_1_idx_x]-=51.2; gemtrkr_1_peak_pos_x[gt_1_idx_x]*=-1.; gemtrkr_1_peak_pos_x[gt_1_idx_x]+=51.2;
          gt_1_idx_x++; Count("gt1_x");
        }
      } else if (gem_peak_plane_name->at(i) == "GEMTR1Y") {
          gemtrkr_1_peak_y_height[gt_1_idx_y] = gem_peak_height->at(i);
          if (gemtrkr_1_peak_y_height[gt_1_idx_y]>TRKR_THRESH) {
            gemtrkr_1_peak_pos_y[gt_1_idx_y] = gem_peak_real_pos->at(i);
            if (gemtrkr_1_peak_pos_y[gt_1_idx_y]<=0) gemtrkr_1_peak_pos_y[gt_1_idx_y]+=51.2; else if (gemtrkr_1_peak_pos_y[gt_1_idx_y]>0) gemtrkr_1_peak_pos_y[gt_1_idx_y]-=51.2;  gemtrkr_1_peak_pos_y[gt_1_idx_y]*=-1.;  gemtrkr_1_peak_pos_y[gt_1_idx_y]+=51.2;
            gt_1_idx_y++; Count("gt1_y");
          }
      } else if (gem_peak_plane_name->at(i) == "GEMTR2X") {
        gemtrkr_2_peak_x_height[gt_2_idx_x] = gem_peak_height->at(i);
        if (gemtrkr_2_peak_x_height[gt_2_idx_x]>TRKR_THRESH) {
          gemtrkr_2_peak_pos_x[gt_2_idx_x] = gem_peak_real_pos->at(i);
          if (gemtrkr_2_peak_pos_x[gt_2_idx_x]<=0) gemtrkr_2_peak_pos_x[gt_2_idx_x]+=51.2; else if (gemtrkr_2_peak_pos_x[gt_2_idx_x]>0) gemtrkr_2_peak_pos_x[gt_2_idx_x]-=51.2;  gemtrkr_2_peak_pos_x[gt_2_idx_x]*=-1.;  gemtrkr_2_peak_pos_x[gt_2_idx_x]+=51.2;
          gt_2_idx_x++; Count("gt2_x");
        }
      } else if (gem_peak_plane_name->at(i) == "GEMTR2Y") {
          gemtrkr_2_peak_y_height[gt_2_idx_y] = gem_peak_height->at(i);
          if (gemtrkr_2_peak_y_height[gt_2_idx_y]>TRKR_THRESH) {
            gemtrkr_2_peak_pos_y[gt_2_idx_y] = gem_peak_real_pos->at(i);
            if (gemtrkr_2_peak_pos_y[gt_2_idx_y]<=0) gemtrkr_2_peak_pos_y[gt_2_idx_y]+=51.2; else if (gemtrkr_2_peak_pos_y[gt_2_idx_y]>0) gemtrkr_2_peak_pos_y[gt_2_idx_y]-=51.2;  gemtrkr_2_peak_pos_y[gt_2_idx_y]*=-1.;  gemtrkr_2_peak_pos_y[gt_2_idx_y]+=51.2;
            gt_2_idx_y++; Count("gt2_y");
        }
      } else if (gem_peak_plane_name->at(i) == "MMG1TRDY") {
        mmg1_peak_y_height[mmg1_idx_y] = gem_peak_height->at(i);
        if (mmg1_peak_y_height[mmg1_idx_y]>TRKR_THRESH+600.) {
          mmg1_peak_pos_y[mmg1_idx_y] = gem_peak_real_pos->at(i);
          if (mmg1_peak_pos_y[mmg1_idx_y]<=0) mmg1_peak_pos_y[mmg1_idx_y]+=51.2; else if (mmg1_peak_pos_y[mmg1_idx_y]>0) mmg1_peak_pos_y[mmg1_idx_y]-=51.2;  mmg1_peak_pos_y[mmg1_idx_y]*=-1.;  mmg1_peak_pos_y[mmg1_idx_y]+=51.2;
          mmg1_idx_y++; Count("mmg1_y");
        }
      } else if (gem_peak_plane_name->at(i) == "VU_GEMTRDY") {
        tgem_peak_y_height[tgem_idx_y] = gem_peak_height->at(i);
        if (tgem_peak_y_height[tgem_idx_y]>TRKR_THRESH+600.) {
          tgem_peak_pos_y[tgem_idx_y] = gem_peak_real_pos->at(i);
          if (tgem_peak_pos_y[tgem_idx_y]<=0) tgem_peak_pos_y[tgem_idx_y]+=51.2; else if (tgem_peak_pos_y[tgem_idx_y]>0) tgem_peak_pos_y[tgem_idx_y]-=51.2;  tgem_peak_pos_y[tgem_idx_y]*=-1.;  tgem_peak_pos_y[tgem_idx_y]+=51.2;
          tgem_idx_y++; Count("tgem_y");
        }
      }
    } //--End SRS peaks loop
    
    for (ULong64_t j=0; j<gt_1_idx_y; j++) {
      if (gemtrkr_1_peak_pos_y[j]>-1.) {
        hgemtrkr_1_peak_y->Fill(gemtrkr_1_peak_pos_y[j]);
        for (ULong64_t k=0; k<gt_1_idx_x; k++) {
          if (gemtrkr_1_peak_pos_x[k]>-1.) {
            hgemtrkr_1_peak_xy->Fill(gemtrkr_1_peak_pos_x[k], gemtrkr_1_peak_pos_y[j]);
            gt1_nhit++;
          }
        }
      }
      hgemtrkr_1_peak_y_height->Fill(gemtrkr_1_peak_y_height[j]);
      if (gemtrkr_1_peak_y_height[j]>gemtrkr1_yamp_max) {
        gemtrkr1_yamp_max=gemtrkr_1_peak_y_height[j];
        gemtrkr1_ych_max=gemtrkr_1_peak_pos_y[j];
      }
    }
    
    for (ULong64_t j=0; j<gt_2_idx_y; j++) {
      if (gemtrkr_2_peak_pos_y[j]>-1.) {
        hgemtrkr_2_peak_y->Fill(gemtrkr_2_peak_pos_y[j]);
        for (ULong64_t k=0; k<gt_2_idx_x; k++) {
          if (gemtrkr_2_peak_pos_x[k]>-1.) {
            hgemtrkr_2_peak_xy->Fill(gemtrkr_2_peak_pos_x[k], gemtrkr_2_peak_pos_y[j]);
            gt2_nhit++;
          }
        }
      }
      hgemtrkr_2_peak_y_height->Fill(gemtrkr_2_peak_y_height[j]);
      if (gemtrkr_2_peak_y_height[j]>gemtrkr2_yamp_max) {
        gemtrkr2_yamp_max=gemtrkr_2_peak_y_height[j];
        gemtrkr2_ych_max=gemtrkr_2_peak_pos_y[j];
      }
    }
    
    for (ULong64_t k=0; k<gt_1_idx_x; k++) {
      if (gemtrkr_1_peak_pos_x[k]>-1.) hgemtrkr_1_peak_x->Fill(gemtrkr_1_peak_pos_x[k]);
      hgemtrkr_1_peak_x_height->Fill(gemtrkr_1_peak_x_height[k]);
      if (gemtrkr_1_peak_x_height[k]>gemtrkr1_xamp_max) {
        gemtrkr1_xamp_max=gemtrkr_1_peak_x_height[k];
        gemtrkr1_xch_max=gemtrkr_1_peak_pos_x[k];
      }
    }
    
    for (ULong64_t k=0; k<gt_2_idx_x; k++) {
      if (gemtrkr_2_peak_pos_x[k]>-1.) hgemtrkr_2_peak_x->Fill(gemtrkr_2_peak_pos_x[k]);
      hgemtrkr_2_peak_x_height->Fill(gemtrkr_2_peak_x_height[k]);
      if (gemtrkr_2_peak_x_height[k]>gemtrkr2_xamp_max) {
        gemtrkr2_xamp_max=gemtrkr_2_peak_x_height[k];
        gemtrkr2_xch_max=gemtrkr_2_peak_pos_x[k];
      }
    }
    
    for (ULong64_t k=0; k<mmg1_idx_y; k++) {
      if (mmg1_peak_pos_y[k]>-1.) mmg1_peak_y->Fill(mmg1_peak_pos_y[k]);
      hmmg1_peak_y_height->Fill(mmg1_peak_y_height[k]);
      if (mmg1_peak_y_height[k]>mmg1_yamp_max) {
        mmg1_yamp_max=mmg1_peak_y_height[k];
        mmg1_ych_max=mmg1_peak_pos_y[k];
      }
    }
    
    for (ULong64_t k=0; k<tgem_idx_y; k++) {
      if (tgem_peak_pos_y[k]>-1.) tgem_peak_y->Fill(tgem_peak_pos_y[k]);
      htgem_peak_y_height->Fill(tgem_peak_y_height[k]);
      if (tgem_peak_y_height[k]>tgem_yamp_max) {
        tgem_yamp_max=tgem_peak_y_height[k];
        tgem_ych_max=tgem_peak_pos_y[k];
      }
    }  
    
    if (tgem_ych_max>0. && gemtrkr1_ych_max>0.) hgemtrkr_1_tgem->Fill(tgem_ych_max, gemtrkr1_ych_max);
    if (mmg1_ych_max>0. && gemtrkr1_ych_max>0.) hgemtrkr_1_mmg1->Fill(mmg1_ych_max, gemtrkr1_ych_max);
    if (gemtrkr1_xch_max>0.) hgemtrkr_1_max_xch->Fill(gemtrkr1_xch_max);
    if (gemtrkr1_xamp_max>0.) { 
      hgemtrkr_1_max_xamp->Fill(gemtrkr1_xamp_max);
      if (gemtrkr2_xamp_max>0.) {
        hgemtrkr_max_xcorr->Fill(gemtrkr1_xch_max, gemtrkr2_xch_max);
        hgemtrkr_1D_xcorr->Fill((gemtrkr2_xch_max-gemtrkr1_xch_max*1.18854)+18.8911);
        if (RunNum>=8228 && RunNum<=8244) { //TRKR1 APV card issue
          if (abs((gemtrkr2_xch_max-gemtrkr1_xch_max*1.18854)+18.8911)<6. && gemtrkr1_xch_max<51.2) { GoodXCorr = true; }
        } else {
          if (abs((gemtrkr2_xch_max-gemtrkr1_xch_max*1.18854)+18.8911)<6.) { GoodXCorr = true; }
        }
      }
    }
    if (gemtrkr2_xch_max>0.) hgemtrkr_2_max_xch->Fill(gemtrkr2_xch_max);
    if (gemtrkr2_xamp_max>0.) hgemtrkr_2_max_xamp->Fill(gemtrkr2_xamp_max);
    if (gemtrkr1_xch_max>0. && gemtrkr1_ych_max>0.) hgemtrkr_1_max_xy->Fill(gemtrkr1_xch_max, gemtrkr1_ych_max);
    if (gemtrkr2_xch_max>0. && gemtrkr2_ych_max>0.) hgemtrkr_2_max_xy->Fill(gemtrkr2_xch_max, gemtrkr2_ych_max);
    if (gemtrkr1_yamp_max>0. && gemtrkr2_yamp_max>0.) {
      hgemtrkr_max_ycorr->Fill(gemtrkr1_ych_max, gemtrkr2_ych_max);
      hgemtrkr_1D_ycorr->Fill((gemtrkr2_ych_max-gemtrkr1_ych_max*1.25149)+5.35196);
      if (abs((gemtrkr2_ych_max-gemtrkr1_ych_max*1.25149)+5.35196)<5.) { GoodYCorr = true; }
    }
    
    if (mmg1_yamp_max>TRKR_THRESH+600. && tgem_yamp_max>TRKR_THRESH+600.) {
      double yCorrection = ymmg1 - ygem; //1.03
      hmmg1_tgem_ydiff->Fill(tgem_ych_max-mmg1_ych_max-yCorrection);
      //if (abs(tgem_ych_max-mmg1_ych_max-yCorrection)<5) GoodYCorr = true;
    }
    
    #endif
    
    //=============== END SRS Data Processing & Correlations ==============
    
    //==================================================================================================
    //                    Process Fa125  Pulse  data
    //==================================================================================================
    
    #if 1
      
      f125_fit->Reset();
      mmg1_f125_fit->Reset();
      urw_f125_fit->Reset();
      #if (USE_PULSE>0)
        hevtk->Reset();
        hevtck->Reset();
      #endif
      
      double chi2cc_tgem=-999., chi2cc_mmg1=-999.,  chi2cc_urw=-999;
      double integral_tgem=0., integral_mmg1=0., integral_urw=0.;
      double tgem_ampmax=-1., mmg1_ampmax=-1., urw_xampmax=-1., urw_yampmax=-1., tgem_xchanmax=-1, mmg1_xchanmax=-1, urw_xchanmax=-1.,  urw_ychanmax=-1.;
      double tgem_ampmax_x=-1., mmg1_ampmax_x=-1., urw_xampmax_x=-1., urw_yampmax_x=-1.;
      int tgem_timemax=0, mmg1_timemax=0, urw_xtimemax=0, urw_ytimemax=0;
      int tgem_x_timemax=-1, mmg1_x_timemax=-1, urw_x_timemax=-1, urw_y_timemax=-1;
      ULong64_t tgem_idx_x=0, mmg1_idx_x=0, urw_idx_x=0, urw_idx_y=0;
      double tgem_pos_x[f125_pulse_count], mmg1_pos_x[f125_pulse_count], urw_pos_x[f125_pulse_count], urw_pos_y[f125_pulse_count], tgem_amp_x[f125_pulse_count], mmg1_amp_x[f125_pulse_count], urw_amp_x[f125_pulse_count], urw_amp_y[f125_pulse_count];
      int tgem_time_x[f125_pulse_count], mmg1_time_x[f125_pulse_count], urw_time_x[f125_pulse_count];
      int gem_trk_hit=0, mmg1_trk_hit=0, urw_trk_hit=0;
      int tgem_channel_max=-1, mmg1_channel_max=-1, urw_channelX_max=-1, urw_channelY_max=-1;
      
      
      //==============================================
      //        First fa125 pulse loop (Find max hit position in TRDs, convert to mm, for X correlations)
      //==============================================
      
      #ifdef GAIN_CALIB
        float mmg1GainAmps[240];
        float gemGainAmps[240];
        float urwXGainAmps[120];
        float urwYGainAmps[120];
        for (ULong64_t i=0; i<240; i++) {
          mmg1GainAmps[i] = -1.;
          gemGainAmps[i] = -1.;
          if (i<120) {
            urwXGainAmps[i] = -1.;
            urwYGainAmps[i] = -1.;
          }
        }
      #endif
      
      for (ULong64_t i=0; i<f125_pulse_count; i++) { //--- Fadc125 Pulse Loop
        
        float peak_amp = f125_pulse_peak_amp->at(i);
        float ped = f125_pulse_pedestal->at(i);
       	if (0 > ped || ped > 200 ) ped = 100;
       	float amp=peak_amp-ped;
       	if (amp<0) amp=0;
       	float time=f125_pulse_peak_time->at(i);
       	int fADCSlot = f125_pulse_slot->at(i);
       	int fADCChan = f125_pulse_channel->at(i);
        if (time>177.) continue;
        
        int tripGemChan = Get3GEMChan(fADCChan, fADCSlot, RunNum);
        float tripGemChan_x = tripGemChan*0.4+3.2; // to [mm]
       	int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
       	float mmg1Chan_x = mmg1Chan*0.4+3.2; // to [mm]
        int urwXChan = GetRWELLXChan(fADCChan, fADCSlot, RunNum);
       	float urwChan_x = urwXChan*0.8+3.2; // to [mm]
        int urwYChan = GetRWELLYChan(fADCChan, fADCSlot, RunNum);
       	float urwChan_y = urwYChan*0.8+3.2; // to [mm]
        
       	if (tripGemChan>-1) {
          amp = Get3GEMCalib(amp, tripGemChan, RunNum);
          if (amp>TGEM_THRESH) {
          #ifdef GAIN_CALIB
            gemGainAmps[tripGemChan] = amp;
          #endif
          f125_el_amp2ds->Fill(time,tripGemChan,amp);
          tgem_tmp_nhit++;
          if (tgem_ampmax_x<amp) {
            tgem_ampmax_x=amp;
            tgem_xchanmax=tripGemChan_x;
            tgem_x_timemax=time;
          }
          }
        }
        else if (mmg1Chan>-1) {
          amp = GetMMGCalib(amp, mmg1Chan, RunNum);
          if (amp>MMG1_THRESH) {
          #ifdef GAIN_CALIB
            mmg1GainAmps[mmg1Chan] = amp;
          #endif
          mmg1_f125_el_amp2ds->Fill(time,mmg1Chan,amp);
          mmg1_tmp_nhit++;
          if (mmg1_ampmax_x<amp) {
            mmg1_ampmax_x=amp;
            mmg1_xchanmax=mmg1Chan_x;
            mmg1_x_timemax=time;
          }
          }
        }
        else if (amp>URW_THRESH && urwXChan>-1) {
          #ifdef GAIN_CALIB
            urwXGainAmps[urwXChan] = amp; 
          #endif
          urw_f125_x_amp2ds->Fill(time,urwXChan,amp);
          urw_tmp_nhit++;
          if (urw_xampmax_x<amp) {
            urw_xampmax_x=amp;
            urw_xchanmax=urwChan_x;
          }
        }
        else if (amp>URW_THRESH && urwYChan>-1) {
          #ifdef GAIN_CALIB
            urwYGainAmps[urwYChan] = amp;
          #endif
          urw_f125_y_amp2ds->Fill(time,urwYChan,amp);
          if (urw_yampmax_x<amp) {
            urw_yampmax_x=amp;
            urw_ychanmax=urwChan_y;
          }
        }
      } //--END first f125 pulse loop
      
      if (tgem_tmp_nhit>0.) htgem_tmp_nhits->Fill(tgem_tmp_nhit);
      if (mmg1_tmp_nhit>0.) hmmg1_tmp_nhits->Fill(mmg1_tmp_nhit);
      if (urw_tmp_nhit>0.) hurw_tmp_nxhits->Fill(urw_tmp_nhit);
      
      float um_slope=1., ug_slope=1., gm_slope=1.;
      float um_offset=0., ug_offset=0., gm_offset=0.;
      //--Initially set the same corrections
      if (RunNum<argonRunStop || RunNum>argonRunStart) {
        um_slope = 0.88487, um_offset = 11.1783;
        ug_slope = 0.93706, ug_offset = 5.2953;
        gm_slope = 0.94538, gm_offset = 6.1718;
      } else { //Xe run
        um_slope = 0.846771, um_offset = 12.9725;
        ug_slope = 0.895908, ug_offset = 7.2492;
        gm_slope = 0.918727, gm_offset = 7.33594;
      }
      
      if (urw_xampmax_x>URW_THRESH && mmg1_ampmax_x>MMG1_THRESH) hmmg1_urw_xdiff->Fill((mmg1_xchanmax-urw_xchanmax*um_slope)-um_offset);
      if (urw_xampmax_x>URW_THRESH && tgem_ampmax_x>TGEM_THRESH) htgem_urw_xdiff->Fill((tgem_xchanmax-urw_xchanmax*ug_slope)-ug_offset);
      if (tgem_ampmax_x>TGEM_THRESH && mmg1_ampmax_x>MMG1_THRESH) hmmg1_tgem_xdiff->Fill((mmg1_xchanmax-tgem_xchanmax*gm_slope)-gm_offset);
      
      #ifdef GAIN_CALIB
        float mmg1AmpSum=0.;
        int mmg1StripSum=0;
        
        if (mmg1_ampmax_x>FE55_THRESH) {
          for (int i=0; i<10; i++) {
            if (mmg1GainAmps[(int)((mmg1_xchanmax-3.2)/0.4) + i]>MMG1_THRESH) {
              mmg1AmpSum+=mmg1GainAmps[(int)((mmg1_xchanmax-3.2)/0.4) + i];
              mmg1StripSum++;
            } else {break;}
          }
          for (int i=1; i<10; i++) {
            if (mmg1GainAmps[(int)((mmg1_xchanmax-3.2)/0.4) - i]>MMG1_THRESH) {
              mmg1AmpSum+=mmg1GainAmps[(int)((mmg1_xchanmax-3.2)/0.4) - i];
              mmg1StripSum++;
            } else {break;}
          }
          if (mmg1AmpSum>0. && mmg1StripSum>1 && mmg1StripSum<10) {
            hmmg1GainSum->Fill(mmg1AmpSum);
            hmmg1GainSpread->Fill(mmg1StripSum);
            hmmg1GainMultiplicity->Fill(mmg1StripSum,mmg1AmpSum);
            hmmg12DMaxGain->Fill(((mmg1_xchanmax-3.2)/0.4),mmg1AmpSum);
            hmmg12DMaxGainSingle->Fill(((mmg1_xchanmax-3.2)/0.4),mmg1_ampmax_x);
            hmmg12DGainLateTime->Fill(mmg1_x_timemax,mmg1_ampmax_x);
            if (mmg1_x_timemax>90) hmmg1GainMaxLate->Fill(mmg1_ampmax_x);
            hmmg12DGainSumLateTime->Fill(mmg1_x_timemax,mmg1AmpSum);
            if (mmg1_x_timemax>90) hmmg1GainSumLate->Fill(mmg1AmpSum);
            for (int i=0; i<240; i++) {
              if (mmg1GainAmps[i]>MMG1_THRESH) { hmmg12DGain->Fill(i,mmg1GainAmps[i]); }
            }
          }
        }
        
        float gemAmpSum=0.;
        int gemStripSum=0;
        
        if (tgem_ampmax_x>FE55_THRESH) {
          for (int i=0; i<10; i++) {
            if (gemGainAmps[(int)((tgem_xchanmax-3.2)/0.4) + i]>TGEM_THRESH) {
              gemAmpSum+=gemGainAmps[(int)((tgem_xchanmax-3.2)/0.4) + i];
              gemStripSum++;
            } else {break;}
          }
          for (int i=1; i<10; i++) {
            if (gemGainAmps[(int)((tgem_xchanmax-3.2)/0.4) - i]>TGEM_THRESH) {
              gemAmpSum+=gemGainAmps[(int)((tgem_xchanmax-3.2)/0.4) - i];
              gemStripSum++;
            } else {break;}
          }
          if (gemAmpSum>0. && gemStripSum>1 && gemStripSum<10) {
            hgemGainSum->Fill(gemAmpSum);
            hgemGainSpread->Fill(gemStripSum);
            hgemGainMultiplicity->Fill(gemStripSum,gemAmpSum);
            hgem2DMaxGain->Fill(((tgem_xchanmax-3.2)/0.4),gemAmpSum);
            hgem2DMaxGainSingle->Fill(((tgem_xchanmax-3.2)/0.4),tgem_ampmax_x);
            hgem2DGainLateTime->Fill(tgem_x_timemax,tgem_ampmax_x);
            if (tgem_x_timemax>90) hgemGainMaxLate->Fill(tgem_ampmax_x);
            hgem2DGainSumLateTime->Fill(tgem_x_timemax,gemAmpSum);
            if (tgem_x_timemax>90) hgemGainSumLate->Fill(gemAmpSum);
            for (int i=0; i<240; i++) {
              if (gemGainAmps[i]>TGEM_THRESH) { hgem2DGain->Fill(i,gemGainAmps[i]); }
            }
          }
        }
        
        float urwXAmpSum=0.;
        int urwXStripSum=0;
        
        if (urw_xampmax_x>FE55_THRESH) {
          for (int i=0; i<8; i++) {
            if (urwXGainAmps[(int)((urw_xchanmax-3.2)/0.8) + i]>URW_THRESH) {
              urwXAmpSum+=urwXGainAmps[(int)((urw_xchanmax-3.2)/0.8) + i];
              urwXStripSum++;
            } else {break;}
          }
          for (int i=1; i<8; i++) {
            if (urwXGainAmps[(int)((urw_xchanmax-3.2)/0.8) - i]>URW_THRESH) {
              urwXAmpSum+=urwXGainAmps[(int)((urw_xchanmax-3.2)/0.8) - i];
              urwXStripSum++;
            } else {break;}
          }
          if (urwXAmpSum>0. && urwXStripSum>1 && urwXStripSum<8) {
            hurwXGainSum->Fill(urwXAmpSum);
            hurwXGainSpread->Fill(urwXStripSum);
            hurwXGainMultiplicity->Fill(urwXStripSum,urwXAmpSum);
            hurwX2DMaxGain->Fill(((urw_xchanmax-3.2)/0.8),urwXAmpSum);
            hurwX2DMaxGainSingle->Fill(((urw_xchanmax-3.2)/0.8),urw_xampmax_x);
            hurwX2DGainLateTime->Fill(urw_x_timemax,urw_xampmax_x);
            if (urw_x_timemax>90) hurwXGainMaxLate->Fill(urw_xampmax_x);
            hurwX2DGainSumLateTime->Fill(urw_x_timemax,urwXAmpSum);
            if (urw_x_timemax>90) hurwXGainSumLate->Fill(urwXAmpSum);
            for (int i=0; i<120; i++) {
              if (urwXGainAmps[i]>URW_THRESH) { hurwX2DGain->Fill(i,urwXGainAmps[i]); }
            }
          }
        }
        
        float urwYAmpSum=0.;
        int urwYStripSum=0;

        if (urw_yampmax_x>FE55_THRESH) {
          for (int i=0; i<8; i++) {
            if (urwYGainAmps[(int)((urw_ychanmax-(3.2))/0.8) + i]>URW_THRESH) {
              urwYAmpSum+=urwYGainAmps[(int)((urw_ychanmax-(3.2))/0.8) + i];
              urwYStripSum++;
            } else {break;}
          }
          for (int i=1; i<8; i++) {
            if (urwYGainAmps[(int)((urw_ychanmax-(3.2))/0.8) - i]>URW_THRESH) {
              urwYAmpSum+=urwYGainAmps[(int)((urw_ychanmax-(3.2))/0.8) - i];
              urwYStripSum++;
            } else {break;}
          }
          if (urwYAmpSum>0. && urwYStripSum>1 && urwYStripSum<8) {
            hurwYGainSum->Fill(urwYAmpSum);
            hurwYGainSpread->Fill(urwYStripSum);
            hurwYGainMultiplicity->Fill(urwYStripSum,urwYAmpSum);
            hurwY2DMaxGain->Fill(((urw_ychanmax-(3.2))/0.8),urwYAmpSum);
            hurwY2DMaxGainSingle->Fill(((urw_ychanmax-(3.2))/0.8),urw_yampmax_x);
            hurwY2DGainLateTime->Fill(urw_y_timemax,urw_yampmax_x);
            if (urw_y_timemax>90) hurwYGainMaxLate->Fill(urw_yampmax_x);
            hurwY2DGainSumLateTime->Fill(urw_y_timemax,urwYAmpSum);
            if (urw_y_timemax>90) hurwYGainSumLate->Fill(urwYAmpSum);
            for (int i=0; i<120; i++) {
              if (urwYGainAmps[i]>URW_THRESH) { hurwY2DGain->Fill(i,urwYGainAmps[i]); }
            }
          }
        }
        
      #endif
      
      
      //========================================================
      //        Second fa125 pulse loop (Build External Track, Calc. Efficiencies)
      //========================================================
      //-- fADC-Based TRD External Tracking
      #if USE_TRD_EXT_TRACK
          float gem_extr = 0.;
          float mmg1_extr = 0.;
          float urw_extr = 0.;
          if (abs((mmg1_xchanmax-urw_xchanmax*um_slope)-um_offset)<4.75) trackFound_gem=true;
          if (abs((mmg1_xchanmax-tgem_xchanmax*gm_slope)-gm_offset)<3.25) trackFound_urw=true;
          if (abs((tgem_xchanmax-urw_xchanmax*ug_slope)-ug_offset)<4.) trackFound_mmg1=true; 
          
          //cout<<"======================"<<endl;
          float x1=mmg1_xchanmax, x2=tgem_xchanmax, x3=urw_xchanmax;
          if (trackFound_urw) {
            float a=(x2-x1)/(zgem-zmmg1);
            float b=((x1*zgem)-(x2*zmmg1))/(zgem-zmmg1);
            urw_extr = a*zurw+b;
            //cout<<"urw_extr ======="<<urw_extr<<endl;
            urw_f125_x_tracker_hits->Fill(urw_extr);
          }
          if (trackFound_gem) {
            float a_1=(x3-x1)/(zurw-zmmg1);
            float b_1=((x1*zurw)-(x3*zmmg1))/(zurw-zmmg1);
            gem_extr = a_1*zgem+b_1;
            //cout<<"gem_extr ======="<<gem_extr<<endl;
            f125_el_tracker_hits->Fill(gem_extr);
          }
          if (trackFound_mmg1) {
            float a_2=(x3-x2)/(zurw-zgem);
            float b_2=((x2*zurw)-(x3*zgem))/(zurw-zgem);
            mmg1_extr = a_2*zmmg1+b_2;
            //cout<<"mmg1_extr ======="<<mmg1_extr<<endl;
            mmg1_f125_el_tracker_hits->Fill(mmg1_extr);
          }
          
          for (ULong64_t i=0; i<f125_pulse_count; i++) { //--- Fadc125 Pulse Loop
          
        	float peak_amp = f125_pulse_peak_amp->at(i);
        	float ped = f125_pulse_pedestal->at(i);
        	if (0 > ped || ped > 200 ) ped = 100;
        	float amp=peak_amp-ped;
        	if (amp<0) amp=0;
        	float time=f125_pulse_peak_time->at(i);
        	int fADCSlot = f125_pulse_slot->at(i);
        	int fADCChan = f125_pulse_channel->at(i);
        	
          if (time>177.) continue;
          
        	int tripGemChan = Get3GEMChan(fADCChan, fADCSlot, RunNum);
        	int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
          int urwXChan = GetRWELLXChan(fADCChan, fADCSlot, RunNum);
          
          if (tripGemChan>-1) {
            amp = Get3GEMCalib(amp, tripGemChan, RunNum);
            if (amp>TGEM_THRESH) {
            float gemChan_x = tripGemChan*0.4+3.2; // to [mm]
            float gem_correction = 0.;
            if (RunNum<argonRunStop || RunNum>argonRunStart) {
              gem_correction = -0.30203 + (gemChan_x)*-0.005875;
            } else {
              gem_correction = -0.669259 + (gemChan_x)*-0.000868769;
            }
            //cout<<"gem_correction ======="<<gem_correction<<endl;
            tgem_residuals->Fill((gemChan_x-gem_extr));
            tgem_residual_ch->Fill(gemChan_x, (gemChan_x-gem_extr));
            gem_trk_hit=0;
            if (trackFound_gem && abs(gemChan_x-gem_extr-gem_correction)<3.85/* && 52.<=time && time<=118.*/) {
              tgem_residualscorr->Fill((gemChan_x-gem_extr)-gem_correction);
              tgem_residual_chcorr->Fill(gemChan_x, (gemChan_x-gem_extr-gem_correction));
              Count ("gem_trk_hit");
              htgem_trdTrackCorr->Fill(gemChan_x,(gem_extr-gem_correction));
              if (!match && tgem_tmp_nhit>=8) {
                f125_el_tracker_eff->Fill(gem_extr);
                match = true;
              }
              gem_trk_hit++;
            } //-- END within Xmm
            }
      	  }
      	  else if (mmg1Chan>-1) {
            amp = GetMMGCalib(amp, mmg1Chan, RunNum);
            if (amp>MMG1_THRESH) {
            float mmg1Chan_x = mmg1Chan*0.4+3.2; //-- to [mm]
            float mmg1_correction = 0.;
            if (RunNum<argonRunStop || RunNum>argonRunStart) {
              mmg1_correction = 0.085491 + (mmg1Chan_x)*0.0190135;
            } else {
              mmg1_correction = 0.351553 + (mmg1Chan_x)*0.000157831;
            }
            //cout<<"mmg1_correction ======="<<mmg1_correction<<endl;
            mmg1_residuals->Fill((mmg1Chan_x-mmg1_extr));
            mmg1_residual_ch->Fill(mmg1Chan_x, (mmg1Chan_x-mmg1_extr));
            mmg1_trk_hit=0;
            if (trackFound_mmg1 && abs(mmg1Chan_x-mmg1_extr-mmg1_correction)<3. /*&& 50.<=time && time<=135*/) {
              mmg1_residualscorr->Fill((mmg1Chan_x-mmg1_extr)-mmg1_correction);
              mmg1_residual_chcorr->Fill(mmg1Chan_x, (mmg1Chan_x-mmg1_extr)-mmg1_correction);
              Count ("mmg1_trk_hit");
              hmmg1_trdTrackCorr->Fill(mmg1Chan_x,(mmg1_extr-mmg1_correction));
              if (!match_mmg1 && mmg1_tmp_nhit>=8) {
                mmg1_f125_el_tracker_eff->Fill(mmg1_extr);
                match_mmg1 = true;
              }
              mmg1_trk_hit++;
            } //-- END within Xmm
            }
      	  }
          else if (amp>URW_THRESH && urwXChan>-1) {
            float urwChan_x = urwXChan*0.8+3.2; //-- to [mm]
            float urw_correction = 0.;
            if (RunNum<argonRunStop || RunNum>argonRunStart) {
              urw_correction = 1.1586 + (urwChan_x)*0.0031787;
            } else {
              urw_correction = 0.695631 + (urwChan_x)*0.00658071;
            }
            //cout<<"urw_correction ======="<<urw_correction<<endl;
            urw_x_residuals->Fill((urwChan_x-urw_extr));
            urw_x_residual_ch->Fill(urwChan_x, (urwChan_x-urw_extr));
            urw_trk_hit=0;
            if (trackFound_urw && abs(urwChan_x-urw_extr-urw_correction)<3.75 /*&& 48.<=time && time<=122*/) {
              urw_x_residualscorr->Fill((urwChan_x-urw_extr)-urw_correction);
              urw_x_residual_chcorr->Fill(urwChan_x, (urwChan_x-urw_extr-urw_correction));
              Count ("urw_trk_hit");
              hurw_trdTrackCorr->Fill(urwChan_x,(urw_extr-urw_correction));
              if (!match_urw && urw_tmp_nhit>=4) {
                urw_f125_x_tracker_eff->Fill(urw_extr);
                match_urw = true;
              }
              urw_trk_hit++;
            } //-- END within Xmm
      	  }
      	} //--- end Fa125 Pulse Loop ---
          
        
      #else //-- SRS-Based GEMTracker External Tracking
      
      //if (gt_2_idx_x>0 && gt_1_idx_x>0) { //--External tracking condition
      if (GoodXCorr && GoodYCorr) { //--External tracking condition
        
        Count("trk_hit");
        float x1=gemtrkr1_xch_max, x2=gemtrkr2_xch_max;
        float a=(x2-x1)/(z2-z1);
        float b=((x1)*z2-(x2)*z1)/(z2-z1);
        float gem_extr = a*zgem+b;
        float mmg1_extr = a*zmmg1+b;
        float urw_extr = a*zurw+b;
        f125_el_tracker_hits->Fill(gem_extr);
        mmg1_f125_el_tracker_hits->Fill(mmg1_extr);
        urw_f125_x_tracker_hits->Fill(urw_extr);
        
        for (ULong64_t i=0; i<f125_pulse_count; i++) { //--- Fadc125 Pulse Loop
          
        	float peak_amp = f125_pulse_peak_amp->at(i);
        	float ped = f125_pulse_pedestal->at(i);
        	if (0 > ped || ped > 200 ) ped = 100;
        	float amp=peak_amp-ped;
        	if (amp<0) amp=0;
        	float time=f125_pulse_peak_time->at(i);
        	int fADCSlot = f125_pulse_slot->at(i);
        	int fADCChan = f125_pulse_channel->at(i);
        	
          if (time>177.) continue;
          
        	int tripGemChan = Get3GEMChan(fADCChan, fADCSlot, RunNum);
        	int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
          int urwXChan = GetRWELLXChan(fADCChan, fADCSlot, RunNum);
        	
        	if (tripGemChan>-1) {
            amp = Get3GEMCalib(amp, tripGemChan, RunNum);
            if (amp>TGEM_THRESH) {
            float gemChan_x = tripGemChan*0.4+3.2; // to [mm]
            float gem_correction = 0.;
            if (RunNum<argonRunStop || RunNum>argonRunStart) {
              gem_correction = -0.30203 + (gemChan_x)*-0.005875;
            } else {
              gem_correction = -0.669259 + (gemChan_x)*-0.000868769;
            }
            tgem_residuals->Fill((gemChan_x-gem_extr));
            tgem_residual_ch->Fill(gemChan_x, (gemChan_x-gem_extr));
            gem_trk_hit=0;
            if (abs(gemChan_x-gem_extr-gem_correction)<3.85) { //within 7.7 mm
              tgem_residualscorr->Fill((gemChan_x-gem_extr)-gem_correction);
              tgem_residual_chcorr->Fill(gemChan_x, (gemChan_x-gem_extr-gem_correction));
              Count ("gem_trk_hit");
              htgem_trdTrackCorr->Fill(gemChan_x,(gem_extr-gem_correction));
              if (!match && tgem_tmp_nhit>=8) {
                f125_el_tracker_eff->Fill(gem_extr);
                match = true;
              }
              gem_trk_hit++;
            } //-- END within 7.7mm
            }
      	  }
      	  else if (mmg1Chan>-1) {
            amp = GetMMGCalib(amp, mmg1Chan, RunNum);
            if (amp>MMG1_THRESH) {
            float mmg1Chan_x = mmg1Chan*0.4+3.2; //-- to [mm]
            float mmg1_correction = 0.;
            if (RunNum<argonRunStop || RunNum>argonRunStart) {
              mmg1_correction = 0.085491 + (mmg1Chan_x)*0.0190135;
            } else {
              mmg1_correction = 0.351553 + (mmg1Chan_x)*0.000157831;
            }
            mmg1_residuals->Fill((mmg1Chan_x-mmg1_extr));
            mmg1_residual_ch->Fill(mmg1Chan_x, (mmg1Chan_x-mmg1_extr));
            mmg1_trk_hit=0;
            if (abs(mmg1Chan_x-mmg1_extr-mmg1_correction)<3.) { //within 6 mm
              mmg1_residualscorr->Fill((mmg1Chan_x-mmg1_extr)-mmg1_correction);
              mmg1_residual_chcorr->Fill(mmg1Chan_x, (mmg1Chan_x-mmg1_extr)-mmg1_correction);
              Count ("mmg1_trk_hit");
              hmmg1_trdTrackCorr->Fill(mmg1Chan_x,(mmg1_extr-mmg1_correction));
              if (!match_mmg1 && mmg1_tmp_nhit>=8) {
                mmg1_f125_el_tracker_eff->Fill(mmg1_extr);
                match_mmg1 = true;
              }
              mmg1_trk_hit++;
            } //-- END within 6mm
            }
      	  }
          else if (amp>URW_THRESH && urwXChan>-1) {
            float urwChan_x = urwXChan*0.8+3.2; //-- to [mm]
            float urw_correction = 0.;
            if (RunNum<argonRunStop || RunNum>argonRunStart) {
              urw_correction = -1.1586 + (urwChan_x)*0.0031787;
            } else {
              urw_correction = -0.695631 + (urwChan_x)*0.00658071;
            }
            urw_x_residuals->Fill((urwChan_x-urw_extr));
            urw_x_residual_ch->Fill(urwChan_x, (urwChan_x-urw_extr));
            urw_trk_hit=0;
            if (abs(urwChan_x-urw_extr-urw_correction)<3.75) { //within 7.5 mm
              urw_x_residualscorr->Fill((urwChan_x-urw_extr)-urw_correction);
              urw_x_residual_chcorr->Fill(urwChan_x, (urwChan_x-urw_extr-urw_correction));
              Count ("urw_trk_hit");
              hurw_trdTrackCorr->Fill(urwChan_x,(urw_extr-urw_correction));
              if (!match_urw && urw_tmp_nhit>=4) {
                urw_f125_x_tracker_eff->Fill(urw_extr);
                match_urw = true;
              }
              urw_trk_hit++;
            } //-- END within 7.5mm
      	  }
          
      	} //--- end second Fa125 Pulse Loop ---
      } //-- End external tracker condition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #endif
      
      
      //==============================================
      //        Third fa125 pulse loop (Check all selections before saving hit info to TTrees)
      //==============================================
      
      int tgem_npulse[240], mmg1_npulse[240], urw_npulse[120];
      int tgem_npulse_sum=0, mmg1_npulse_sum=0, urw_npulse_sum=0;
      int tgem_npulse_weighted=0., mmg1_npulse_weighted=0., urw_npulse_weighted=0.;
      int tgem_npulse_avgChan=-1, mmg1_npulse_avgChan=-1, urw_npulse_avgChan=-1;
      for (ULong64_t i=0; i<240; i++) {
          tgem_npulse[i] = 0;
          mmg1_npulse[i] = 0;
       }
       for (ULong64_t i=0; i<120; i++) {
          urw_npulse[i] = 0;
       }
      
      for (ULong64_t i=0; i<f125_pulse_count; i++) {
        
      	float peak_amp = f125_pulse_peak_amp->at(i);
      	float ped = f125_pulse_pedestal->at(i);
      	if (0 > ped || ped > 200 ) ped = 100;
      	float amp=peak_amp-ped;
      	float time=f125_pulse_peak_time->at(i);
      	int fADCSlot = f125_pulse_slot->at(i);
      	int fADCChan = f125_pulse_channel->at(i);
      	#ifdef GAIN_CALIB
        if (time>177.) continue;
        #else
        if ((!match && !match_mmg1 && !match_urw) || time>177.) continue;
        #endif
        
      	int tripGemChan = Get3GEMChan(fADCChan, fADCSlot, RunNum);
      	double tripGemChan_x = tripGemChan*0.4 + 3.2;
      	int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
      	double mmg1Chan_x = mmg1Chan*0.4 + 3.2;
        int urwXChan = GetRWELLXChan(fADCChan, fADCSlot, RunNum);
       	double urwChan_x = urwXChan*0.8+3.2; // to [mm]
        int urwYChan = GetRWELLYChan(fADCChan, fADCSlot, RunNum);
       	double urwChan_y = urwYChan*0.8+3.2; // to [mm]
        
      	if (amp<0) amp=0;
      	
      	if (tripGemChan>-1) {
          amp = Get3GEMCalib(amp, tripGemChan, RunNum);
          if (amp>TGEM_THRESH) {
            f125_fit->Fill(time,tripGemChan,amp);
            if (/*52.<=time && time<=120. &&*/ (abs((mmg1_xchanmax-(tripGemChan_x*gm_slope))-gm_offset))<3.25 && abs((tripGemChan_x-urw_xchanmax*ug_slope)-ug_offset)<4.) {
              tgem_pos_x[tgem_idx_x] = tripGemChan_x;
              tgem_amp_x[tgem_idx_x] = amp;
              tgem_time_x[tgem_idx_x] = time;
              tgem_idx_x++;
            }
            #ifdef GAIN_CALIB
        
            #else
            if (match) { //tgem hit matched with external track
            #endif
            if (tgem_ampmax<amp) {
              tgem_ampmax=amp;
              tgem_timemax=time;
              tgem_channel_max=tripGemChan;
            }
            if (mmg1_xchanmax>0.) hgemPulseDiff_el->Fill(mmg1_xchanmax-tripGemChan_x);//*0.908525)-9.22087));
      	    f125_el->Fill(amp);
            f125_el_amp2d->Fill(time,tripGemChan,amp);
            f125_xVSamp->Fill(tripGemChan,amp);
            f125_timeVSamp->Fill(time,amp);
      	    tgem_xpos.push_back(tripGemChan);
      	    tgem_dedx.push_back(amp);
      	    tgem_zpos.push_back(time);
      	    tgem_nhit++;
            tgem_npulse[tripGemChan]++;
            tgem_npulse_sum++;
      	    tgem_zHist->Fill(time,amp);
            bool foundRun = false;
            for (int j=0; j<listSize; j++) {
              if (RunNum==noRadList[j]) {foundRun=true; break;}
            }
            if (foundRun==true) {tgem_parID.push_back(0);}
            else {tgem_parID.push_back(1);}
            #ifdef GAIN_CALIB
        
            #else
            }
            #endif
          }
    	  }
    	  else if (mmg1Chan>-1) {
          amp = GetMMGCalib(amp, mmg1Chan, RunNum);
          if (amp>MMG1_THRESH) {
            mmg1_f125_fit->Fill(time,mmg1Chan,amp);
            if (/*52.<=time && time<=135. && */abs((mmg1Chan_x-tgem_xchanmax*gm_slope)-gm_offset)<3.25 && abs((mmg1Chan_x-urw_xchanmax*um_slope)-um_offset)<4.75) {
              mmg1_pos_x[mmg1_idx_x] = mmg1Chan_x;
              mmg1_amp_x[mmg1_idx_x] = amp;
              mmg1_time_x[mmg1_idx_x] = time;
              mmg1_idx_x++;
            }
            #ifdef GAIN_CALIB
            
            #else
            if (match_mmg1) { //mmg1 hit matched with external track
            #endif
            if (mmg1_ampmax<amp) {
              mmg1_ampmax=amp;
              mmg1_timemax=time;
              mmg1_channel_max=mmg1Chan;
            }
            if (tgem_xchanmax>0.) hmmg1PulseDiff_el->Fill(mmg1Chan_x-tgem_xchanmax);//*0.908525-mmg1Chan_x-9.22087);
      	    mmg1_f125_el_amp2d->Fill(time,mmg1Chan,amp);
            mmg1_f125_xVSamp->Fill(mmg1Chan,amp);
            mmg1_f125_timeVSamp->Fill(time,amp);
      	    mmg1_f125_el->Fill(amp);
      	    mmg1_xpos.push_back(mmg1Chan);
      	    mmg1_dedx.push_back(amp);
      	    mmg1_zpos.push_back(time);
      	    mmg1_nhit++;
            mmg1_npulse[mmg1Chan]++;
            mmg1_npulse_sum++;
      	    mmg1_zHist->Fill(time,amp);
            bool foundRun = false;
            for (int j=0; j<listSize; j++) {
              if (RunNum==noRadList[j]) {foundRun=true; break;}
            }
            if (foundRun==true) {mmg1_parID.push_back(0);}
            else {mmg1_parID.push_back(1);}
            #ifdef GAIN_CALIB
        
            #else
            }
            #endif
          }
        }
        else if (urwXChan>-1) {
          if (amp>URW_THRESH) {
            urw_f125_fit->Fill(time,urwXChan,amp);
            if (/*48.<=time && time<=125. && */abs((tgem_xchanmax-urwChan_x*ug_slope)-ug_offset)<4. && abs((mmg1_xchanmax-urwChan_x*um_slope)-um_offset)<4.75) {
              urw_pos_x[urw_idx_x] = urwChan_x;
              urw_amp_x[urw_idx_x] = amp; 
              urw_time_x[urw_idx_x] = time;
              urw_idx_x++;
            }
            #ifdef GAIN_CALIB

            #else
            if (match_urw) {
            #endif
            if (urw_xampmax<amp) {
              urw_xampmax=amp;
              urw_xtimemax=time;
              urw_channelX_max=urwXChan;
            }
            if (tgem_xchanmax>0.) hurwPulseDiff_el->Fill(tgem_xchanmax-urwChan_x);//*0.881493-7.886-tgem_xchanmax);
            if (mmg1_xchanmax>0.) hurwPulseDiff_mmg->Fill(mmg1_xchanmax-urwChan_x);//*0.829002-mmg1_xchanmax-18.4591);
            urw_f125_x_amp2d->Fill(time,urwXChan,amp);
            urw_f125_xVSamp->Fill(urwXChan,amp);
            urw_f125_x_timeVSamp->Fill(time,amp);
            urw_f125_el_x->Fill(amp);
            urw_xpos.push_back(urwXChan);
            urw_dedx.push_back(amp);
            urw_zpos.push_back(time);
            urw_nhit++;
            urw_npulse[urwXChan]++;
            urw_npulse_sum++;
            urw_zHist->Fill(time,amp);
            bool foundRun = false;
            for (int j=0; j<listSize; j++) {
              if (RunNum==noRadList[j]) {foundRun=true; break;}
            }
            if (foundRun==true) {urw_parID.push_back(0);}
              else {urw_parID.push_back(1);}
            #ifdef GAIN_CALIB

            #else
            }
            #endif
          }
        }
        else if (amp>URW_THRESH && urwYChan>-1) {
            //if ( /* 48.<=time && time<=125. && */ (abs((tgem_xchanmax-urwChan_x*ug_slope)-ug_offset)<5. || abs((mmg1_xchanmax-urwChan_x*um_slope)-um_offset)<5.)) { //2.5
            if (match_urw) {
              urw_pos_y[urw_idx_y] = urwChan_y;
              urw_amp_y[urw_idx_y] = amp;
              urw_idx_y++;
            }
            #ifdef GAIN_CALIB

            #else
            if (match_urw) {
            #endif
            if (urw_yampmax<amp) {
              urw_yampmax=amp;
              urw_ytimemax=time;
              urw_channelY_max=urwYChan;
            }
            urw_f125_y_amp2d->Fill(time,urwYChan,amp);
            urw_f125_yVSamp->Fill(urwYChan,amp);
            urw_f125_el_y->Fill(amp);
            urw_f125_y_timeVSamp->Fill(time,amp);
            
            //urw_ypos.push_back(urwYChan);
            //urw_dedx.push_back(amp);
            //urw_zpos.push_back(time);
            urw_nyhit++;
            //urw_zHist->Fill(time,amp);
            #ifdef GAIN_CALIB

            #else
            }
            #endif
        }
    	} //--- end Fa125 Pulse Loop ---
      
      for (int j=0; j<240; j++) {
        tgem_npulse_weighted+=j*tgem_npulse[j];
        if (tgem_npulse[j]>0) htgem_2DPulseMultiplicity->Fill(j,tgem_npulse[j]);
        mmg1_npulse_weighted+=j*mmg1_npulse[j];
        if (mmg1_npulse[j]>0) hmmg1_2DPulseMultiplicity->Fill(j,mmg1_npulse[j]);
      }
      for (int j=0; j<120; j++) {
        urw_npulse_weighted+=j*urw_npulse[j];
        if (urw_npulse[j]>0) hurw_2DPulseMultiplicity->Fill(j,urw_npulse[j]);
      }
      
      if (tgem_npulse_sum>0) {
        tgem_npulse_avgChan=tgem_npulse_weighted/tgem_npulse_sum;
        htgem_2DPulseVsChan->Fill(tgem_npulse_avgChan, tgem_npulse_sum);
      }
      if (mmg1_npulse_sum>0) {
        mmg1_npulse_avgChan=mmg1_npulse_weighted/mmg1_npulse_sum;
        hmmg1_2DPulseVsChan->Fill(mmg1_npulse_avgChan, mmg1_npulse_sum);
      }
      if (urw_npulse_sum>0) {
        urw_npulse_avgChan=urw_npulse_weighted/urw_npulse_sum;
        hurw_2DPulseVsChan->Fill(urw_npulse_avgChan, urw_npulse_sum);
      }
      //===================================================================
      //                    Chi^2 Fit Calculation
      //===================================================================
      
      char f125Title[80]; sprintf(f125Title,"GEM-TRD: Event=%lld Run=%d; z pos [time *8ns]; y pos [ch #]",jentry,RunNum);
      f125_fit->SetTitle(f125Title);
      if (f125_fit->GetEntries()!=0) {
        std::pair<Double_t, Double_t> fitResult  = TrkFit(f125_fit,fx1,"fx1",1);
        chi2cc_tgem = fitResult.first;
        integral_tgem = fitResult.second;
        a0 = fx1.GetParameter(0);
        a1 = fx1.GetParameter(1);
      }
      
      char mmg1f125Title[80]; sprintf(mmg1f125Title,"MMG1-TRD: Event=%lld Run=%d; z pos [time *8ns]; y pos [ch #]",jentry,RunNum);
      mmg1_f125_fit->SetTitle(mmg1f125Title);
      if (mmg1_f125_fit->GetEntries()!=0) {
        std::pair<Double_t, Double_t> fitResult  = TrkFit(mmg1_f125_fit,fx2,"fx2",1);
        chi2cc_mmg1 = fitResult.first;
        integral_mmg1 = fitResult.second;
        a0_mmg1 = fx2.GetParameter(0);
        a1_mmg1 = fx2.GetParameter(1);
      }
      
      char urwf125Title[80]; sprintf(urwf125Title,"uRWELL-TRD:  Event=%lld Run=%d; z pos [time *8ns]; y pos [ch #]",jentry,RunNum);
      urw_f125_fit->SetTitle(urwf125Title);
      if (urw_f125_fit->GetEntries()!=0) {
        std::pair<Double_t, Double_t> fitResult  = TrkFit(urw_f125_fit,fx2,"fx2",1);
        chi2cc_urw = fitResult.first;
        integral_urw = fitResult.second;
        a0_urw = fx2.GetParameter(0);
        a1_urw = fx2.GetParameter(1);
      }
      
      //==================== Max Amplitude histos ============================
        if (tgem_ampmax>TGEM_THRESH) {
          f125_el_max->Fill(tgem_ampmax);
          f125_timeVSamp_max->Fill(tgem_timemax,tgem_ampmax);
          if (tgem_timemax>112) f125_el_max_late->Fill(tgem_ampmax);
          f125_el_amp2d_max->Fill(tgem_timemax,(tgem_channel_max*0.4+3.2),tgem_ampmax);
          f125_xVSamp_max->Fill(tgem_channel_max,tgem_ampmax);
          if (mmg1_ampmax>MMG1_THRESH) tgem_mmg1_max_xcorr->Fill((tgem_channel_max*0.4+3.2),(mmg1_channel_max*0.4+3.2));
          if (urw_xampmax>URW_THRESH) tgem_urw_max_xcorr->Fill((urw_channelX_max*0.8+3.2),(tgem_channel_max*0.4+3.2));
          if (tgem_yamp_max>TRKR_THRESH+600.) htgem_max_xy->Fill((tgem_channel_max*0.4+3.2),tgem_ych_max);
        }
        if (mmg1_ampmax>MMG1_THRESH) {
          mmg1_f125_el_max->Fill(mmg1_ampmax);
          mmg1_f125_timeVSamp_max->Fill(mmg1_timemax,mmg1_ampmax);
          if (mmg1_timemax>125) mmg1_f125_el_max_late->Fill(mmg1_ampmax);
          mmg1_f125_el_amp2d_max->Fill(mmg1_timemax,(mmg1_channel_max*0.4+3.2),mmg1_ampmax);
          mmg1_f125_xVSamp_max->Fill(mmg1_channel_max,mmg1_ampmax);
          if (urw_xampmax>URW_THRESH) urw_mmg1_max_xcorr->Fill((urw_channelX_max*0.8+3.2),(mmg1_channel_max*0.4+3.2));
          if (mmg1_yamp_max>TRKR_THRESH+600.) hmmg1_max_xy->Fill((mmg1_channel_max*0.4+3.2),mmg1_ych_max);
        }
        if (urw_xampmax>URW_THRESH) {
          urw_f125_el_xmax->Fill(urw_xampmax);
          urw_f125_x_timeVSamp_max->Fill(urw_xtimemax,urw_xampmax);
          if (urw_xtimemax>116) urw_f125_el_xmax_late->Fill(urw_xampmax);
          urw_f125_x_amp2d_max->Fill(urw_xtimemax,(urw_channelX_max*0.8+3.2),urw_xampmax);
          urw_f125_xVSamp_max->Fill(urw_channelX_max,urw_xampmax);
          if (urw_yampmax>URW_THRESH) hurw_max_xy->Fill((urw_channelX_max*0.8+3.2),urw_ychanmax);
        }
        if (urw_yampmax>URW_THRESH) {
          urw_f125_el_ymax->Fill(urw_yampmax);
          urw_f125_y_amp2d_max->Fill(urw_ytimemax,(urw_channelY_max*0.8+3.2),urw_yampmax);
          urw_f125_y_timeVSamp_max->Fill(urw_ytimemax,urw_yampmax);
          urw_f125_yVSamp_max->Fill(urw_channelY_max,urw_yampmax);
        }
      
      tgem_amp_max=tgem_ampmax;
      tgem_time_max=tgem_timemax;
      tgem_xch_max=tgem_channel_max;
      tgem_chi2cc.push_back(chi2cc_tgem);
      tgem_integral.push_back(integral_tgem);
      mmg1_amp_max=mmg1_ampmax;
      mmg1_time_max=mmg1_timemax;
      mmg1_xch_max=mmg1_channel_max;
      mmg1_chi2cc.push_back(chi2cc_mmg1);
      mmg1_integral.push_back(integral_mmg1);
      urw_amp_max=urw_xampmax;
      urw_time_max=urw_xtimemax;
      urw_xch_max=urw_xchanmax;
      urw_chi2cc.push_back(chi2cc_urw);
      urw_integral.push_back(integral_urw);
      
      //if (abs(tgem_el_chan_max-(mmg1_xchanmax-1.35))<100. || abs(mmg1_el_chan_max-(tgem_xchanmax-1.35))<100.) gem_mmg1_max_xcorr->Fill(tgem_xchanmax, mmg1_xchanmax);
      
      for (int i=1; i<21; i++) {
        tgem_zHist_vect.push_back(tgem_zHist->GetBinContent(i));
        mmg1_zHist_vect.push_back(mmg1_zHist->GetBinContent(i));
        urw_zHist_vect.push_back(urw_zHist->GetBinContent(i));
      }
      
      if (tgem_nhit>0.) htgem_nhits->Fill(tgem_nhit);
      if (mmg1_nhit>0.) hmmg1_nhits->Fill(mmg1_nhit);
      if (urw_nhit>0.) hurw_nxhits->Fill(urw_nhit);
      if (urw_nyhit>0.) hurw_nyhits->Fill(urw_nyhit);
      if (gt1_nhit>0.) hgt1_nhits->Fill(gt1_nhit);
      if (gt2_nhit>0.) hgt2_nhits->Fill(gt2_nhit);
      
        if (gemtrkr1_xamp_max>TRKR_THRESH && tgem_ampmax>TGEM_THRESH) tgem_gt1_xcorr->Fill(tgem_xchanmax, gemtrkr1_xch_max);
        if (gemtrkr1_xamp_max>TRKR_THRESH && urw_xampmax>URW_THRESH) urw_gt1_xcorr->Fill(urw_xchanmax, gemtrkr1_xch_max);
        if (gemtrkr1_xamp_max>TRKR_THRESH && mmg1_ampmax>MMG1_THRESH) mmg1_gt1_xcorr->Fill(mmg1_xchanmax, gemtrkr1_xch_max);
        if (gemtrkr2_xamp_max>TRKR_THRESH && tgem_ampmax>TGEM_THRESH) tgem_gt2_xcorr->Fill(tgem_xchanmax, gemtrkr2_xch_max);
        if (gemtrkr2_xamp_max>TRKR_THRESH && urw_xampmax>URW_THRESH) urw_gt2_xcorr->Fill(urw_xchanmax, gemtrkr2_xch_max);
        if (gemtrkr2_xamp_max>TRKR_THRESH && mmg1_ampmax>MMG1_THRESH) mmg1_gt2_xcorr->Fill(mmg1_xchanmax, gemtrkr2_xch_max);
        if (gemtrkr1_xamp_max>TRKR_THRESH && gemtrkr2_xamp_max>TRKR_THRESH) hTrackDiff->Fill(gemtrkr1_xch_max - gemtrkr2_xch_max);
      
      for (ULong64_t j=0; j<urw_idx_x; j++) {
        if (urw_xtimemax>0 && urw_xtimemax!=urw_time_x[j]) hurw_timeDiff->Fill((urw_pos_x[j]-3.2)/0.8,(urw_xtimemax - urw_time_x[j]));
          for (ULong64_t i=0; i<urw_idx_y; i++) {
            if (urw_amp_y[i]>URW_THRESH && urw_amp_x[j]>URW_THRESH) hurw_xy->Fill(urw_pos_x[j], urw_pos_y[i]);
          }
          for (ULong64_t i=0; i<tgem_idx_x; i++) {
            if (tgem_amp_x[i]>TGEM_THRESH && urw_amp_x[j]>URW_THRESH) urw_tgem_xcorr->Fill(urw_pos_x[j], tgem_pos_x[i]);
          }
          for (ULong64_t i=0; i<mmg1_idx_x; i++) {
            if (mmg1_amp_x[i]>MMG1_THRESH && urw_amp_x[j]>URW_THRESH) urw_mmg1_xcorr->Fill(urw_pos_x[j], mmg1_pos_x[i]);
          }
      }
      for (ULong64_t j=0; j<tgem_idx_x; j++) {
          if (tgem_timemax>0 && tgem_timemax!=tgem_time_x[j]) htgem_timeDiff->Fill((tgem_pos_x[j]-3.2)/0.4,(tgem_timemax - tgem_time_x[j]));
          for (ULong64_t i=0; i<tgem_idx_y; i++) {
            if (tgem_peak_y_height[i]>TRKR_THRESH+600. && tgem_amp_x[j]>TGEM_THRESH) htgem_xy->Fill(tgem_pos_x[j], tgem_peak_pos_y[i]);
          }
          for (ULong64_t i=0; i<mmg1_idx_x; i++) {
            if (mmg1_amp_x[i]>MMG1_THRESH && tgem_amp_x[j]>TGEM_THRESH) tgem_mmg1_xcorr->Fill(tgem_pos_x[j], mmg1_pos_x[i]);
          }
      }
      for (ULong64_t i=0; i<tgem_idx_y; i++) {
          for (ULong64_t j=0; j<mmg1_idx_y; j++) {
            if (mmg1_peak_y_height[j]>TRKR_THRESH+600. && tgem_peak_y_height[i]>TRKR_THRESH+600.) tgem_mmg1_ycorr->Fill(tgem_peak_pos_y[i], mmg1_peak_pos_y[j]);
          }
      }
      for (ULong64_t i=0; i<mmg1_idx_x; i++) {
        if (mmg1_timemax>0 && mmg1_timemax!=mmg1_time_x[i]) hmmg1_timeDiff->Fill((mmg1_pos_x[i]-3.2)/0.4,(mmg1_timemax - mmg1_time_x[i]));
        for (ULong64_t j=0; j<mmg1_idx_y; j++) {
          if (mmg1_amp_x[i]>MMG1_THRESH && mmg1_peak_y_height[j]>TRKR_THRESH+600.) hmmg1_xy->Fill(mmg1_pos_x[i], mmg1_peak_pos_y[j]);
          }
      }
      
      
    #endif //if 1 (compilation condition)
    //======================= End Process Fa125 Pulse data ================================
    
    //==================================================================================================
    //                    Process Fa125  RAW data
    //==================================================================================================
    
    #ifdef USE_125_RAW
      #ifdef VERBOSE
        if (jentry<MAX_PRINT) printf("------------------ Fadc125  wraw_count = %llu ---------\n", f125_wraw_count);
      #endif
      mhevt->Reset();
      mhevtc->Reset();
      mhevtf->Reset();
      uhevt->Reset();
      uhevtc->Reset();
      uhevtf->Reset();
      hevt->Reset();
      hevtc->Reset();
      hevtf->Reset();
      
      for (ULong64_t i=0; i<f125_wraw_count; i++) { // --- fadc125 channels loop
        
        if (!match && !match_mmg1 && !match_urw) continue;
        //cout<<"******** EVT# "<<jentry<<" RAW LOOP COND.: match=="<<match<<",  match_mmg1=="<<match_mmg1<<" **********"<<endl;
        int fadc_window = f125_wraw_samples_count->at(i);
        int fADCSlot = f125_wraw_slot->at(i);
        int fADCChan = f125_wraw_channel->at(i);
        int tripGemChan = Get3GEMChan(fADCChan, fADCSlot, RunNum);
        int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
        int urwXChan = GetRWELLXChan(fADCChan, fADCSlot, RunNum);
        double DEDX_THR = TGEM_THRESH, mDEDX_THR = MMG1_THRESH, uDEDX_THR = URW_THRESH;
        int TimeWindowStart = 60;
        int TimeWindowStart_m = 50;
        int TimeWindowStart_u = 45;
        int TimeMin = 0;
        int TimeMax = 140;
        
        //--ped calculation for f125 raw data
        int nped = 0, ped = 100;
        double ped_sum = 0.;
        for (int si=TimeWindowStart-15; si<TimeWindowStart; si++) {
          int ped_samp = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si);
          ped_sum += ped_samp;
          nped++;
        }
        ped = ped_sum / nped;
        if (0. > ped || ped > 200 ) ped = 100;
        //--FOR MMG
        int nped_m = 0, ped_m = 100;
        double ped_m_sum = 0.;
        for (int si=TimeWindowStart_m-15; si<TimeWindowStart_m; si++) {
          int ped_m_samp = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si);
          ped_m_sum += ped_m_samp;
          nped_m++;
        }
        ped_m = ped_m_sum / nped_m;
        if (0. > ped_m || ped_m > 200 ) ped_m = 100;
        //--FOR URW
        int nped_u = 0, ped_u = 100;
        double ped_u_sum = 0.;
        for (int si=TimeWindowStart_u-15; si<TimeWindowStart_u; si++) {
          int ped_u_samp = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si);
          ped_u_sum += ped_u_samp;
          nped_u++;
        }
        ped_u = ped_u_sum / nped_u;
        if (0. > ped_u || ped_u > 200 ) ped_u = 100;
        
        for (int si=0; si<fadc_window; si++) {
          int time=si;
          int adc = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si);
          if (tripGemChan>-1) {
            adc = adc - ped;
            if (adc<0) adc=0;
            adc = Get3GEMCalib(adc, tripGemChan, RunNum);
            //if (adc>4090) printf("!!!!!!!!!!!!!!!!!!!!!! ADC 125 overflow: %d \n",adc);
            if (adc>DEDX_THR) {
              time-=TimeWindowStart;
              ///////////////if ( TimeMin > time || time > TimeMax ) continue; // --- drop early and late hits ---
              //if ((100-time)<1 || (100-time)>hevt->GetNbinsX() || (tripGemChan+1)<1 || (tripGemChan+1)>hevt->GetNbinsY()) {
                //cout<<"Warning:: hevtc Histo bin out of range: time="<<time<<", chan="<<tripGemChan<<endl;
              //}
              hevtc->SetBinContent(100-time,tripGemChan+1,adc/100.);
              hevt->SetBinContent(100-time,tripGemChan+1,adc/100.);
            }
          }
          if (mmg1Chan>-1) {
            adc = adc - ped_m;
            if (adc<0) adc=0;
            adc = GetMMGCalib(adc, mmg1Chan, RunNum);
            if (adc > mDEDX_THR) {
              time-=(TimeWindowStart_m); //+40
              mhevtc->SetBinContent(140-time,mmg1Chan+1,adc/100.); //120-
              mhevt->SetBinContent(140-time,mmg1Chan+1,adc/100.); //120-
            }
          }
          if (urwXChan>-1) {
            adc = adc - ped_u;
            if (adc<0) adc=0;
            if (adc > uDEDX_THR) {
              time-=(TimeWindowStart_u); //+40
              uhevtc->SetBinContent(100-time,urwXChan+1,adc/100.);
              uhevt->SetBinContent(100-time,urwXChan+1,adc/100.);
            }
          }
        } // --  end of samples loop
      } // -- end of fadc125 raw channels loop
      
      #ifdef SHOW_EVT_DISPLAY
          #if (USE_PULSE>0)
            c2->cd(1); hevtk->Draw("box");
          #else
            c2->cd(1); hevt->Draw("colz");
            c2->cd(6); mhevt->Draw("colz");
            c2->cd(11); uhevt->Draw("colz");
          #endif
          c2->cd(2);   hevtf->Draw("text");
          c2->cd(7);   mhevtf->Draw("text");
          c2->cd(12);  uhevtf->Draw("text");
          c2->Modified(); c2->Update();
      #endif
      
      //==================================================================================================
      //            Begin NN Clustering & Track Fitting
      //==================================================================================================
      #if (USE_CLUST>0)
        // -------------------------------   hist dist clustering         ------------------------
        //--GEM-TRD
        float clust_Xmax[MAX_CLUST];
        float clust_Zmax[MAX_CLUST];
        float clust_Emax[MAX_CLUST];
        float clust_Xpos[MAX_CLUST];
        float clust_Zpos[MAX_CLUST];
        float clust_dEdx[MAX_CLUST];
        float clust_Size[MAX_CLUST];
        float clust_Width[MAX_CLUST][3];  // y1, y2, dy ; strips
        float clust_Length[MAX_CLUST][3]; // x1, x2, dx ; time
        float hits_Xpos[500];
        float hits_Zpos[500];
        float hits_dEdx[500];
        float hits_Size[MAX_CLUST];
        float hits_Width[MAX_CLUST];  // y1, y2, dy ; strips
        float hits_Length[MAX_CLUST]; // x1, x2, dx ; time
        //--MMG1TRD
        float mmg1_clust_Xmax[MAX_CLUST];
        float mmg1_clust_Zmax[MAX_CLUST];
        float mmg1_clust_Emax[MAX_CLUST];
        float mmg1_clust_Xpos[MAX_CLUST];
        float mmg1_clust_Zpos[MAX_CLUST];
        float mmg1_clust_dEdx[MAX_CLUST];
        float mmg1_clust_Size[MAX_CLUST];
        float mmg1_clust_Width[MAX_CLUST][3];  // y1, y2, dy ; strips
        float mmg1_clust_Length[MAX_CLUST][3]; // x1, x2, dx ; time
        float mmg1_hits_Xpos[500];
        float mmg1_hits_Zpos[500];
        float mmg1_hits_dEdx[500];
        float mmg1_hits_Size[MAX_CLUST];
        float mmg1_hits_Width[MAX_CLUST];  // y1, y2, dy ; strips
        float mmg1_hits_Length[MAX_CLUST]; // x1, x2, dx ; time
        //--uRWELL-TRD
        float urw_clust_Xmax[MAX_CLUST];
        float urw_clust_Zmax[MAX_CLUST];
        float urw_clust_Emax[MAX_CLUST];
        float urw_clust_Xpos[MAX_CLUST];
        float urw_clust_Zpos[MAX_CLUST];
        float urw_clust_dEdx[MAX_CLUST];
        float urw_clust_Size[MAX_CLUST];
        float urw_clust_Width[MAX_CLUST][3];  // y1, y2, dy ; strips
        float urw_clust_Length[MAX_CLUST][3]; // x1, x2, dx ; time
        float urw_hits_Xpos[500];
        float urw_hits_Zpos[500];
        float urw_hits_dEdx[500];
        float urw_hits_Size[MAX_CLUST];
        float urw_hits_Width[MAX_CLUST];  // y1, y2, dy ; strips
        float urw_hits_Length[MAX_CLUST]; // x1, x2, dx ; time
        
        for (int k=0; k<MAX_CLUST; k++) {
	        clust_Xpos[k]=0; clust_Zpos[k]=0; clust_dEdx[k]=0;  clust_Size[k]=0;
	        clust_Xmax[k]=0; clust_Zmax[k]=0; clust_Emax[k]=0;
          clust_Width[k][0]=999999;   	clust_Width[k][1]=-999999;   	clust_Width[k][2]=0;
          clust_Length[k][0]=999999;  	clust_Length[k][1]=-999999;  	clust_Length[k][2]=0;
          
          mmg1_clust_Xpos[k]=0; mmg1_clust_Zpos[k]=0; mmg1_clust_dEdx[k]=0;  mmg1_clust_Size[k]=0;
          mmg1_clust_Xmax[k]=0; mmg1_clust_Zmax[k]=0; mmg1_clust_Emax[k]=0;
          mmg1_clust_Width[k][0]=999999;     mmg1_clust_Width[k][1]=-999999;    mmg1_clust_Width[k][2]=0;
          mmg1_clust_Length[k][0]=999999;    mmg1_clust_Length[k][1]=-999999;   mmg1_clust_Length[k][2]=0;
          
          urw_clust_Xpos[k]=0; urw_clust_Zpos[k]=0; urw_clust_dEdx[k]=0;  urw_clust_Size[k]=0;
	        urw_clust_Xmax[k]=0; urw_clust_Zmax[k]=0; urw_clust_Emax[k]=0;
          urw_clust_Width[k][0]=999999;   	urw_clust_Width[k][1]=-999999;   	urw_clust_Width[k][2]=0;
          urw_clust_Length[k][0]=999999;  	urw_clust_Length[k][1]=-999999;  	urw_clust_Length[k][2]=0;
        }
        int nclust=0, mmg1_nclust=0, urw_nclust=0;
        #if (USE_PULSE>0)
          TH2F* hp = hevtk; // -- hevtk and hevtck should be same bin size
          TH2F* hpc = hevtck;
        #else
          TH2F* hp = hevt; // -- hevt and hevtc should be same bin size
          TH2F* hpc = hevtc;
          TH2F* hmp = mhevt;
          TH2F* hmpc = mhevtc;
          TH2F* hup = uhevt;
          TH2F* hupc = uhevtc;
        #endif
        //--GEM
        int nx=hp->GetNbinsX();    int ny=hp->GetNbinsY();
        double xmi=hp->GetXaxis()->GetBinLowEdge(1);     double xma=hp->GetXaxis()->GetBinUpEdge(nx);
        double ymi=hp->GetYaxis()->GetBinLowEdge(1);     double yma=hp->GetYaxis()->GetBinUpEdge(ny);
        double binx = (xma-xmi)/nx;      double biny = (yma-ymi)/ny;
        //--MMG1
        int nmx=hmp->GetNbinsX();    int nmy=hmp->GetNbinsY();
        double xmmi=hmp->GetXaxis()->GetBinLowEdge(1);   double xmma=hmp->GetXaxis()->GetBinUpEdge(nmx);
        double ymmi=hmp->GetYaxis()->GetBinLowEdge(1);   double ymma=hmp->GetYaxis()->GetBinUpEdge(nmy);
        double binmx = (xmma-xmmi)/nmx;      double binmy = (ymma-ymmi)/nmy;
        //--URW
        int nux=hup->GetNbinsX();    int nuy=hup->GetNbinsY();
        double xumi=hup->GetXaxis()->GetBinLowEdge(1);     double xuma=hup->GetXaxis()->GetBinUpEdge(nux);
        double yumi=hup->GetYaxis()->GetBinLowEdge(1);     double yuma=hup->GetYaxis()->GetBinUpEdge(nuy);
        double binux = (xuma-xumi)/nux;      double binuy = (yuma-yumi)/nuy;
        #ifdef VERBOSE
          printf("nx=%d,ny=%d,xmi=%f,xma=%f,ymi=%f,yma=%f\n",nx,ny,xmi,xma,ymi,yma);
        #endif
        #if (USE_PULSE>0)
          float CL_DIST=3.3; // mm
          double THR2 = 0.01;
        #else
          float CL_DIST=3.; //2.9; // mm
          double THR2 = 2.; //1.2
        #endif
        
        for (int iy=0; iy<ny; iy++) {  //-------------------- Clustering Loop (GEMTRD) ------------------------------------
          for (int ix=0; ix<nx; ix++) {
            double c1 = hpc->GetBinContent(ix+1,iy+1);                    // energy
            double x1=double(ix)/double(nx)*(xma-xmi)+xmi+binx/2.;    // drift time
            double y1=double(iy)/double(ny)*(yma-ymi)+ymi+biny/2.;    // X strip
            ////if (c1>0) cout<<"      *** GEM-TRD EVENT="<<event_num<<" CL ENERGY="<<c1<<" FOR ix,iy="<<ix<<","<<iy<<" AND x1,y1="<<x1<<","<<y1<<endl;
            if (c1<THR2) continue;
            if (nclust==0) {
	            clust_Xpos[nclust]=y1; clust_Zpos[nclust]=x1;  clust_dEdx[nclust]=c1;  clust_Size[nclust]=1;
	            clust_Xmax[nclust]=y1; clust_Zmax[nclust]=x1;  clust_Emax[nclust]=c1;
              clust_Width[nclust][0]=y1;   	clust_Width[nclust][1]=y1;   	clust_Width[nclust][2]=0;
              clust_Length[nclust][0]=x1;  	clust_Length[nclust][1]=x1;  	clust_Length[nclust][2]=0;
              nclust++; continue;
            }
            int added=0;
            for (int k=0; k<nclust; k++) {
              double dist=sqrt(pow((y1-clust_Xpos[k]),2.)+pow((x1-clust_Zpos[k]),2.)); //--- dist hit to clusters
              ////if (dist>0) cout<<"      ****** GEM-TRD EVENT="<<event_num<<" CL_DIST="<<dist<<" FOR NCLUST="<<nclust<<" k="<<k<<", AND x1,y1="<<x1<<","<<y1<<endl;
              if (dist<CL_DIST) {
                clust_Xpos[k]=(y1*c1+clust_Xpos[k]*clust_dEdx[k])/(c1+clust_dEdx[k]);  //--  new X pos
                clust_Zpos[k]=(x1*c1+clust_Zpos[k]*clust_dEdx[k])/(c1+clust_dEdx[k]);  //--  new Z pos
	              if (c1>clust_Emax[k]) {
		              clust_Xmax[k]=y1;
		              clust_Zmax[k]=x1;
		              clust_Emax[k]=c1;
	              }
                clust_dEdx[k]=c1+clust_dEdx[k];  // new dEdx
                clust_Size[k]=1+clust_Size[k];  // clust size in pixels
                if (y1<clust_Width[k][0]) clust_Width[k][0]=y1; if (y1>clust_Width[k][1]) clust_Width[k][1]=y1; clust_Width[k][2]=clust_Width[k][1]-clust_Width[k][0];
                if (x1<clust_Length[k][0]) clust_Length[k][0]=x1;if (x1>clust_Length[k][1]) clust_Length[k][1]=x1;clust_Length[k][2]=clust_Length[k][1]-clust_Length[k][0];
                hpc->SetBinContent(ix,iy,k+1.);
                added=1; break;
              }
            }
            if (added==0) {
              if (nclust+1>=MAX_CLUST) continue;
	            clust_Xpos[nclust]=y1; clust_Zpos[nclust]=x1;  clust_dEdx[nclust]=c1;  clust_Size[nclust]=1;
	            clust_Xmax[nclust]=y1; clust_Zmax[nclust]=x1;  clust_Emax[nclust]=c1;
              clust_Width[nclust][0]=y1;   	clust_Width[nclust][1]=y1;   	clust_Width[nclust][2]=0;
              clust_Length[nclust][0]=x1;  	clust_Length[nclust][1]=x1;  	clust_Length[nclust][2]=0;
              nclust++;
            }
          }
        } //---------------------- End Clustering Loop (GEMTRD) ------------------------------
        
        for (int iy=0; iy<nmy; iy++) {  //-------------------- Clustering Loop (MMG1TRD) ------------------------------------
          for (int ix=0; ix<nmx; ix++) {
            double c1 = hmpc->GetBinContent(ix+1,iy+1);                       // energy
            //if (c1>0) cout<<"      *** MMG-TRD CL ENERGY="<<c1<<" FOR ix,iy="<<ix<<endl;
            double x1=double(ix)/double(nmx)*(xmma-xmmi)+xmmi+binmx/2.;   // drift time
            //cout<<"***MMG1*** nmx="<<nmx<<" xmma="<<xmma<<" xmmi="<<xmmi<<" binmx="<<binmx<<endl;
            double y1=double(iy)/double(nmy)*(ymma-ymmi)+ymmi+binmy/2.;        // X strip
            if (c1<THR2) continue;
            if (mmg1_nclust==0) {
              mmg1_clust_Xpos[mmg1_nclust]=y1; mmg1_clust_Zpos[mmg1_nclust]=x1;  mmg1_clust_dEdx[mmg1_nclust]=c1;  mmg1_clust_Size[mmg1_nclust]=1;
              mmg1_clust_Xmax[mmg1_nclust]=y1; mmg1_clust_Zmax[mmg1_nclust]=x1;  mmg1_clust_Emax[mmg1_nclust]=c1;
              mmg1_clust_Width[mmg1_nclust][0]=y1;    mmg1_clust_Width[mmg1_nclust][1]=y1;    mmg1_clust_Width[mmg1_nclust][2]=0;
              mmg1_clust_Length[mmg1_nclust][0]=x1;   mmg1_clust_Length[mmg1_nclust][1]=x1;   mmg1_clust_Length[mmg1_nclust][2]=0;
              mmg1_nclust++; continue;
            }
            int mmg1_added=0;
            for (int k=0; k<mmg1_nclust; k++) {
              double dist=sqrt(pow((y1-mmg1_clust_Xpos[k]),2.)+pow((x1-mmg1_clust_Zpos[k]),2.)); //--- dist hit to clusters
              //if (dist>0) cout<<"      *** MMG-TRD DIST="<<dist<<" FOR NCLUST="<<mmg1_nclust<<" k="<<k<<endl;
              if (dist<CL_DIST) {
                mmg1_clust_Xpos[k]=(y1*c1+mmg1_clust_Xpos[k]*mmg1_clust_dEdx[k])/(c1+mmg1_clust_dEdx[k]);  //--  new X pos
                mmg1_clust_Zpos[k]=(x1*c1+mmg1_clust_Zpos[k]*mmg1_clust_dEdx[k])/(c1+mmg1_clust_dEdx[k]);  //--  new Z pos
                if (c1>mmg1_clust_Emax[k]) {
		              mmg1_clust_Xmax[k]=y1;
		              mmg1_clust_Zmax[k]=x1;
		              mmg1_clust_Emax[k]=c1;
	              }
                mmg1_clust_dEdx[k]=c1+mmg1_clust_dEdx[k];  // new dEdx
                mmg1_clust_Size[k]=1+mmg1_clust_Size[k];  // clust size in pixels
                if (y1<mmg1_clust_Width[k][0]) mmg1_clust_Width[k][0]=y1; if (y1>mmg1_clust_Width[k][1]) mmg1_clust_Width[k][1]=y1; mmg1_clust_Width[k][2]=mmg1_clust_Width[k][1]-mmg1_clust_Width[k][0];
                if (x1<mmg1_clust_Length[k][0]) mmg1_clust_Length[k][0]=x1;if (x1>mmg1_clust_Length[k][1]) mmg1_clust_Length[k][1]=x1; mmg1_clust_Length[k][2]=mmg1_clust_Length[k][1]-mmg1_clust_Length[k][0];
                hmpc->SetBinContent(ix,iy,k+1.);
                mmg1_added=1; break;
              }
            }
            if (mmg1_added==0) {
              if (mmg1_nclust+1>=MAX_CLUST) continue;
              mmg1_clust_Xpos[mmg1_nclust]=y1; mmg1_clust_Zpos[mmg1_nclust]=x1;  mmg1_clust_dEdx[mmg1_nclust]=c1;  mmg1_clust_Size[mmg1_nclust]=1;
              mmg1_clust_Xmax[mmg1_nclust]=y1; mmg1_clust_Zmax[mmg1_nclust]=x1;  mmg1_clust_Emax[mmg1_nclust]=c1;
              mmg1_clust_Width[mmg1_nclust][0]=y1;    mmg1_clust_Width[mmg1_nclust][1]=y1;    mmg1_clust_Width[mmg1_nclust][2]=0;
              mmg1_clust_Length[mmg1_nclust][0]=x1;   mmg1_clust_Length[mmg1_nclust][1]=x1;   mmg1_clust_Length[mmg1_nclust][2]=0;
              mmg1_nclust++;
            }
          }
        } //---------------------- End Clustering Loop (MMG1TRD) ------------------------------
        
        for (int iy=0; iy<nuy; iy++) {  //-------------------- Clustering Loop (uRWELL-TRD) ------------------------------------
          for (int ix=0; ix<nux; ix++) {
            double c1 = hupc->GetBinContent(ix+1,iy+1);                       // energy
            //if (c1>0) cout<<"      *** URW-TRD CL ENERGY="<<c1<<" FOR ix,iy="<<ix<<endl;
            double x1=double(ix)/double(nux)*(xuma-xumi)+xumi+binux/2.;   // drift time
            //cout<<"***URW*** nux="<<nux<<" xuma="<<xuma<<" xumi="<<xumi<<" binux="<<binux<<endl;
            double y1=double(iy)/double(nuy)*(yuma-yumi)+yumi+binuy/2.;        // X strip
            if (c1<(THR2-1.0)) continue;
            if (urw_nclust==0) {
              urw_clust_Xpos[urw_nclust]=y1; urw_clust_Zpos[urw_nclust]=x1;  urw_clust_dEdx[urw_nclust]=c1;  urw_clust_Size[urw_nclust]=1;
              urw_clust_Xmax[urw_nclust]=y1; urw_clust_Zmax[urw_nclust]=x1;  urw_clust_Emax[urw_nclust]=c1;
              urw_clust_Width[urw_nclust][0]=y1;    urw_clust_Width[urw_nclust][1]=y1;    urw_clust_Width[urw_nclust][2]=0;
              urw_clust_Length[urw_nclust][0]=x1;   urw_clust_Length[urw_nclust][1]=x1;   urw_clust_Length[urw_nclust][2]=0;
              urw_nclust++; continue;
            }
            int urw_added=0;
            for (int k=0; k<urw_nclust; k++) {
              double dist=sqrt(pow((y1-urw_clust_Xpos[k]),2.)+pow((x1-urw_clust_Zpos[k]),2.)); //--- dist hit to clusters
              //if (dist>0) cout<<"      *** URW-TRD DIST="<<dist<<" FOR NCLUST="<<urw_nclust<<" k="<<k<<endl;
              if (dist<CL_DIST) {
                urw_clust_Xpos[k]=(y1*c1+urw_clust_Xpos[k]*urw_clust_dEdx[k])/(c1+urw_clust_dEdx[k]);  //--  new X pos
                urw_clust_Zpos[k]=(x1*c1+urw_clust_Zpos[k]*urw_clust_dEdx[k])/(c1+urw_clust_dEdx[k]);  //--  new Z pos
                if (c1>urw_clust_Emax[k]) {
		              urw_clust_Xmax[k]=y1;
		              urw_clust_Zmax[k]=x1;
		              urw_clust_Emax[k]=c1;
	              }
                urw_clust_dEdx[k]=c1+urw_clust_dEdx[k];  // new dEdx
                urw_clust_Size[k]=1+urw_clust_Size[k];  // clust size in pixels
                if (y1<urw_clust_Width[k][0]) urw_clust_Width[k][0]=y1; if (y1>urw_clust_Width[k][1]) urw_clust_Width[k][1]=y1; urw_clust_Width[k][2]=urw_clust_Width[k][1]-urw_clust_Width[k][0];
                if (x1<urw_clust_Length[k][0]) urw_clust_Length[k][0]=x1; if (x1>urw_clust_Length[k][1]) urw_clust_Length[k][1]=x1; urw_clust_Length[k][2]=urw_clust_Length[k][1]-urw_clust_Length[k][0];
                hupc->SetBinContent(ix,iy,k+1.);
                urw_added=1; break;
              }
            }
            if (urw_added==0) {
              if (urw_nclust+1>=MAX_CLUST) continue;
              urw_clust_Xpos[urw_nclust]=y1; urw_clust_Zpos[urw_nclust]=x1;  urw_clust_dEdx[urw_nclust]=c1;  urw_clust_Size[urw_nclust]=1;
              urw_clust_Xmax[urw_nclust]=y1; urw_clust_Zmax[urw_nclust]=x1;  urw_clust_Emax[urw_nclust]=c1;
              urw_clust_Width[urw_nclust][0]=y1;    urw_clust_Width[urw_nclust][1]=y1;    urw_clust_Width[urw_nclust][2]=0;
              urw_clust_Length[urw_nclust][0]=x1;   urw_clust_Length[urw_nclust][1]=x1;   urw_clust_Length[urw_nclust][2]=0;
              urw_nclust++;
            }
          }
        } //---------------------- End Clustering Loop (uRWELL-TRD) ------------------------------
        
        #if (USE_PULSE>0)
          int MinClustSize=1;
          double MinClustWidth=0.00;
          double MinClustLength=0.0;
          double MaxClustLength=5.;
          double zStart =  0.; // mm
          double zEnd   = 29.; // mm
        #else
          int MinClustSize=5;//1; // pixels
          double MinClustWidth=0.05;//0.001;
          double MinClustLength=0.01;
          double MaxClustLength=3.5; //4
          double zStart =  0.; // mm
          double zEnd   = 30.; // mm
        #endif
        
        double maxClust_dEdx=0., totalClust_dEdx=0.;
        int ii=0;
        #ifdef VERBOSE
          printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
          printf("                Xpos   Ypos   Zpos       E    Width  Length   Size \n");
        #endif
        for (int k=0; k<nclust; k++) {
          #ifdef VERBOSE
            if (k<30) printf("%2d Clust(%2d): %6.1f %6.1f %8.1f %6.2f %6.2f %8.1f  ",k,k+1,clust_Xpos[k],clust_Zpos[k],clust_dEdx[k],clust_Width[k][2],clust_Length[k][2],clust_Size[k]);
          #endif
          //-------------  Cluster Filter (GEMTRD) -----------------
          if ((clust_Size[k]>=MinClustSize && zStart<clust_Zpos[k] && clust_Zpos[k]<zEnd && clust_Width[k][2]>MinClustWidth) || clust_Length[k][2]<MaxClustLength) {
            if (match) {
        	    #if (USE_MAXPOS>0)
	              hgemClusterDiff_el->Fill(mmg1_xchanmax-((clust_Xmax[k]+3.2)*0.908525-9.22087));
                hits_Xpos[ii]=clust_Xmax[k];
	              hits_Zpos[ii]=clust_Zmax[k];
                //cout<<"  *** GEM-TRD HITS_ZPOSITION="<<clust_Zmax[k]<<" FOR k="<<k<<endl;
              #else
                hgemClusterDiff_el->Fill(mmg1_xchanmax-((clust_Xpos[k]+3.2)*0.908525-9.22087));
                hits_Xpos[ii]=clust_Xpos[k];
        	      hits_Zpos[ii]=clust_Zpos[k];
                //cout<<"  *** GEM-TRD HITS_ZPOSITION="<<clust_Zpos[k]<<" FOR k="<<k<<endl;
              #endif
        	    hits_dEdx[ii]=clust_dEdx[k];
              hits_Width[ii]=clust_Width[k][2];
              hits_Length[ii]=clust_Length[k][2];
        	    ii++;
        	    if (clust_dEdx[k]>maxClust_dEdx) maxClust_dEdx=clust_dEdx[k];
              totalClust_dEdx+=clust_dEdx[k];
              #ifdef VERBOSE
                if (k<30) printf("\n");
              #endif
            }
        	} else {
          #ifdef VERBOSE
            if (k<30) printf(" <--- skip \n");
          #endif
          }
        }
        int nhits=ii;
        clu_dedx_tot=totalClust_dEdx;
        
        double maxClust_m_dEdx=0., totalClust_m_dEdx=0.;
        int mmg1_ii=0;
        for (int k=0; k<mmg1_nclust; k++) {
          //-------------  Cluster Filter (MMG1TRD) -----------------
          if ((mmg1_clust_Size[k]>=MinClustSize && zStart<mmg1_clust_Zpos[k] && mmg1_clust_Zpos[k]<zEnd && mmg1_clust_Width[k][2]>MinClustWidth) || mmg1_clust_Length[k][2]<MaxClustLength+2.) {
            if (match_mmg1) {
              #if (USE_MAXPOS>0)
                hmmg1ClusterDiff_el->Fill(tgem_xchanmax*0.908525-(mmg1_clust_Xmax[k]+3.2)-9.22087);
                mmg1_hits_Xpos[mmg1_ii]=mmg1_clust_Xmax[k];
                mmg1_hits_Zpos[mmg1_ii]=mmg1_clust_Zmax[k];
              #else
                hmmg1ClusterDiff_el->Fill(tgem_xchanmax*0.908525-(mmg1_clust_Xpos[k]+3.2)-9.22087);
                mmg1_hits_Xpos[mmg1_ii]=mmg1_clust_Xpos[k];
                mmg1_hits_Zpos[mmg1_ii]=mmg1_clust_Zpos[k];
              #endif
              mmg1_hits_dEdx[mmg1_ii]=mmg1_clust_dEdx[k];
              mmg1_hits_Width[mmg1_ii]=mmg1_clust_Width[k][2];
              mmg1_hits_Length[mmg1_ii]=mmg1_clust_Length[k][2];
              mmg1_ii++;
              if (mmg1_clust_dEdx[k]>maxClust_m_dEdx) maxClust_m_dEdx=mmg1_clust_dEdx[k];
              totalClust_m_dEdx+=mmg1_clust_dEdx[k];
            }
          }
        }
        int mmg1_nhits=mmg1_ii;
        mmg1_clu_dedx_tot=totalClust_m_dEdx;
        
        double maxClust_u_dEdx=0., totalClust_u_dEdx=0.;
        int urw_ii=0;
        for (int k=0; k<urw_nclust; k++) {
          //-------------  Cluster Filter (uRWELL-TRD) -----------------
          if ((urw_clust_Size[k]>=(MinClustSize-2) && zStart<urw_clust_Zpos[k] && urw_clust_Zpos[k]<zEnd && urw_clust_Width[k][2]>MinClustWidth) || urw_clust_Length[k][2]<MaxClustLength) {
            if (match_urw) {
              #if (USE_MAXPOS>0)
                hurwClusterDiff_el->Fill((urw_clust_Xmax[k]+3.2)*0.881493-tgem_xchanmax-7.886);
                urw_hits_Xpos[urw_ii]=urw_clust_Xmax[k];
                urw_hits_Zpos[urw_ii]=urw_clust_Zmax[k];
              #else
                hurwClusterDiff_el->Fill((urw_clust_Xpos[k]+3.2)*0.881493-tgem_xchanmax-7.886);
                urw_hits_Xpos[urw_ii]=urw_clust_Xpos[k];
                urw_hits_Zpos[urw_ii]=urw_clust_Zpos[k];
              #endif
              urw_hits_dEdx[urw_ii]=urw_clust_dEdx[k];
              urw_hits_Width[urw_ii]=urw_clust_Width[k][2];
              urw_hits_Length[urw_ii]=urw_clust_Length[k][2];
              urw_ii++;
              if (urw_clust_dEdx[k]>maxClust_u_dEdx) maxClust_u_dEdx=urw_clust_dEdx[k];
              totalClust_u_dEdx+=urw_clust_dEdx[k];
            }
          }
        }
        int urw_nhits=urw_ii;
        urw_clu_dedx_tot=totalClust_u_dEdx;
        // ----------------------- end hist dist clustering ---------------------------------
        
        //=================================== Draw HITS and CLUST  ============================================
        #ifdef SHOW_EVTbyEVT
          char hevtTitle[80]; sprintf(hevtTitle,"GEM-TRD: Event=%lld Run=%d; z pos [mm]; y pos [mm]",jentry,RunNum);
          hevt->SetTitle(hevtTitle);
          #if (USE_PULSE>0)
            hevtk->SetTitle(hevtTitle);
          #endif
          char mhevtTitle[80]; sprintf(mhevtTitle,"MMG1-TRD: Event=%lld Run=%d; z pos [mm]; y pos [mm]",jentry,RunNum);
          mhevt->SetTitle(mhevtTitle);
          char uhevtTitle[80]; sprintf(uhevtTitle,"URW-TRD: Event=%lld Run=%d; z pos [mm]; y pos [mm]",jentry,RunNum);
          uhevt->SetTitle(uhevtTitle);
          #ifdef VERBOSE
            printf("hits_SIZE=%d  Clust size = %d \n",nhits,nclust);
          #endif
          	c2->cd(1); gPad->Modified(); gPad->Update();
          	int COLMAP[]={1,2,3,4,6,5};
          	int pmt=22, pmt0 = 20; // PM type
            int max2draw = nclust;
            for (int i=0; i<max2draw; i++) {
              #if (USE_MAXPOS>0)
	              TMarker m = TMarker(clust_Zmax[i],clust_Xmax[i],pmt);
              #else
          	    TMarker m = TMarker(clust_Zpos[i], clust_Xpos[i], pmt);
              #endif
          	  int tcol=2;
          	  if (clust_Size[i]<MinClustSize) pmt=22; else pmt=pmt0;
          	  int mcol = COLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(pmt);
          	  m.SetMarkerSize(0.7+clust_dEdx[i]/300);
          	  m.DrawClone();
          	}
            gPad->Modified(); gPad->Update();
            
            //--MMG1TRD
            c2->cd(6); gPad->Modified(); gPad->Update();
            //int mCOLMAP[]={1,2,3,4,6,5};
            //int mpmt=22, mpmt0 = 20; // PM type
            int mmax2draw = mmg1_nclust;
            for (int i=0; i<mmax2draw; i++) {
              #if (USE_MAXPOS>0)
	              TMarker m = TMarker(mmg1_clust_Zmax[i], mmg1_clust_Xmax[i], pmt);
              #else
                TMarker m = TMarker(mmg1_clust_Zpos[i], mmg1_clust_Xpos[i], pmt);
              #endif
              int tcol=2;
              if (mmg1_clust_Size[i]<MinClustSize) pmt=22; else pmt=pmt0;
              int mcol = COLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(pmt);
              m.SetMarkerSize(0.7+mmg1_clust_dEdx[i]/300);
              m.DrawClone();
            }
            gPad->Modified(); gPad->Update();
            
            //--uRWELL-TRD
            c2->cd(11); gPad->Modified(); gPad->Update();
            //int uCOLMAP[]={1,2,3,4,6,5};
            //int upmt=22, upmt0 = 20; // PM type
            int umax2draw = urw_nclust;
            for (int i=0; i<umax2draw; i++) {
              #if (USE_MAXPOS>0)
	              TMarker m = TMarker(urw_clust_Zmax[i], urw_clust_Xmax[i], pmt);
              #else
                TMarker m = TMarker(urw_clust_Zpos[i], urw_clust_Xpos[i], pmt);
              #endif
              int tcol=2;
              if (urw_clust_Size[i]<MinClustSize) pmt=22; else pmt=pmt0;
              int mcol = COLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(pmt);
              m.SetMarkerSize(0.7+urw_clust_dEdx[i]/300);
              m.DrawClone();
            }
            gPad->Modified(); gPad->Update();
        #endif
        
        #if (USE_GNN==1)   // GNN MC
          //----------------------------------------------------------
          //--   Send to Model simulation
          //----------------------------------------------------------
          
          #ifdef VERBOSE
            printf("**> Start Model simulation nclust=%d nhits=%d \n",nclust,nhits);
          #endif
          std::vector<int> tracks(nhits, 0);
          std::vector<float> Xcl;
          std::vector<float> Zcl;
          Xcl.clear();
          Zcl.clear();
          for (int n=0; n<nhits; n++) {
          	Xcl.push_back(hits_Xpos[n]);
          	Zcl.push_back(hits_Zpos[n]);
          }
          doPattern(Xcl, Zcl, tracks);  //---- call GNN ---
          //-- MMG1
          std::vector<int> mmg1_tracks(mmg1_nhits, 0);
          std::vector<float> mmg1_Xcl;
          std::vector<float> mmg1_Zcl;
          mmg1_Xcl.clear();
          mmg1_Zcl.clear();
          for (int n=0; n<mmg1_nhits; n++) {
            mmg1_Xcl.push_back(mmg1_hits_Xpos[n]);
            mmg1_Zcl.push_back(mmg1_hits_Zpos[n]);
          }
          doPattern(mmg1_Xcl, mmg1_Zcl, mmg1_tracks);  //---- call GNN ---
          //-- URW
          std::vector<int> urw_tracks(urw_nhits, 0);
          std::vector<float> urw_Xcl;
          std::vector<float> urw_Zcl;
          urw_Xcl.clear();
          urw_Zcl.clear();
          for (int n=0; n<urw_nhits; n++) {
            urw_Xcl.push_back(urw_hits_Xpos[n]);
            urw_Zcl.push_back(urw_hits_Zpos[n]);
          }
          doPattern(urw_Xcl, urw_Zcl, urw_tracks);  //---- call GNN ---
          #ifdef VERBOSE
            printf("**> End Model simulation \n"); //===================================================
          #endif
          #ifdef SHOW_EVTbyEVT
              c2->cd(2); gPad->Modified(); gPad->Update();
              int COLMAP2[]={1,2,3,4,6,5};
              for(ULong64_t i=0; i<tracks.size(); i++) {
                #ifdef VERBOSE
                  if (i<30) printf("i=%d trk=%d |  %8.2f,%8.2f\n",i, tracks[i], Xcl[i], Zcl[i]);
              	#endif
                TMarker m = TMarker(hits_Zpos[i], hits_Xpos[i], 24);
              	int tcol = min(tracks[i], 6);
              	int mcol = COLMAP2[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(41);     m.SetMarkerSize(1.5);
            	  m.DrawClone();  gPad->Modified(); gPad->Update();
              }
              //-- MMG1
              c2->cd(7); gPad->Modified(); gPad->Update();
              //int mCOLMAP2[]={1,2,3,4,6,5};
              for(ULong64_t i=0; i<mmg1_tracks.size(); i++) {
                TMarker m = TMarker(mmg1_hits_Zpos[i], mmg1_hits_Xpos[i], 24);
                int tcol = min(mmg1_tracks[i], 6);
                int mcol = COLMAP2[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(41);     m.SetMarkerSize(1.5);
                m.DrawClone();  gPad->Modified(); gPad->Update();
              }
              //-- uRWELL
              c2->cd(12); gPad->Modified(); gPad->Update();
              //int uCOLMAP2[]={1,2,3,4,6,5};
              for(ULong64_t i=0; i<urw_tracks.size(); i++) {
                TMarker m = TMarker(urw_hits_Zpos[i], urw_hits_Xpos[i], 24);
                int tcol = min(urw_tracks[i], 6);
                int mcol = COLMAP2[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(41);     m.SetMarkerSize(1.5);
                m.DrawClone();  gPad->Modified(); gPad->Update();
              }
              #ifdef VERBOSE
                printf("\n\n");
                printf("**> End Cluster Plot \n");
              #endif
          #endif
          //--------------------------------------------------
          //----           Track fitting                 -----
          //--------------------------------------------------
          
          #ifdef VERBOSE
            printf("==> GNN: tracks sort  : trk_siz=%ld \r\n", tracks.size());
          #endif
          //-----------------   tracks sorting -------------
          std::vector<std::vector<float>> TRACKS;
          TRACKS.resize(nhits);
          std::vector<int>  TRACKS_N(nhits, 0);
          for (int i=0; i<nhits; i++)  { TRACKS_N[i] = 0;  }
          for (int i2=0; i2<nhits; i2++) {
          	int num =  tracks[i2];
          	int num2 = std::max(0, std::min(num, nhits - 1));
            #ifdef VERBOSE
              if (i2<20) printf("==> lstm3:track sort i=%d  : num=%d(%d) x=%f z=%f \n", i2, num, num2,  Xcl[i2],Zcl[i2]);
          	#endif
          	TRACKS[num2].push_back(Xcl[i2]);
            TRACKS[num2].push_back(Zcl[i2]);
          	TRACKS_N[num2]++;
          }
          //-- MMG1
          std::vector<std::vector<float>> mmg1_TRACKS;
          mmg1_TRACKS.resize(mmg1_nhits);
          std::vector<int>  mmg1_TRACKS_N(mmg1_nhits, 0);
          for (int i=0; i<mmg1_nhits; i++)  { mmg1_TRACKS_N[i] = 0;  }
          for (int i2=0; i2<mmg1_nhits; i2++) {
            int num =  mmg1_tracks[i2];
            int num2 = std::max(0, std::min(num, mmg1_nhits - 1));
            mmg1_TRACKS[num2].push_back(mmg1_Xcl[i2]);
            mmg1_TRACKS[num2].push_back(mmg1_Zcl[i2]);
            mmg1_TRACKS_N[num2]++;
          }
          //-- uRWELL
          std::vector<std::vector<float>> urw_TRACKS;
          urw_TRACKS.resize(urw_nhits);
          std::vector<int>  urw_TRACKS_N(urw_nhits, 0);
          for (int i=0; i<urw_nhits; i++)  { urw_TRACKS_N[i] = 0;  }
          for (int i2=0; i2<urw_nhits; i2++) {
            int num =  urw_tracks[i2];
            int num2 = std::max(0, std::min(num, urw_nhits - 1));
            urw_TRACKS[num2].push_back(urw_Xcl[i2]);
            urw_TRACKS[num2].push_back(urw_Zcl[i2]);
            urw_TRACKS_N[num2]++;
          }
          #if (DEBUG > 1)
            for (int i2 = 0; i2 < nhits; i2++) {
              printf(" trdID=%d n_hits=%d v_size=%d \n",i2,TRACKS_N[i2],TRACKS[i2].size());
              for (ULong64_t i3 = 0; i3 < TRACKS[i2].size(); i3+=2) {
                printf(" trkID=%d  hit=%d x=%f z=%f \n",i2,i3/2,TRACKS[i2].at(i3),TRACKS[i2].at(i3+1));
              }
              if ( TRACKS_N[i2]>0) printf("\n");
            }
          #endif
          //------------------ end tracks sorting --------------------
          
          #if (USE_FIT==1)
            //-----------------------------------
            //---       linear fitting        ---
            //-----------------------------------
            static TMultiGraph *mg;
            if (mg != NULL ) delete mg;
            mg = new TMultiGraph();
            int NTRACKS=0;
            int MIN_HITS=2;
            Double_t p0, p1;
            
            for (int i2=1; i2<nhits; i2++) {  //-- GEM tracks loop; zero track -> noise
              
             if (TRACKS_N[i2]<MIN_HITS) continue;   //---- select 2 (x,z) and more hits on track ----
            	#ifdef VERBOSE
                printf("==> fit: start trk: %d \r\n", i2);
              #endif
            	std::vector<Double_t> x;
            	std::vector<Double_t> y;
            	for (int i3=0; i3<(int)TRACKS[i2].size(); i3+=2) {
            	  #ifdef VERBOSE
                  printf(" trkID=%d  hit=%d x=%f z=%f \n",i2,i3/2,TRACKS[i2].at(i3),TRACKS[i2].at(i3+1));
            	  #endif
                x.push_back(TRACKS[i2].at(i3+1));
            	  y.push_back(TRACKS[i2].at(i3));
            	}
              #ifdef SHOW_EVTbyEVT
              	gErrorIgnoreLevel = kBreak; // Suppress warning messages from empty fit data
              	TGraph *g = new TGraph(TRACKS_N[i2], &x[0], &y[0]);  g->SetMarkerStyle(21); g->SetMarkerColor(i2);
              	TF1 *f = new TF1("f", "[1] * x + [0]");
              	g->Fit(f,"Q");
                //  --- get fit parameters ---
                TF1 *ffunc=g->GetFunction("f");
                p0=ffunc->GetParameter(0);
                p1=ffunc->GetParameter(1);
                Double_t chi2x_nn = ffunc->GetChisquare();
                Double_t Ndfx_nn = ffunc->GetNDF();
                double chi2nn=chi2x_nn/Ndfx_nn;
                #ifdef VERBOSE
                  printf("+++++>>  Track = %d fit: p0=%f p1=%f (%f deg) ff(15)=%f chi2nn=%f \n",i2,p0,p1,p1/3.1415*180.,ffunc->Eval(15.),chi2nn);
                #endif
            	  mg->Add(g,"p");
              #endif
              NTRACKS++;
            }  //-- end GEM tracks loop --
            
            //-- MMG1
            static TMultiGraph *mmg1_mg;
            if (mmg1_mg != NULL ) delete mmg1_mg;
            mmg1_mg = new TMultiGraph();
            int mmg1_NTRACKS=0;
            int mmg1_MIN_HITS=2;
            Double_t mp0, mp1;
            
            for (int i2=1; i2<mmg1_nhits; i2++) {  //-- MMG1 tracks loop; zero track -> noise
              
              if (mmg1_TRACKS_N[i2]<mmg1_MIN_HITS) continue;   //---- select 2 (x,z) or more hits on track ----
              std::vector<Double_t> mmg1_x;
              std::vector<Double_t> mmg1_y;
              for (int i3=0; i3<(int)mmg1_TRACKS[i2].size(); i3+=2) {
                mmg1_x.push_back(mmg1_TRACKS[i2].at(i3+1));
                mmg1_y.push_back(mmg1_TRACKS[i2].at(i3));
              }
              #ifdef SHOW_EVTbyEVT
                gErrorIgnoreLevel = kBreak; // Suppress warning messages from empty fit data
                TGraph *mmg1_g = new TGraph(mmg1_TRACKS_N[i2], &mmg1_x[0], &mmg1_y[0]);  mmg1_g->SetMarkerStyle(21); mmg1_g->SetMarkerColor(i2);
                TF1 *mmg1_f = new TF1("mmg1_f", "[1] * x + [0]");
                mmg1_g->Fit(mmg1_f,"Q");
                //  --- get fit parameters ---
                TF1 *mmg1_ffunc = mmg1_g->GetFunction("mmg1_f");
                mp0 = mmg1_ffunc->GetParameter(0);
                mp1 = mmg1_ffunc->GetParameter(1);
                Double_t mmg1_chi2x_nn = mmg1_ffunc->GetChisquare();
                Double_t mmg1_Ndfx_nn = mmg1_ffunc->GetNDF();
                double mmg1_chi2nn = mmg1_chi2x_nn/mmg1_Ndfx_nn;
                mmg1_mg->Add(mmg1_g, "p");
              #endif
              mmg1_NTRACKS++;
            }  //-- end MMG1 tracks loop --
            
            //-- uRWELL
            static TMultiGraph *urw_mg;
            if (urw_mg != NULL ) delete urw_mg;
            urw_mg = new TMultiGraph();
            int urw_NTRACKS=0;
            int urw_MIN_HITS=2;
            Double_t up0, up1;
            
            for (int i2=1; i2<urw_nhits; i2++) {  //-- uRWELL tracks loop; zero track -> noise
              
              if (urw_TRACKS_N[i2]<urw_MIN_HITS) continue;   //---- select 2 (x,z) or more hits on track ----
              std::vector<Double_t> urw_x;
              std::vector<Double_t> urw_y;
              for (int i3=0; i3<(int)urw_TRACKS[i2].size(); i3+=2) {
                urw_x.push_back(urw_TRACKS[i2].at(i3+1));
                urw_y.push_back(urw_TRACKS[i2].at(i3));
              }
              #ifdef SHOW_EVTbyEVT
                gErrorIgnoreLevel = kBreak; // Suppress warning messages from empty fit data
                TGraph *urw_g = new TGraph(urw_TRACKS_N[i2], &urw_x[0], &urw_y[0]);  urw_g->SetMarkerStyle(21); urw_g->SetMarkerColor(i2);
                TF1 *urw_f = new TF1("urw_f", "[1] * x + [0]");
                urw_g->Fit(urw_f,"Q");
                //  --- get fit parameters ---
                TF1 *urw_ffunc = urw_g->GetFunction("urw_f");
                up0 = urw_ffunc->GetParameter(0);
                up1 = urw_ffunc->GetParameter(1);
                Double_t urw_chi2x_nn = urw_ffunc->GetChisquare();
                Double_t urw_Ndfx_nn = urw_ffunc->GetNDF();
                double urw_chi2nn = urw_chi2x_nn/urw_Ndfx_nn;
                urw_mg->Add(urw_g, "p");
              #endif
              urw_NTRACKS++;
            }  //-- end uRWELL tracks loop --
            
            #ifdef SHOW_EVTbyEVT
              if (NTRACKS<1 && mmg1_NTRACKS<1 && urw_NTRACKS<1) continue;  // --- skip event ----
              //if (nhits<3)    continue;  // --- skip event ----
              ///////////if (gem_trk_hit<1) continue;
                char mgTitle[80]; sprintf(mgTitle,"GEM ML-FPGA response, #Tracks=%d; z pos [mm]; y pos [mm]",NTRACKS);
                mg->SetTitle(mgTitle);
                c2->cd(3); mg->Draw("AP");
                mg->GetXaxis()->SetLimits(Xmin,Xmax);
                mg->SetMinimum(Ymin);
                mg->SetMaximum(Ymax);
                gPad->Modified(); gPad->Update();
                //-- MMG1
                char mmg1_mgTitle[80]; sprintf(mmg1_mgTitle,"MMG1 ML-FPGA response, #Tracks=%d; z pos [mm]; y pos [mm]",mmg1_NTRACKS);
                mmg1_mg->SetTitle(mmg1_mgTitle);
                c2->cd(8); mmg1_mg->Draw("AP");
                mmg1_mg->GetXaxis()->SetLimits(Xmin,Xmax); //XMAX NEEDS FIXED!!!!
                mmg1_mg->SetMinimum(Ymin);
                mmg1_mg->SetMaximum(Ymax);
                gPad->Modified(); gPad->Update();
                //-- uRWELL
                char urw_mgTitle[80]; sprintf(urw_mgTitle,"URW ML-FPGA response, #Tracks=%d; z pos [mm]; y pos [mm]",urw_NTRACKS);
                urw_mg->SetTitle(urw_mgTitle);
                c2->cd(13); urw_mg->Draw("AP");
                urw_mg->GetXaxis()->SetLimits(Xmin,Xmax); //XMAX NEEDS FIXED!!!! (?)
                urw_mg->SetMinimum(Ymin);
                urw_mg->SetMaximum(uYmax);
                gPad->Modified(); gPad->Update();
            #endif
          #endif // USE_FIT
        #endif // USE_GNN MC
        
        //******************************************************************************
        #ifdef SHOW_EVTbyEVT
            cout<<"Event#="<<event_num<<" #ofGEMTracks="<<NTRACKS<<" #ofMMGTracks="<<mmg1_NTRACKS<<endl;
            cout<<"GEM_External_Track_Match="<<match<<" MMG_External_Track_Match="<<match_mmg1<<" URW_External_Track_Match="<<match_urw<<endl;
            #ifdef WRITE_CSV
              WriteToCSV(csvFile,event_num,NTRACKS,chi2cc_gem);
            #endif
            //c3->cd(1);  hCal_sum->Draw();         gPad->Modified();   gPad->Update();
            //c3->cd(2);  hCal_pulse->Draw("hist"); gPad->Modified();   gPad->Update();
            //c3->cd(3);  hCher_pulse->Draw("hist"); gPad->Modified();  gPad->Update();
            //c3->cd(4);  hPresh_pulse->Draw("hist"); gPad->Modified(); gPad->Update();
            //c3->cd(5);  hMult_pulse->Draw("hist"); gPad->Modified(); gPad->Update();
            
            c2->cd(4);  f125_fit->Draw("box");    gPad->Modified(); gPad->Update();
            c2->cd(5);  hevt->Draw("colz");       gPad->Modified();   gPad->Update();
            c2->cd(9);  mmg1_f125_fit->Draw("box"); gPad->Modified(); gPad->Update();
            c2->cd(10);  mhevt->Draw("colz");       gPad->Modified(); gPad->Update();
            c2->cd(14);  urw_f125_fit->Draw("box");    gPad->Modified(); gPad->Update();
            c2->cd(15);  uhevt->Draw("colz");       gPad->Modified();   gPad->Update();
            printf("All done, click middle of canvas ...\n");
            #ifdef VERBOSE
              printf(" a0=%f a1=%f (%f deg)  fx1(150)=%f chi2cc_gem=%f  \n",a0,a1,a1/3.1415*180.,fx1.Eval(150.),chi2cc_gem);
            #endif
            if (match || match_mmg1 || match_urw/*NTRACKS>0 && mmg1_NTRACKS>0*/) c2->cd(2); gPad->WaitPrimitive();
        #endif
        if (maxClust_dEdx!=0.) hClusterMaxdEdx_el->Fill(maxClust_dEdx);
        if (totalClust_dEdx!=0.) hClusterTotaldEdx_el->Fill(totalClust_dEdx);
        if (maxClust_m_dEdx!=0.) hmmg1ClusterMaxdEdx_el->Fill(maxClust_m_dEdx);
        if (totalClust_m_dEdx!=0.) hmmg1ClusterTotaldEdx_el->Fill(totalClust_m_dEdx);
        if (maxClust_u_dEdx!=0.) hurwClusterMaxdEdx_el->Fill(maxClust_u_dEdx);
        if (totalClust_u_dEdx!=0.) hurwClusterTotaldEdx_el->Fill(totalClust_u_dEdx);
      #endif   // --- End if USE_CLUST>0 ---
    #endif   //=======================  End Fa125 RAW process Loop  =====================================
    
    //============ END GEMTRD Pattern Recognition Tracking ==================
    
    //=====================================================================================
    //===                Fill Root TTree Hits                                            ===
    //=====================================================================================
    
    #ifdef SAVE_TRACK_HITS
      #if USE_CLUST
      tgem_nclu=nhits;
      tgem_ntracks=NTRACKS;
      for (int n=0; n<nhits; n++) {
        clu_xpos.push_back(hits_Xpos[n]);
        clu_zpos.push_back(hits_Zpos[n]);
        clu_dedx.push_back(hits_dEdx[n]);
        clu_width.push_back(hits_Width[n]);
        if (hits_dEdx[n] > clu_dedx_max) {
          clu_dedx_max=hits_dEdx[n];
          clu_xpos_max=hits_Xpos[n];
          clu_zpos_max=hits_Zpos[n];
          clu_width_max=hits_Width[n];
        }
        //clu_length.push_back(hits_Length[n]);
        //f125_el_clu2d->Fill(hits_Zpos[n],hits_Xpos[n],hits_dEdx[n]);
      }
      
      mmg1_nclu=mmg1_nhits;
      mmg1_ntracks=mmg1_NTRACKS;
      for (int n=0; n<mmg1_nhits; n++) {
        mmg1_clu_xpos.push_back(mmg1_hits_Xpos[n]);
        mmg1_clu_zpos.push_back(mmg1_hits_Zpos[n]);
        mmg1_clu_dedx.push_back(mmg1_hits_dEdx[n]);
        mmg1_clu_width.push_back(mmg1_hits_Width[n]);
        if (mmg1_hits_dEdx[n] > mmg1_clu_dedx_max) {
          mmg1_clu_dedx_max=mmg1_hits_dEdx[n];
          mmg1_clu_xpos_max=mmg1_hits_Xpos[n];
          mmg1_clu_zpos_max=mmg1_hits_Zpos[n];
          mmg1_clu_width_max=mmg1_hits_Width[n];
        }
      }
      
      urw_nclu=urw_nhits;
      urw_ntracks=urw_NTRACKS;
      for (int n=0; n<urw_nhits; n++) {
        urw_clu_xpos.push_back(urw_hits_Xpos[n]);
        urw_clu_zpos.push_back(urw_hits_Zpos[n]);
        urw_clu_dedx.push_back(urw_hits_dEdx[n]);
        urw_clu_width.push_back(urw_hits_Width[n]);
        if (urw_hits_dEdx[n] > urw_clu_dedx_max) {
          urw_clu_dedx_max=urw_hits_dEdx[n];
          urw_clu_xpos_max=urw_hits_Xpos[n];
          urw_clu_zpos_max=urw_hits_Zpos[n];
          urw_clu_width_max=urw_hits_Width[n];
        }
      }
      #endif
      if (tgem_nhit>0) EVENT_VECT_GEM->Fill();
      if (mmg1_nhit>0) EVENT_VECT_MMG1->Fill();
      if (urw_nhit>0 /*|| urw_nyhit>0.*/) EVENT_VECT_URW->Fill();
    #endif
  } // ------------------------ END of event loop  ------------------------------
  
  timer.Stop();
  cout<<"***>>> End Event Loop, Elapsed Time:"<<endl; timer.Print();
  #ifdef WRITE_CSV
    csvFile.close();
  #endif
  cout<<" Total events = "<<(MaxEvt - FirstEvt)<<endl;
  /*
  //---Drift time distribution plot ---
  TH1D *f125_drift = f125_el_amp2d->ProjectionX("f125_drift",100,135);
  TH1D *f125_drift_c = new TH1D(*f125_drift); f125_drift_c->SetStats(0);
  double tgemDriftScale = 1./f125_drift_c->GetEntries();
  f125_drift_c->Scale(tgemDriftScale);
  TH1D *mmg1_f125_drift = mmg1_f125_el_amp2d->ProjectionX("mmg1_f125_drift",95,130);
  TH1D *mmg1_f125_drift_c = new TH1D(*mmg1_f125_drift); mmg1_f125_drift_c->SetStats(0);
  double mmg1DriftScale = 1./mmg1_f125_drift_c->GetEntries();
  mmg1_f125_drift_c->Scale(mmg1DriftScale);
  TH1D *urw_f125_xdrift = urw_f125_x_amp2d->ProjectionX("urw_f125_xdrift",25,50);
  TH1D *urw_f125_xdrift_c = new TH1D(*urw_f125_xdrift); urw_f125_xdrift_c->SetStats(0);
  double urwxDriftScale = 1./urw_f125_xdrift_c->GetEntries();
  urw_f125_xdrift_c->Scale(urwxDriftScale);
  f125_drift_c->SetLineColor(4);  f125_drift_c->SetLineWidth(2);  HistList->Add(f125_drift_c);
  mmg1_f125_drift_c->SetLineColor(2); mmg1_f125_drift_c->SetLineWidth(2);  HistList->Add(mmg1_f125_drift_c);
  urw_f125_xdrift_c->SetLineColor(3); urw_f125_xdrift_c->SetLineWidth(2);  HistList->Add(urw_f125_xdrift_c);
  f125_drift_c->GetXaxis()->SetTitle("Time Response (*8ns)");
  f125_drift_c->GetYaxis()->SetTitle("1. / nEntries");
  f125_drift_c->SetTitle("Drift Time Distribution");
  TLegend *l1 = new TLegend(0.75,0.65,0.9,0.9);
  l1->SetNColumns(2);
  l1->SetTextSize(0.025);
  l1->AddEntry(f125_drift_c, "TRIP-GEM", "l");
  l1->AddEntry(mmg1_f125_drift_c, "MMG1", "l");
  l1->AddEntry(urw_f125_xdrift_c, "uRW X", "l");
  
  TCanvas *c4 = new TCanvas("c4","Drift Time Distribution", 1200, 800);
  c4->cd();
  f125_drift_c->Draw();
  mmg1_f125_drift_c->Draw("same");
  urw_f125_xdrift_c->Draw("same");
  l1->Draw();
  */
  if(TEfficiency::CheckConsistency(*f125_el_tracker_eff,*f125_el_tracker_hits)) {
    p_tgem_eff = new TEfficiency(*f125_el_tracker_eff,*f125_el_tracker_hits);
    p_tgem_eff->SetTitle("Triple GEM-TRD Tracking Efficiency");
  }
  if(TEfficiency::CheckConsistency(*mmg1_f125_el_tracker_eff,*mmg1_f125_el_tracker_hits)) {
    p_mmg1_eff = new TEfficiency(*mmg1_f125_el_tracker_eff,*mmg1_f125_el_tracker_hits);
    p_mmg1_eff->SetTitle("MMG1-TRD Tracking Efficiency");
  }
  if(TEfficiency::CheckConsistency(*urw_f125_x_tracker_eff,*urw_f125_x_tracker_hits)) {
    p_urw_eff = new TEfficiency(*urw_f125_x_tracker_eff,*urw_f125_x_tracker_hits);
    p_urw_eff->SetTitle("uRWELL-TRD Tracking Efficiency");
  }
  
  //=====================================================================================
  //===                 S A V E   H I S T O G R A M S                                ====
  //=====================================================================================
  TFile* fOut;
  #if ANALYZE_MERGED
    char rootFileName[256]; sprintf(rootFileName, "RootOutput/ps26/merged/Run_%06d_%0dEntries_Output.root",RunNum,nEntries);
  #else
    char rootFileName[256]; sprintf(rootFileName, "RootOutput/ps26/Run_%06d_Output.root",RunNum);
  #endif
  fOut = new TFile(rootFileName, "RECREATE");
  fOut->cd();
  cout<<"Writing Output File: "<<rootFileName<<endl;
  HistList->Write("HistDQM", TObject::kSingleKey);
  //c4->Write();
  fOut->Close();
  delete fOut;
  
  //=====================================================================================
  //===                 S A V E   T R A C K   H I T   T T R E E S                    ====
  //=====================================================================================
  #ifdef SAVE_TRACK_HITS
    printf("Writing TTree Hit Info File... \n");
    fHits->cd();
    EVENT_VECT_GEM->Write();
    EVENT_VECT_MMG1->Write();
    EVENT_VECT_URW->Write();
    fHits->Close();
    printf("TTree File Written & Closed OK \n");
  #endif

  //=====================================================================================
  //===                 P L O T     H I S T O G R A M S                               ===
  //=====================================================================================
  #if ANALYZE_MERGED
    const char *OutputDir="RootOutput/ps26/merged";
  #else
    const char *OutputDir="RootOutput/ps26";
  #endif
  #ifdef SAVE_PDF
    char ctit[120];
    #if ANALYZE_MERGED
      sprintf(G_DIR,"%s/Run_%06d_%06dEntries",OutputDir,RunNum,nEntries);
    #else
      sprintf(G_DIR,"%s/Run_%06d",OutputDir,RunNum);
    #endif
    sprintf(ctit,"File=%s",G_DIR);
    bool COMPACT=false;
    TCanvas *cc;
    int nxd=3;
    int nyd=5;
    char pdfname[120];  sprintf(pdfname,"%s_evdisp.pdf",G_DIR);  //c0->Print(pdfname);
    
    auto AddFitText = [](TF1* fitFunc) -> TPaveText* {
      TPaveText* pt = new TPaveText(0.2,0.15,0.6,0.3,"NDC");
      pt->SetFillColor(0);
      pt->SetTextColor(kRed);
      pt->SetTextSize(0.035);
      pt->SetBorderSize(1);
      pt->AddText(Form("Fit constant = %.3f #pm %.4f",fitFunc->GetParameter(0), fitFunc->GetParError(0)));
      pt->AddText(Form("#chi^{2} = %.2f    NDF = %d",fitFunc->GetChisquare(), fitFunc->GetNDF()));
      return pt;
    };
    
    //---------------------  page 1 --------------------
    htitle(" Count ");   // if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hcount->Draw();
    cc=NextPlot(nxd,nyd);  hurw_nyhits->Draw();
    cc=NextPlot(nxd,nyd);  hurw_tmp_nxhits->Draw();
    cc=NextPlot(nxd,nyd);  hurw_nxhits->Draw();
    cc=NextPlot(nxd,nyd);  htgem_tmp_nhits->Draw();
    cc=NextPlot(nxd,nyd);  htgem_nhits->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1_tmp_nhits->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1_nhits->Draw();
    
    htitle(" TRD Correlations  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  hgt1_nhits->Draw();
    cc=NextPlot(nxd,nyd);  hgt2_nhits->Draw();
    cc=NextPlot(nxd,nyd);  tgem_mmg1_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  tgem_mmg1_max_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  urw_mmg1_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  urw_mmg1_max_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  urw_tgem_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  tgem_urw_max_xcorr->Draw("colz");
    
    htitle(" TRD Correlations  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  hgemtrkr_max_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  hgemtrkr_max_ycorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  hgemtrkr_1D_xcorr->Draw("");
    cc=NextPlot(nxd,nyd);  hgemtrkr_1D_ycorr->Draw("");
    
    cc=NextPlot(nxd,nyd);  htgem_trdTrackCorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  hmmg1_trdTrackCorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  hurw_trdTrackCorr->Draw("colz");
    
    //---------------------  page 2a --------------------
    htitle(" TRD Correlations  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  htgem_xy->Draw("colz");
    cc=NextPlot(nxd,nyd);  htgem_max_xy->Draw("colz");
    cc=NextPlot(nxd,nyd);  hmmg1_xy->Draw("colz");
    cc=NextPlot(nxd,nyd);  hmmg1_max_xy->Draw("colz");
    cc=NextPlot(nxd,nyd);  hurw_xy->Draw("colz");
    cc=NextPlot(nxd,nyd);  hurw_max_xy->Draw("colz");
    #if USE_TRD_EXT_TRACK
    //Skip SRS info
    #else
    htitle(" X Correlations  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  tgem_gt1_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  tgem_gt2_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  mmg1_gt1_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  mmg1_gt2_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  urw_gt1_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  urw_gt2_xcorr->Draw("colz");  
  
    //---------------------  page 2a --------------------
    htitle(" SRS GEM-TRKR 1  ");   if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_1_max_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_x->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_x_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_y_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_max_xch->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_max_xamp->Draw();

    htitle(" SRS GEM-TRKR 2  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_2_max_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_x->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_x_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_y_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_max_xch->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_max_xamp->Draw();

    htitle(" SRS GEM / MMG1  ");   if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=3;
    cc=NextPlot(nxd,nyd); mmg1_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hmmg1_peak_y_height->Draw();
    cc=NextPlot(nxd,nyd); tgem_peak_y->Draw();
    cc=NextPlot(nxd,nyd); htgem_peak_y_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_tgem->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_1_mmg1->Draw("colz");
    cc=NextPlot(nxd,nyd); tgem_mmg1_ycorr->Draw("colz");
    
    #endif
    
    htitle(" Plane Diffs  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd); hmmg1_urw_xdiff->Draw();
    cc=NextPlot(nxd,nyd); htgem_urw_xdiff->Draw();
    cc=NextPlot(nxd,nyd); hmmg1_tgem_xdiff->Draw();
    cc=NextPlot(nxd,nyd); hmmg1_tgem_ydiff->Draw();
    cc=NextPlot(nxd,nyd); hTrackDiff->Draw();
    cc=NextPlot(nxd,nyd); htgem_timeDiff->Draw("colz");
    cc=NextPlot(nxd,nyd); hmmg1_timeDiff->Draw("colz");
    cc=NextPlot(nxd,nyd); hurw_timeDiff->Draw("colz");
    
    htitle(" TRD 2D Pulse Differences ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd); htgem_2DPulseMultiplicity->Draw("colz");
    cc=NextPlot(nxd,nyd); htgem_2DPulseVsChan->Draw("colz");
    cc=NextPlot(nxd,nyd); hmmg1_2DPulseMultiplicity->Draw("colz");
    cc=NextPlot(nxd,nyd); hmmg1_2DPulseVsChan->Draw("colz");
    cc=NextPlot(nxd,nyd); hurw_2DPulseMultiplicity->Draw("colz");
    cc=NextPlot(nxd,nyd); hurw_2DPulseVsChan->Draw("colz");
    
    htitle(" fADC125 Raw (Clustering) Track Differences ");   if (!COMPACT) cc=NextPlot(0,0);
    #if (USE_CLUST>0)
    cc=NextPlot(nxd,nyd);  hgemClusterDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1ClusterDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hClusterMaxdEdx_el->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1ClusterMaxdEdx_el->Draw();
    cc=NextPlot(nxd,nyd);  hClusterTotaldEdx_el->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1ClusterTotaldEdx_el->Draw();
    
    //htitle(" fADC125 Raw and Pulse Track Differences ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  hurwClusterDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hurwClusterMaxdEdx_el->Draw();
    cc=NextPlot(nxd,nyd);  hurwClusterTotaldEdx_el->Draw();
    #endif

   //---------------------  page 3 --------------------
    htitle("  TRD (fa125) Amp Distributions ");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  f125_el->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  mmg1_f125_el->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  urw_f125_el_x->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  urw_f125_el_y->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  f125_el_max->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  mmg1_f125_el_max->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  urw_f125_el_xmax->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  urw_f125_el_ymax->Draw();

   //---------------------  page 3a --------------------
    htitle("  GEM-TRD (fa125) Amp 2D");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   f125_el_amp2ds->Draw("colz");
    cc=NextPlot(nxd,nyd);   f125_el_amp2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_amp2ds->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_amp2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_x_amp2ds->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_x_amp2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_y_amp2ds->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_y_amp2d->Draw("colz");

    htitle("  Amplitude vs Channel, 2D");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   f125_xVSamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   f125_xVSamp_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_xVSamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_xVSamp_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_xVSamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_xVSamp_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_yVSamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_yVSamp_max->Draw("colz");

    htitle("  Amplitude vs Time, 2D");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   f125_timeVSamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   f125_timeVSamp_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_timeVSamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_timeVSamp_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_x_timeVSamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_x_timeVSamp_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_y_timeVSamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_y_timeVSamp_max->Draw("colz");

    htitle("  fa125 Max Amp 2D");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   f125_el_amp2d_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_amp2d_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_x_amp2d_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_y_amp2d_max->Draw("colz");
    cc=NextPlot(nxd,nyd);  hgemPulseDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1PulseDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hurwPulseDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hurwPulseDiff_mmg->Draw();
    
    #ifdef GAIN_CALIB
    htitle("  Gain Calib");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   hmmg1GainSum->Draw("");
    cc=NextPlot(nxd,nyd);   hmmg1GainSpread->Draw("");
    cc=NextPlot(nxd,nyd);   hmmg1GainMultiplicity->Draw("colz");
    cc=NextPlot(nxd,nyd);   hmmg12DGain->Draw("colz");
    cc=NextPlot(nxd,nyd);   hmmg12DMaxGain->Draw("colz");
    cc=NextPlot(nxd,nyd);   hmmg12DMaxGainSingle->Draw("colz");
    
    htitle("  Gain Calib");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   hgemGainSum->Draw("");
    cc=NextPlot(nxd,nyd);   hgemGainSpread->Draw("");
    cc=NextPlot(nxd,nyd);   hgemGainMultiplicity->Draw("colz");
    cc=NextPlot(nxd,nyd);   hgem2DGain->Draw("colz");
    cc=NextPlot(nxd,nyd);   hgem2DMaxGain->Draw("colz");
    cc=NextPlot(nxd,nyd);   hgem2DMaxGainSingle->Draw("colz");
    
    htitle("  Gain Calib");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   hurwXGainSum->Draw("");
    cc=NextPlot(nxd,nyd);   hurwXGainSpread->Draw("");
    cc=NextPlot(nxd,nyd);   hurwXGainMultiplicity->Draw("colz");
    cc=NextPlot(nxd,nyd);   hurwX2DGain->Draw("colz");
    cc=NextPlot(nxd,nyd);   hurwX2DMaxGain->Draw("colz");
    cc=NextPlot(nxd,nyd);   hurwX2DMaxGainSingle->Draw("colz");
    
    htitle("  Gain Calib");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   hurwYGainSum->Draw("");
    cc=NextPlot(nxd,nyd);   hurwYGainSpread->Draw("");
    cc=NextPlot(nxd,nyd);   hurwYGainMultiplicity->Draw("colz");
    cc=NextPlot(nxd,nyd);   hurwY2DGain->Draw("colz");
    cc=NextPlot(nxd,nyd);   hurwY2DMaxGain->Draw("colz");
    cc=NextPlot(nxd,nyd);   hurwY2DMaxGainSingle->Draw("colz");
    
    htitle("  Gain Calib");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   hgem2DGainLateTime->Draw("colz");
    cc=NextPlot(nxd,nyd);   hgemGainMaxLate->Draw("");
    cc=NextPlot(nxd,nyd);   hgem2DGainSumLateTime->Draw("colz");
    cc=NextPlot(nxd,nyd);   hgemGainSumLate->Draw("");
    cc=NextPlot(nxd,nyd);   hmmg12DGainLateTime->Draw("colz");
    cc=NextPlot(nxd,nyd);   hmmg1GainMaxLate->Draw("");
    cc=NextPlot(nxd,nyd);   hmmg12DGainSumLateTime->Draw("colz");
    cc=NextPlot(nxd,nyd);   hmmg1GainSumLate->Draw("");
    
    htitle("  Gain Calib");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   hurwX2DGainLateTime->Draw("colz");
    cc=NextPlot(nxd,nyd);   hurwXGainMaxLate->Draw("");
    cc=NextPlot(nxd,nyd);   hurwX2DGainSumLateTime->Draw("colz");
    cc=NextPlot(nxd,nyd);   hurwXGainSumLate->Draw("");
    cc=NextPlot(nxd,nyd);   hurwY2DGainLateTime->Draw("colz");
    cc=NextPlot(nxd,nyd);   hurwYGainMaxLate->Draw("");
    cc=NextPlot(nxd,nyd);   hurwY2DGainSumLateTime->Draw("colz");
    cc=NextPlot(nxd,nyd);   hurwYGainSumLate->Draw("");
    
    #else

    htitle("  External Tracking");    if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=3;
    cc=NextPlot(nxd,nyd); if (p_tgem_eff!=0) {
      p_tgem_eff->Draw("AP");
      TGraphAsymmErrors *tg = p_tgem_eff->CreateGraph();
      tg->Draw("P same");
      #if USE_TRD_EXT_TRACK
      tg->Fit("pol0","Q","",25.,80.);
      #else
      if (RunNum>=8228 && RunNum<=8244) {
        tg->Fit("pol0","Q","",14.,47.);
      } else {
        tg->Fit("pol0","Q","",15.,85.);
      }
      #endif
      TF1* fitFunc=tg->GetFunction("pol0");
      if (fitFunc) {
        fitFunc->SetLineColor(kRed);
        fitFunc->Draw("same");
        AddFitText(fitFunc)->Draw("same");
      }
    }
    cc=NextPlot(nxd,nyd); if (p_mmg1_eff!=0) {
      p_mmg1_eff->Draw("AP");
      TGraphAsymmErrors *mg = p_mmg1_eff->CreateGraph();
      mg->Draw("P same");
      #if USE_TRD_EXT_TRACK
      mg->Fit("pol0","Q","",25.,80.);
      #else
      if (RunNum>=8228 && RunNum<=8244) {
        mg->Fit("pol0","Q","",15.,47.);
      } else {
        mg->Fit("pol0","Q","",20.,85.);
      }
      #endif
      TF1* fitFunc=mg->GetFunction("pol0");
      if (fitFunc) {
        fitFunc->SetLineColor(kRed);
        fitFunc->Draw("same");
        AddFitText(fitFunc)->Draw("same");
      }
    }
    cc=NextPlot(nxd,nyd); if (p_urw_eff!=0) {
      p_urw_eff->Draw("AP");
      TGraphAsymmErrors *ug = p_urw_eff->CreateGraph();
      ug->Draw("P same");
      #if USE_TRD_EXT_TRACK
      ug->Fit("pol0","Q","",25.,80.);
      #else
      if (RunNum>=8228 && RunNum<=8244) {
        ug->Fit("pol0","Q","",9.,45.);
      } else {
        ug->Fit("pol0","Q","",15.,90.);
      }
      #endif
      TF1* fitFunc=ug->GetFunction("pol0");
      if (fitFunc) {
        fitFunc->SetLineColor(kRed);
        fitFunc->Draw("same");
        AddFitText(fitFunc)->Draw("same");
      }
    }

    htitle("  External Tracking");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   f125_el_tracker_hits->Draw();
    cc=NextPlot(nxd,nyd);   f125_el_tracker_eff->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_tracker_hits->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_tracker_eff->Draw();
    cc=NextPlot(nxd,nyd);   urw_f125_x_tracker_hits->Draw();
    cc=NextPlot(nxd,nyd);   urw_f125_x_tracker_eff->Draw();
/*
    htitle("  External Tracking");    if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=3;
    cc=NextPlot(nxd,nyd);   if (f125_el_tracker_hits->GetEntries()>0 && f125_el_tracker_eff->GetEntries()>0 ) {  f125_el_tracker_eff->Divide(f125_el_tracker_hits);  f125_el_tracker_eff->SetMaximum(2.);  f125_el_tracker_eff->Draw();  f125_el_tracker_eff->Fit("pol0","Q","",33.,70.);  TF1* fitFunc=f125_el_tracker_eff->GetFunction("pol0");  fitFunc->SetLineColor(kRed);  fitFunc->Draw("same");  AddFitText(fitFunc)->Draw("same"); }
    cc=NextPlot(nxd,nyd);   if (mmg1_f125_el_tracker_hits->GetEntries()>0 && mmg1_f125_el_tracker_eff->GetEntries()>0 ) {  mmg1_f125_el_tracker_eff->Divide(mmg1_f125_el_tracker_hits);  mmg1_f125_el_tracker_eff->SetMaximum(2.);  mmg1_f125_el_tracker_eff->Draw();  mmg1_f125_el_tracker_eff->Fit("pol0","Q","",15.,54.);  TF1* fitFunc=mmg1_f125_el_tracker_eff->GetFunction("pol0");  fitFunc->SetLineColor(kRed);  fitFunc->Draw("same");  AddFitText(fitFunc)->Draw("same"); }
    cc=NextPlot(nxd,nyd);   if (urw_f125_x_tracker_hits->GetEntries()>0 && urw_f125_x_tracker_eff->GetEntries()>0 ) {  urw_f125_x_tracker_eff->Divide(urw_f125_x_tracker_hits);  urw_f125_x_tracker_eff->SetMaximum(2.);  urw_f125_x_tracker_eff->Draw();  urw_f125_x_tracker_eff->Fit("pol0","Q","",42.,85.);  TF1* fitFunc=urw_f125_x_tracker_eff->GetFunction("pol0");  fitFunc->SetLineColor(kRed);  fitFunc->Draw("same");  AddFitText(fitFunc)->Draw("same"); }
*/
    htitle("  External Tracking");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   tgem_residuals->Draw();
    cc=NextPlot(nxd,nyd);   tgem_residualscorr->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_residuals->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_residualscorr->Draw();
    cc=NextPlot(nxd,nyd);   urw_x_residuals->Draw();
    cc=NextPlot(nxd,nyd);   urw_x_residualscorr->Draw();

    htitle("  External Tracking");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   tgem_residual_ch->Draw("colz");
    cc=NextPlot(nxd,nyd);   tgem_residual_chcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_residual_ch->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_residual_chcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_x_residual_ch->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_x_residual_chcorr->Draw("colz");
    
#endif
   //------------- MAX COMPARISONS ---------------
   
    htitle("  TRD (fa125) Max Amp Comparisons ");    if (!COMPACT) cc=NextPlot(0,0);
    if (f125_el_max_late->GetEntries()>0.) { f125_el_max_late->Scale(1./(f125_el_max_late->GetEntries()));  f125_el_max_late->SetLineColor(2); f125_el_max_late->SetStats(0);
    cc=NextPlot(nxd,nyd);  f125_el_max_late->Draw("same"); }
    if (mmg1_f125_el_max_late->GetEntries()>0.) { mmg1_f125_el_max_late->Scale(1./(mmg1_f125_el_max_late->GetEntries()));  mmg1_f125_el_max_late->SetLineColor(2); mmg1_f125_el_max_late->SetStats(0);
    cc=NextPlot(nxd,nyd);  mmg1_f125_el_max_late->Draw("same"); }
    if (urw_f125_el_xmax_late->GetEntries()>0.) { urw_f125_el_xmax_late->Scale(1./(urw_f125_el_xmax_late->GetEntries()));  urw_f125_el_xmax_late->SetLineColor(2); urw_f125_el_xmax_late->SetStats(0);
    cc=NextPlot(nxd,nyd);  urw_f125_el_xmax_late->Draw("same"); }

    //--- close PDF file ----
    cc=NextPlot(-1,-1);
  #endif
  cout<<"========== END OF RUN "<<RunNum<<" =========="<<endl;
}
