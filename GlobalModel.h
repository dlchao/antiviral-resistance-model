/* class GlobalModel
 *
 * Dennis Chao
 * 6/2010
 */

#ifndef __GLOBALMODEL_H
#define __GLOBALMODEL_H

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include "CityModel.h"

class GlobalModelParams {
 public:
  GlobalModelParams() {
    mutation = mutationTreatment = globalgamma = 0.0;
    R0min = R0max = 1.6;
    alpha = 1.0;
    phi = 1.0;
    AVES = 0.3;
    AVEI = 0.5;
    nStartDay = 0;
    fDrift = 0.1;
    fCrossStrain = 0.6;
    fProphylaxis = 0.0;
    fInfectedTravel=1.0;
    fResistanceCost = 0.0;
    bResistance = true;
    fReassortFactor = 0.0;
  }
  ~GlobalModelParams() {}
  double R0max;        // transmissibility
  double R0min;        // transmissibility
  double mutation;     // daily mutation probability
  double mutationTreatment; // daily mutation probability when under treatment
  bool bResistance;
  double alpha;        // transmission of new strains wrt strain 0
  double phi;          // transmission of resistant strain wrt wildtype in treated
  double fResistanceCost; // transmission cost of resistance
  double globalgamma;  // daily treatment probability
  double AVES;         // antiviral efficacy for susceptibility
  double AVEI;         // antiviral efficacy for infectiousness
  double fDrift;       // how much immunity loss per year
  double fCrossStrain; // how much VE_S-like protection between strains
  double fProphylaxis; // fraction of susceptibles with antiviral prophylaxis
  double fInfectedTravel; // probability multiplier for travel for infected people (set to 1.0 for no difference)
  int nStartDay;
  double fReassortFactor; // multiplier for reassortment probability
};

class GlobalModel {
 public:
  GlobalModel() {
    nDay = 0;
    nYear = 0;
    R0min = R0max = 1.6;
    mutation = 0.0;
    mutationTreatment = 0.0;
    bResistance = true;
    alpha = 1.0;
    phi = 1.0;
    globalgamma = 0.0;
    AVES = 0.3;
    AVEI = 0.5;
    fDrift = 0.05;
    fProphylaxis = 0.0;
    fInfectedTravel=1.0;
    fTransport = NULL;
    fReassortFactor = 0.0;
  }
  GlobalModel(const GlobalModelParams &params) {
    nYear = 0;
    R0min = params.R0min;
    R0max = params.R0max;
    nDay = params.nStartDay;
    mutation = params.mutation;
    mutationTreatment = params.mutationTreatment;
    bResistance = params.bResistance;
    alpha = params.alpha;
    fResistanceCost = params.fResistanceCost;
    phi = params.phi;
    globalgamma = params.globalgamma;
    AVES = params.AVES;
    AVEI = params.AVEI;
    fCrossStrain = params.fCrossStrain;
    fDrift = params.fDrift;
    fProphylaxis = params.fProphylaxis;
    fInfectedTravel=params.fInfectedTravel;
    fReassortFactor = params.fReassortFactor;
    double sum = 0.0;
    for (int i=0; i<nMaxDaysInfected; i++)
      sum+=fViralLoad[i];
    p = R0max/sum;
    fTransport = NULL;
  }
  ~GlobalModel() {
    if (fTransport)
      delete fTransport;
  }
  inline unsigned int getDay() { return nDay; }
  inline unsigned int getYear() { return nYear; }
  inline double getp() { return p; }
  inline double getR0min() { return R0min; }
  inline double getR0max() { return R0max; }
  inline double getMutation() { return mutation; }
  inline double getMutationTreatment() { return mutationTreatment; }
  inline void setMutation(double x) { mutation=x; }
  inline void setMutationTreatment(double x) { mutationTreatment=x; }
  inline double getAlpha() { return alpha; }
  inline bool useResistance() { return bResistance; }
  inline double getResistanceCost() { return fResistanceCost; }
  inline void setReassortFactor(double x) { fReassortFactor = x; }
  inline double getReassortFactor() { return fReassortFactor; }
  inline double getPhi() { return phi; }
  inline void setPhi(double x) { phi=x; }
  inline double getGlobalGamma() { return globalgamma; }
  inline void setGlobalGamma(double x) { globalgamma=x; }
  inline double getAVES() { return AVES; }
  inline double getAVEI() { return AVEI; }
  inline double getDrift() { return fDrift; }
  inline double getCrossStrain() { return fCrossStrain; }
  inline double getProphylaxis() { return fProphylaxis; }
  inline void setProphylaxis(double x) { fProphylaxis = x; }
  void tick(gsl_rng *rng);
  int readCityFile(const char *cityfilename);
  bool readTransportFile(const char *transportfilename);
  bool readSeasonalityFile(const char *seasonalityfilename);
  vector<CityModel> city;

 protected:
  double *fTransport;
  double R0min;        // out-of-season transmissibility
  double R0max;        // in-season transmissibility
  double p;            // transmissibility
  double mutation;     // daily mutation probability
  double mutationTreatment; // daily mutation probability when under treatment
  bool bResistance;
  double alpha;        // transmission of new strains wrt strain 0
  double fResistanceCost; // transmission cost of resistance
  double phi;          // transmission of resistant strain wrt wildtype in treated
  double globalgamma;  // daily treatment probability
  double AVES;         // antiviral efficacy for susceptibility
  double AVEI;         // antiviral efficacy for infectiousness
  double fDrift;       // how much immunity loss per year
  double fCrossStrain; // how much VE_S-like protection between strains
  double fProphylaxis; // number of susceptibles to prophylax per infected
  double fInfectedTravel; // probability multiplier for travel for infected people (set to 1.0 for no difference)
  unsigned int nDay;   // simulation day (0,365,730...=January 1)
  unsigned int nYear;  // simulation year (starts at 0)
  double fReassortFactor; // multiplier for reassortment probability
};
#endif
