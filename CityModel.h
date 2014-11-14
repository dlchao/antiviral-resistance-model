/* class CityModel
 *
 * Dennis Chao
 * 7/2010
 */

#ifndef __CITYMODEL_H
#define __CITYMODEL_H

#include <vector>
using namespace std;
class GlobalModel;

const int nMaxDaysInfected=6; // length of infection
const int nTreatmentDelay=2; // how many days before treatment can possibly start
const int nMonthsImmune=3;    // minimum number of months that an infected individual is totally immune to new infection
const double fViralLoad[nMaxDaysInfected] = {0.2241997, 0.9722359, 1.0000000, 0.8530784, 0.6263502, 0.1214552};

extern int nMaxStrains; // number of possible strains (set by user)
const int NMAXSTRAINS=2; // number of possible strains (hard limit)
const int nMaxMonth=12*31+1; // maximum number of simulation months
const int nDaysInYear=365; // length of year
const int nMonthIndex[nDaysInYear] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // January
  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,       // February
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, // March
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,   // April
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, // May
  5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,   // June
  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6, // July
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, // August
  8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,   // September
  9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, // October
  10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,   // November
  11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11}; // December

class Recovered {
 public:
  Recovered() {
    nNum = 0;
    nLastInfectionMonth = -1;
    for (int i=0; i<nMaxStrains; i++)
      nRecoveryMonth[i] = -1;
    memset(nRecoverQueue, 0, nMaxDaysInfected*sizeof(int));
  }
  int nNum;   // number of people
  int nRecoveryMonth[NMAXSTRAINS]; // month in which individuals recovered from strain n nMaxStrains
  int nRecoverQueue[nMaxDaysInfected]; // circular queue of infected people about to enter this recovered group
  int nLastInfectionMonth; // month of the last infection
};

class CityModel {
 public:
  CityModel(GlobalModel *g, int nNumSusceptibles);
  ~CityModel() {}
  void tick(gsl_rng *rng);
  void addInfected(int wt, int r); 
  int getNumInfectedTotal() {
    int total = 0;
    for (int strain=0; strain<nMaxStrains; strain++)
      for (int res=0; res<2; res++)
	for (int treat=0; treat<2; treat++)
	  for (int d=nMaxDaysInfected-1; d>0; d--)
	    total += nNumInfected[strain][res][treat][d];
    return total;
  }
  int getNumInfected(int strain, int resistant) {
    int total = 0;
    for (int treat=0; treat<2; treat++)
      for (int d=nMaxDaysInfected-1; d>0; d--)
	total += nNumInfected[strain][resistant][treat][d];
    return total;
  }
  void setCityGamma(double x) { fCityGamma=x; }

  int nID;            // unique ID of city
  string szName;      // name of city
  string szNation;    // name of country
  double fLatitude;   // latitude of city
  double fCityGamma;  // gamma for this city (daily probability of treatment of cases)
  unsigned int nNumSusceptiblesStart; // initial number of susceptible people (possible proxy for city size)
  int nNumInfected[NMAXSTRAINS][2][2][nMaxDaysInfected]; // number infected by strain, antiviral resistance status (0/1), untreated/treated (0/1), and number of days infected
  vector<Recovered *> vRecovered; // susceptible and recovered individuals
  double p00,p01,p10,p11; // transmission probabilities based on treatment of susceptible or infectious person
  GlobalModel *gm;
  double fSeasonality[365]; // relative transmissibility because of seasonality (0..1)
  static int nNextID; // the next city ID to assign
};
#endif
