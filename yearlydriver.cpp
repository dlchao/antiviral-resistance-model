/* yearlydriver
 *
 * creates a GlobalModel and outputs yearly proportion of resistant cases.
 *
 * Takes command line options:
 *   -alpha n : sets alpha (fitness of resistant strain) to "n"
 *   -mu n : sets mu (mutation rate in untreated) to "n"
 *   -muT n : sets mu_T (mutation rate in treated) to "n"
 *   -randomseed n : sets the random number generator seed to "n"
 *   -prophy n : sets the number of prophylaxsis to "n" times the number of treateds
 *   -o filename : outputs results to "filename" (default is "temp.csv")
 */

#include <cstdlib>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string.h> 
#include <math.h>
#include "GlobalModel.h"
#include "CityModel.h"

using namespace std;

int main(int argc, char *argv[]) {
  //  gsl_rng * rng = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus2);
  struct GlobalModelParams params;
  params.R0min = 0.8;
  params.R0max = 1.1;
  params.mutation = 0.0;
  params.mutationTreatment = 0.0;
  params.mutation = 0.0;
  params.mutationTreatment = 0.0;
  params.nStartDay=150;
  params.fCrossStrain = 1.0;
  params.fDrift = 0.35;
  params.fProphylaxis = 0.0;
  params.fInfectedTravel = 1.0; //0.75;
  params.globalgamma = 0.0;
  params.alpha=1.0;
  params.fResistanceCost = 0.0;
  params.AVES = 0.63;
  params.AVEI = 0.15;
  params.bResistance = true;

  int randomseed = 5489;
  double muT = 0.01;
  double mu = 0.000001;
  double phi = 1.0;
  double prophy = 0.0;
  double fTreat = -1.0; // default fraction of cases treated
  double fReassort = 0.0;
  string filename = "temp.csv";
  int runlength = 30; // length of run in years

  if (argc>1) {
    for (int i=1; i<argc; i++) {
      char **end = NULL;
      if (strcmp(argv[i], "-alpha")==0) {
	params.alpha = strtof(argv[i+1],end);
	cerr << "alpha = " << params.alpha << endl;
	i++;
      } else if (strcmp(argv[i], "-c")==0) {
	params.fResistanceCost = strtof(argv[i+1],end);
	cerr << "c = " << params.fResistanceCost << endl;
	i++;
      } else if (strcmp(argv[i], "-R0")==0 || strcmp(argv[i], "-r0")==0) {
	params.R0max = strtof(argv[i+1],end);
	cerr << "R0max = " << params.R0max << endl;
	i++;
      } else if (strcmp(argv[i], "-strains")==0) {
	nMaxStrains = strtof(argv[i+1],end);
	cerr << "strains = " << nMaxStrains << endl;
	i++;
      } else if (strcmp(argv[i], "-drift")==0) {
	params.fDrift = strtof(argv[i+1],end);
	cerr << "drift = " << params.fDrift << endl;
	i++;
      } else if (strcmp(argv[i], "-cross")==0) {
	params.fCrossStrain = strtof(argv[i+1],end);
	cerr << "cross = " << params.fCrossStrain << endl;
	i++;
      } else if (strcmp(argv[i], "-muT")==0) {
	muT = strtof(argv[i+1],end);
	cerr << "muT = " << muT << endl;
	i++;
      } else if (strcmp(argv[i], "-mu")==0) {
	mu = strtof(argv[i+1],end);
	cerr << "mu = " << mu << endl;
	i++;
      } else if (strcmp(argv[i], "-prophy")==0) {
	prophy = strtof(argv[i+1],end);
	cerr << "prophy = " << prophy << endl;
	i++;
      } else if (strcmp(argv[i], "-phi")==0) {
	phi = strtof(argv[i+1],end);
	cerr << "phi = " << phi << endl;
	i++;
      } else if (strcmp(argv[i], "-reassort")==0) {
	fReassort = strtof(argv[i+1],end);
	cerr << "reassort = " << fReassort << endl;
	i++;
      } else if (strcmp(argv[i], "-AVES")==0) {
	params.AVES = strtof(argv[i+1],end);
	cerr << "AVES = " << params.AVES << endl;
	i++;
      } else if (strcmp(argv[i], "-AVEI")==0) {
	params.AVEI = strtof(argv[i+1],end);
	cerr << "AVEI = " << params.AVEI << endl;
	i++;
      } else if (strcmp(argv[i], "-randomseed")==0) {
	randomseed = strtol(argv[i+1],end,10);
	cerr << "seed = " << randomseed << endl;
	i++;
      } else if (strcmp(argv[i], "-runlength")==0) {
	runlength = strtol(argv[i+1],end,10);
	cerr << "runlength = " << runlength << endl;
	i++;
      } else if (strcmp(argv[i], "-notreat")==0) {
	fTreat = 0.0;
	cerr << "no treatment" << endl;
      } else if (strcmp(argv[i], "-noresistance")==0) {
	params.bResistance = false;
	cerr << "no resistance" << endl;
      } else if (strcmp(argv[i], "-treat")==0) {
	fTreat = strtof(argv[i+1],end);
	cerr << "treatment = " << fTreat << endl;
	cerr << "gamma = " << (1.0-exp(log(1-fTreat)/6.0)) << endl;
	i++;
      } else if (strcmp(argv[i], "-o")==0) {
	filename = argv[i+1];
	cerr << "file = " << filename << endl;
	i++;
      } else {
	cerr << "Unknown option: " << argv[i] << endl;
      }
    }
  }

  ofstream outfile;
  outfile.open(filename.c_str());
  if(outfile.fail()) {
    cerr << "Could not create " << filename << endl;
    return -1;
  }

  gsl_rng_set(rng, randomseed);
  GlobalModel g(params);

  // create 321 cities based on data files
  if (!g.readCityFile("DataFiles/population_321_age.txt")) {
    cerr << "Could not find DataFiles/population_321_age.txt" << endl;
    return -1;
  }
  if (!g.readTransportFile("DataFiles/travel_321.txt")) {
    cerr << "Could not find DataFiles/travel_321.txt" << endl;
    return -1;
  }
  if (!g.readSeasonalityFile("DataFiles/seasonality_321.csv")) {
    cerr << "Could not find DataFiles/seasonality_321.csv" << endl;
    return -1;
  }
  g.city[114].addInfected(1000,0);

  unsigned int nSumCases = 0;
  unsigned int nSumWildtypeCases = 0;
  unsigned int nSumResistant = 0;
  unsigned int nSumResistant0 = 0;
  
  for (int t=params.nStartDay; t<runlength*365; t++) {
    if (g.getDay()==3*365) {
      g.setMutationTreatment(muT); // mutation in treated infecteds
      g.setMutation(mu); // mutation in untreated infecteds
      g.setProphylaxis(prophy);
      g.setPhi(phi);
      g.setReassortFactor(fReassort);
      if (fTreat>0.0) {
	g.setGlobalGamma(1.0-exp(log(1-fTreat)/(nMaxDaysInfected-nTreatmentDelay)));
      } else if (fTreat<0.0)
	for (unsigned int i=0; i<g.city.size(); i++) {
	  CityModel &source = g.city[i];
	  if (source.szNation.compare("United_States")==0) {
	    source.setCityGamma(1.0-exp(log(1-0.057)/(nMaxDaysInfected-nTreatmentDelay)));
	  } else if (source.szNation.compare("Japan")==0) {
	    source.setCityGamma(1.0-exp(log(1-0.375)/(nMaxDaysInfected-nTreatmentDelay)));
	  } else if (source.szNation.compare("Austria")==0) { // 2005 data
	    source.setCityGamma(1.0-exp(log(1-0.0032)/(nMaxDaysInfected-nTreatmentDelay)));
	  } else if (source.szNation.compare("Belgium")==0) {
	    source.setCityGamma(1.0-exp(log(1-0.0104)/(nMaxDaysInfected-nTreatmentDelay)));
	  } else if (source.szNation.compare("Germany")==0) {
	    source.setCityGamma(1.0-exp(log(1-0.0129)/(nMaxDaysInfected-nTreatmentDelay)));
	  } else if (source.szNation.compare("Denmark")==0) {
	    source.setCityGamma(1.0-exp(log(1-0.0043)/(nMaxDaysInfected-nTreatmentDelay)));
	  } else if (source.szNation.compare("Greece")==0) {
	    source.setCityGamma(1.0-exp(log(1-0.0024)/(nMaxDaysInfected-nTreatmentDelay)));
	  } else if (source.szNation.compare("Finland")==0) {
	    source.setCityGamma(1.0-exp(log(1-0.0082)/(nMaxDaysInfected-nTreatmentDelay)));
	  } else if (source.szNation.compare("France")==0) {
	    source.setCityGamma(1.0-exp(log(1-0.0096)/(nMaxDaysInfected-nTreatmentDelay)));
	  } else if (source.szNation.compare("Norway")==0) {
	    source.setCityGamma(1.0-exp(log(1-0.0074)/(nMaxDaysInfected-nTreatmentDelay)));
	  }
	}
    }
    g.tick(rng);

    for (unsigned int i=0; i<g.city.size(); i++) 
      for (int treat=0; treat<2; treat++) {
	nSumWildtypeCases += g.city[i].nNumInfected[0][0][treat][0];
	nSumResistant0 += g.city[i].nNumInfected[0][1][treat][0];
	for (int strain=0; strain<nMaxStrains; strain++) {
	  nSumResistant += g.city[i].nNumInfected[strain][1][treat][0];
	  nSumCases += g.city[i].nNumInfected[strain][0][treat][0] +
	    g.city[i].nNumInfected[strain][1][treat][0];;
	}
      }
    
    // output
    if (g.getDay()%365==0) {
      outfile << g.getAlpha() << "," << g.getResistanceCost() << "," << g.getMutationTreatment() << ","  << g.getMutation() << ","  << g.getProphylaxis() << "," << g.getPhi() << "," << g.getReassortFactor() << ","  << randomseed << "," << g.getDay()/365 << "," << nSumCases << "," << nSumWildtypeCases << "," << nSumResistant0 << "," << nSumResistant << endl;
      if (nSumResistant==nSumCases && nSumCases>100000)
	break;
      nSumCases = 0;
      nSumWildtypeCases = 0;
      nSumResistant = 0;
      nSumResistant0 = 0;
    }
  }
  outfile.close();
  gsl_rng_free(rng);
  return 0;
}
