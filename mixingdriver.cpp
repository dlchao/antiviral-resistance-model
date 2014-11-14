/* mixingdriver
 *
 * creates a "well-mixed" GlobalModel and outputs daily number of infecteds.
 * each day, infected people are teleported randomly around the world.
 */

#include <cstdlib>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string.h> 
#include "GlobalModel.h"
#include "CityModel.h"

using namespace std;

int main(int argc, char *argv[]) {
  struct GlobalModelParams params;
  params.R0min = 0.8;
  params.R0max = 1.6;
  params.mutation = 0.0;
  params.mutationTreatment = 0.0;
  params.nStartDay=150;
  params.fDrift = 0.9;
  params.fProphylaxis = 0.0;
  params.fInfectedTravel = 0.75;
  params.globalgamma = 0.0;
  params.alpha=1.0;
  params.R0max = 1.3;
  int randomseed = 5489;
  double muT = 0.001;
  double mu = 0.00001;

  if (argc>1) {
    for (int i=1; i<argc; i++) {
      char **end=NULL;
      if (strcmp(argv[i], "-alpha")==0) {
	params.alpha = strtof(argv[i+1],end);
	cerr << "alpha = " << params.alpha << endl;
	i++;
      } else if (strcmp(argv[i], "-muT")==0) {
	muT = strtof(argv[i+1],end);
	cerr << "muT = " << muT << endl;
	i++;
      } else if (strcmp(argv[i], "-mu")==0) {
	mu = strtof(argv[i+1],end);
	cerr << "mu = " << mu << endl;
	i++;
      } else if (strcmp(argv[i], "-R0")==0) {
	params.R0max = strtof(argv[i+1],end);
	cerr << "R0 (max) = " << params.R0max << endl;
	i++;
      } else if (strcmp(argv[i], "-randomseed")==0) {
	randomseed = strtol(argv[i+1],end,10);
	cerr << "seed = " << randomseed << endl;
	i++;
      }
    }
  }

  gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus2);
  gsl_rng_set(rng, randomseed);
  GlobalModel g(params);
  if (true) {
    // create 321 cities based on data files
    if (!g.readCityFile("DataFiles/population_321_age.txt")) {
      cerr << "Could not find DataFiles/population_321_age.txt" << endl;
      return -1;
    }
    /*    if (!g.readTransportFile("DataFiles/travel_321.txt")) {
      cerr << "Could not find DataFiles/travel_321.txt" << endl;
      return -1;
      } */ // no transportation!
    if (!g.readSeasonalityFile("DataFiles/seasonality_321.csv")) {
      cerr << "Could not find DataFiles/seasonality_321.csv" << endl;
      return -1;
    }
    g.city[114].addInfected(1000,0);
  } else if (false) {
    g.readCityFile("smallpop.txt");
    g.readTransportFile("smalltravel.txt");
    g.readSeasonalityFile("smallseason.txt");
    g.city[0].addInfected(100,0);
  } else {
    // create a single city with a small infected population
    CityModel *c = new CityModel(&g, 10000);
    c->addInfected(10,0);
    g.city.push_back(*c);
  }

  double travelprobs[1000];
  assert(g.city.size()<1000);
  unsigned int total = 0;
  for (unsigned int j=0; j<g.city.size(); j++)
    total += g.city[j].nNumSusceptiblesStart;
  for (unsigned int j=0; j<g.city.size(); j++)
    travelprobs[j] = g.city[j].nNumSusceptiblesStart/(double)total;

  // main loop
  for (int t=0; t<20*365; t++) {
    if (g.getDay()==4*365) {
      g.setMutationTreatment(muT); // mutation in treated infecteds
      g.setMutation(mu); // mutation in untreated infecteds
      //      g.setGlobalGamma(0.008512445); // 5% receive treatment
      //      g.setPhi(0.25);
      for (unsigned int i=0; i<g.city.size(); i++) {
	CityModel &source = g.city[i];
	if (source.szNation.compare("United_States")==0) {
	  source.setCityGamma(0.003173215);
	} else if (source.szNation.compare("Japan")==0) {
	  source.setCityGamma(0.02721289);
	} else if (source.szNation.compare("Austria")==0) {
	  source.setCityGamma(0.0001790275);
	} else if (source.szNation.compare("Belgium")==0) {
	  source.setCityGamma(0.0005794352);
	} else if (source.szNation.compare("Germany")==0) {
	  source.setCityGamma(0.0007200669);
	} else if (source.szNation.compare("Denmark")==0) {
	  source.setCityGamma(0.0002417250);
	} else if (source.szNation.compare("Greece")==0) {
	  source.setCityGamma(0.0001342556);
	} else if (source.szNation.compare("Finland")==0) {
	  source.setCityGamma(0.0004568372);
	} else if (source.szNation.compare("France")==0) {
	  source.setCityGamma(0.0005345736);
	} else if (source.szNation.compare("Norway")==0) {
	  source.setCityGamma(0.0004120031);
	}
      }
    }
    g.tick(rng);

    // count global infecteds
    int nNumInfected[nMaxStrains][2][2][nMaxDaysInfected]; 
    memset(nNumInfected, 0, sizeof(nNumInfected));
    for (unsigned int citynum=0; citynum<g.city.size(); citynum++) {
      CityModel &source = g.city[citynum];
      for (int strain=0; strain<nMaxStrains; strain++)
	for (int res=0; res<2; res++)
	  for (int treat=0; treat<2; treat++)
	    for (int day=0; day<nMaxDaysInfected; day++)
	      nNumInfected[strain][res][treat][day] += source.nNumInfected[strain][res][treat][day];
    }

    // distribute infecteds randomly around globe
    for (int strain=0; strain<nMaxStrains; strain++)
      for (int res=0; res<2; res++)
	for (int treat=0; treat<2; treat++)
	  for (int day=0; day<nMaxDaysInfected; day++) 
	    if (nNumInfected[strain][res][treat][day]>0) {
	      unsigned int travel[1000];
	      gsl_ran_multinomial(rng, g.city.size(), nNumInfected[strain][res][treat][day],
				  travelprobs, travel);
	      for (unsigned int citynum=0; citynum<g.city.size(); citynum++)
		g.city[citynum].nNumInfected[strain][res][treat][day] = travel[citynum];
	    }
    
    // daily output to stdout
    cout << g.getDay() << ",T";
    for (unsigned int i=0; i<g.city.size(); i++) {
      CityModel &source = g.city[i];
      cout << "," << source.getNumInfectedTotal();
    }
    cout << endl;
    cout << g.getDay() << ",R";
    for (unsigned int i=0; i<g.city.size(); i++) {
      CityModel &source = g.city[i];
      cout << "," << (source.getNumInfectedResistantUntreated()+source.getNumInfectedResistantTreated());
    }
    cout << endl;
  }
  gsl_rng_free(rng);
  return 0;
}
