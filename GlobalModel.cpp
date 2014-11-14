/* class GlobalModel implementation
 *
 * Links a network of CityModels into a global network of cities
 */

#include <assert.h>
#include <string.h> 
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "GlobalModel.h"
#include "CityModel.h"

const int nVersionMajor = 1;
const int nVersionMinor = 0;

// GlobalModel::tick
//   returns number of people infected during this time step
void GlobalModel::tick(gsl_rng *rng) {
  // infection within city
  for (unsigned int i=0; i<city.size(); i++)
    city[i].tick(rng);

  // transport
  if (fTransport)
    for (unsigned int i=0; i<city.size(); i++) {
      CityModel &source = city[i];
      for (int strain=0; strain<nMaxStrains; strain++) {
	for (int res=0; res<2; res++) {
	  int nNumInfected = 0;
	  // assume infected people under treatment do not travel
	  for (int day=0; day<nMaxDaysInfected; day++) 
	    nNumInfected += source.nNumInfected[strain][res][0][day];
	  if (nNumInfected>0) {
	    double travelprobs[1000];
	    unsigned int travelresult[1000];
	    assert(city.size()<1000);
	    int popsize = source.nNumSusceptiblesStart;
	    double sum=0.0;
	    for (unsigned int j=0; j<city.size(); j++)
	      if (j!=i) {
		travelprobs[j] = fTransport[city.size()*i+j]/popsize*fInfectedTravel;
		sum += travelprobs[j];
	      }
	    travelprobs[i] = 1.0-sum;
	    gsl_ran_multinomial(rng, city.size(), nNumInfected,
				travelprobs, travelresult);
	    for (unsigned int j=0; j<city.size(); j++) {
	      if (i!=j && travelresult[j]>0) {
		CityModel &dest = city[j];
		unsigned totalinfected = nNumInfected;
		assert(travelresult[j]<=totalinfected);
		for (int day=0; day<nMaxDaysInfected; day++)
		  if (source.nNumInfected[strain][res][0][day]>0 && travelresult[j]>0) {
		    totalinfected -= source.nNumInfected[strain][res][0][day];
		    unsigned int result = gsl_ran_hypergeometric(rng, source.nNumInfected[strain][res][0][day], totalinfected, travelresult[j]);
		//		cerr << k << "?" <<  (source.nNumInfectedWildtypeUntreated[k]) << "," << (totalinfected) << "," << travelresult[j] << endl;
		    source.nNumInfected[strain][res][0][day]-=result;
		    dest.nNumInfected[strain][res][0][day]+=result;
		    nNumInfected-=result;
		    travelresult[j]-=result;
		  }
	      }
	    }
	  }
	}
      }
    }

  // increment time
  nDay++;
  if (nDay%nDaysInYear==0)
    nYear++;
}

int GlobalModel::readCityFile(const char *cityfilename) {
  ostringstream oss;
  if (cityfilename)
    oss.str(cityfilename);
  else
    return -1;
  ifstream iss(oss.str().c_str());
  if (!iss) {
    cerr << "ERROR: " << cityfilename << " not found." << endl;
    return -1;
  }
  string line;
  getline(iss, line); // throw away header line
  int nNumLinesRead = 0;
  while (getline(iss, line)) {
    istringstream iss;
    iss.str(line);
    string name;
    iss >> name;
    string temp;
    iss >> temp; // number
    iss >> temp; // zone
    int pop;
    iss >> pop;
    iss >> temp; // longitude
    double latitude;
    iss >> latitude;
    iss >> temp; // country
    CityModel *c = new CityModel(this, pop); //round(pop*fProphylaxis));
    c->szName = name;
    c->fLatitude = latitude;
    c->szNation = temp;
    city.push_back(*c);
    nNumLinesRead++;
  }
  return nNumLinesRead;
}

bool GlobalModel::readTransportFile(const char *transportfilename) {
  ostringstream oss;
  if (transportfilename)
    oss.str(transportfilename);
  else
    return false;
  ifstream iss(oss.str().c_str());
  if (!iss) {
    cerr << "ERROR: " << transportfilename << " not found." << endl;
    return false;
  }
  if (fTransport)
    delete fTransport;
  fTransport = new double[city.size()*city.size()];
  string line;
  getline(iss, line); // throw away header line
  unsigned int nNumLinesRead = 0;
  while (getline(iss, line)) {
    istringstream iss;
    iss.str(line);
    string temp;
    iss >> temp; // city name
    iss >> temp; // number
    for (unsigned int i=0; i<city.size(); i++)
      iss >> fTransport[city.size()*nNumLinesRead+i];
    nNumLinesRead++;
  }
  if (nNumLinesRead!=city.size()) {
    cerr << "Travel file: Read " << nNumLinesRead << ", wanted " << city.size() << endl;
    return false;
  }
  return true;
}

bool GlobalModel::readSeasonalityFile(const char *seasonalityfilename) {
  ostringstream oss;
  if (seasonalityfilename)
    oss.str(seasonalityfilename);
  else
    return false;
  ifstream iss(oss.str().c_str());
  if (!iss) {
    cerr << "ERROR: " << seasonalityfilename << " not found." << endl;
    return false;
  }
  string line;
  unsigned int nNumLinesRead = 0;
  while (getline(iss, line)) {
    istringstream iss;
    iss.str(line);
    char *p = (char *)line.c_str();
    while (*(p++) != ',')
      ;
    for (unsigned int i=0; i<365; i++) { 
      char *next;
      city[nNumLinesRead].fSeasonality[i] =  strtod ( p, &next);
      p = next+1;
    }
    nNumLinesRead++;
  }
  if (nNumLinesRead!=city.size()) {
    cerr << "Seasonality: Read " << nNumLinesRead << ", wanted " << city.size() << endl;
    return false;
  }
  return true;
}
