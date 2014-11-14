/* class CityModel
 *
 * A stochastic influenza epidemic simulation of a single city.
 * Dennis Chao
 * 7/2010
 */

#include <string>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "CityModel.h"
#include "GlobalModel.h"

using namespace std;

int CityModel::nNextID=0;
int nMaxStrains=2; // number of possible strains

CityModel::CityModel(GlobalModel *g, int nNumSusceptibles) {
  gm = g;
  nNumSusceptiblesStart = nNumSusceptibles;
  p00 = (*gm).getp()/nNumSusceptibles;
  p10 = p00 * (1.0-(*gm).getAVES());
  p01 = p00 * (1.0-(*gm).getAVEI());
  p11 = p00 * (1.0-(*gm).getAVES()) * (1.0-(*gm).getAVEI());
  nID = nNextID++;
  for (int i=0; i<365; i++)
    fSeasonality[i] = 1.0;

  fCityGamma = -1.0; // not defined

  memset(nNumInfected, 0, sizeof(int)*nMaxStrains*2*2*nMaxDaysInfected);
  Recovered *r = new Recovered();
  r->nNum = nNumSusceptibles;
  vRecovered.push_back(r);
}

void CityModel::addInfected(int wt, int r) {
  nNumInfected[0][0][0][0] = wt;
  nNumInfected[0][1][0][0] = r;
}

void CityModel::tick(gsl_rng *rng) {
  int nDayOfYear = (*gm).getDay()%365;
  int nMonth = (*gm).getYear()*12 + nMonthIndex[nDayOfYear] + 1;
  int nRecoverQueueDay = (*gm).getDay() % nMaxDaysInfected;
  int nRecoverQueueYesterday = (nRecoverQueueDay+nMaxDaysInfected-1) % nMaxDaysInfected;
  assert(nMonth<nMaxMonth);

  // update infected pop clock
  int nOldNumInfected[NMAXSTRAINS][2][2][nMaxDaysInfected];
  //  memcpy(nOldNumInfected, nNumInfected, sizeof(nOldNumInfected));
  for (int strain=0; strain<nMaxStrains; strain++)
    for (int res=0; res<2; res++)
      for (int treat=0; treat<2; treat++)
	for (int d=nMaxDaysInfected-1; d>=0; d--)
	  nOldNumInfected[strain][res][treat][d] = nNumInfected[strain][res][treat][d];

  for (vector<Recovered *>::iterator it = vRecovered.begin(); it!=vRecovered.end(); ++it) {
    Recovered *susceptible = (*it);
    susceptible->nNum += susceptible->nRecoverQueue[nRecoverQueueDay];
    susceptible->nRecoverQueue[nRecoverQueueDay]=0;
  }

  // consolidate recovered classes
  {
    Recovered *susceptible0 = vRecovered.front();
    vector<Recovered *>::iterator it = vRecovered.begin(); 
    ++it;
    while (it!=vRecovered.end()) {
      Recovered *susceptible = *it;
      if (nMonth-susceptible->nLastInfectionMonth > nMonthsImmune) {
	susceptible0->nNum += susceptible->nNum;
	for (int i=0; i<nMaxDaysInfected; i++)
	  susceptible0->nRecoverQueue[i]+=susceptible->nRecoverQueue[i];
	vRecovered.erase(it);
	it = vRecovered.begin(); // start over
      }
      ++it;
    }
  }

  // advance infected people by one day
  for (int strain=0; strain<nMaxStrains; strain++)
    for (int res=0; res<2; res++)
      for (int treat=0; treat<2; treat++) {
	for (int d=nMaxDaysInfected-1; d>0; d--)
	  nNumInfected[strain][res][treat][d] = nNumInfected[strain][res][treat][d-1];
	nNumInfected[strain][res][treat][0] = 0;
      }

  // what proportion of susceptibles are prophylaxed
  double fProphy = (*gm).getProphylaxis();
  if (fProphy>0.0) {
    int totaltreated = 0;
    int totalsusceptibles = 0;
    for (int strain=0; strain<nMaxStrains; strain++)
      for (int res=0; res<2; res++)
	for (int d=nMaxDaysInfected-1; d>0; d--)
	  totaltreated+=nNumInfected[strain][res][1][d];
    for (vector<Recovered *>::iterator it = vRecovered.begin(); it!=vRecovered.end(); ++it)
      totalsusceptibles += (*it)->nNum;
    fProphy *= totaltreated/(double)totalsusceptibles;
  }

  // infect people
  for (int rnum=vRecovered.size()-1; rnum>=0; rnum--) {
    Recovered *susceptible = vRecovered[rnum];
    if (susceptible->nNum>0 && 
	(susceptible->nLastInfectionMonth<0 || 
	 nMonth-susceptible->nLastInfectionMonth>nMonthsImmune)) {
      for (int strain=0; strain<nMaxStrains && susceptible->nNum>0; strain++) {
	double susceptibility = ((*gm).getR0min() + fSeasonality[nDayOfYear] * ((*gm).getR0max()-(*gm).getR0min()))/(*gm).getR0max();  // susceptibility scaled using R0 and seasonality

	/* // no waning immunity in this model
	// figure out prior immunity
	if (susceptible->nLastInfectionMonth>0 && 
	    nMonth-susceptible->nLastInfectionMonth>nMonthsImmune) {
	  if (susceptible->nRecoveryMonth[strain]>=0 &&
	      nMonth-susceptible->nRecoveryMonth[strain]>nMonthsImmune) {
	    //	    susceptibility *= 1.0-pow(1.0-(*gm).getDrift()/12.0,nMonth-susceptible->nRecoveryMonth[strain]-nMonthsImmune); // immunity from this strain
	    // no waning immunity in this code!!!!!!!!!
	  } else { // immunity from other strains
	    for (int strain2=0; strain2<nMaxStrains; strain2++)
	      if (strain!=strain2 && 
		  susceptible->nRecoveryMonth[strain2]>=0 &&
		  nMonth-susceptible->nRecoveryMonth[strain2]>nMonthsImmune) {
		//		susceptibility *= 1.0-(*gm).getCrossStrain()*pow(1.0-(*gm).getDrift()/12.0,nMonth-susceptible->nRecoveryMonth[strain2]-nMonthsImmune); // cross-immunity from strain2
		// total cross immunity in this code!!!!!!!!!
		break; 
		// fix this - choose the most recent infection????
	      }
	  }
	}
*/
	// probabilities of escaping infection
	double qU=1.0,  // untreated susceptibles with wildtype
	  qT=1.0,       // treated susceptibles with wildtype
	  qR=1.0;       // any susceptible with resistant
	for (int day=0; day<nMaxDaysInfected; day++) {
	  double mult = susceptibility*fViralLoad[day];
	  if (strain>0)
	    mult *= (*gm).getAlpha();
	  if (nOldNumInfected[strain][0][0][day]>0) { //nNumWU[day]>0) {
	    qU *= pow(1.0-mult*p00,nOldNumInfected[strain][0][0][day]);
	    qT *= pow(1.0-mult*p10,nOldNumInfected[strain][0][0][day]);
	  }
	  if (nOldNumInfected[strain][0][1][day]>0) { //nNumWT[day]>0) {
	    qU *= pow(1.0-(*gm).getPhi()*mult*p01,nOldNumInfected[strain][0][1][day]);
	    qT *= pow(1.0-(*gm).getPhi()*mult*p11,nOldNumInfected[strain][0][1][day]);
	  }
	  if (nOldNumInfected[strain][1][0][day]>0) //nNumRU[day]>0)
	    qR *= pow(1.0-mult*(1.0-(*gm).getResistanceCost())*p00,nOldNumInfected[strain][1][0][day]);
	  if (nOldNumInfected[strain][1][1][day]>0) //nNumRT[day]>0)
	    qR *= pow(1.0-mult*(1.0-(*gm).getResistanceCost())*p00*(*gm).getPhi(),nOldNumInfected[strain][1][1][day]);
	}

	int infectedU = 0;
	int infectedUR = 0;
	int infectedT = 0;
	int infectedTR = 0;
	  assert(qU<=1.0 && qU>=0.0);
	  assert(qR<=1.0 && qR>=0.0);
	  int nNumProphy = round(susceptible->nNum * fProphy);
	  int nNumNoProphy = susceptible->nNum-nNumProphy;
	  if (nNumNoProphy>0) {
	    if (qU<1.0)
	      infectedU = gsl_ran_binomial(rng, 
					   1.0-qU,
					   nNumNoProphy);
	    susceptible->nNum -= infectedU;
	    if (qR<1.0)
	      infectedUR = gsl_ran_binomial(rng,  
					    1.0-qR,
					    nNumNoProphy);
	    susceptible->nNum -= infectedUR;
	  }
	  if (nNumProphy>0) {
	    assert(qT<=1.0 && qT>=0.0);
	    assert(qR<=1.0 && qR>=0.0);
	    if (qT<1.0)
	      infectedT = gsl_ran_binomial(rng, 
					   1.0-qT,
					   nNumProphy);
	    susceptible->nNum -= infectedT;
	    if (qR<1.0)
	      infectedTR = gsl_ran_binomial(rng, 
					    1.0-qR,
					    nNumProphy);
	    susceptible->nNum -= infectedTR;
	  }
	
	nNumInfected[strain][0][0][0] += infectedU;
	nNumInfected[strain][0][1][0] += infectedT;
	nNumInfected[strain][1][0][0] += infectedUR;
	nNumInfected[strain][1][1][0] += infectedTR;

	// recover in a few days
	if (infectedU>0 || infectedT>0 || infectedUR>0 || infectedTR>0) {
	  bool bFound=false;
	  int temp[NMAXSTRAINS];
	  memcpy(temp, susceptible->nRecoveryMonth, nMaxStrains*sizeof(int));
	  temp[strain] = nMonth;
	  for (unsigned int i=rnum+1; i<vRecovered.size(); i++) {
	    if (vRecovered[i]->nRecoveryMonth[strain]==nMonth &&
		memcmp(temp, vRecovered[i]->nRecoveryMonth, nMaxStrains*sizeof(int))==0) {
	      vRecovered[i]->nRecoverQueue[nRecoverQueueYesterday]+=infectedU+infectedT+infectedUR+infectedTR;
	      bFound=true;
	      break;
	    }
	  }
	  if (!bFound) { // create new recovered class
	    Recovered *r = new Recovered();
	    r->nLastInfectionMonth = nMonth;
	    memcpy(r->nRecoveryMonth, temp, nMaxStrains*sizeof(int));
	    r->nRecoverQueue[nRecoverQueueYesterday]=infectedU+infectedT+infectedUR+infectedTR;
	    vRecovered.push_back(r);
	  }
	}
      }
    }
  }

  // reassortment?
  if (nMaxStrains>1 && (*gm).useResistance() && (*gm).getReassortFactor()>0.0) {
    double fReassortFactor = (*gm).getReassortFactor();
    for (int strain=0; strain<nMaxStrains; strain++)
      for (int res=0; res<2; res++)
	for (int treat=0; treat<2; treat++)
	  if (nNumInfected[strain][res][treat][0]>0) {
	    double fPropInfected1 = nNumInfected[strain][res][treat][0]/(double)nNumSusceptiblesStart;
	    for (int strain2=strain; strain2<nMaxStrains; strain2++)
	      for (int res2=res; res2<2; res2++)
		for (int treat2=treat; treat2<2; treat2++)
		  if ((strain!=strain2 || res!=res2 || treat!=treat2) &&
		      nNumInfected[strain2][res2][treat2][0]>0) {
		    //		    double fPropInfected2 = nNumInfected[strain2][res2][treat2][0]/(double)nNumSusceptiblesStart;
		    int nNumReassort = gsl_ran_binomial(rng, 
							fReassortFactor * fPropInfected1,
							nNumInfected[strain2][res2][treat2][0]);
		    if (nNumReassort>0) { // slow!
		      nNumInfected[strain][res][treat][0]-=nNumReassort;
		      nNumInfected[strain2][res2][treat2][0]-=nNumReassort;
		      assert(nNumInfected[strain][res][treat][0]>=0);
		      assert(nNumInfected[strain2][res2][treat2][0]>=0);
		      while (nNumReassort-- > 0) {
			int which = gsl_rng_uniform_int(rng, 8);
			nNumInfected[(which&0x01)?strain:strain2][(which&0x02)?res:res2][(which&0x04)?treat:treat2][0]++;
		      }
		    }
		  }
	  }
  }

  // strain mutation and resistance mutation in untreated
  if ((*gm).getMutation()>0.0) {
    for (int strain=0; strain<nMaxStrains; strain++) {
      for (int day=0; day<nMaxDaysInfected; day++) {
	if (strain+1<nMaxStrains)
	  for (int res=0; res<2; res++)
	    for (int treat=0; treat<2; treat++)
	      if (nNumInfected[strain][res][treat][day]>0) { // change strains
		int nNumMutants = gsl_ran_binomial(rng, 
						   (*gm).getMutation(),
						   nNumInfected[strain][res][treat][day]);
		if (nNumMutants>0) {
		  assert(strain+1<nMaxStrains);
		  nNumInfected[strain][res][treat][day]-=nNumMutants;
		  nNumInfected[strain+1][res][treat][day]+=nNumMutants;
		}
	      }

	if ((*gm).useResistance()) {
	  int nNumMutants=0;
	  if (nNumInfected[strain][0][0][day]>0) // gain resistance
	    nNumMutants = gsl_ran_binomial(rng, 
					   (*gm).getMutation(),
					   nNumInfected[strain][0][0][day]);
	  if (nNumInfected[strain][1][0][day]>0) // lose resistance
	    nNumMutants -= gsl_ran_binomial(rng, 
					    (*gm).getMutation(),
					    nNumInfected[strain][1][0][day]);
	  if (nNumMutants!=0) {
	    nNumInfected[strain][0][0][day]-=nNumMutants;
	    nNumInfected[strain][1][0][day]+=nNumMutants;
	  }
	}
      }
    }
  }

  // resistance mutation in treated
  if ((*gm).getMutationTreatment()>0.0 && (*gm).useResistance()) {
    for (int strain=0; strain<nMaxStrains; strain++) {
      for (int day=0; day<nMaxDaysInfected; day++)
	if (nNumInfected[strain][0][1][day]>0) { // gain resistance
	  int nNumMutants = gsl_ran_binomial(rng, 
					     (*gm).getMutationTreatment(), 
					     nNumInfected[strain][0][1][day]);
	  nNumInfected[strain][0][1][day]-=nNumMutants;
	  nNumInfected[strain][1][1][day]+=nNumMutants;
	}
    }
  }

  // start treatment
  double gamma = (fCityGamma>=0.0?fCityGamma:(*gm).getGlobalGamma());
  if (gamma>0.0)
    for (int strain=0; strain<nMaxStrains; strain++) {
      for (int day=nTreatmentDelay; day<nMaxDaysInfected; day++) {
	int n = gsl_ran_binomial(rng, gamma, nNumInfected[strain][0][0][day]);
	if (n>0) {
	  nNumInfected[strain][0][0][day] -= n;
	  nNumInfected[strain][0][1][day] += n;
	}
	n=gsl_ran_binomial(rng, gamma, nNumInfected[strain][1][0][day]);
	if (n>0) {
	  nNumInfected[strain][1][0][day] -= n;
	  nNumInfected[strain][1][1][day] += n;
	}
      }
    }
}
