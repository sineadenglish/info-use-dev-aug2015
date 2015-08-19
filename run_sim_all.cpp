/***************************************************************************************** 
 GROWTH MODEL SIMULATIONS
 
 run_sim_all.cpp
 
 Runs forward simulations based on parameter combinations and decision array produced in growth_info.cpp
 *All file - all for mortality/not, whether have experiment etc.
 
 Created by:  Sinead English 
 Modified:		Jun 2015
 
 ******************************************************************************************/ 

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <vector> 
#include <iostream>
#include "growthfun.h"

using namespace::std;

double decArray[SIZESTEPS][PSTEPS];
double infoArray[SIZESTEPS][PSTEPS];
double binArray[SIZESTEPS][PSTEPS];
int i,j;
double PHI, KAPPA_GOOD, KAPPA_BAD, CONST_GOOD, CONST_BAD, MU0, MU1_GOOD, MU1_BAD, GAMMA, PFOOD_GOOD, PFOOD_BAD; // save parameter values for that particular input

// function to read in data from decision array and parameters files 
void readData(int parno)
{
	
	char decname[50];
	char parname[50];
	char infoname[50];
	
	sprintf(decname, "Decision%1d.txt",parno);
	sprintf(parname, "Params%1d.txt",parno);
	sprintf(infoname, "InfoVal%1d.txt",parno);
	
	// read in decision array
	FILE *decfile;
  double myvariable;
  decfile=fopen(decname, "r");
  for(i = 0; i < SIZESTEPS; i++)
  {
    for (j = 0 ; j < PSTEPS; j++)
    {
      fscanf(decfile,"%lf",&myvariable);
			decArray[i][j] = myvariable;
    }
  }
	fclose(decfile);
	
	// read in value of information array
	FILE *infofile;
  double myvariable2;
  infofile=fopen(infoname, "r");
  for(i = 0; i < SIZESTEPS; i++)
  {
    for (j = 0 ; j < PSTEPS; j++)
    {
      fscanf(infofile,"%lf",&myvariable2);
			infoArray[i][j] = myvariable2;
    }
  }
	fclose(infofile);
	
	
	// read in parameter values
	FILE *parfile;		
  parfile=fopen(parname, "r");
	fscanf(parfile,"%lf",&PHI);
	fscanf(parfile,"%lf",&KAPPA_GOOD);
	fscanf(parfile,"%lf",&KAPPA_BAD);
	fscanf(parfile,"%lf",&CONST_GOOD);
	fscanf(parfile,"%lf",&CONST_BAD);
	fscanf(parfile,"%lf",&MU0);
	fscanf(parfile,"%lf",&MU1_GOOD);
	fscanf(parfile,"%lf",&MU1_BAD);
	fscanf(parfile,"%lf",&GAMMA);
	fscanf(parfile,"%lf",&PFOOD_GOOD);
	fscanf(parfile,"%lf",&PFOOD_BAD);
	fclose(parfile);
	
}

// convert decision array to binary array (for interpolating at boundary)
void convertBin()
{
	for(i = 0; i < SIZESTEPS; i++)
  {
    for (j = 0 ; j < PSTEPS; j++)
    {
      if(decArray[i][j] > 0.0) binArray[i][j] = 1.0; else binArray[i][j] = 0.0;
    }
  }	
}

// Monte-Carlo forward simulation for one individual: run simulations and output to file
void simulate(FILE* simfile, FILE* simfile2, double startPrior, double pfoodExpt, int AGE1, int AGE2, int env_type, int mort, int simno, int it, vector<int> age, vector<double> mass, vector<double> belief, vector<double> forage, vector<int> env, vector<int> state, vector<double> info)   
{	
	int thisAge, thisEnv, thisState;			
	double thisMass, thisBelief, thisAct, thisInfo, mortGood, mortBad;
	double pFindFood = 0.0; 
	
	// initialise values
	thisAge = 1;
	thisMass = 1.0; 
	thisBelief = startPrior;
	thisAct = 0.0;
	thisEnv = env_type;   
	thisState = 1;
	thisInfo = 0.0;
	mortGood = 0.0;
	mortBad = 0.0;
	
	// place values into vector
	age.push_back(thisAge);
	mass.push_back(thisMass);
	belief.push_back(thisBelief); 
	forage.push_back(thisAct);
	env.push_back(thisEnv);
	state.push_back(thisState);
	info.push_back(thisInfo);
	
	while (thisState == 1) // continue iterations while individual still growing (not dead or metamorphosed) 
	{		
		// work out if forage or metamorphose (using binary array)
		double probGrow = interpolate(binArray, thisMass, thisBelief);  
		
		// if forages, find out whether finds food, dies and next size and belief
		if ((rand()*1.0)/(RAND_MAX*1.0) < probGrow) 
		{
			thisState = 1;
			thisAct =	interpolate(decArray, thisMass, thisBelief);
			if(isnan(thisAct)) thisAct = nearestval(decArray, thisMass, thisBelief);
			thisInfo = interpolate(infoArray, thisMass, thisBelief);

			// probability of dying in each environment
			mortGood = MU0*(1.0 + MU1_GOOD*pow(thisAct,GAMMA)); 
			mortBad = MU0*(1.0 + MU1_BAD*pow(thisAct,GAMMA));
			
			//// whether finds food 
			
			// if in experimental window
			if(thisAge >= AGE1 && thisAge < AGE2) pFindFood = pfoodExpt*thisAct;
			
			// otherwise evolved environment
			else 
			{
				if(thisEnv==0) pFindFood = PFOOD_GOOD*thisAct;
				if(thisEnv==1) pFindFood = PFOOD_BAD*thisAct;
			}

			//// update mass and belief
			
			// if finds food
			if((rand()*1.0)/(RAND_MAX*1.0) < pFindFood) 
			{
				thisBelief = thisBelief*PFOOD_GOOD*(1.0-mortGood)/(thisBelief*PFOOD_GOOD*(1.0-mortGood) + (1.0-thisBelief)*PFOOD_BAD*(1.0-mortBad));				
				thisMass += 1.0; // increment mass by 1 unit
			}
			
			// if does not find food
			else
			{
				thisBelief = thisBelief*(1-PFOOD_GOOD*thisAct)*(1.0-mortGood)/(thisBelief*(1.0-PFOOD_GOOD*thisAct)*(1.0-mortGood) + (1.0-thisBelief)*(1.0-PFOOD_BAD*thisAct)*(1.0-mortBad));			
				thisMass = thisMass; // mass remains the same
			}
			
			// delimit mass and belief
			if(thisMass > SIZEMAX) thisMass = SIZEMAX;
			if(thisBelief < 10e-6) thisBelief = 10e-6;
			if(thisBelief > 1-10e-6) thisBelief = 1-10e-6;
			
			// mortality
			if(mort==1)
			{
				if(thisEnv==0 && (rand()*1.0)/(RAND_MAX*1.0) < mortGood) thisState = 0;
				if(thisEnv==1 && (rand()*1.0)/(RAND_MAX*1.0) < mortBad) thisState = 0; 
				if(thisEnv==2 && (rand()*1.0)/(RAND_MAX*1.0) < (mortGood + mortBad)/2.0) thisState = 0;
			}
		}
		
		// if metamorphose, set next state=2, next act=NAN
		else 
		{
			thisState = 2;
			thisAct = NAN;
		}
		
		age.push_back(thisAge += 1);
		mass.push_back(thisMass);
		belief.push_back(thisBelief); 
		forage.push_back(thisAct);
		env.push_back(thisEnv);
		state.push_back(thisState);
		info.push_back(thisInfo);
	}
	
	// store data by age for first 100 individuals for each line of simulation run (with corresponding iteration no.) in file
	if (it <= 100)
	{
		for (unsigned i=0; i<mass.size(); i++)
		{
			fprintf(simfile, "%4d %4d %4d %5.2f %5.3f %5.3f %4d %4d %5.15f", simno, it, age.at(i), mass.at(i), belief.at(i), forage.at(i), env.at(i), state.at(i), info.at(i));
			fprintf(simfile,"\n");
		}
	}
	
	// store end results for all individuals	
	fprintf(simfile2, "%4d %4d %4d %5.2f %5.3f %4d %4d", simno, it, thisAge, thisMass, thisBelief, thisEnv, thisState); 
	fprintf(simfile2,"\n");
	
}


int main(int argc, char* argv[])
{	
	// read in values from argv character string
	int parno = atoi(argv[2]);				// parameter number
	int	nsims = atoi(argv[4]);				// number of simulations to run
	int nind = atoi(argv[6]); 				// number of individuals in each environment to simulate
	int mort = atoi(argv[8]);					// whether allow mortality (1) or not (0)
	double startPrior = atof(argv[10]); 
	double pfoodExpt = atof(argv[12]);		// experimental pfood (rather than that in parameter values) to simulate (if no experiment, set to 0.0)
	int AGE1 = atoi(argv[14]); // expt period start
	int AGE2 = atoi(argv[15]); // expt period end
	
	readData(parno);
	convertBin();

	FILE *simfile, *simfile2; // files to store sample of data by age and final status data for all individuals 
	
	char simname[50];
	char simname2[50];
	
	sprintf(simname, "SimOutputInd.%1d.%1d.%0.1f.%0.1f.%1d.txt",parno,mort,startPrior,pfoodExpt,AGE1);  
	sprintf(simname2, "SimOutputAll.%1d.%1d.%0.1f.%0.1f.%1d.txt",parno,mort,startPrior,pfoodExpt,AGE1);  
	
	if ((simfile = fopen(simname,"wt")) == NULL)
	{
		printf("%s\n", "Error opening file");
		exit(1);
	}
	
	if ((simfile2 = fopen(simname2,"wt")) == NULL)
	{
		printf("%s\n", "Error opening file");
		exit(1);
	}
	
	fprintf(simfile, "  simno it age mass belief forage env state info\n"); 
	fprintf(simfile2, "  simno it age mass belief env state\n"); 

	int ltime = time(NULL);              			
	int utime = (unsigned int) ltime/2;
	srand(utime);

	// create empty vectors
	vector<int> age;
	vector<double> mass;
	vector<double> belief;
	vector<double> forage;
	vector<int> env;
	vector<int> state;
	vector<double> info;
	
	for(int simno = 1; simno <= nsims; simno++)
	{

		// run simulations
			for (int env_type=0; env_type <= 1; env_type++) 
			{
				for (int it=1; it<=nind; it++) simulate(simfile, simfile2, startPrior, pfoodExpt, AGE1, AGE2, env_type, mort, simno, it, age, mass, belief, forage, env, state, info);			
			}
	
		cout << "Finished running simulation number: " << simno << endl;		
	}

	fclose(simfile);
	fclose(simfile2);
	
	return 0;	
}