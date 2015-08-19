/***************************************************************************************** 
 GROWTH MODEL
 
 Calculates optimal strategy in terms of how much effort to put into foraging and when to mature; 
 and decision array for information
 
 Created by:  Sinead English 
 Date:				Dec 2012
 Updated:			Jun 2015
 
 ******************************************************************************************/

// Define all constants and functions

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <vector> 
#include <iostream>
#include "growthfun.h"

using namespace std;

// ************ declare constants ************

const double	DIFFCRIT = 10e-7; // convergence criterion

// create fitness value and decision arrays (each cell describes size and belief state)
double old_value[SIZESTEPS][PSTEPS];  // store previous value of fitness
double new_value[SIZESTEPS][PSTEPS];  // store new value of fitness (after one interation)  
double decArray[SIZESTEPS][PSTEPS];   // the optimal decision array (decision = alpha (0-1) or mature (NAN)) 
double fitArray[SIZESTEPS][PSTEPS];		// the reproductive value array (assuming optimal alpha)
double fitAlphaArray[SIZESTEPS][PSTEPS][ALPHASTEPS];	// the reproductive value array for all levels of alpha
double infoArray[SIZESTEPS][PSTEPS];									// the value of information array 

// declare parameters
int i, j, k, p_val, size_val, alpha_i;

double alpha, mu_good, mu_bad;
double pgood, size;
double exp_value_success, exp_value_failure;
double current_value, best_value;
double pgood_given_success, pgood_given_failure;
double best_alpha, diff;

// create vector of differences between arrays from one iteration to the next (for convergence)
std::vector<double> diffvec (1,diff); 

// save results to file
FILE *file1, *file2, *file3, *file4, *file5, *file0;				
/* 
 file0 = information array (for size x belief)
 file1 = reproductive value array (for size x belief states) 
 file2 = optimal alpha array (for size x belief states)
 file3 = vector of convergence values
 file4 = parameters used in model
 file5 = reproductive value array (for size x belief x alpha states) => used for calculating fitness value of information
*/


// ************ define functions ************

// NOTE TO ANDY: may be better way of doing functions (not separate for each environment) but stuck like this for now as reading in parameters as unique string of values

// find fitness at maturity in good environment
double mature_good(double size, double &PHI, double &KAPPA_GOOD, double &CONST_GOOD)
{
	return 1.0/(CONST_GOOD + exp((-1.0*PHI)*(size - KAPPA_GOOD)));      
}    

// find fitness at maturity in bad environment
double mature_bad(double size, double &PHI, double &KAPPA_BAD, double &CONST_BAD)
{
	return 1.0/(CONST_BAD + exp((-1.0*PHI)*(size - KAPPA_BAD)));      
}      

// increase in size on finding food 
double grow(double size, double inc)
{
	if((size+inc) > SIZEMAX) return SIZEMAX;
	else if((size+inc) < SIZEMIN) return SIZEMIN;
	else return size+inc;
}

// mortality as function of alpha and size in good environment
double pdead_good(double alpha, double size, double &MU0, double &MU1_GOOD, double &GAMMA)
{
	return MU0*(1.0 + MU1_GOOD*pow(alpha,GAMMA));
}

// mortality as function of alpha and size in bad environment
double pdead_bad(double alpha, double size, double &MU0, double &MU1_BAD, double &GAMMA)
{
	return MU0*(1.0 + MU1_BAD*pow(alpha,GAMMA));
}


// Initialise values [NOTE TO ANDY: *at random, zero or fitness if mature? => shouldn't make difference, for now done at random] 

void Initialise() 
{
	for (p_val = 0; p_val < PSTEPS; p_val++)
	{ 
		for (size_val = 0; size_val < SIZESTEPS; size_val++)
		{ 
			old_value[size_val][p_val] = (rand()*1.0)/(RAND_MAX*1.0);
      new_value[size_val][p_val] = 0.0;    
			decArray[size_val][p_val] = 0.0;
			infoArray[size_val][p_val] = NAN;

			for (alpha_i = 0; alpha_i < ALPHASTEPS; alpha_i++)
			{
				fitAlphaArray[size_val][p_val][alpha_i] = 0.0/0;  // so NaN [Not sure why didn't just do this as NAN?]
			}
		} 
	} 
}

// First find the optimal decisions (alpha, mature) at 

void Iterate_Ends(double &PHI, double &KAPPA_GOOD, double &KAPPA_BAD, double &CONST_GOOD, double &CONST_BAD, double &MU0, double &MU1_GOOD, double &MU1_BAD, double &GAMMA, double &PFOOD_GOOD, double &PFOOD_BAD) 
{
	diff = 10.0; // arbitrary starting difference between old and new value arrays
	i = 0;	// starting iteration
	
	while(diff > DIFFCRIT) // continue iterating while maximum difference between equivalent cells in new and old fitness arrays exceeds convergence criterion
	{  
		i = i+1; 
		
		// cycle through levels of belief (here p=0 and p=1) and size (all levels)
		for (p_val = 0; p_val < PSTEPS; p_val+=(PSTEPS-1)) 
		{
			// convert integer placer value to absolute value (probably the same?)
			pgood = (double)p_val/((double)(PSTEPS-1));
			
			for (size_val = 0; size_val < SIZESTEPS; size_val++)
			{ 
        // Calc optimal activity level, alpha_hat
        best_value = 0.0;
        best_alpha = 0.0;
				std::vector<double> bestAlpha; // vector to store ties
        				
        for (alpha_i = 0; alpha_i < ALPHASTEPS; alpha_i++)
				{					
					size = (double)size_val*(SIZEMAX-SIZEMIN)/((double)(SIZESTEPS-1)) + SIZEMIN;
					
					alpha = (double)alpha_i/((double)(ALPHASTEPS - 1));      
					
					// probability of survival given particular alpha and size value for each environment
					mu_good = (1.0 - pdead_good(alpha, size, MU0, MU1_GOOD, GAMMA));  
					mu_bad = (1.0 - pdead_bad(alpha, size, MU0, MU1_BAD, GAMMA));  
						
					// find expected value for alpha depending on belief in environment (*not updated here* - individuals at p=0 or 1) and whether find food
					exp_value_success = interpolate(old_value, grow(size, 1.0), pgood);
					exp_value_failure = interpolate(old_value, size, pgood);
					
					current_value	= pgood*mu_good*(PFOOD_GOOD*alpha*exp_value_success + (1.0 - PFOOD_GOOD*alpha)*exp_value_failure) 
					+ (1.0 - pgood)*mu_bad*(PFOOD_BAD*alpha*exp_value_success + (1.0 - PFOOD_BAD*alpha)*exp_value_failure);
					
					// dealing with ties
					if (current_value == best_value) 
					{
						bestAlpha.push_back(alpha);
					}
					
					if (current_value > best_value)   
					{
						best_value = current_value;         
						bestAlpha.resize(1); 
						bestAlpha[0] = alpha;
					}
				} // for each alpha
	
				// get best_alpha in case of ties (is equivalent to alpha if not) 
				random_shuffle(bestAlpha.begin(), bestAlpha.end());
				best_alpha = bestAlpha[0];
        
				// compare value if forage to value if mature (depends on probability of being in good or bad environment) 
        if (best_value > pgood*mature_good(size, PHI, KAPPA_GOOD, CONST_GOOD) + (1.0 - pgood)*mature_bad(size, PHI, KAPPA_BAD, CONST_BAD))
				{
          new_value[size_val][p_val] = best_value;
          decArray[size_val][p_val] = best_alpha;
				}
        else  // mature
				{
          new_value[size_val][p_val] = pgood*mature_good(size, PHI, KAPPA_GOOD, CONST_GOOD) + (1.0 - pgood)*mature_bad(size, PHI, KAPPA_BAD, CONST_BAD);
					decArray[size_val][p_val] = NAN;
				} 
			}  // for each p value
		}  // for each size
		
		cout << endl;
		
		// calculate absolute difference between old_value and new_value (all values)
		double diff_value[SIZESTEPS][PSTEPS];
		
		for(size_val = 0; size_val < SIZESTEPS; size_val++)   
		{
			for(p_val = 0; p_val < PSTEPS; p_val++) 
			{
				diff_value[size_val][p_val] = absol(new_value[size_val][p_val]-old_value[size_val][p_val]);
			}
		}
		
		// find maximum value for diff_value array
		diff = maxval(*diff_value, SIZESTEPS*PSTEPS);
		
		// add difference to vector so can save it
		diffvec.insert (diffvec.end(), diff); 
		
		// now old_value is new_value for next iteration
		for (size_val = 0; size_val< SIZESTEPS; size_val++)
		{  
			for (p_val = 0; p_val < PSTEPS; p_val++)
			{  
				old_value[size_val][p_val]= new_value[size_val][p_val];
			} 
		}  
		
		// print i
		cout << endl << "iteration: " << i;
		cout << "   " << "difference: "  << diff << endl;
		
	}  // end of convergence criterion
}

// and do iterations again for 0 < p_e1 < 1
void Iterate_Mid(double &PHI, double &KAPPA_GOOD, double &KAPPA_BAD, double &CONST_GOOD, double &CONST_BAD, double &MU0, double &MU1_GOOD, double &MU1_BAD, double &GAMMA, double &PFOOD_GOOD, double &PFOOD_BAD) 
{	
	diff = 10.0; // arbitrary starting difference between old and new value arrays
	
	std::vector<double> diffvec (1,diff);  // vector to store difference values
	
	i = 0; // reset i
	
	while(diff > DIFFCRIT) // if reaching convergence criterion
	// for(int it = 0; it < ITERATIONS; it++) // if doing fixed no. iterations
	{  
		
		i = i+1;
		
		for (p_val = 1; p_val < (PSTEPS-1); p_val++)
		{  
			pgood = (double)(p_val)/((double)(PSTEPS-1)); // convert integer placer value to absolute value (btween 0 and 1) 
			
			for (size_val = 0; size_val < SIZESTEPS; size_val++)
			{ 
        // Calc optimal activity level, alpha_hat
        best_value = 0.0;
        best_alpha = 0.0;
				std::vector<double> bestAlpha; // vector to store ties

				// convert integer placer value to absolute value (probably the same?)
				size = (double)size_val*(SIZEMAX-SIZEMIN)/((double)(SIZESTEPS-1)) + SIZEMIN;
				
				for (alpha_i = 0; alpha_i < ALPHASTEPS; alpha_i++)
				{
					
					size = (double)size_val*(SIZEMAX-SIZEMIN)/((double)(SIZESTEPS-1)) + SIZEMIN;
					
					alpha = (double)alpha_i/((double)(ALPHASTEPS - 1));      
					
					// probability of death given that alpha and size in each environment 
					mu_good = (1.0 - pdead_good(alpha, size, MU0, MU1_GOOD, GAMMA));  
					mu_bad = (1.0 - pdead_bad(alpha, size, MU0, MU1_BAD, GAMMA));  
	
					// update belief if finds food or not (for a given foraging effort, and knowing that it survived)
					pgood_given_success = pgood*PFOOD_GOOD*mu_good/(pgood*PFOOD_GOOD*mu_good + (1.0-pgood)*PFOOD_BAD*mu_bad);				
					pgood_given_failure = pgood*(1.0-PFOOD_GOOD*alpha)*mu_good/(pgood*(1.0-PFOOD_GOOD*alpha)*mu_good + (1.0-pgood)*(1.0-PFOOD_BAD*alpha)*mu_bad); 
					
					// find expected value for alpha depending on belief in environment and whether find food
					exp_value_success = interpolate(old_value, grow(size, 1.0), pgood_given_success);
					exp_value_failure = interpolate(old_value, size, pgood_given_failure);
					
					current_value	= pgood*mu_good*(PFOOD_GOOD*alpha*exp_value_success + (1.0 - PFOOD_GOOD*alpha)*exp_value_failure) 
					+ (1.0-pgood)*mu_bad*(PFOOD_BAD*alpha*exp_value_success + (1.0 - PFOOD_BAD*alpha)*exp_value_failure);
					
					// dealing with ties
					if (current_value == best_value) 
					{
						bestAlpha.push_back(alpha);
					}
					
					if (current_value > best_value)   
					{
						best_value = current_value;         
						bestAlpha.resize(1);
						bestAlpha[0] = alpha;
					}
				} // for each alpha
				
				// get best_alpha in case of ties (is equivalent to alpha if not) 				
				random_shuffle(bestAlpha.begin(), bestAlpha.end());
				best_alpha = bestAlpha[0];

        if (best_value > pgood*mature_good(size, PHI, KAPPA_GOOD, CONST_GOOD) + (1.0 - pgood)*mature_bad(size, PHI, KAPPA_BAD, CONST_BAD))
				{
          new_value[size_val][p_val] = best_value;
          decArray[size_val][p_val] = best_alpha;
				}
        else // mature
				{
          new_value[size_val][p_val] = pgood*mature_good(size, PHI, KAPPA_GOOD, CONST_GOOD) + (1.0 - pgood)*mature_bad(size, PHI, KAPPA_BAD, CONST_BAD);
          decArray[size_val][p_val] = 0.0;
				} 
        
			}  // for each p value
		}  // for each size
		
		
		// calculate absolute difference between old_value and new_value (all values)
		double diff_value[SIZESTEPS][PSTEPS];
		
		for(size_val = 0; size_val < SIZESTEPS; size_val++)   
		{
			for (p_val = 0; p_val < PSTEPS; p_val++)
			{
				diff_value[size_val][p_val] = absol(new_value[size_val][p_val]-old_value[size_val][p_val]);
			}
		}
		
		// find maximum value for diff_value array
		diff = maxval(*diff_value, SIZESTEPS*PSTEPS);
		
		// add difference to vector so can plot it
		diffvec.insert (diffvec.end(), diff); 
		
		// now old_value is new_value for next iteration
		for (size_val = 0; size_val < SIZESTEPS; size_val++)
		{  
			for (p_val = 0; p_val < PSTEPS; p_val++)
			{  
				old_value[size_val][p_val]= new_value[size_val][p_val];
			} 
		}  
		
		
		// print i
		cout << endl << "iteration: " << i;
		cout << "   " << "difference: "  << diff << endl;
		
	}  // end of convergence criterion or iteration steps
}


// After convergence, run one last iteration to get RV values for all levels of alpha at all {S,P} combinations (for information equation) 
void Iterate_Final(double &PHI, double &KAPPA_GOOD, double &KAPPA_BAD, double &CONST_GOOD, double &CONST_BAD, double &MU0, double &MU1_GOOD, double &MU1_BAD, double &GAMMA, double &PFOOD_GOOD, double &PFOOD_BAD) 
{
	
	cout << "Start final iteration (to save RV array for all combinations of size, belief and alpha)"  << endl;
	
	for (p_val = 0; p_val < PSTEPS; p_val++)
	{  
		pgood = (double)(p_val)/((double)(PSTEPS-1)); // convert integer placer value to absolute value (btween 0 and 1) 
		
		for (size_val = 0; size_val < SIZESTEPS; size_val++)
		{
			// Calc optimal activity level, alpha_hat
			best_value = 0.0;
			best_alpha = 0.0;
			std::vector<double> bestAlpha; // vector to store ties
						
			// convert integer placer value to absolute value (probably the same?)
			size = (double)size_val*(SIZEMAX-SIZEMIN)/((double)(SIZESTEPS-1)) + SIZEMIN;
			
			for (alpha_i = 0; alpha_i < ALPHASTEPS; alpha_i++)
			{
				
				size = (double)size_val*(SIZEMAX-SIZEMIN)/((double)(SIZESTEPS-1)) + SIZEMIN;
				
				alpha = (double)alpha_i/((double)(ALPHASTEPS - 1));      
				
				// probability of death given that alpha and size in each environment 
				mu_good = (1.0 - pdead_good(alpha, size, MU0, MU1_GOOD, GAMMA));  
				mu_bad = (1.0 - pdead_bad(alpha, size, MU0, MU1_BAD, GAMMA));  
						
				// update belief if finds food or not (for a given foraging effort, and knowing that it survived)
				pgood_given_success = pgood*PFOOD_GOOD*mu_good/(pgood*PFOOD_GOOD*mu_good + (1.0-pgood)*PFOOD_BAD*mu_bad);				
				pgood_given_failure = pgood*(1.0-PFOOD_GOOD*alpha)*mu_good/(pgood*(1.0-PFOOD_GOOD*alpha)*mu_good + (1.0-pgood)*(1.0-PFOOD_BAD*alpha)*mu_bad); 
				
				// find expected value for alpha depending on belief in environment and whether find food
				exp_value_success = interpolate(old_value, grow(size, 1.0), pgood_given_success);
				exp_value_failure = interpolate(old_value, size, pgood_given_failure);
				
				current_value	= pgood*mu_good*(PFOOD_GOOD*alpha*exp_value_success + (1.0 - PFOOD_GOOD*alpha)*exp_value_failure) 
				+ (1.0-pgood)*mu_bad*(PFOOD_BAD*alpha*exp_value_success + (1.0 - PFOOD_BAD*alpha)*exp_value_failure);
				
				// store in 3D array for info equation
				fitAlphaArray[size_val][p_val][alpha_i] = current_value;
				
				// dealing with ties
				if (current_value == best_value) 
				{
					bestAlpha.push_back(alpha);
				}
				
				if (current_value > best_value)   
				{
					best_value = current_value;         
					bestAlpha.resize(1);
					bestAlpha[0] = alpha;
				}				
			} // for each alpha
			
			// get best_alpha in case of ties (is equivalent to alpha if not) 				
			random_shuffle(bestAlpha.begin(), bestAlpha.end());
			best_alpha = bestAlpha[0];
			
			if (best_value > pgood*mature_good(size, PHI, KAPPA_GOOD, CONST_GOOD) + (1.0 - pgood)*mature_bad(size, PHI, KAPPA_BAD, CONST_BAD))
			{
				new_value[size_val][p_val] = best_value;
				decArray[size_val][p_val] = best_alpha;
			}
			else // mature
			{
				new_value[size_val][p_val] = pgood*mature_good(size, PHI, KAPPA_GOOD, CONST_GOOD) + (1.0 - pgood)*mature_bad(size, PHI, KAPPA_BAD, CONST_BAD);
				decArray[size_val][p_val] = NAN;
			} 
			
		}  // for each p value
	}  // for each size
	
	// calculate absolute difference between old_value and new_value (all values)
	double diff_value[SIZESTEPS][PSTEPS];
	
	for(size_val = 0; size_val < SIZESTEPS; size_val++)   
	{
		for (p_val = 0; p_val < PSTEPS; p_val++)
		{
			diff_value[size_val][p_val] = absol(new_value[size_val][p_val]-old_value[size_val][p_val]);
		}
	}
	
	// find maximum value for diff_value array
	diff = maxval(*diff_value, SIZESTEPS*PSTEPS);
	
	// add difference to vector so can plot it
	diffvec.insert (diffvec.end(), diff); 
	
	// now old_value is new_value for next iteration
	for (size_val = 0; size_val < SIZESTEPS; size_val++)
	{  
		for (p_val = 0; p_val < PSTEPS; p_val++)
		{  
			old_value[size_val][p_val]= new_value[size_val][p_val];
		} 
	}  
	
	// save fitness values into fitArray
	for (size_val = 0; size_val < SIZESTEPS; size_val++)
	{  
		for (p_val = 0; p_val < PSTEPS; p_val++)
		{  
			fitArray[size_val][p_val]= old_value[size_val][p_val];
		} 
	}
	
	cout << "End final iteration (to save RV array for all combinations of size, belief and alpha)"  << endl;

}

// function to calculate, for all size/belief combinations, fitness value of information
void Val_Info_Grow(double &PHI, double &KAPPA_GOOD, double &KAPPA_BAD, double &CONST_GOOD, double &CONST_BAD, double &MU0, double &MU1_GOOD, double &MU1_BAD, double &GAMMA, double &PFOOD_GOOD, double &PFOOD_BAD)
{
	cout << "Begin value of information/growth calculation"  << endl;

	double thisvalInfo, pgood0, size0, alpha_size0_pgood0, alpha_sizeSucc_pgoodSucc, alpha_size0_pgoodFail, alpha_size0_pgoodSucc, alpha_sizeSucc_pgood0, survGood, survBad, pgoodSucc, pgoodFail, probSurvSucc, probSurvFail, RVSuccBoth, RVFailBoth, RVSuccSize, RVFailSize, RVSuccInfo, valInfo;
	
	for(size_val=0; size_val<SIZESTEPS; size_val++)
	{
		size0 = (double)size_val*(SIZEMAX-SIZEMIN)/((double)(SIZESTEPS-1)) + SIZEMIN;
		
		for(p_val=0;p_val<PSTEPS; p_val++)
		{
			pgood0 = (double)(p_val)/((double)(PSTEPS-1));
			
			alpha_size0_pgood0 = interpolate(decArray, size0, pgood0);
			thisvalInfo = NAN;
			
			// work out pr(env_good, succeed) and pr(env_good, fail); know next size (S' = S+1 or S' = S)
			
			if(!isnan(alpha_size0_pgood0))
			{
				// probability of survival given that alpha and size in each environment 
				survGood = (1.0 - pdead_good(alpha_size0_pgood0, size0, MU0, MU1_GOOD, GAMMA));
				survBad = (1.0 - pdead_bad(alpha_size0_pgood0, size0, MU0, MU1_BAD, GAMMA));
				
				// update belief if finds food or not (for a given foraging effort, and knowing that it survived)
				pgoodSucc = pgood0*PFOOD_GOOD*survGood/(pgood0*PFOOD_GOOD*survGood + (1.0-pgood0)*PFOOD_BAD*survBad);				
				pgoodFail = pgood0*(1.0-(PFOOD_GOOD*alpha_size0_pgood0))*survGood/(pgood0*(1.0-(PFOOD_GOOD*alpha_size0_pgood0))*survGood + (1.0-pgood0)*(1.0-(PFOOD_BAD*alpha_size0_pgood0))*survBad); 
				
				// find alpha at size/belief state depending on whether finds food or not
				alpha_size0_pgoodFail = interpolate(decArray, size0, pgoodFail);
				alpha_size0_pgoodSucc = interpolate(decArray, size0, pgoodSucc);
				
				// can only calculate value of information for when not at maturation threshold (i.e. value of alpha is a number not NaN) 
				if(!isnan(alpha_size0_pgoodSucc) && !isnan(alpha_size0_pgoodFail))
				{
					
					// to interpolate fitness by alpha array, need to get P-by-alpha arrays for S0
					double RV_size0[PSTEPS][ALPHASTEPS];
					
					for (j = 0; j < PSTEPS; j++)
					{
						for (k = 0; k < ALPHASTEPS; k++)
						{
							RV_size0[j][k] = fitAlphaArray[size_val][j][k];
						}
					}
					
					// calculate probability individual survives and finds food or fails to find food
					probSurvSucc = pgood0*survGood*PFOOD_GOOD*alpha_size0_pgood0 + (1-pgood0)*survBad*PFOOD_BAD*alpha_size0_pgood0;		// P(+) in pete notation
					probSurvFail = pgood0*survGood*(1-(PFOOD_GOOD*alpha_size0_pgood0)) + (1-pgood0)*survBad*(1-(PFOOD_BAD*alpha_size0_pgood0));	// P(-) in pete notation
					
					// calculate value of information and value of growth					
					thisvalInfo = probSurvSucc*(interpolate2(RV_size0, pgoodSucc, alpha_size0_pgoodSucc) - interpolate2(RV_size0, pgoodSucc, alpha_size0_pgood0)) + probSurvFail*(interpolate2(RV_size0, pgoodFail, alpha_size0_pgoodFail) - interpolate2(RV_size0, pgoodFail, alpha_size0_pgood0));
				}
			}
			
			infoArray[size_val][p_val] = thisvalInfo;

		}
	}

	cout << "End value of information/growth calculation"  << endl;

}


void Output(double &PHI, double &KAPPA_GOOD, double &KAPPA_BAD, double &CONST_GOOD, double &CONST_BAD, double &MU0, double &MU1_GOOD, double &MU1_BAD, double &GAMMA, double &PFOOD_GOOD, double &PFOOD_BAD) 
{
	// save fitness, decision and information array
	for (size_val = 0; size_val < SIZESTEPS; size_val++)
	{  
		for (p_val = 0; p_val < PSTEPS; p_val++)
		{  
			fprintf(file1, "%6.7f ", old_value[size_val][p_val]);	
			fprintf(file2, "%6.7f ", decArray[size_val][p_val]);		
			fprintf(file0, "%10.10f ", infoArray[size_val][p_val]);	
		}
		fprintf(file1, "\n");
		fprintf(file2, "\n");
		fprintf(file0, "\n");
	}        
	
	// save convergence vector
	for (vector<double>::iterator vec=diffvec.begin(); vec<diffvec.end(); vec++)
		fprintf(file3, "%12.10f\n", *vec); 
	
	// save parameters
	fprintf(file4, "%f %f %f %f %f %f %f %f %f %f %f\n", PHI, KAPPA_GOOD, KAPPA_BAD, CONST_GOOD, CONST_BAD, MU0, MU1_GOOD, MU1_BAD, GAMMA, PFOOD_GOOD, PFOOD_BAD);
	
	// save fitness-by-alpha array
	for (size_val = 0; size_val < SIZESTEPS; size_val++)
	{  		
		for (p_val = 0; p_val < PSTEPS; p_val++)
		{
			for (alpha_i = 0; alpha_i < ALPHASTEPS; alpha_i++)
			{
				fprintf(file5, "%6.7f ", fitAlphaArray[size_val][p_val][alpha_i]);	
			}
			fprintf(file5,"\n");
		}
	}        
}

void PressEnterToContinue()
{
  int c;
  printf( "Press ENTER to continue... " );
  fflush( stdout );
  do c = getchar(); while ((c != '\n') && (c != EOF));
}

string usage_message() 
{
	return "Usage: -parno value -PHI value -KAPPA_GOOD value -KAPPA_BAD value -CONST_GOOD value -CONST_BAD value -MU0 value -MU1_GOOD value -MU1_BAD value -GAMMA value -PFOOD_GOOD value -PFOOD_BAD value.\n";
}


// Import parameters from the run code

//simple params values struct
struct starting_values {
  double parno;				// integer to relate to set of parameter combinations (for naming files)
  double PHI;	
	double KAPPA_GOOD;
	double KAPPA_BAD;
	double CONST_GOOD;
	double CONST_BAD;
	double MU0;
	double MU1_GOOD;
	double MU1_BAD;
	double GAMMA;
	double PFOOD_GOOD;
	double PFOOD_BAD;	
  //add new params here
} values; //declare a struct member

void get_params_from_args(int &argc, char* argv[], starting_values &values)
{
	// Check the value of argc. If not enough parameters have been passed, inform user and exit.
	// argv[0] = path and name of program
	// argv[1]..argv[n] are the params.
	// e.g. c:\prog.exe -em 1 -epmax 2 -another_param 3
	// argv[0] = c:\prog.exe
	// argv[1] = -em
	// argv[2] = 1
	// argv[3] = -epmax, etc
	
	// argc = number of parameters including name of program.
	if (argc != 25) {
		cout << usage_message(); // Inform the user of how to use the program
		PressEnterToContinue();
		exit(0);
	} else {
		//capture args from command prompt - ensure the correct switches are used
		double arg_parno = 0.0;				// integer to relate to set of parameter combinations (for naming files)
		double arg_PHI = 0.0;	
		double arg_KAPPA_GOOD = 0.0;
		double arg_KAPPA_BAD = 0.0;
		double arg_CONST_GOOD = 0.0;
		double arg_CONST_BAD = 0.0;
		double arg_MU0 = 0.0;
		double arg_MU1_GOOD = 0.0;
		double arg_MU1_BAD = 0.0;
		double arg_GAMMA = 0.0;
		double arg_PFOOD_GOOD = 0.0;
		double arg_PFOOD_BAD = 0.0;	
		
		// loop over all args and skip the values just check the flags
		for (int i = 1; i < argc; i+=2) {
			// Check that we haven't finished parsing already
			if (i + 1 != argc) {
				//convert the argv vars into their correct types, switches (-"par name") are strings, switch values are doubles.
				string curr_arg = string(argv[i]);
				double curr_arg_val = atof(argv[i+1]);
				//ensure that only known switches are used (fail if switch doesn't match these)
				if (curr_arg == "-parno") {
					values.parno = curr_arg_val;
				}				
				else if (curr_arg == "-PHI") {
					values.PHI = curr_arg_val;
				} 
				else if (curr_arg == "-KAPPA_GOOD") {
					values.KAPPA_GOOD = curr_arg_val;
				} 
				else if (curr_arg == "-KAPPA_BAD") {
					values.KAPPA_BAD = curr_arg_val;
				} 
				else if (curr_arg == "-CONST_GOOD") {
					values.CONST_GOOD = curr_arg_val;
				} 
				else if (curr_arg == "-CONST_BAD") {
					values.CONST_BAD = curr_arg_val;
				} 								
				else if (curr_arg == "-MU0") {
					values.MU0 = curr_arg_val;
				} 
				else if (curr_arg == "-MU1_GOOD") {
					values.MU1_GOOD = curr_arg_val;
				} 
				else if (curr_arg == "-MU1_BAD") {
					values.MU1_BAD = curr_arg_val;
				} 
				else if (curr_arg == "-GAMMA") {
					values.GAMMA = curr_arg_val;
				} 				
				else if (curr_arg == "-PFOOD_GOOD") {
					values.PFOOD_GOOD = curr_arg_val;
				} 
				else if (curr_arg == "-PFOOD_BAD") {
					values.PFOOD_BAD = curr_arg_val;
				} else {
					cout << "Not enough or invalid arguments." + usage_message();
					PressEnterToContinue();
					exit(0);
				}
			}
		}
	}
}



/***************************************************************************************** 
 //  main
 //
 //  Dummy return value.
 ******************************************************************************************/
int main(int argc, char* argv[])
{
	// set seed for random numbers
	int ltime = time(NULL);              			
	int utime = (unsigned int) ltime/2;
	srand(utime);
		
	//extract and assign the vars using a struct ptr to a ref of my starting_values struct called values.
	starting_values* param_values = &values;
	get_params_from_args(argc, argv, *param_values);
	
	// create file names	
 	char fullname1[50];
	char fullname2[50];
	char fullname3[50];
	char fullname4[50];
	char fullname5[50];
	char fullname0[50];

	sprintf(fullname1, "FitnessValue%1.0f.txt",param_values->parno);
	sprintf(fullname2, "Decision%1.0f.txt",param_values->parno);
	sprintf(fullname3, "ConvergeVector%1.0f.txt",param_values->parno);
	sprintf(fullname4, "Params%1.0f.txt",param_values->parno);	
	sprintf(fullname5, "FitnessAlpha%1.0f.txt",param_values->parno);	
	sprintf(fullname0, "InfoVal%1.0f.txt",param_values->parno);

	file1 = fopen(fullname1,"w");
	file2 = fopen(fullname2,"w");
	file3 = fopen(fullname3,"w");
	file4 = fopen(fullname4,"w");
	file5 = fopen(fullname5,"w");
	file0 = fopen(fullname0,"w");

	Initialise();
	Iterate_Ends(param_values->PHI, param_values->KAPPA_GOOD, param_values->KAPPA_BAD, param_values->CONST_GOOD, param_values->CONST_BAD, param_values->MU0, param_values->MU1_GOOD, param_values->MU1_BAD, 
							 param_values->GAMMA, param_values->PFOOD_GOOD, param_values->PFOOD_BAD);
	Iterate_Mid(param_values->PHI, param_values->KAPPA_GOOD, param_values->KAPPA_BAD, param_values->CONST_GOOD, param_values->CONST_BAD, param_values->MU0, param_values->MU1_GOOD, param_values->MU1_BAD, 
							param_values->GAMMA, param_values->PFOOD_GOOD, param_values->PFOOD_BAD);
	Iterate_Final(param_values->PHI, param_values->KAPPA_GOOD, param_values->KAPPA_BAD, param_values->CONST_GOOD, param_values->CONST_BAD, param_values->MU0, param_values->MU1_GOOD, param_values->MU1_BAD, 
							param_values->GAMMA, param_values->PFOOD_GOOD, param_values->PFOOD_BAD);
	Val_Info_Grow(param_values->PHI, param_values->KAPPA_GOOD, param_values->KAPPA_BAD, param_values->CONST_GOOD, param_values->CONST_BAD, param_values->MU0, param_values->MU1_GOOD, param_values->MU1_BAD, param_values->GAMMA, param_values->PFOOD_GOOD, param_values->PFOOD_BAD);	
	Output(param_values->PHI, param_values->KAPPA_GOOD, param_values->KAPPA_BAD, param_values->CONST_GOOD, param_values->CONST_BAD, param_values->MU0, param_values->MU1_GOOD, param_values->MU1_BAD, 
				 param_values->GAMMA, param_values->PFOOD_GOOD, param_values->PFOOD_BAD);
	
	fclose(file1);
	fclose(file2);
	fclose(file3);
	fclose(file4);
	fclose(file5);
	fclose(file0);

	return 0;
}

