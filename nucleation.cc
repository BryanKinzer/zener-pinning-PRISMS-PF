#include <iostream>
using namespace std;
double cum_prob {0.0};
// =================================================================================
// NUCLEATION FUNCTIONS
// =================================================================================

// =================================================================================
// Nucleation probability
// =================================================================================
template <int dim, int degree>
double customPDE<dim,degree>::getNucleationProbability(variableValueContainer variable_value, double dV, dealii::Point<dim> p, unsigned int variable_index) const
{
      
    //unsigned int num_ODS = 16;
    double x_diff = 0.1;
    double y_diff = 0.1;
    double d2_diff = 0.1;
    double buffer = rad + 0.01;
    unsigned int close_nuc = 0;
    //2.0e-4 is the grid size from parameters.prm, adjust accordingly
    for (unsigned int i = 0; i < num_ODS; i++){
	buffer = rad + 0.01;
	x_diff = std::abs(x_coord[i] - (p[0]/(userInputs.domain_size[0])));
	y_diff = std::abs(y_coord[i] - (p[1]/(userInputs.domain_size[1])));
	d2_diff = x_diff*x_diff + y_diff*y_diff;
	if (d2_diff > (rad*rad) && d2_diff < (buffer*buffer)){
		close_nuc = 1;
	}
    }

	// Calculate the nucleation rate
	double retProb;
	//double T = Tliquidus - 0.1 - Tdiff*(this->currentTime/userInputs.finalTime); //second term is slope
	double T = Tliquidus - 0.1 - Tdiff*(this->currentTime/userInputs.finalTime) - (p[0])*T_grad; 
	if (T < 300){
		T = 300;
	}
	double Del_T = Tliquidus - T;
	double delG_V = L_m*Del_T/Tliquidus;
	double sigma_p = sigma_g;
	//double DelG_V_hom = (3.1416)*(sigma_p*sigma_p)*(40e-9)/(delG_V);  //for a 40 nm tall cylinder
	double DelG_V_hom = (16*3.1416/3)*(sigma_p*sigma_p*sigma_p)/(delG_V *delG_V); //for 3D, what I have used before.
	double f_theta = 4.447e-6; //for 4 degrees
	double k = 1.381e-23;
	double h = 6.626e-34;
	double Qd = Qg;
	double J_Prob = T*N_sites*(k/h)*exp(-Qd/(R*T))*exp(-DelG_V_hom*f_theta/(k*T));

	
    //heteregenous nucleation can occur within 0.01 of the mold wall or within 1 particle radius of a nanoparticle
    if (p[0] > 0.99*userInputs.domain_size[0] || close_nuc == 1){
	p[0] = userInputs.domain_size[0];
	//cout << "dV: " << userInputs.refine_factor << endl;
	retProb=1.0-exp(-J_Prob*userInputs.dtValue*((double)userInputs.steps_between_nucleation_attempts)*(5.9605e-23)); //Poission Distribution Approach for messh size of 5*10^-6/2^7
	//retProb=1.0-exp(-J_Prob*userInputs.dtValue*((double)userInputs.steps_between_nucleation_attempts)*(7.4506e-24)); //Poission Distribution Approach for mesh size of 5*10^-6/2^8
	//(5.9605e-23) is hard coded for 3D element size with refine factor of 7 and 5 um x 5 um domain uniform. dV only gives 2D, but we are assuming 3D nucleation of a slice of a sphere a single layer tall.
    }else{
	retProb = 0.0;
    }
    cum_prob = cum_prob + retProb;
    if (cum_prob > 200.0){
	retProb= 0.0;
    }else{
    }

    return retProb;
}
