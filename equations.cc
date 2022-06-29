#include <iostream>
using namespace std;
#include <fstream>
#include <vector>
// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for each
// function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

void variableAttributeLoader::loadVariableAttributes(){

    for (unsigned int var_index=0; var_index<17; var_index++){
        std::string var_name = "n";
        var_name.append(std::to_string(var_index));

        set_variable_name				(var_index,var_name);
    	set_variable_type				(var_index,SCALAR);
    	set_variable_equation_type		(var_index,EXPLICIT_TIME_DEPENDENT);

        set_dependencies_value_term_RHS(var_index, "n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16");
        set_dependencies_gradient_term_RHS(var_index, "grad(n0), grad(n1), grad(n2), grad(n3), grad(n4), grad(n5), grad(n6), grad(n7), grad(n8), grad(n9), grad(n10), grad(n11), grad(n12), grad(n13), grad(n14), grad(n15), grad(n16)");
	if (var_index < 16){
		set_allowed_to_nucleate			(var_index, true);
		set_need_value_nucleation		(var_index, true);
	}else{
		set_allowed_to_nucleate			(var_index, false);
		set_need_value_nucleation		(var_index, false);
	}
    }

	//set_variable_name				(8,"T");
	//set_variable_type				(8,SCALAR);
	//set_variable_equation_type		(8,EXPLICIT_TIME_DEPENDENT);

    	//set_dependencies_value_term_RHS(8, "T");
    	//set_dependencies_gradient_term_RHS(8, "grad(T)");
	//set_allowed_to_nucleate			(8, false);
        //set_need_value_nucleation		(8, false);

}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// --- Getting the values and derivatives of the model variables ---

dealii::VectorizedArray<double> fnV = constV(0.0);
scalarvalueType ni, nj;
scalargradType nix;

//scalarvalueType T = variable_list.get_scalar_value(8);
//scalargradType Tx = variable_list.get_scalar_gradient(8);

// In this application, create temporary variables for the residual terms. We cannot
// call 'set_scalar_value_residual_term' and 'set_scalar_gradient_residual_term'in the
// for loop below because those functions write over the scalar value and scalar gradient
// internal variables in 'variable_list' (for performance reasons). Therefore, we wait
// to set the residual terms until all the residuals have been calculated.

std::vector<scalarvalueType> value_terms;
value_terms.resize(userInputs.number_of_variables);
std::vector<scalargradType> gradient_terms;
gradient_terms.resize(userInputs.number_of_variables);

double T = Tliquidus - 0.1 - Tdiff*(this->currentTime/userInputs.finalTime) - (q_point_loc[0][0])*T_grad; //second term is slope
if (T < 300){
	//Mob_frac = 0.1;
	T = 300;
}

//cout << "lg: " << 1/lg << endl;
double Kg = ak*sigma_g*lg;
double mg = 0.75*sigma_g*(1/Delta_fg)*(1/lg);
double Lg = (1/(ak*lg))*D0*(exp(-Qg/(R*T)));

for (unsigned int i=0; i<userInputs.number_of_variables; i++){

	ni = variable_list.get_scalar_value(i);
	nix = variable_list.get_scalar_gradient(i);
	fnV = - ni + ni*ni*ni;
	//check if  -1 variables induces pinning.
	for (unsigned int j=0; j<userInputs.number_of_variables; j++){
		if (i != j){
			nj = variable_list.get_scalar_value(j);
			fnV += constV(2.0*gamma) * ni * nj*nj;
		}
	}

	if (i == 16){
		//manually set values for ODS Particles, no gradient so no evolution of ODS
		value_terms[16] = variable_list.get_scalar_value(16);
		gradient_terms[16] = constV(0.0);
	}else{
		//Lg = 0.5 Lg
		value_terms[i] = ni-constV(userInputs.dtValue*Lg*mg)*fnV;
		gradient_terms[i] = constV(-userInputs.dtValue*Kg*Lg)*nix;
		//if (nix[1][i] > 0.1){
		//	cout << "grad: " << nix[1][i] << endl;
		//}
	}

}

// --- Submitting the terms for the governing equations ---

// ---------------------------------------------------------------------
// Nucleation section to seed nuclei and freeze the Allen-Cahn mobility
// ---------------------------------------------------------------------
//probably a way to reduce all source_term doubles into single variable
dealii::AlignedVector<dealii::VectorizedArray<double > > source_terms(16,constV(0.0));
dealii::VectorizedArray<double> lambda = constV(1.0);
seedNucleus(q_point_loc,source_terms,lambda);
// ---------------------------------------------------------------------


for (unsigned int i=0; i<userInputs.number_of_variables; i++){
	if (i == 16){
		variable_list.set_scalar_value_term_RHS(i,value_terms[i]);
	}else{
		variable_list.set_scalar_value_term_RHS(i,value_terms[i]+source_terms[i]);
	}
	variable_list.set_scalar_gradient_term_RHS(i,gradient_terms[i]);

}

//scalarvalueType eq_T=T;
//scalargradType eqx_T=(constV(-D*userInputs.dtValue)*Tx);

//variable_list.set_scalar_value_term_RHS(8,eq_T);
//variable_list.set_scalar_gradient_term_RHS(8,eqx_T);

}

// =================================================================================
// seedNucleus: a function particular to this app
// =================================================================================
template <int dim,int degree>
void customPDE<dim,degree>::seedNucleus(const dealii::Point<dim, dealii::VectorizedArray<double> > & q_point_loc,
	dealii::AlignedVector<dealii::VectorizedArray<double> > & source_terms,
	dealii::VectorizedArray<double> & lambda) const {

		for (typename std::vector<nucleus<dim> >::const_iterator thisNucleus=this->nuclei.begin(); thisNucleus!=this->nuclei.end(); ++thisNucleus){
			if (thisNucleus->seededTime + thisNucleus->seedingTime > this->currentTime){
				// Calculate the weighted distance function to the order parameter freeze boundary (weighted_dist = 1.0 on that boundary)
				dealii::VectorizedArray<double> weighted_dist = this->weightedDistanceFromNucleusCenter(thisNucleus->center, userInputs.get_nucleus_freeze_semiaxes(thisNucleus->orderParameterIndex), q_point_loc, thisNucleus->orderParameterIndex);

				for (unsigned i=0; i<lambda.n_array_elements;i++){
					if (weighted_dist[i] <= 1.0){
						lambda[i] = 0.0;

						// Seed a nucleus if it was added to the list of nuclei this time step
						if (thisNucleus->seedingTimestep == this->currentIncrement){
							// Find the weighted distance to the outer edge of the nucleus and use it to calculate the order parameter source term
							dealii::Point<dim,double> q_point_loc_element;
							for (unsigned int j=0; j<dim; j++){
								q_point_loc_element(j) = q_point_loc(j)[i];
							}
							double r = this->weightedDistanceFromNucleusCenter(thisNucleus->center, userInputs.get_nucleus_semiaxes(thisNucleus->orderParameterIndex), q_point_loc_element, thisNucleus->orderParameterIndex);

							double avg_semiaxis = 0.0;
							for (unsigned int j=0; j<dim; j++){
								avg_semiaxis += thisNucleus->semiaxes[j];
							}
							avg_semiaxis /= dim;

							if (thisNucleus->orderParameterIndex == 0){
								//was interface_coeff instead of Lp
								source_terms[0][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 1){
								source_terms[1][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 2){
								source_terms[2][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 3){
								source_terms[3][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 4){
								source_terms[4][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 5){
								source_terms[5][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 6){
								source_terms[6][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 7){
								source_terms[7][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 8){
								source_terms[8][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 9){
								source_terms[9][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 10){
								source_terms[10][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 11){
								source_terms[11][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 12){
								source_terms[12][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 13){
								source_terms[13][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 14){
								source_terms[14][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else if (thisNucleus->orderParameterIndex == 15){
								source_terms[15][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							else {
								//source_terms[2][i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/Lp));
							}
							//source_term[i] =0.5*(1.0-std::tanh(avg_semiaxis*(r-1.0)/interface_coeff));

						}
					}
				}
			}
		}
	}



// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are not
// explicit time-dependent equations. It takes "variable_list" as an input, which is
// a list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one proportional
// to the test function and one proportional to the gradient of the test function. The
// index for each variable in this list corresponds to the index given at the top of
// this file.

template <int dim, int degree>
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of that
// quadrature point is given by "q_point_loc". The function outputs two terms to
// variable_list -- one proportional to the test function and one proportional to the
// gradient of the test function -- for the left-hand-side of the equation. The index
// for each variable in this list corresponds to the index given at the top of this
// file. If there are multiple elliptic equations, conditional statements should be
// sed to ensure that the correct residual is being submitted. The index of the field
// being solved can be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {
}
