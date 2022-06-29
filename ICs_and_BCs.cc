#include <iostream>
using namespace std;
#include <fstream>
#include <vector>
double x_coord [1280];
double y_coord [1280];
// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

int first_time {1};

template <int dim, int degree>
void customPDE<dim,degree>::setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC){
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index

	  // The initial condition is a set of overlapping circles/spheres defined
	  // by a hyperbolic tangent function. The center of each circle/sphere is
	  // given by "center" and its radius is given by "radius".
    if (index == 16){
	std::vector<dealii::Point<dim>> center;
	          // The big grains 	//0.2, 0.15 and 0.25, 0.7

	std::ifstream fin("10p0_uniform.csv");
	//std::ifstream fin("10p0_mod-------.csv");
    	std::vector<double> numbers;
	

    	double num;
    	numbers.push_back(num);  // store the number
    	while (fin >> num)           // read a number
    	{
		numbers.push_back(num);  // store the number 
    	}
    	// Separate into even and odd
    	double temp_Coord = 0;
    	for (unsigned int i = 0; i < numbers.size(); i++) {
        	temp_Coord = numbers[i+1];
		if (i % 2 == 0) {
            		x_coord[(i/2)] = temp_Coord;
        	} else {
			y_coord[((i-1)/2)] = temp_Coord;
        	}
        	//cout << "Coord: " << temp_Coord << endl;
    	}
	//x_coord[0] = numbers[1];
	//manually set single particle in center of domain
	x_coord[0] = 0.5;
	y_coord[0] = 0.5;
	//unsigned int num_ODS = 16;
	for (unsigned int i = 0; i < num_ODS; i++){
		{dealii::Point<dim> p(x_coord[i], y_coord[i]); center.push_back(p);}
		//{dealii::Point<dim> p(0.015, 0.015); center.push_back(p);}
		//cout << "x" << i << ": " << x_coord[i] << endl;
		//cout << "y" << i << ": " << y_coord[i] << endl;
	}
	
	//double rad = 0.01;
	double dist = 0.0;
    	scalar_IC = 0;
	double init_factor = theta;
//1e6 = init_factor
	//double init_factor = 1.0e5; makes grains larger when init_factor decreases
	for (unsigned int i = 0; i < num_ODS; i++){
		dist = 0.0;
		for (unsigned int dir = 0; dir < dim; dir++){
    			dist += (p[dir]-center[i][dir]*userInputs.domain_size[dir])*(p[dir]-center[i][dir]*userInputs.domain_size[dir]);
    	  	}
    	  	dist = std::sqrt(dist);
    	  	     scalar_IC +=	0.5*(1.0-std::tanh(init_factor*(dist-rad*userInputs.domain_size[0])/0.5));
	}



	
    }

	  // --------------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC)
{
    // --------------------------------------------------------------------------
    // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
    // --------------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the boundary condition for each variable
    // according to its variable index. This function can be left blank if there
    // are no non-uniform Dirichlet boundary conditions. For BCs that change in
    // time, you can access the current time through the variable "time". The
    // boundary index can be accessed via the variable "direction", which starts
    // at zero and uses the same order as the BC specification in parameters.in
    // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).


    // -------------------------------------------------------------------------

}
