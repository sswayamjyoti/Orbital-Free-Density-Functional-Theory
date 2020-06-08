double NegativeCharge; 
double Ncore;
int int_order;
int* boundary_nodes;
int total_count;
double penalty;
double tol_penalty;
double tol_lagrange;
double pi = 3.14159265358979323846;
double lam_kin;
double gam = -0.1471;
double beta1 = 1.1581;
double beta2 = 0.3446;
double A = 0.0311;
double B = -0.048;
double C = 0.0014;
double D = -0.0108;

struct ATOM
	{
		Dune::FieldVector<double,3> pos;
		int charge;
	};

std::list<ATOM> atoms;
