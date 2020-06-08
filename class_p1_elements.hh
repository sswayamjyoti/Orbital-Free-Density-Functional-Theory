template<class GV>
class P1Elements
{
public:
    static const int dim = GV::dimension;

    typedef typename GV::ctype ctype;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Matrix;
		//typedef Dune::FieldMatrix<double,1,1>  Matrix;
    typedef Dune::BlockVector<Dune::FieldVector<double,1>> ScalarField;


public:
    Matrix PoissonMatrix;
		Matrix poisson_diagonal;
		Matrix ThomasFermiMatrix;
		Matrix MassMatrix;
		Matrix TF_diagonal;
    ScalarField b;
		ScalarField b1;
		ScalarField u;
		ScalarField u1;
		ScalarField gradLag;
    std::vector< std::set<int> > adjacencyPattern;
		std::vector< std::set<int> > adjacencyPattern_u_i;

P1Elements(const GV& gv_, Dune::BlockVector<Dune::FieldVector<double,1> >& u_i_, Dune::BlockVector<Dune::FieldVector<double,1> >& phi_i_, Dune::BlockVector<Dune::FieldVector<double,1>>& atom_node_, double& integral_u_, Dune::BlockVector<Dune::FieldVector<double,1> >& gradL_, double& lagrange_, double& penalty_) : gv(gv_), u_i(u_i_), phi_i(phi_i_), atom_node(atom_node_), integral_u(integral_u_), gradL(gradL_), lagrange(lagrange_), penalty(penalty_){}


    // store adjacency information in a vector of sets
    void determineAdjacencyPattern();

		void determineAdjacencyPattern_u_i();

    // assemble stiffness matrix A and right side b
    void assemble_poisson();

		void assemble_poisson_rhs();

		void assemble_mass();

		void assemble_u_i();

		void assemble_TF_diagonal();

    // finally solve Au = b for u
    void solve();

		void solve_u_i();


private:
    typedef typename GV::template Codim<0>::Iterator LeafIterator;
    typedef typename GV::IntersectionIterator IntersectionIterator;
    typedef typename GV::IndexSet LeafIndexSet;

    const GV& gv;
		//Dune::BlockVector<Dune::FieldVector<double,1>>& u_i;
		Dune::BlockVector<Dune::FieldVector<double,1>>& u_i;
		Dune::BlockVector<Dune::FieldVector<double,1>>& phi_i;
		//const int& atom_node;
		Dune::BlockVector<Dune::FieldVector<double,1>>& atom_node;
		//Dune::BlockVector<Dune::FieldVector<double,1>>& gradL;
		double& integral_u;
		Dune::BlockVector<Dune::FieldVector<double,1>>& gradL;
		double& lagrange;
		double& penalty;
   // const F& f;
};

