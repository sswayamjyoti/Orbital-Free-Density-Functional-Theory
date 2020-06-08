template<class GV>
void P1Elements<GV>::assemble_mass()
{
    const int N = gv.size(dim);
		double sum, u;

    const LeafIndexSet& set = gv.indexSet();
    const LeafIterator itend = gv.template end<0>();

    // set sizes of A and b
    MassMatrix.setSize(N, N, (N + 2*gv.size(dim-1)));
    MassMatrix.setBuildMode(Matrix::random);

		for (int i = 0; i < N; i++)
        MassMatrix.setrowsize(i,adjacencyPattern[i].size());
    MassMatrix.endrowsizes();

    // set sparsity pattern of A with the information gained in determineAdjacencyPattern 
    for (int i = 0; i < N; i++) 
    {
        std::template set<int>::iterator setend = adjacencyPattern[i].end();
        for (std::template set<int>::iterator setit = adjacencyPattern[i].begin();setit != setend; ++setit)
            MassMatrix.addindex(i,*setit);
    }

    MassMatrix.endindices(); 
		
    // initialize A and b
    MassMatrix = 0.0;

    // get a set of P1 shape functions
    const P1ShapeFunctionSet<ctype,ctype,dim>& basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

    for (LeafIterator it = gv.template begin<0>(); it != itend; ++it) 
    {
        // determine geometry type of the current element and get the matching reference element
        Dune::GeometryType gt = it->type();
        const Dune::template GenericReferenceElement<double,dim> &ref =
        Dune::GenericReferenceElements<double,dim>::general(gt);
        int vertexsize = ref.size(dim);

        // get a quadrature rule of order one for the given geometry type
        const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,int_order);
        for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
             r != rule.end() ; ++r)
        {
            // get the weight at the current quadrature point
            double weight = r->weight();

            // compute Jacobian determinant for the transformation formula
            double detjac = it->geometry().integrationElement(r->position());
            for (int i = 0; i < vertexsize; i++)
            {
                for (int j = 0; j < vertexsize; j++)
                {
                    // gain global inidices of vertices i and j and update associated matrix entry
                    MassMatrix[set.subIndex(*it,i,dim)][set.subIndex(*it,j,dim)] += (basis[i].evaluateFunction(r->position()))*(basis[j].evaluateFunction(r->position()))*weight*detjac; 
										
                }
            }
        } 
    }
 }

