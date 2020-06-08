
template<class GV>
void P1Elements<GV>::assemble_poisson()
{
    const int N = gv.size(dim);
		double sum, u;

    const LeafIndexSet& set = gv.indexSet();
    const LeafIterator itend = gv.template end<0>();

    // set sizes of A and b
    PoissonMatrix.setSize(N, N, (N + 2*gv.size(dim-1)));
    PoissonMatrix.setBuildMode(Matrix::random);
    b.resize(N, false);

		poisson_diagonal.setSize(N, N, (N + 2*gv.size(dim-1)));
    poisson_diagonal.setBuildMode(Matrix::random);

		for (int i = 0; i < N; i++)
        PoissonMatrix.setrowsize(i,adjacencyPattern[i].size());
    PoissonMatrix.endrowsizes();

    // set sparsity pattern of A with the information gained in determineAdjacencyPattern 
    for (int i = 0; i < N; i++) 
    {
        std::template set<int>::iterator setend = adjacencyPattern[i].end();
        for (std::template set<int>::iterator setit = adjacencyPattern[i].begin();setit != setend; ++setit)
            PoissonMatrix.addindex(i,*setit);
    }

    PoissonMatrix.endindices(); 
		
    PoissonMatrix = 0.0;

    // Poisson Diagonal (Preconditioner)

		for (int i = 0; i < N; i++)
        poisson_diagonal.setrowsize(i,adjacencyPattern[i].size());
    poisson_diagonal.endrowsizes();

    // set sparsity pattern of A with the information gained in determineAdjacencyPattern 
    for (int i = 0; i < N; i++) 
    {
        std::template set<int>::iterator setend = adjacencyPattern[i].end();
        for (std::template set<int>::iterator setit = adjacencyPattern[i].begin();setit != setend; ++setit)
            poisson_diagonal.addindex(i,*setit);
    }

    poisson_diagonal.endindices(); 
		
    


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
            // compute the jacobian inverse transposed to transform the gradients
            Dune::FieldMatrix<double,dim,dim> jacInvTra = it->geometry().jacobianInverseTransposed(r->position());

            // get the weight at the current quadrature point
            double weight = r->weight();

            // compute Jacobian determinant for the transformation formula
            double detjac = it->geometry().integrationElement(r->position());
            for (int i = 0; i < vertexsize; i++)
            {
                // compute transformed gradients
                Dune::FieldVector<double,dim> grad1;
                jacInvTra.mv(basis[i].evaluateGradient(r->position()),grad1);
                for (int j = 0; j < vertexsize; j++)
                {
                    Dune::FieldVector<double,dim> grad2;
                    jacInvTra.mv(basis[j].evaluateGradient(r->position()),grad2);

                    // gain global inidices of vertices i and j and update associated matrix entry
                    PoissonMatrix[set.subIndex(*it,i,dim)][set.subIndex(*it,j,dim)] += (-1/(4*pi))*(grad1*grad2)* weight*detjac; 
										
                }
            }
        } 

        // get a quadrature rule of order two for the given geometry type
     /*  const Dune::QuadratureRule<ctype,dim>& rule2 = Dune::QuadratureRules<ctype,dim>::rule(gt,int_order);
        for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule2.begin();
             r != rule2.end() ; ++r)
        {
            double weight = r->weight();
            double detjac = it->geometry().integrationElement(r->position());
						
						u = 0.0;
						for (int i = 0 ; i<vertexsize; i++)
						{
							sum = basis[i].evaluateFunction(r->position())*u_i[set.subIndex(*it,i,dim)];
			    		u  += sum;
						}
						
            for (int i = 0 ; i<vertexsize; i++)
            {
                // evaluate the integrand of the right side
                double fval = basis[i].evaluateFunction(r->position()) * pow(u,2); 
                b[set.subIndex(*it,i,dim)] += fval * weight * detjac; 
            }
        } */
    }

	  //	b[atom_node] += -Ncore; 
			
		//	b += atom_node;

    // Dirichlet boundary conditions:
    // replace lines in A related to Dirichlet vertices by trivial lines
    for ( LeafIterator it = gv.template begin<0>() ; it != itend ; ++it) 
    {
        const IntersectionIterator isend = gv.iend(*it);
        for (IntersectionIterator is = gv.ibegin(*it) ; is != isend ; ++is)
        {
            // determine geometry type of the current element and get the matching reference element
            Dune::GeometryType gt = it->type();
            const Dune::template GenericReferenceElement<double,dim> &ref =
            Dune::GenericReferenceElements<double,dim>::general(gt);

            // check whether current intersection is on the boundary
            if (is->boundary())
            {
                // traverse all vertices the intersection consists of
                for (int i=0; i < ref.size(is->indexInInside(),1,dim); i++)
                {
                    // and replace the associated line of A and b with a trivial one
                    int indexi = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,i,dim),dim);

                    PoissonMatrix[indexi] = 0.0; 
                    PoissonMatrix[indexi][indexi] = 1.0;
                   // b[indexi] = 0.0;
                }
            }
        }
    }  

		poisson_diagonal = 0.0;

		 for (int i=0;i<N;i++)
		 {
				poisson_diagonal[i][i] = 1/(PoissonMatrix[i][i]);
		 }
}
