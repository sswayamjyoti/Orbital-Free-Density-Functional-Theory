template<class GV>
void penalty_gradLphi(const GV& gv, Dune::BlockVector<Dune::FieldVector<double,1>>& gradLphi1, Dune::BlockVector<Dune::FieldVector<double,1>>& u_i, Dune::BlockVector<Dune::FieldVector<double,1>>& phi_i, Dune::BlockVector<Dune::FieldVector<double,1>>& atom_node, double side_length, double penalty, double integral_u)
{
	static const int dim = 3;
	typedef Dune::UGGrid<dim>::ctype ctype;
	Dune::BlockVector<Dune::FieldVector<double,1>> gradN(dim); 
	Dune::BlockVector<Dune::FieldVector<double,1>> gradu(dim);
	Dune::BlockVector<Dune::FieldVector<double,1>> gradphi(dim);
	Dune::BlockVector<Dune::FieldVector<double,1>> x(gv.size(dim));
	Dune::BlockVector<Dune::FieldVector<double,1>> gradL(gv.size(dim));

	double C_F, Ts1, Ts2, Ts, L1, L2, ES1, ES2, Exc, u, phi, sum, sum1, eps_x, eps_c, rs, lagrange;
	const P1ShapeFunctionSet<ctype,ctype,dim> &basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

	typedef typename GV::template Codim<dim>::Iterator VertexLeafIterator;
	typedef typename GV::template Codim<0>::Iterator ElementLeafIterator;
	typedef typename GV::IntersectionIterator IntersectionIterator;
	typedef typename GV::IndexSet LeafIndexSet;
	const LeafIndexSet &set = gv.indexSet();
	const ElementLeafIterator elend = gv.template end<0>();
	const ElementLeafIterator itend = gv.template end<0>();
	gradLphi1 = 0.0;

	C_F = 0.3*pow((3*pi*pi),0.66666666666);
			
	// Iteration over elements to compute the gradLphi1
	Dune::Timer watch_GridIter;
	watch_GridIter.reset();
		
	//	gradL = 1.0;
		penalty = 1.0; //set arbitrarily
		P1Elements<GV> p1(gv, u_i, phi_i, atom_node, integral_u, gradL, lagrange, penalty);
    p1.determineAdjacencyPattern();
    p1.assemble_poisson(); 
		p1.assemble_poisson_rhs();
	  (p1.PoissonMatrix).mv(phi_i,x);		

/*	for (ElementLeafIterator el = gv.template begin<0>(); el != elend; ++el)
	{
		// Determine geometry type of the current element and get the matching reference element
		Dune::GeometryType gt = el->type();
		const Dune::template GenericReferenceElement<double,dim> &ref =  Dune::GenericReferenceElements<double,dim>::general(gt);
		int vertexsize = ref.size(dim);
		 
		// Get a quadrature rule of order one for the given geometry type
		const Dune::QuadratureRule<ctype,dim>& rule2 = Dune::QuadratureRules<ctype,dim>::rule(gt,int_order);
		for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule2.begin();r != rule2.end() ; ++r)
		{
			Dune::FieldMatrix<ctype,dim,dim> jacInvTra = el->geometry().jacobianInverseTransposed(r->position());		    
			double weight = r->weight();		   
			double detjac = el->geometry().integrationElement(r->position());
			u = 0.0;
			phi = 0.0;
			gradu = 0.0;
			gradphi = 0.0;

			for (int i = 0 ; i<vertexsize; i++)
			{
				sum = basis[i].evaluateFunction(r->position())*u_i[set.subIndex(*el,i,dim)];
				sum1 = basis[i].evaluateFunction(r->position())*phi_i[set.subIndex(*el,i,dim)];
				u  += sum;
				phi += sum1;

				// Calculation of grad u
				jacInvTra.mv(basis[i].evaluateGradient(r->position()),gradN);
				gradN *= u_i[set.subIndex(*el,i,dim)];
				gradu += gradN;

			        // Calculation of gradphi

			        jacInvTra.mv(basis[i].evaluateGradient(r->position()),gradN);
			        gradN *= phi_i[set.subIndex(*el,i,dim)];
			        gradphi += gradN;						
			}		     
			 			       		    
			for (int i = 0 ; i<vertexsize; i++)
			{
			jacInvTra.mv(basis[i].evaluateGradient(r->position()),gradN);		     			
			gradLphi1[set.subIndex(*el,i,dim)] += ( (-1/(4*pi))*(gradphi*gradN)*weight*detjac + (pow(u,2)*basis[i].evaluateFunction(r->position()))*weight*detjac);
			}

				rs = pow((3/(4*pi*u*u)),0.33333333333333);
   			eps_x = (-3./4.)*pow((3./pi),0.3333333)*(pow(pow(u,2),0.3333333333)); 
                       
				if(rs > 1.0 || rs == 1.0)
				{
					eps_c = gam/(1 + (beta1*pow(rs,0.5)) + (beta2*rs));
				}
				else
				{
					eps_c = A*log(rs) + B + C*rs*log(rs) + D*rs;
				}

				if (u < pow(10, -10))
				{
					eps_c = 0.0;
				}
		}		
	} */ 
	
 // gradLphi1[atom_node] += -Ncore;
		
	gradLphi1  = x;
	gradLphi1 += p1.b;
	//gradLphi1 += atom_node;

	// Set gradient to zero in the boundary	
	watch_GridIter.reset();
	for (int count = 0; count < total_count; count++)
	{
		gradLphi1[boundary_nodes[count]] = 0.0;
	}

	  /*for ( ElementLeafIterator it = gv.template begin<0>() ; it != itend ; ++it) 
    {
        const IntersectionIterator isend = gv.iend(*it);
        for (IntersectionIterator is = gv.ibegin(*it) ; is != isend ; ++is)
        {
            // determine geometry type of the current element and get the matching reference element
            Dune::GeometryType gt = it->type();
            const Dune::template GenericReferenceElement<double,dim> &ref =
                Dune::GenericReferenceElements<double,dim>::general(gt);

            // check whether current intersection is on the boundary
            if ( is->boundary() )
            {
                // traverse all vertices the intersection consists of
                for (int i=0; i < ref.size(is->indexInInside(),1,dim); i++)
                {
                    // and replace the associated line of A and b with a trivial one
                    int indexi = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,i,dim),dim);

                    gradLphi1[indexi] = 0.0; 
                    
                }
            }
        }
    } */ 

}
