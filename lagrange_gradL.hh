template<class GV>
void lagrange_gradL(const GV& gv, Dune::BlockVector<Dune::FieldVector<double,1>>& gradL, Dune::BlockVector<Dune::FieldVector<double,1>>& u_i, Dune::BlockVector<Dune::FieldVector<double,1>>& phi_i, Dune::BlockVector<Dune::FieldVector<double,1>>& atom_node, double side_length, double lagrange, double integral_u)
{
	static const int dim = 3;
	typedef Dune::UGGrid<dim>::ctype ctype;
	Dune::FieldVector<double,3> qp;
	
	const int N = gv.size(dim);
	Dune::BlockVector<Dune::FieldVector<double,1>> gradN(dim), gradNi(dim), gradu(dim), gradphi(dim);	
	double C_F, Ts1, Ts2, Ts, L1, L2, ES1, ES2, Exc, u, phi, eps_x, gradeps_x, sum, sum1, rs, deriv_rs, eps_c, gradeps_c, ep_interaction;

	const P1ShapeFunctionSet<ctype,ctype,dim> &basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

	typedef typename GV::template Codim<dim>::Iterator VertexLeafIterator;
	typedef typename GV::template Codim<0>::Iterator ElementLeafIterator;
	typedef typename GV::IntersectionIterator IntersectionIterator;
	

	typedef typename GV::IndexSet LeafIndexSet;
	const LeafIndexSet &set = gv.indexSet();
	const ElementLeafIterator elend = gv.template end<0>();
	const ElementLeafIterator itend = gv.template end<0>();
	

	C_F = 0.3*pow((3*pi*pi),0.66666666666);						
	gradL = 0.0;

	// Iteration over elements
	for (ElementLeafIterator el = gv.template begin<0>(); el != elend; ++el)
	{				
		// determine geometry type of the current element and get the matching reference element
		Dune::GeometryType gt = el->type();
		const Dune::template GenericReferenceElement<double,dim> &ref =  Dune::GenericReferenceElements<double,dim>::general(gt);
		int vertexsize = ref.size(dim);

		// get a quadrature rule of order one for the given geometry type
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
				// Calculation of u and phi
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

				rs = pow((3/(4*pi*u*u)),0.33333333333333);
				deriv_rs = (pow((3/(4*pi)),0.3333333333333333))*(-2./3.)*((pow(u,-3))*(pow(pow(u,4),0.333333333333)));	  
   			eps_x = (-3./4.)*pow((3./pi),0.3333333)*(pow(pow(u,2),0.3333333333)); 
                       
				if(rs > 1.0 || rs == 1.0)
				{
					eps_c = gam/(1 + (beta1*pow(rs,0.5)) + (beta2*rs));
				}
				else
				{
					eps_c = A*log(rs) + B + C*rs*log(rs) + D*rs;
				}
        // Gradient of the Exchange-Coorelation Term
  			gradeps_x = (-0.5)*pow((3./pi),0.33333333)*u*(pow(pow(u,-4),0.3333333333));
				if(rs > 1.0 || rs == 1.0)
				{
					gradeps_c = (-gam/(pow((1 + beta1*pow(rs,0.5) + beta2*rs),2)))*((beta1/(2*pow(rs,0.5))) + beta2)*(deriv_rs);
				}
				else
				{
					gradeps_c = ((A/rs) + C + C*log(rs) + D)*(deriv_rs);
				}

					if (u < pow(10, -10))
				{
					eps_c = 0.0;
					gradeps_x = 0.0;
					gradeps_c = 0.0;
				}

			qp = el->geometry().global(r->position());
			  ep_interaction = POT_EL_ATOMS(atoms, qp);
							    
			for (int i = 0 ; i<vertexsize; i++)
			{
				// Compute gradient of kinetic energy terms
				jacInvTra.mv(basis[i].evaluateGradient(r->position()),gradNi);	  
   			
				gradL[set.subIndex(*el,i,dim)] += (C_F*(10./3.)*(((pow(pow(u,4),0.3333333333))*u)*basis[i].evaluateFunction(r->position()))
																		 + lam_kin*(gradu*gradNi)
																		 + ((gradeps_x)*u + (gradeps_c)*pow(u,2) + 2*(eps_x + eps_c)*u)*basis[i].evaluateFunction(r->position())
																		 + (2.0*u*phi*basis[i].evaluateFunction(r->position()))
																		 /*+ ((phi - ep_interaction)*u*basis[i].evaluateFunction(r->position())) */
																		 + 2.0*lagrange*u*basis[i].evaluateFunction(r->position()))*weight*detjac; 
			}
		}	
	}
	
	for (int count = 0; count < total_count; count++)
	{
		gradL[boundary_nodes[count]] = 0.0;	
	} 

	 /* for ( ElementLeafIterator it = gv.template begin<0>() ; it != itend ; ++it) 
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

                    gradL[indexi] = 0.0; 
                    
                }
            }
        }
    } */ 
			
	//return(gradL);
}
