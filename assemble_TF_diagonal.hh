template<class GV>
void P1Elements<GV>::assemble_TF_diagonal()
{
    const int N = gv.size(dim);
		double sum, sum1, u, phi, C_F, eps_x, gradeps_x, grad2eps_x, rs, deriv_rs, deriv2_rs, eps_c, gradeps_c, grad2eps_c, ep_interaction;
		C_F = 0.3*pow((3*pi*pi),0.66666666666);
			Dune::FieldVector<double,3> qp;

    const LeafIndexSet& set = gv.indexSet();
    const LeafIterator itend = gv.template end<0>();

		Dune::BlockVector<Dune::FieldVector<double,1>> mass_u_i(N);
		
		TF_diagonal.setSize(N, N, (N + 2*gv.size(dim-1)));
    TF_diagonal.setBuildMode(Matrix::random);
   
		for (int i = 0; i < N; i++)
        TF_diagonal.setrowsize(i,adjacencyPattern[i].size());
    TF_diagonal.endrowsizes();

    // set sparsity pattern of A with the information gained in determineAdjacencyPattern 
    for (int i = 0; i < N; i++) 
    {
        std::template set<int>::iterator setend = adjacencyPattern[i].end();
        for (std::template set<int>::iterator setit = adjacencyPattern[i].begin();
             setit != setend; ++setit)
            TF_diagonal.addindex(i,*setit);
    }

    TF_diagonal.endindices(); 
		
    // initialize A and b
    TF_diagonal = 0.0;

		MassMatrix.mv(u_i,mass_u_i);	
		double integral_u = u_i*mass_u_i;

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
        const Dune::QuadratureRule<ctype,dim>& rule = Dune::QuadratureRules<ctype,dim>::rule(gt,1);
        for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule.begin();
             r != rule.end() ; ++r)
        {
            // compute the jacobian inverse transposed to transform the gradients
            Dune::FieldMatrix<double,dim,dim> jacInvTra =
                it->geometry().jacobianInverseTransposed(r->position());

            // get the weight at the current quadrature point
            double weight = r->weight();

            // compute Jacobian determinant for the transformation formula
            double detjac = it->geometry().integrationElement(r->position());

						u = 0.0;
						phi = 0.0;
						for (int i = 0 ; i<vertexsize; i++)
						{
							sum = basis[i].evaluateFunction(r->position())*u_i[set.subIndex(*it,i,dim)];
							sum1 = basis[i].evaluateFunction(r->position())*phi_i[set.subIndex(*it,i,dim)];
			    		u  += sum;
							phi += sum1;
						}

						rs = pow((3/(4*pi*u*u)),0.33333333333333);
				    deriv_rs = (pow((3/(4*pi)),0.3333333333333333))*(-2./3.)*pow(std::abs(u),-1.666667);//((pow(u,-3))*(pow(pow(u,4),0.333333333333)));	  
						deriv2_rs = (pow((3/(4*pi)),0.3333333333333333))*(10./9.)*pow(std::abs(u),1.33333333);//*(pow(pow(u,4),0.333333333333));
						
						eps_x = (-3./4.)*pow((3./pi),0.3333333)*pow(std::abs(u),0.66666667);//*(pow(pow(u,2),0.3333333333));                        
  					gradeps_x = (-0.5)*pow((3./pi),0.33333333)*pow(std::abs(u),-0.33333333);//*u*(pow(pow(u,-4),0.3333333333));
						grad2eps_x = (1./6.)*pow((3./pi),0.33333333)*pow(std::abs(u),-1.3333333);//(pow(pow(u,-4),0.3333333333));

						if(rs > 1.0 || rs == 1.0)
						{
							eps_c = gam/(1 + (beta1*pow(rs,0.5)) + (beta2*rs));
						}
						else
						{
							eps_c = A*log(rs) + B + C*rs*log(rs) + D*rs;
						}

						if(rs > 1.0 || rs == 1.0)
						{
							gradeps_c = (-gam/(pow((1 + beta1*pow(rs,0.5) + beta2*rs),2)))*((beta1/(2*pow(rs,0.5))) + beta2)*(deriv_rs);
						}
						else
						{
							gradeps_c = ((A/rs) + C + C*log(rs) + D)*(deriv_rs);
						}

						if(rs > 1.0 || rs == 1.0)
						{
							grad2eps_c = (-gam/(pow((1 + beta1*pow(rs,0.5) + beta2*rs),2)))*((beta1/(2*pow(rs,0.5))) + beta2)*(deriv2_rs)
													 + (-gam/(pow((1 + beta1*pow(rs,0.5) + beta2*rs),2)))*(-beta1/(4*pow(rs,1.5)))*pow(deriv_rs,2)
													 + (2*gam/(pow((1 + beta1*pow(rs,0.5) + beta2*rs),3)))*(pow(((beta1/(2*pow(rs,0.5))) + beta2),2))*(pow(deriv_rs,2));
						}
						else
						{
							grad2eps_c = ((A/rs) + C + C*log(rs) + D)*(deriv2_rs) + ((-A/(pow(rs,2))) + (C/rs))*pow(deriv_rs,2);
						}

						if (u < pow(10, -10))
						{
							eps_c = 0.0;
							eps_x = 0.0;
							gradeps_c = 0.0;
							gradeps_x = 0.0;
							grad2eps_x = 0.0;
							grad2eps_c = 0.0;
						}

							qp = it->geometry().global(r->position());
			 			  ep_interaction = POT_EL_ATOMS(atoms, qp);
											
            for (int i = 0; i < vertexsize; i++)
            {
                // compute transformed gradients
                Dune::FieldVector<double,dim> grad1;
                jacInvTra.mv(basis[i].evaluateGradient(r->position()),grad1);
               
                    Dune::FieldVector<double,dim> grad2;
                    jacInvTra.mv(basis[i].evaluateGradient(r->position()),grad2);

                    // gain global inidices of vertices i and j and update associated matrix entry
                   TF_diagonal[set.subIndex(*it,i,dim)][set.subIndex(*it,i,dim)] 
+= (C_F*(10./3.)*(7./3.)*((pow(std::abs(u),1.333333))*basis[i].evaluateFunction(r->position())*basis[i].evaluateFunction(r->position())))*weight*detjac 
+ lam_kin*(grad1*grad2)*weight*detjac 
+ (2.0*phi*basis[i].evaluateFunction(r->position())*basis[i].evaluateFunction(r->position()))*weight*detjac 
/* + ((phi - ep_interaction)*(basis[i].evaluateFunction(r->position()))*(basis[i].evaluateFunction(r->position())))*weight*detjac */
+ ((grad2eps_x + grad2eps_c)*u*u + (gradeps_x + gradeps_c)*4.0*u + 2.0*(eps_x + eps_c))*(basis[i].evaluateFunction(r->position()))*(basis[i].evaluateFunction(r->position()))*weight*detjac; 
//+ 2.0*lagrange*(basis[i].evaluateFunction(r->position()))*(basis[i].evaluateFunction(r->position()))*weight*detjac;               

            }
        }
			}
	
			 for (int i=0;i<N;i++)
		 {
				TF_diagonal[i][i] += 8*penalty*(mass_u_i[i])*(mass_u_i[i]) + 4*penalty*(integral_u - NegativeCharge)*(MassMatrix[i][i]);
		 } 


		 for (int i=0;i<N;i++)
		 {
				TF_diagonal[i][i] = 1/(TF_diagonal[i][i]);
		 }
}
