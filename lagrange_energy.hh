template<class GV>
void lagrange_energy(const GV& gv, Dune::BlockVector<Dune::FieldVector<double,1>>& u_i, Dune::BlockVector<Dune::FieldVector<double,1>>& phi_i, Dune::BlockVector<Dune::FieldVector<double,1>>& atom_node, double penalty, double lagrange, double integral_u, double& E, double& E1, double& E_penalty, double& E_lagrange, double& Enucnuc, const std::list<ATOM> atoms, bool energybreakup)
{
	Dune::Timer watch_energy;
	static const int dim = 3;
	typedef Dune::UGGrid<3>::ctype ctype;
	Dune::FieldVector<double,3> qp;

	Dune::BlockVector<Dune::FieldVector<double,1>> gradN(dim), gradu(dim), gradphi(dim);

	double C_F, Ts1, Ts2, Ts, L1, L2, ES1, ES2, u, phi, Exc, eps_x, eps_c, rs, sum, sum1, ep_interaction;
	const P1ShapeFunctionSet<ctype,ctype,dim> &basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

	typedef typename GV::template Codim<dim>::Iterator VertexLeafIterator;
	typedef typename GV::template Codim<0>::Iterator ElementLeafIterator;
	typedef typename GV::IndexSet LeafIndexSet;
	const LeafIndexSet &set = gv.indexSet();
	const ElementLeafIterator elend = gv.template begin<0>();
 	Dune::GeometryType gt = elend->type();

	watch_energy.reset();

	Ts1 = 0.0;
	Ts2 = 0.0;		
	Ts = 0.0;
	Exc = 0.0;
	ES1 = 0.0;
	ES2 = 0.0;
	L1 = 0.0;
  L2 = 0.0;
  C_F = 0.3*pow((3*3.14*3.14),0.66666666666);
			
  const Dune::QuadratureRule<ctype,dim>& rule2 = Dune::QuadratureRules<ctype,dim>::rule(gt,int_order);
 	const Dune::template GenericReferenceElement<double,dim> &ref =  Dune::GenericReferenceElements<double,dim>::general(gt);
  int vertexsize = ref.size(dim);
	for (ElementLeafIterator el = gv.template begin<0>(); el != gv.template end<0>(); ++el)
	{		
		// determine geometry type of the current element and get the matching reference element
		// get a quadrature rule of order one for the given geometry type
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

			// Evaluate energy terms

		  	rs = pow((3/(4*pi*u*u)),0.33333333333333);
   			eps_x = (-3./4.)*pow((3./pi),0.3333333)*pow(std::abs(u),0.666666666667);//*(pow(pow(u,2),0.3333333333)); 
                       
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

			qp = el->geometry().global(r->position());
			ep_interaction = POT_EL_ATOMS(atoms, qp);
				
		//Ts1 += C_F*((pow(pow(u,4),0.3333333333))*pow(u,2))*weight*detjac; 
			Ts1 += C_F*(pow(std::abs(u),3.3333333333))*weight*detjac; 
			Ts2 += (gradu*gradu)*weight*detjac;
			//L1 += 0.5*(phi + Ncore/(el->geometry().global(r->position()).two_norm()))*(pow(u,2))*weight*detjac;
			//L2 += (-Ncore/(el->geometry().global(r->position()).two_norm()))*(pow(u,2))*weight*detjac;	
			L1 += 0.5*(phi + ep_interaction)*u*u*weight*detjac;
			L2 += (-ep_interaction)*u*u*weight*detjac;				
			Exc += (eps_x + eps_c)*u*u*weight*detjac;
			ES1 += -(1/(8*3.14))*(gradphi*gradphi)*weight*detjac;
			ES2 += (u*u*phi)*weight*detjac;
			}
	}
			//ES2 += -Ncore*phi_i[atom_node];
			ES2 += atom_node*phi_i;
		  Ts = Ts1 + 0.5*lam_kin*Ts2; 
			E_lagrange = lagrange*(integral_u - NegativeCharge);
			E_penalty = penalty*(integral_u - NegativeCharge)*(integral_u - NegativeCharge);	
			E = Ts + Exc + ES1 + ES2 + Enucnuc;  
		  E1 = Ts + Exc + L1 + L2 + E_lagrange + Enucnuc;  

		if (energybreakup)
		{
		std::cout << "Total Energy  = " << E1 << "\n" << "E_kinetic  = " << Ts << "\n" << "E_xc  = " << Exc << "\n" << "E_coul  = " << L1 << "\n" << "E_ext  = " << L2 << "\n" << "E_nucnuc  = " << Enucnuc << "\n";
		}
} 
