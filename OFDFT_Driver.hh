template<class G>
void OFDFT_Driver (const G& grid, double side_length, double alpha)
{
	static const int dim = 3;
	typedef Dune::UGGrid<dim> GridType;
	typedef GridType::ctype ctype;

	typedef Dune::UGGrid<dim>::LeafGridView GV;
	const GV& gv=grid.leafView();		

	Dune::Timer watch1;
	
	//Declaring arrays for storing coefficients
	const int N = gv.size(dim);
	Dune::BlockVector<Dune::FieldVector<double,1>> u_i(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> mass_u_i(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> phi_i(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> rho_i(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> g(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> h(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> normalized_h(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> temp_h(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> temp_u_i(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> old_u_i(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> temp_phi_i(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> old_phi_i(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> gold(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> hold(N);

	Dune::BlockVector<Dune::FieldVector<double,1>> gradLold(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> gradLnew(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> gradLphi1(N);
	Dune::BlockVector<Dune::FieldVector<double,1>> gradL1(N);
	

	Dune::BlockVector<Dune::FieldVector<double,1>> gradu(dim); 
	Dune::BlockVector<Dune::FieldVector<double,1>> gradN(dim); 
	Dune::BlockVector<Dune::FieldVector<double,1>> gradphi(dim); 
	
Dune::BlockVector<Dune::FieldVector<double,1>> gradLphi(N), gradL(N), ubar_i(N), x(N), y(N), b_temp(N), b1_temp(N), old_gradL(N), new_gradL(N), h_temp(N), phi_i_old(N), atom_node(N), temp_g(N), z(N), zold(N);

Dune::BlockVector<Dune::FieldVector<double,1>> r(N), d(N), s(N), temp_d(N);
double i, k, imax, delta_0, delta_new, epsilon, eta_new, eta_prev, j, jmax, delta_d, beta, delta_old, delta_mid;

	double mu;
	double Ts;
	double Ts1;
	double Ts2;
	double Exc;
	double LT,LBT;
	double L,E,Eold, delta_u_i;
	double C_F;
	//int atom_node,counter;
	int counter;
	double penalty;
	//double integral_u;
	double factor1,factor2,rs,gradrs,eps_x,gradeps_x,eps_c,gradeps_c;
  double A,B,C,D,beta1,beta2,gamma_1;
	double L1,L2;
	double deriv_f_old,deriv_f, h_norm;
	double lambda, lambda_old, lambda_new;
  double a,b;
	double E_penalty, E1, lagrange, E_lagrange, old_lagrange;
	double integral_u;
	bool energybreakup = false;

	typedef typename GV::template Codim<dim>::Iterator VertexLeafIterator;
	typedef typename GV::template Codim<0>::Iterator ElementLeafIterator;

	//Get a set of P1 shape functions

	const P1ShapeFunctionSet<ctype,ctype,dim> &basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

 	typedef typename GV::IndexSet LeafIndexSet;
  const LeafIndexSet& set = gv.indexSet();
 	const ElementLeafIterator elend = gv.template end<0>();
	const ElementLeafIterator itend = gv.template end<0>();
	typedef typename GV::IntersectionIterator IntersectionIterator;

  //u_i *= pow((Ncore/integral_u_initial),0.5); 
	//u_i = 0.089;
	phi_i = -0.089;
	mu=1.0;
	gradu = 0.0;
	gradphi = 0.0;
	double phi = 0.0;
	Ts = 0.0;
	Ts1 = 0.0;
	Ts2 = 0.0;
	g = 0.0;
	h = 1.0;
  penalty = 10.0; 
	E = 0.0;
	E_penalty = 0.0;
	E1 = 0.0;
	lagrange = 0.6;
	E_lagrange = 0.0;
	gradLphi1 = 0.0;
	gradL1 = 0.0;
	atom_node = 0.0;

  //ATOM at;
	//at.pos[0] = 0; at.pos[1] = 1.5; at.pos[2] = 1.5; at.charge =	8; atoms.push_back(at);
	//at.pos[0] = 0; at.pos[1] = -1.5; at.pos[2] = -1.5; at.charge = 8; atoms.push_back(at); 
	//at.pos[0] = -0.5; at.pos[1] = 0; at.pos[2] = 0; at.charge = 6; atoms.push_back(at);
	//at.pos[0] = 0; at.pos[1] = 0; at.pos[2] = 0; at.charge = 10; atoms.push_back(at);

 int vertexsize = 4;
	for (ElementLeafIterator it = gv.template begin<0>(); it != elend; ++it)
	{
		//Dune::GeometryType gt = it->type();
		//const Dune::template GenericReferenceElement<double,dim> &ref =  Dune::GenericReferenceElements<double,dim>::general(gt);
		//int vertexsize = ref.size(dim);


		
		for (int i = 0 ; i<vertexsize; i++)
		{
			for (typename std::list<ATOM>::iterator list_iter = atoms.begin(); list_iter != atoms.end(); list_iter++)
			 {
       		if (it->geometry().corner(i)[0] == list_iter->pos[0] && it->geometry().corner(i)[1] == list_iter->pos[1] && it->geometry().corner(i)[2] == list_iter->pos[2])
					{
						atom_node[set.subIndex(*it,i,dim)] = -(list_iter->charge);
				    //std::cout<<it->geometry().corner(i) << "Node Number =" << set.subIndex(*it,i,dim) << std::endl;	
					} 	 	
			 }
	 
		}
	}

		Ncore = 0.0;
	 	for (typename std::list<ATOM>::iterator list_iter = atoms.begin(); list_iter != atoms.end(); list_iter++)
		{
			Ncore += (list_iter->charge);
		}

	 NegativeCharge = Ncore;
	 std::cout << "Ncore =" << Ncore << std::endl;	


		//Initializing coefficients and terms
	u_i = 0.089;

  for (ElementLeafIterator it = gv.template begin<0>(); it != elend; ++it)
  {
     Dune::GeometryType gt = it->type();
     const Dune::template GenericReferenceElement<double,dim> &ref =  Dune::GenericReferenceElements<double,dim>::general(gt);
     int vertexsize = ref.size(dim);
     for (int i = 0 ; i<vertexsize; i++)
     {

			 for (typename std::list<ATOM>::iterator list_iter = atoms.begin(); list_iter != atoms.end(); list_iter++)
			 {
       		
						 
						u_i[set.subIndex(*it,i,dim)] += exp(-alpha*pow((pow((it->geometry().corner(i)[0] - list_iter->pos[0]),2) + pow((it->geometry().corner(i)[1] - list_iter->pos[1]),2) + pow((it->geometry().corner(i)[2] - list_iter->pos[2]),2)), 0.5));
				

					 	 	
			 }
        
     }
   }    

    // Rescaling
  double integral_u_initial = 0.0;
  // Iteration over elements
  double u,sum;
  for (ElementLeafIterator el = gv.template begin<0>(); el != elend; ++el)
  {
     Dune::GeometryType gt = el->type();
     const Dune::template GenericReferenceElement<double,dim> &ref =  Dune::GenericReferenceElements<double,dim>::general(gt);
     int vertexsize = ref.size(dim);
    
     // get a quadrature rule of order one for the given geometry type
     const Dune::QuadratureRule<ctype,dim>& rule2 = Dune::QuadratureRules<ctype,dim>::rule(gt,int_order);
     for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r = rule2.begin();r != rule2.end() ; ++r)
     {    
            u = 0.0;
            // Iteration over vertices of an element
	          for (int i = 0 ; i<vertexsize; i++)
            {
                    // Calculation of u and phi
                    sum = basis[i].evaluateFunction(r->position())*u_i[set.subIndex(*el,i,dim)];
                    u  += sum;
                    
             }
             double weight = r->weight();
             double detjac = el->geometry().integrationElement(r->position());
             integral_u_initial += weight*detjac*u*u;
     }
  }   


	u_i *= pow((Ncore/integral_u_initial),0.5); 

	double Enucnuc(0.0);
	for (auto it1 = atoms.begin(); it1 != atoms.end();it1++)
	{
       Dune::FieldVector<double,3> posA = it1->pos;
       //std::cout << "posA = " << posA << std::endl;
       auto ithelp = it1;
       for (auto it2 = ++ithelp; it2 != atoms.end(); it2++)
       {
               Dune::FieldVector<double,3> posB = it2->pos;
               //std::cout << "posB = " << posB << std::endl;
               posB -=posA;
               Enucnuc += (it1->charge * it2->charge)/(posB.two_norm());
       }
}	

std::cout << "Enucnuc" << Enucnuc << std::endl;	
		
/*inline double POT_EL_ATOMS(const std::list<ATOM> atoms, const Dune::FieldVector<double,3>& pos_qp)
{
       double result = 0.0;
       Dune::FieldVector<double,3> dist;
       for (auto ait= atoms.begin(); ait != atoms.end(); ++ait)
       {
         dist = ait->pos;
         dist -= pos_qp;
         result += (ait->charge / dist.two_norm());
       }
       return result;
}; */



	/*for (ElementLeafIterator it = gv.template begin<0>(); it != elend; ++it)
	{
		Dune::GeometryType gt = it->type();
		const Dune::template GenericReferenceElement<double,dim> &ref =  Dune::GenericReferenceElements<double,dim>::general(gt);
		int vertexsize = ref.size(dim);

		for (int i = 0 ; i<vertexsize; i++)
		{								
			if (it->geometry().corner(i)[0] == -side_length/2.0 || it->geometry().corner(i)[0] == side_length/2.0 || it->geometry().corner(i)[1] == -side_length/2.0 || it->geometry().corner(i)[1] == side_length/2.0 || it->geometry().corner(i)[2] == -side_length/2.0 || it->geometry().corner(i)[2] == side_length/2.0) 
			{
				u_i[set.subIndex(*it,i,dim)] = 0.0;
				phi_i[set.subIndex(*it,i,dim)] = 0.0;
			}
		}
	} */
		
	

	  for ( ElementLeafIterator it = gv.template begin<0>() ; it != itend ; ++it) 
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

                    u_i[indexi] = 0.0; 
										phi_i[indexi] = 0.0; 
                 
                }
            }
        }
    } 

	int count = 0;
	for (int k = 0; k < gv.size(dim); k++)
	{
		if (u_i[k] == 0.0)
		{ 
        count += 1;
		}
	}
	total_count = count;
	std::cout << "Amount of boundary nodes =   " << total_count <<"\n";
	//Dune::DynamicVector<int> bc(total_count);
	boundary_nodes = new int[total_count];
	count = 0;

	for (int k = 0; k < gv.size(dim); k++)
	{
		if (u_i[k] == 0.0)
		{ 
			boundary_nodes[count] = k;
      count += 1;
		}
	}


	
	int iter = 0;
	//for (int iter = 0; iter<14; iter++)
		P1Elements<GV> p3(gv, u_i, phi_i, atom_node, integral_u, gradL, lagrange, penalty);
    p3.determineAdjacencyPattern();
    p3.assemble_mass();
		(p3.MassMatrix).mv(u_i,mass_u_i);	
		integral_u = u_i*mass_u_i;
	//	integral_u = penalty_integral_rho(gv,u_i);
  //	while((compute_gradLphi1(gv, u_i, phi_i, atom_node, side_length, penalty, E, E_penalty, E1, integral_u).two_norm() > pow(10,-4)) ||    (compute_gradL1(gv, u_i, phi_i, atom_node, side_length, penalty, E, E_penalty, E1, integral_u).two_norm() > pow(10,-4)))

		P1Elements<GV> p4(gv, u_i, phi_i, atom_node, integral_u, gradL, lagrange, penalty);
    p4.determineAdjacencyPattern();
		p4.determineAdjacencyPattern_u_i();
		p4.assemble_mass();
		p4.assemble_u_i();
		p4.assemble_TF_diagonal();
    p4.assemble_poisson();

		/*typedef typename p4.TF_diagonal :: ConstRowIterator RowI ;
		typedef typename p4.TF_diagonal :: ConstColIterator ColI ;
		for ( RowI row = matrix . begin (); row != matrix . end (); ++ row )
		{
			std :: cout << " row " << row . index () < < " : "
			for ( ColI col = row - > begin (); col != row - > end (); ++ col )
				std :: cout < < col . index () < < " " ;
			std :: cout < < std :: endl ;
		}*/


//	for (int r = 1; r < 8; r++) 
	gradL1 = 1.0;
	gradLphi1 = 1.0;
	while(gradL1.two_norm() > tol_penalty || gradLphi1.two_norm() > tol_penalty)
	{	
		 
		std::cout << "********************************************" << "\n";
		std::cout << "BIG ITERATION   " << iter <<  "\n" ;
		std::cout << "********************************************" << "\n";
		
		std::cout << "********************************************" << "\n";
		std::cout << "Maximize w.r.t phi_i" <<"    Print values after each CG iteration" << "\n" ; 
		std::cout << "********************************************" << "\n";

		# ifdef debug
		std::cout << "||gradLphi||  " <<std::endl;
		# endif		

		#ifdef debug2
		std::cout << "E(pen+self+)       " << "E(pen+self-)    " << "E_penalty  "<< "       ||gradL||    "  <<  "  ||gradLphi||  " <<std::endl;
		#endif
		
	

		/*P1Elements<GV> p4(gv, u_i, phi_i, atom_node, integral_u, gradL, lagrange, penalty);
    p4.determineAdjacencyPattern();
    p4.assemble_poisson();
		p4.assemble_poisson_rhs();*/

	  int n = 0;
    penalty_gradLphi(gv, gradLphi1, u_i, phi_i, atom_node, side_length, penalty, integral_u);
		g = gradLphi1;
	  g *= -1; 
		phi_i_old = phi_i;
		while(g.two_norm() > 1e-4)
		{
			// Set the values of gold and hold from the previous iteration
			// Negative sign not included because it's a maximization problem

			(p4.poisson_diagonal).mv(g,z);
			
			h = z; 	
			if (n>0)
			{
		 		//double gamma =  (g*g)/(gold*gold);
				double gamma =  (z*g)/(zold*gold);
				hold *= gamma;
		 		h += hold;		 		
			}
			hold = h;

			(p4.PoissonMatrix).mv(h,h_temp);
			lambda = (z*g)/(h*h_temp);

			
		  h *= lambda;
		  phi_i +=  h; 
			
			gold = g;
			zold = z;
		 //	penalty_gradLphi(gv, gradLphi1, u_i, phi_i, atom_node, side_length, penalty, integral_u);
		 // g = gradLphi1;
			
			h_temp *= lambda;
			g -= h_temp;
			
			# ifdef debug
			penalty_gradLphi(gv, gradLphi1, u_i, phi_i, atom_node, side_length, penalty, integral_u);
			std::cout << std::fixed  << std::scientific << std::showpos << gradLphi1.two_norm() << std::endl;

			# ifdef debug2
			penalty_gradL(gv, gradL1, u_i, phi_i, atom_node, side_length, penalty, integral_u);
			penalty_energy(gv, u_i, phi_i, atom_node, penalty, lagrange, integral_u, E, E1, E_penalty, E_lagrange, Enucnuc, atoms);			
			std::cout <<  std::fixed  << std::scientific << std::showpos << E  << "\t""\t""\t" << E1 << "\t""\t""\t" << E_penalty << "\t""\t""\t"<< gradL1.two_norm() << "\t""\t""\t" << gradLphi1.two_norm() << std::endl;
			# endif
			# endif
	
			n = n+1;			
		} 

	penalty_energy(gv, u_i, phi_i, atom_node, penalty, lagrange, integral_u, E, E1, E_penalty, E_lagrange, Enucnuc, atoms);	
	std::cout << "Energy after solving poisson =   " << E1 <<"\n";
		
	/*	if (iter > 0)
		{
		std::cout << " Potential Mixing Performed" << std::endl;
		phi_i *= 0.5;
		phi_i_old *= 0.5;		
		phi_i += phi_i_old;
		} */

		std::cout << "********************************************" << "\n";
		std::cout << "Minimize w.r.t u_i" <<  "\n" ;
		std::cout << "********************************************" << "\n";
	# ifdef debug
		std::cout	<< "lambda        "<<"E(pen+self+)  " <<"E(p+/self-)    "<<"E_penalty     " << "||gradL||     "	<< "||gradLphi||  "<< "h*gradL       " << "integral_u    " <<std::endl;
	#else
		std::cout	<< " lambda       " << "     ||gradL||        "	<< "      deriv_f        " << "       integral_u       " << std::endl;

  # endif

	// Minimization w.r.t u_i
		double a = pow(0.2,iter); 
		double cond = std::min(a,pow(10,-5));
		n = 0;
		penalty_gradL(gv, gradL1, u_i, phi_i, atom_node, side_length, penalty, integral_u);
		g = gradL1;
		g *= -1.0;
		
		bool restart = false;
		while(g.two_norm() > cond && n < 15)
		{	
			// Set the values of gold and hold from the previous iteration

			p4.assemble_TF_diagonal();
			(p4.TF_diagonal).mv(g,z);
		

			h = z; //h = g;
	  	if (n>0) 
	   	{
				temp_g = g;
				temp_g -= gold;
			 	double gamma =  (z*temp_g)/(zold*gold); 
				hold *= gamma;
				h += hold;
	   	}		      
			
		/*	if (g*h < 1) 
	   	{
				restart = true;
	   	}		 */
						
			 if (restart)
			{
			h = g;
			std::cout << "RESTART CG PERFORMED" << std::endl;
			restart = false;
			} 

			
    	
			h_norm = h.two_norm();
			normalized_h = h;
			normalized_h /= h_norm;
			
	   	// Line search
			lambda_old = 0.0;
	   	lambda = pow(10,-10);
	   	temp_h = normalized_h;
			temp_h *= lambda;
			temp_u_i = u_i; 
			temp_u_i += temp_h;
			deriv_f_old = -(normalized_h*g); //see above
			(p3.MassMatrix).mv(temp_u_i,mass_u_i);	
  		integral_u = temp_u_i*mass_u_i;
			//integral_u = penalty_integral_rho(gv,temp_u_i);
			penalty_gradL(gv, gradL1, temp_u_i, phi_i, atom_node, side_length, penalty, integral_u);
			gradLnew = gradL1;
			deriv_f = (normalized_h*gradLnew); 
	   	double min_lambda_line =0.0;
			double min_deriv_f_line = 1e6;
			double min_E_line = 1e6;
	   	int n_line = 0;
	   	while( (deriv_f > 1e-5 || deriv_f < -1e-5) && n_line < 25)
	   	//START LINE SEARCH
	 	 	{
				lambda_new = lambda - ((lambda - lambda_old)*deriv_f)/(deriv_f - deriv_f_old);
				lambda_old = lambda;
				lambda = lambda_new;
				temp_h = normalized_h;
				temp_h *= lambda;
				temp_u_i = u_i;
        temp_u_i += temp_h; //u_i + lambda h* 
				
				(p3.MassMatrix).mv(temp_u_i,mass_u_i);	
    		integral_u = temp_u_i*mass_u_i;
				//integral_u = penalty_integral_rho(gv,temp_u_i);		
				penalty_gradL(gv, gradL1, temp_u_i, phi_i, atom_node, side_length, penalty, integral_u);
				gradLnew = gradL1;
				deriv_f_old = deriv_f;	
				deriv_f = (normalized_h * gradLnew);
			
			/*	penalty_energy(gv, u_i, phi_i, atom_node, penalty, lagrange, integral_u, E, E1, E_penalty, E_lagrange, Enucnuc, atoms);	

				if (min_E_line > E1)
				{
					min_E_line = E1;
					min_lambda_line = lambda;
				} */

				# ifdef debug
				penalty_gradLphi(gv, gradLphi1, temp_u_i, phi_i, atom_node, side_length, penalty, integral_u);
				penalty_energy(gv, u_i, phi_i, atom_node, penalty, lagrange, integral_u, E, E1, E_penalty, E_lagrange, Enucnuc, atoms);	
				# endif

				
				# ifdef debug
	  		std::cout	<< std::fixed  << std::scientific << std::showpos << lambda << "\t" << E << "\t" << E1<<"\t"<< E_penalty << "\t" << gradLnew.two_norm() << "\t" << gradLphi1.two_norm() << "\t" 
	  							<< deriv_f << "\t" << integral_u << std::endl;
				#else	
				std::cout	<< std::fixed  << std::scientific << std::showpos << lambda << "\t" << gradLnew.two_norm() << "\t" << deriv_f << "\t" << integral_u << std::endl;
				
				# endif 

	  		n_line++;			
				if (n_line == 25) 
				{
					restart = true;		
 					//lambda = min_lambda_line;
					lambda = 0.0;
					temp_h = normalized_h;
					temp_h *= lambda;
          temp_u_i = u_i; 
					temp_u_i += temp_h;
				}								
			}// End of line search
		  	
		  //UPDATE	
			old_u_i = u_i;
			gold = g;	
			zold = z;
		  u_i =  temp_u_i; 
			(p3.MassMatrix).mv(u_i,mass_u_i);	
  		integral_u = u_i*mass_u_i;
			//integral_u = penalty_integral_rho(gv,u_i);
		 	penalty_gradL(gv, gradL1, u_i, phi_i, atom_node, side_length, penalty, integral_u);
			g = gradL1;
	  	g *= -1.0;
	  	hold = h; //not normalised  
			
			# ifdef debug
	  	  std::cout << "\n" << "NEXT LINE" << "\n";
				std::cout	<< "lambda        "<<"E(pen+self+)  " <<"E(p+/self-)   "<<"E_penalty     " << "||gradL||     "	<< "||gradLphi||  "<< "h*gradL       " << "integral_u    " <<std::endl;
		  #else
				std::cout	<< " lambda       " << "     ||gradL||        "	<< "      deriv_f        " << "       integral_u       " << std::endl;
		  # endif

	  	n = n+1;
	 	} //while u_i min


		penalty_gradLphi(gv, gradLphi1, temp_u_i, phi_i, atom_node, side_length, penalty, integral_u);
		penalty_energy(gv, u_i, phi_i, atom_node, penalty, lagrange, integral_u, E, E1, E_penalty, E_lagrange, Enucnuc, atoms);	
		//penalty_gradL(gv, gradL1, u_i, phi_i, atom_node, side_length, penalty, integral_u);


			# ifdef debug
	    #else
		  std::cout	<<"E(pen+self+)       " <<" E(p+/self-)         "<<"E_penalty          " << " ||gradL||         "	<< " ||gradLphi||         " << "integral_u     " << std::endl;
			std::cout	<< std::fixed  << std::scientific << std::showpos << E << "\t" << E1<<"\t"<< E_penalty << "\t" << gradL1.two_norm() << "\t" << gradLphi1.two_norm() << "\t" << integral_u << std::endl;
      # endif

		
		iter = iter + 1;
		if (iter < 50)
		{
			penalty *= 1.2;
		}

			
			
		
	} //END OF MINIMISATION PART WITH GRADIENT DESCENT
	// PASS VALUES TO NEWTON METHODS

	(p3.MassMatrix).mv(u_i,mass_u_i);	
	integral_u = u_i*mass_u_i;
	lagrange = 2*penalty*(integral_u - NegativeCharge);

	lagrange_gradLphi(gv, gradLphi, u_i, phi_i, atom_node, side_length, lagrange, integral_u);
	lagrange_gradL(gv, gradL, u_i, phi_i, atom_node, side_length, lagrange, integral_u);

		
  std::cout << "Initial Lagrange parameter  =  " <<lagrange<< "\n";

	// Set Up Poisson Matrix
		P1Elements<GV> p1(gv, u_i, phi_i, atom_node, integral_u, gradL, lagrange, penalty);
   	p1.determineAdjacencyPattern();
   	p1.assemble_poisson();

// Set Up Thomas-Fermi Matrix

			P1Elements<GV> p2(gv, u_i, phi_i, atom_node, integral_u, gradL, lagrange, penalty);
      p2.determineAdjacencyPattern_u_i();
	double Eprev = 0;
 	while((gradLphi.two_norm() > tol_lagrange) || (gradL.two_norm() > tol_lagrange) || std::abs(E1 - Eprev) > 0.0001)
	{	
		std::cout << "********************************************" << "\n";
		std::cout << "BIG NEW ITERATION   " << iter <<  "\n" ;
		std::cout << "********************************************" << "\n";
		
		lagrange_energy(gv, u_i, phi_i, atom_node, penalty, lagrange, integral_u, E, E1, E_penalty, E_lagrange, Enucnuc, atoms, energybreakup);
		Eprev = E1;
		
		std::cout << "********************************************" << "\n";
		std::cout << "Minimize w.r.t u_i" <<  "\n" ;
		std::cout << "********************************************" << "\n";
		
		# ifdef debug
		std::cout <<"E (w lagrange & self)        " <<"   E (w/o self)  "<<"     E_lagrange    " << "     ||gradL||   "	<< "     ||gradLphi||   "<< "     integral_u  " << std::endl;
		# else
		std::cout <<"E (w lagrange & self)" <<"   E (w/o self)  "<<"   E_lagrange  " << "   ||gradL||   "	<< "  integral_u " <<std::endl;
		# endif		

		Dune::Timer watch_solve;

		// Minimization w.r.t u_i
		//double a = pow(0.1,iter); 
		//double cond = std::min(a,pow(10,-5));
//	n = 0;
		lagrange_gradL(gv, gradL, u_i, phi_i, atom_node, side_length, lagrange, integral_u);
		delta_u_i = 1.0;
		double tol_inner;
		tol_inner = tol_lagrange;
		if (gradL.two_norm() < tol_inner)
		{ tol_inner = 1e-10;}
		while(gradL.two_norm() > tol_inner)
		{	
		//u_i_old = u_i;
		//P1Elements<GV> p2(gv, u_i, phi_i, atom_node, integral_u, gradL, lagrange, penalty);
    //p2.determineAdjacencyPattern_u_i();
    p2.assemble_u_i();
		b1_temp = p2.b1;
		watch_solve.reset();
    p2.solve_u_i();
	  // std::cout << "Time in solving the system = " << watch_solve.elapsed() << std::endl;	
    ubar_i = p2.u1; 
		
		// TEST std::cout << "ubar_i.two_norm() =  " << ubar_i.two_norm() << "\n";

			old_gradL = gradL;
			old_gradL *= 1.2;
			new_gradL = gradL;
			
			beta = 0.5;

			old_u_i = u_i;
			old_lagrange = lagrange;

			for (int i=0;i<N;i++)
			{
				 u_i[i] += ubar_i[i];
			}
			lagrange += ubar_i[N];

			lagrange_gradL(gv, gradL, u_i, phi_i, atom_node, side_length, lagrange, integral_u);
			new_gradL = gradL;

			double old_norm = old_gradL.two_norm();
			double new_norm = new_gradL.two_norm();
			
			if (old_norm < new_norm)
			{
			while (old_norm < new_norm)
			{
			
			u_i = old_u_i;
			lagrange = old_lagrange; 
			
			ubar_i *= beta;			
		
			for (int i=0;i<N;i++)
			{
				 u_i[i] += ubar_i[i];
			}
			lagrange += ubar_i[N]; 
			lagrange_gradL(gv, gradL, u_i, phi_i, atom_node, side_length, lagrange, integral_u);
			new_gradL = gradL;
			new_norm = gradL.two_norm();
			// TEST std::cout << "    old_gradL =    " <<    old_gradL.two_norm()  << "     new_gradL =    " << new_gradL.two_norm() << "     beta =    " << beta <<"\n";
			}
			}

		/*	for (int i=0;i<N;i++)
			{
				 u_i[i] += ubar_i[i];
			}
			lagrange += ubar_i[N]; */ 
			
				(p3.MassMatrix).mv(u_i,mass_u_i);	
      	integral_u = u_i*mass_u_i;
	  	//integral_u = lagrange_integral_rho(gv,u_i);

	//	std::cout << "lagrange =  " <<lagrange<< "\n"; //TEST

		//	u_i *= pow((Ncore/integral_u),0.5);

		//	(p2.ThomasFermiMatrix).mv(p2.u1,y);
		//	y += b1_temp;
			lagrange_gradL(gv, gradL, u_i, phi_i, atom_node, side_length, lagrange, integral_u);

			# ifdef debug
			lagrange_gradLphi(gv, gradLphi, u_i, phi_i, atom_node, side_length, lagrange, integral_u);
			# endif

			lagrange_energy(gv, u_i, phi_i, atom_node, penalty, lagrange, integral_u, E, E1, E_penalty, E_lagrange, Enucnuc, atoms, energybreakup);

				# ifdef debug
	  	std::cout	 << std::fixed  << std::scientific << std::showpos << E << "\t" << E1<<"\t" << E_lagrange << "\t" <<gradL.two_norm() << "\t" << gradLphi.two_norm() << "\t" << integral_u << "\t" << std::endl;

				#else	
				std::cout	 << std::fixed  << std::scientific << std::showpos << E << "\t" << E1 << "\t" << E_lagrange << "\t" <<gradL.two_norm() << "\t" << integral_u << "\t" << std::endl;

				# endif
		
			delta_u_i = ubar_i.two_norm();
	 	}

		std::cout << "********************************************" << "\n";
		std::cout << "Maximize w.r.t phi_i" <<"    Print values after each CG iteration" << "\n" ; 
		std::cout << "********************************************" << "\n";
		
		# ifdef debug
		std::cout << "E (w lagrange & self)  " << "     E (w/o self)    " << "       E_lagrange    "<< "     ||gradL||    "  <<  "     ||gradLphi||  " << "       integral_u "<<std::endl;
		# else
		std::cout << "E (w lagrange & self)  " << "   E (w/o self)    " << "     E_lagrange    "<<  "     ||gradLphi||  " << "     integral_u " <<std::endl;
		# endif
		Dune::Timer watch_phi;
		watch_phi.reset();

		//Maximization w.r.t \phi_i
		//for (int n = 0; n<10; n++)
		//int n = 0;
		phi_i_old = phi_i;
		//lagrange_gradL(gv, gradL, u_i, phi_i, atom_node, side_length, lagrange, integral_u);
	  //P1Elements<GV> p1(gv, u_i, phi_i, atom_node, integral_u, gradL, lagrange, penalty);
    //p1.determineAdjacencyPattern();
    //p1.assemble_poisson();
		p1.assemble_poisson_rhs();
		b_temp = p1.b;
		p1.solve();
    phi_i = p1.u;

		phi_i *= 0.5;
		phi_i_old *= 0.5;		
		phi_i += phi_i_old;
		
   // (p1.PoissonMatrix).mv(p1.u,x);
	 //	x += b_temp;
		
			(p3.MassMatrix).mv(u_i,mass_u_i);	
			integral_u = u_i*mass_u_i;
			//	integral_u = lagrange_integral_rho(gv,u_i);
			
			# ifdef debug
			lagrange_gradL(gv, gradL, u_i, phi_i, atom_node, side_length, lagrange, integral_u);
			# endif

			lagrange_gradLphi(gv, gradLphi, u_i, phi_i, atom_node, side_length, lagrange, integral_u);
		
    	// TEST std::cout << " residual = "<< x.two_norm() << std::endl;
			
			lagrange_energy(gv, u_i, phi_i, atom_node, penalty, lagrange, integral_u, E, E1, E_penalty, E_lagrange, Enucnuc, atoms, energybreakup);

			# ifdef debug
			std::cout << std::fixed  << std::scientific << std::showpos << E  << "\t" << E1 << "\t" << E_lagrange << "\t"<< gradL.two_norm() << "\t" << gradLphi.two_norm() << "\t" << integral_u << std::endl;
			#else	
			std::cout	 << std::fixed  << std::scientific << std::showpos << E << "\t" << E1<<"\t" << E_lagrange << "\t" << "\t" << gradLphi.two_norm() << "\t" << integral_u << std::endl;
			# endif
					
		iter = iter + 1;	
	}
	
	energybreakup = true;
	lagrange_energy(gv, u_i, phi_i, atom_node, penalty, lagrange, integral_u, E, E1, E_penalty, E_lagrange, Enucnuc, atoms, energybreakup);

	for (int k = 0; k < gv.size(dim); k++)
	{
		rho_i[k] = u_i[k] * u_i[k]; //pow(u_i[k],2);
	}
	 
	Dune:: VTKWriter<GridType::LeafGridView> vtkwriter(grid.leafView());
	vtkwriter.addVertexData(rho_i,"rho_i");
	vtkwriter.write("j_grid");

	//Dune:: VTKWriter<GridType::LeafGridView> vtkwriter(grid.leafView());
	vtkwriter.addVertexData(phi_i,"phi_i");
  vtkwriter.write("j1_grid"); 	
	
						
}
