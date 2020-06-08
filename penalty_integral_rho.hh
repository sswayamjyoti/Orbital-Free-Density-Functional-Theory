template<class GV>
double penalty_integral_rho(const GV& gv, Dune::BlockVector<Dune::FieldVector<double,1>>& u_i)
{
	static const int dim = 3;
	typedef Dune::UGGrid<3>::ctype ctype;
	double u,sum, integral_u;
	const P1ShapeFunctionSet<ctype,ctype,dim> &basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

	typedef typename GV::template Codim<dim>::Iterator VertexLeafIterator;
	typedef typename GV::template Codim<0>::Iterator ElementLeafIterator;
	typedef typename GV::IndexSet LeafIndexSet;
	const LeafIndexSet &set = gv.indexSet();
	const ElementLeafIterator elend = gv.template end<0>();
	
	integral_u = 0.0;

	// Iteration over elements
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
			integral_u += weight*detjac*u*u;
		}
	}	
return(integral_u);
}
