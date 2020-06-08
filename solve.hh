template<class GV>
void P1Elements<GV>::solve()
{
   // make linear operator from A
    Dune::MatrixAdapter<Matrix,ScalarField,ScalarField> op(PoissonMatrix);

    // initialize preconditioner
    Dune::SeqSOR<Matrix,ScalarField,ScalarField> sor1(PoissonMatrix, 6, 0.92);

    // the inverse operator
    Dune::BiCGSTABSolver<ScalarField> bcgs(op, sor1, 1e-15, 5000, 0);	
		//Dune::CGSolver<ScalarField> cg(op, sor1, 1e-15, 5000, 0);
    Dune::InverseOperatorResult r;

    // initialize u to some arbitrary value to avoid u being the exact
    // solution
    u.resize(b.N(), false);
    u = 0.0;

		b *= -1.0;
    // finally solve the system
    bcgs.apply(u, b, r);
		//cg.apply(u, b, r);

	//	A.solve(u,b);

} 

