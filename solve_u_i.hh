template<class GV>
void P1Elements<GV>::solve_u_i()
{
   // make linear operator from A
    Dune::MatrixAdapter<Matrix,ScalarField,ScalarField> op(ThomasFermiMatrix);

    // initialize preconditioner
    Dune::SeqILUn<Matrix,ScalarField,ScalarField> ilu2(ThomasFermiMatrix, 1, 0.92);
		//Dune::SeqSSOR<Matrix,ScalarField,ScalarField> ilu2(ThomasFermiMatrix, 1, 0.9);
		

    // the inverse operator
		 Dune::BiCGSTABSolver<ScalarField> bcgs(op, ilu2, 1e-15, 5000, 0);	
		//Dune::CGSolver<ScalarField> cg(op, ilu2, 1e-15, 5000, 0);
    Dune::InverseOperatorResult r1;

    // initialue u to some arbitrary value to avoid u being the exact 
    // solution
    u1.resize(b1.N(), false);
    u1 = 2.0;

		b1 *= -1.0;
    // finally solve the system
    bcgs.apply(u1, b1, r1);
		//cg.apply(u1, b1, r1);
} 

