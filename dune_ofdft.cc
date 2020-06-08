#define HEADERCHECK
#include"config.h"     
#include<iostream>
#include <fstream>
#include <iomanip>
#include<dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/parser.hh>
#include <dune/grid/common/gridfactory.hh>
#include<dune/common/mpihelper.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/common/quadraturerules.hh>
#include<dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include<dune/common/dynvector.hh>
#include<vector>
#include<map>
#include<set>
#include <cmath>
#include<dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include<dune/common/timer.hh>
#include <fstream>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <dune/grid/uggrid/uggridfactory.hh>
#include<dune/grid/common/gridinfo.hh>
#include"shapefunctions.hh"

//#define debug
//#define debug2

#include"global.hh"
#include"pot_el_atoms.hh"
#include"class_p1_elements.hh"
#include"determineAdjacencyPattern.hh"
#include"determineAdjacencyPattern_u_i.hh"
#include"assemble_poisson.hh"
#include"assemble_poisson_rhs.hh"
#include"assemble_mass.hh"
#include"assemble_u_i.hh"
#include"assemble_TF_diagonal.hh"
#include"solve.hh"
#include"solve_u_i.hh"
#include"lagrange_gradLphi.hh"
#include"lagrange_gradL.hh"
#include"penalty_integral_rho.hh"
#include"lagrange_integral_rho.hh"
#include"penalty_energy.hh"
#include"lagrange_energy.hh"
#include"penalty_gradLphi.hh"
#include"penalty_gradL.hh"
#include"OFDFT_Driver.hh"

int main(int argc, char **argv) 
{
	Dune::MPIHelper::instance(argc,argv);
	Dune::Timer watch;

	try
	{ 
		watch.reset();
		const int dim=3;   
	  std::cout.precision(10);	

		int level;
		int side_length;
		double alpha;
		int alpha10;
		std::string grid_file;
		
		#include"input.hh"

		alpha = alpha10/10.0;
		std::cout << "alpha = " << alpha << std::endl; 
 
	   Dune::UGGrid<dim>::setDefaultHeapSize(5000);
     
		 typedef Dune::UGGrid<dim> Grid;
     //Dune::GridPtr<Grid> gridptr("grids/testtetgen.dgf");
		 Dune::GridPtr<Grid> gridptr(grid_file.data());
     Grid &grid = *gridptr;
     grid.globalRefine(level);
     typedef Dune::UGGrid<dim>::LeafGridView GV;
     const GV& gv=grid.leafView();

	/* typedef Dune::UGGrid<3> GridType;
   Dune::UGGrid<3> grid;
   Dune::GmshReader<Dune::UGGrid<3> >::read(grid,"2atom.msh");

		typedef Dune::UGGrid<dim>::LeafGridView GV;
	  const GV& gv=grid.leafView();	*/

		std::cout << "INITIAL CONDITIONS:"<< "\n"<< "\n";	 
		std::cout << "Shapefunctions P1" << std::endl;	
		std::cout << "Z =  " << Ncore << std::endl;	 
		std::cout << "Domain Size =  " << side_length << std::endl;	 
		std::cout << "No. of nodes =  " << gv.size(dim) << std::endl;
		std::cout << "Lambda_kin =  " << lam_kin << "\n"<< "\n";	
		Dune::gridinfo(grid);
		OFDFT_Driver(grid, side_length, alpha);

/*		Dune:: VTKWriter<GridType::LeafGridView> vtkwriter(grid.leafView());
	  	vtkwriter.write("j7_grid");	 */
		double time = watch.elapsed();	
		std::cout << "time =  " << time << std::endl;	    		
  	}  

 /* try
	{ 
		watch.reset();
		const int dim=3;   
		double side_length = 10.0; 
	  std::cout.precision(10);	

		double alpha = 0.6;
	   
	  	Dune::UGGrid<dim>::setDefaultHeapSize(1000);                    
      typedef Dune::UGGrid<dim> GridType;

	  	Dune::FieldVector<double,dim> lowerLeft(-side_length/2);
    		Dune::FieldVector<double,dim> upperRight(side_length/2);
    		Dune::array<unsigned int, dim> elements;
    		std::fill(elements.begin(), elements.end(), 1);

    		Dune::StructuredGridFactory<GridType> factory;  
    		Dune::shared_ptr<GridType> grid = factory.createSimplexGrid(lowerLeft, upperRight, elements);
		
		grid->globalRefine(3);	

		typedef Dune::UGGrid<dim>::LeafGridView GV;
	  const GV& gv=grid->leafView();	

		std::cout << "INITIAL CONDITIONS:"<< "\n"<< "\n";	 
		std::cout << "Shapefunctions P1" << std::endl;	
		std::cout << "Z =  " << Ncore << std::endl;	 
		std::cout << "Domain Size =  " << side_length << std::endl;	 
		std::cout << "No. of nodes =  " << gv.size(dim) << std::endl;
		std::cout << "Lambda_kin =  " << lam_kin << "\n"<< "\n";	
		OFDFT_Driver(*grid, side_length, alpha);

		//Dune:: VTKWriter<GridType::LeafGridView> vtkwriter(grid.leafView());
		//vtkwriter.write("j_grid");	
		double time = watch.elapsed();	
		std::cout << "time =  " << time << std::endl;	    		
  	} */ 
		
  	catch (std::exception & e)
	{ 
    		std::cout << "STL ERROR: " << e.what() << std::endl;
    		return 1;
  	}
  	catch (Dune::Exception & e) 
	{
    		std::cout << "DUNE ERROR: " << e.what() << std::endl;
        	return 1;
  	}
  	catch (...)
	{
    		std::cout << "Unknown ERROR" << std::endl;
    		return 1;
  	}

  	// done
  	return 0;
} 
