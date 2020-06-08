  std::string line;
  std::vector<bool> CORRECT_INPUT(9,false);
  std::ifstream myfile; 
 	myfile.open(argv[1]); 

  if (myfile) 
  { 
     getline (myfile,line);
     while (line != "FINITO") //whole line must be identical!
    { 
	    if (line == "ALPHA")
  		{
  			 	getline (myfile,line);
 			 	  std::istringstream isstr(line);
					isstr>>alpha10;
					//if(helper.rank()==0) cout << "alpha10 =  " << alpha10 << endl;
					CORRECT_INPUT[0] = true;
			}

			 if (line == "GLOBAL_REFINEMENT")
  		{
  			 	getline (myfile,line);
 			 	  std::istringstream isstr(line);
					isstr>>level;
					//if(helper.rank()==0) cout << "level =  " << level << endl;
					CORRECT_INPUT[1] = true;
			}

			 if (line == "DOMAIN_SIZE")
  		{
  			 	getline (myfile,line);
 			 	  std::istringstream isstr(line);
					isstr>>side_length;
					//if(helper.rank()==0) cout << "side_length =  " << side_length << endl;
					CORRECT_INPUT[2] = true;
			}
		 
			if (line == "LAMBDA_KINETIC")
  		{
  			 	getline (myfile,line);
 			 	  std::istringstream isstr(line);
					isstr>>lam_kin;
					//if(helper.rank()==0) cout << "side_length =  " << side_length << endl;
					CORRECT_INPUT[3] = true;
			}

			if (line == "QUADRATURE_RULE")
  		{
  			 	getline (myfile,line);
 			 	  std::istringstream isstr(line);
					isstr>>int_order;
					//if(helper.rank()==0) cout << "side_length =  " << side_length << endl;
					CORRECT_INPUT[4] = true;
			}
			
				if (line == "CONVERGENCE_CRITERIA_FOR_PENALTY_METHOD")
  		{
  			 	getline (myfile,line);
 			 	  std::istringstream isstr(line);
					isstr>>tol_penalty;
					//if(helper.rank()==0) cout << "side_length =  " << side_length << endl;
					CORRECT_INPUT[5] = true;
			}

				if (line == "CONVERGENCE_CRITERIA_FOR_LAGRANGE_METHOD")
  		{
  			 	getline (myfile,line);
 			 	  std::istringstream isstr(line);
					isstr>>tol_lagrange;
					//if(helper.rank()==0) cout << "side_length =  " << side_length << endl;
					CORRECT_INPUT[6] = true;
			}
			
			if (line == "GRIDFILE")
			{
			  getline (myfile,grid_file);	 	  	
 	 	  	//if(helper.rank()==0){std::cout << "Gridfile = " << grid_file << std::endl;}
				CORRECT_INPUT[7] = true;
			}

			if (line == "ATOM_POSITIONS CHARGE")
			{
			  getline (myfile,line);
 			 	int Nofatoms;
		 	  std::istringstream isstr(line);
				isstr>>Nofatoms;
				for (int i=0;i<Nofatoms;i++)
				{
				  getline (myfile,line);
			 	  std::istringstream isstr2(line);
			 	  ATOM at;
    			isstr2 >> at.pos[0];
    			isstr2 >> at.pos[1];
    			isstr2 >> at.pos[2];
    			isstr2 >> at.charge;
    			atoms.push_back(at);					
    			//TotalCharge += at.charge;
				}					
				CORRECT_INPUT[8] = true;
			}		
		
	   getline (myfile,line);
   } 
   myfile.close(); 
  } //if
  //else cout << "Unable to open file"; 
  //CHECK if input file was complete and correct
  for (std::vector<bool>::iterator it = CORRECT_INPUT.begin(); it < CORRECT_INPUT.end();it++)
  {
  	if (*it == false){ std::cout << "Error in input file \n"; return 1;}
	} 
 
  

