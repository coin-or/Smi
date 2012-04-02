/*
	Scenario-tree classes for the FlopCpp-Smi modelling framework
	authors: Michal Kaut and Alan King

	Part of the bundle described in paper "A C++ Modelling Environment for
	Stochastic Programming" by Michal Kaut, Alan King and Tim Hultberg,
	IBM Technical report RC24662, http://domino.watson.ibm.com/library/
	cyberdig.nsf/papers/3E80629707DD1782852574E300592E33

	This file includes implementations of misc. methods used for the examples,
	such as creating the SMI core models or solving the deterministic equivalent.
*/

#include "flop-smi_methods.hpp"

#include <vector>

using namespace FlopSmiEx;
using std::string;


SmiCoreData * FlopSmiEx::create_smi_core (vector<CoreNodeBase *> & coreNodes,
                                          string problemName)
{
	int i, j, t;
	int nStages = static_cast<int>(coreNodes.size());

	MP_model & coreModel = MP_model::getDefaultModel();
	coreModel.setSolver(new OSI_SOLVER_INTERFACE);
	#ifndef NDEBUG
		coreModel.verbose(); // more output from FlopC++
	#else
		coreModel.silent(); // less output from FlopC++
	#endif

	CoreNodeBase * p2root = coreNodes[0];

	coreModel.setObjective(p2root->get_obj_func()); // set the objective func.
	coreModel.attach();                             // create the OSI model

	// Get number of variables and constraints from the OSI model
	// Remember that the "->" operator on an MP_model returns the OSI model
	int nmbCoreCols = coreModel->getNumCols();
	int nmbCoreRows = coreModel->getNumRows();

	#ifndef NDEBUG
		// Write an MPS file + print some info
		coreModel->writeMps((problemName + ".core").c_str());
		cout << endl << "The core (deterministic) model has " << nmbCoreCols
			  << " variables and " << nmbCoreRows << " constraints." << endl;
	#endif

	// Get the stage numbers for all variables and constraints
	// Note that this can be done first after we have attached the model!
	int * colStages = new int[nmbCoreCols];
	for (t = 0; t < nStages; t++) {
		for (j = 0; j < (int) coreNodes[t]->all_variables.size(); j++) {
			int colIndx = coreNodes[t]->all_variables[j]->getColumn();
			#ifndef NDEBUG
				cout << "stage " << t << ": var no. " << j+1 << " is in column "
					  << colIndx << endl;
			#endif
			colStages[colIndx] = t;
		}
	}

	// The same for constraints
	int * rowStages = new int[nmbCoreRows];
	for (t = 0; t < nStages; t++) {
		for (i = 0; i < (int) coreNodes[t]->all_constraints.size(); i++) {
			int rowIndx = coreNodes[t]->all_constraints[i]->row_number()
			              + coreNodes[t]->constr_row_offsets[i];
			#ifndef NDEBUG
				cout << "stage " << t << ": constraint no. " << i+1
				     << " is in row " << rowIndx << endl;
			#endif
			rowStages[rowIndx] = t;
		}
	}

	// Build the CORE problem, i.e. the deterministic version, as an SMI object
	SmiCoreData * p2stochCore
		= new SmiCoreData(coreModel.Solver, nStages, colStages, rowStages);

	// cleaning
	delete[] colStages;
	delete[] rowStages;
	delete & coreModel; // also deletes coreModel.Solver

	return p2stochCore;
}


double FlopSmiEx::solve_det_equiv(SmiScnModel & stochModel,
                                  MP_model::MP_direction const minOrMax,
                                  double const * & p2solVector,
                                  string const problemName) {
	stochModel.setOsiSolverHandle(new OSI_SOLVER_INTERFACE());

	// 'loadOsiSolverData()' loads the deterministic equivalent into an
	// internal OSI data structure and returns a handle (pointer) to it.
	OsiSolverInterface * p2detEqModel = stochModel.loadOsiSolverData();
	p2detEqModel->setObjSense(minOrMax);
	p2detEqModel->writeMps((problemName + ".det-equiv").c_str());
	cout << endl << "The deterministic equivalent model has "
	     << p2detEqModel->getNumCols() << " variables and "
	     << p2detEqModel->getNumRows() << " constraints." << endl;

	cout << endl << "Solving the deterministic equivalent:" << endl;
	p2detEqModel->initialSolve();   // solve the LP problem/relaxation
	p2detEqModel->branchAndBound(); // if needed, run B&B
	cout << "Optimal objective value = " << p2detEqModel->getObjValue()
	     << endl << endl;

	p2solVector = p2detEqModel->getColSolution();
	return p2detEqModel->getObjValue();
}


double FlopSmiEx::minimize_det_equiv(SmiScnModel & stochModel,
                                     double const * & p2solVector,
                                     string const problemName)
{
	return solve_det_equiv(stochModel, MP_model::MINIMIZE,
	                       p2solVector, problemName);
}


double FlopSmiEx::maximize_det_equiv(SmiScnModel & stochModel,
                                     double const * & p2solVector,
                                     string const problemName)
{
	return solve_det_equiv(stochModel, MP_model::MAXIMIZE,
	                       p2solVector, problemName);
}
