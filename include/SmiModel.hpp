// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef SmiModel_HPP
#define SmiModel_HPP

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "OsiSolverInterface.hpp"
#include "CoinPackedVector.hpp"
#include "scn_c_api.h"

typedef int SmiCoreIndex;
typedef int SmiScenarioIndex;
typedef int SmiStageIndex;

typedef struct
{
	OsiSolverInterface *osi;
	SmiStageIndex nstag;
	SmiStageIndex *cstag,*rstag;
	int *rowIndex,*colIndex;
	double *rowSol,*colSol;
	int nrow,ncol,numints;
} SmiCoreModel;

class SmiModel
{
public:

	// set core model from Osi data
	SmiCoreIndex setCore(OsiSolverInterface *osi, int nstage, 
				SmiStageIndex *cstage, SmiStageIndex *rstage);

	// generate scenario 
	SmiScenarioIndex genScenarioReplaceCoreValues(SmiCoreIndex core, 
				CoinPackedMatrix *matrix,
				CoinPackedVector *drlo, CoinPackedVector *drup,
				CoinPackedVector *dobj,
				CoinPackedVector *dclo, CoinPackedVector *dcup,
				SmiStageIndex branch, SmiScenarioIndex anc, double prob);
		
	// get scenario problem data
	SmiScenarioIndex getNumScenarios(){ return scen_;}

	// getXXX by scenario
	double getObjectiveValue(SmiScenarioIndex ns);
	const double *getColSolution(SmiScenarioIndex ns);
	const double *getRowSolution(SmiScenarioIndex ns);

	// methods implented by OSLSE
	void decompSolve();
	void writeSMPS(const char *c,const char *t, const char* s);
	void readSMPS(const char *c,const char *t, const char* s);

	// OsiSolverInterface
	void setOsiSolverHandle(OsiSolverInterface *osi)
	{
		osiStoch_ = osi->clone(false);
		osiStoch_->reset();
	}
	OsiSolverInterface * getOsiSolverInterface();

	// load det equiv model into osi and return handle
	OsiSolverInterface * loadOsiSolverData();

	// constructor -- OSLSE requires max number of scenarios
	SmiModel(int maxScen): 
		scen_(0),solve_synch_(false)
	{ gutsOfConstructor(maxScen); }

	// destructor
	~SmiModel(){ gutsOfDelete();}

private:
	void gutsOfDelete();
	void gutsOfConstructor(int ns);

private:
	EKKContext *env_;
	EKKStoch *ekkStoch_;
	OsiSolverInterface * osiStoch_;
	int scen_;
	int minrow_;
	bool solve_synch_;
	std::vector<SmiCoreModel *> core_vec_;
};


#endif //#define SmiModel_HPP
