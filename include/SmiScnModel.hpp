// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef SmiScnModel_HPP
#define SmiScnModel_HPP

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "SmiScenarioTree.hpp"
#include "SmiScnData.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinPackedVector.hpp"

#include <map>
#include <vector>

using namespace std;

class SmiScnNode;
//#############################################################################

/** SmiScnModel: COIN-SMI Scenario Model Class
	
	Concrete class for generating scenario stochastic linear programs.

	This class implements the Scenarios format of the Stochastic MPS
	modeling system (TODO: web pointer).  Core data and Scenarios data
	can be passed using COIN/OSI structures, or can be read from SMPS
	formatted files.

	Typical driver fragment looks like this
	\code
	SmiScnModel smi;
	smi.readSmps("app0110R");
	smi.setOsiSolverHandle(OsiClpSolverInterface());
	OsiSolverInterface *osiStoch = smi.loadOsiSolverData();
	osiStoch->initialSolve();
	\encode

	The setOsiSolverHandle method allows the user to pass in any OSI
	compatible solver.

  */
class SmiScnModel
{
	
public:

    /**@name Read SMPS files.

		There should be three files: {name}.[core, time, stoch].
		If you have different extension conventions, then you can
		hack the method yourself.
		The files can be compressed. The object that reads the files is
		derived from CoinMpsIO.
	*/
	int readSmps(const char *name);


	/**@name Direct methods.
		
		Direct methods require the user to create instances of Core files
		and Scenarios.
		Currently, the dimension of the core file nodes determines the
		dimension of the scenario nodes, but this is something that 
		could easily be changed.
	*/

//@{
	/// set core model from Osi data
	SmiCoreIndex setCore(OsiSolverInterface *osi, int nstage, 
				SmiStageIndex *cstage, SmiStageIndex *rstage);
	/// roll your own core model
	SmiCoreIndex setCore(SmiCoreData *s)
	{core_vec_.push_back(s);	return core_vec_.size()-1; }

	/// get core model 
	SmiCoreData *getCore(SmiCoreIndex i){ return core_vec_[i]; }

	/** generate scenario replacing core values
		
		Core argument must be supplied. 
		Data values replace corresponding core values, if found, 
		or creates them if not. 
		
		Scenario nodes need to have same dimensions as core nodes.
		
		Data field arguments can be NULL, or empty. 
		
		branch, anc, arguments must be supplied.  These
		identify the branching node according to the Stochastic MPS
		standard.
		
	*/
	SmiScenarioIndex genScenarioReplaceCoreValues(SmiCoreIndex core, 
				CoinPackedMatrix *matrix,
				CoinPackedVector *dclo, CoinPackedVector *dcup,
				CoinPackedVector *dobj,
				CoinPackedVector *drlo, CoinPackedVector *drup,
				SmiStageIndex branch, SmiScenarioIndex anc, double prob);
//@}
		
	// get scenario problem data
	SmiScenarioIndex getNumScenarios(){ return scen_;}
	double getScenarioProb(SmiScenarioIndex ns);
	SmiScnNode * getLeafNode(SmiScenarioIndex i){ return smiTree_.getLeaf(i)->getDataPtr(); }
	SmiScnNode * getRootNode(){ return smiTree_.getRoot()->getDataPtr(); }


	// getXXX by scenario
	double getObjectiveValue(SmiScenarioIndex ns);
	const double *getColSolution(SmiScenarioIndex ns);
	const double *getRowSolution(SmiScenarioIndex ns);

	// OsiSolverInterface
	void setOsiSolverHandle(OsiSolverInterface &osi)
	{
		osiStoch_ = osi.clone(false);
	}
	OsiSolverInterface * getOsiSolverInterface();

	// load det equiv model into osi and return handle
	OsiSolverInterface * loadOsiSolverData();

	// constructor 
	SmiScnModel(): 
		scen_(-1),solve_synch_(false),totalProb_(0)
	{ }

	// destructor
	~SmiScnModel();

	void addNode(SmiScnNode *node);

private:
	OsiSolverInterface * osiStoch_;
	int nrow_;
	int ncol_;
	int nels_;
	double *drlo_; 
	double *drup_;
	double *dobj_;
	double *dclo_; 
	double *dcup_;
	CoinPackedMatrix *matrix_;
	int scen_;
	int minrow_;
	bool solve_synch_;
	double totalProb_;
	std::vector<SmiCoreData *> core_vec_;
	SmiScenarioTree<SmiScnNode *> smiTree_;
};

class SmiScnNode
{
	friend class SmiScnModel;
public:
	int getCoreColIndex(int i);
	int getCoreRowIndex(int i);
	inline int  getColStart() {return coffset_;}	
	inline int  getRowStart() {return roffset_;}
	inline int getNumCols(){ return node_->getCore()->getNumCols(node_->getStage());}
	inline int getNumRows(){ return node_->getCore()->getNumRows(node_->getStage());}
	inline double getModelProb(){return mdl_prob_;}
	inline SmiScnNode * getParent(){ return parent_;}

	~SmiScnNode(){}

private:
	inline void setRowOffset(int r) {roffset_ = r;}
	inline void setParent(SmiScnNode * g) {parent_=g;}
	inline void setColOffset(int c) {coffset_ = c;}
	inline double addProb(double prob){ return prob_+=prob;}
	inline double getProb(){return prob_;}
	inline void setProb(double p){prob_=p;}
	inline void setModelProb(double p){mdl_prob_=p;}
	inline SmiNodeData *getNode() {return node_;}
	SmiScnNode(SmiNodeData *&node)	{node_=node;prob_=0;parent_=NULL;}

private:
	SmiNodeData *node_;
	SmiScnNode *parent_;  
	double prob_;
	double mdl_prob_;
	int coffset_;
	int roffset_;
};

// function object for addnode loop
class SmiScnModelAddNode{
public:
	void operator() (SmiScnNode *node)
	{
		s_->addNode(node);
	}

	SmiScnModelAddNode(SmiScnModel *s) { s_ = s;}
private:
	SmiScnModel *s_;


};


#endif //#define SmiScnModel_HPP
