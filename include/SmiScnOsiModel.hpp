// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef SmiScnModel_HPP
#define SmiScnModel_HPP

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "SmiScenarioTree.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinPackedVector.hpp"

#include <map>

using namespace std;


typedef int SmiCoreIndex;
typedef int SmiScenarioIndex;
typedef int SmiStageIndex;

class SmiScnOsiCoreModel;

class SmiScnOsiNode
{
public:
	typedef map<int,CoinPackedVector *> SmiRowMap;
	inline  CoinPackedVector * getRowLower(){return drlo_;}
	inline  CoinPackedVector * getRowUpper(){return drup_;}
	inline  CoinPackedVector * getColLower(){return dclo_;}
	inline  CoinPackedVector * getColUpper(){return dcup_;}
	inline  CoinPackedVector * getObjCoefficients(){return dobj_;}
	CoinPackedVector * getRow(int i) { 
		SmiRowMap::iterator r=rowMap.find(i); 
		if (r!=rowMap.end()) return r->second; 
		else return NULL;}
	inline SmiScnOsiCoreModel * getCore() { return core_;}
	inline int getStage() { return stg_;}
	inline  int getNumElements(){ return nels_; }
	SmiScnOsiNode(SmiStageIndex stg, SmiScnOsiCoreModel *core,
				 const CoinPackedMatrix *const matrix,
				 CoinPackedVector *dclo, 
				 CoinPackedVector *dcup,
				 CoinPackedVector *dobj,
				 CoinPackedVector *drlo, 
				 CoinPackedVector *drup);
	~SmiScnOsiNode();
private:
	SmiStageIndex stg_;
	SmiRowMap rowMap;
	int nels_;
	CoinPackedVector *drlo_; 
	CoinPackedVector *drup_;
	CoinPackedVector *dobj_;
	CoinPackedVector *dclo_; 
	CoinPackedVector *dcup_;
	SmiScnOsiCoreModel *core_;
};


class SmiScnOsiCoreModel
{
public:
	inline int getNumCols(){ return ncol_;}
	inline int getNumRows(){ return nrow_;}
	inline int getNumStages(){ return nstag_;}
	inline int getNumCols(SmiStageIndex t){ return nColInStage_[t];}
	inline int getNumRows(SmiStageIndex t){ return nRowInStage_[t];}
	inline int getColStart(SmiStageIndex t){ return stageColPtr_[t];}
	inline int getRowStart(SmiStageIndex t){ return stageRowPtr_[t];}
	inline int getColStage(int i){ return colStage_[i];}
	inline int getRowStage(int i){ return rowStage_[i];}
	inline int getRowInternalIndex(int i){ return rowEx2In_[i];}
	inline int getColInternalIndex(int i){ return colEx2In_[i];}
	inline int getRowExternalIndex(int i){ return rowIn2Ex_[i];}
	inline int getColExternalIndex(int i){ return colIn2Ex_[i];}
	inline CoinPackedVector * getMatrixRow(SmiStageIndex t, int i){ return nodes_[t]->getRow(i);}
	inline CoinPackedVector *getRowLower(SmiStageIndex t){return nodes_[t]->getRowLower();}
	inline CoinPackedVector *getRowUpper(SmiStageIndex t){return nodes_[t]->getRowUpper();}
	inline CoinPackedVector *getColLower(SmiStageIndex t){return nodes_[t]->getColLower();}
	inline CoinPackedVector *getColUpper(SmiStageIndex t){return nodes_[t]->getColUpper();}
	inline CoinPackedVector *getObjCoefficients(SmiStageIndex t){return nodes_[t]->getObjCoefficients();}
	inline SmiScnOsiNode * getNode(SmiStageIndex t){return nodes_[t];}
	SmiScnOsiCoreModel(OsiSolverInterface *osi, int nstag, int *cstag, int *rstag);
	~SmiScnOsiCoreModel(){}
private:
	int nrow_;
	int ncol_;
	SmiStageIndex nstag_;
	int *nColInStage_;
	int *nRowInStage_;
	int *stageColPtr_;
	int *stageRowPtr_;
	int *colStage_;
	int *rowStage_;
	int *colEx2In_;
	int *rowEx2In_;
	int *colIn2Ex_;
	int *rowIn2Ex_;
	vector<SmiScnOsiNode*> nodes_;
};


class SmiScnOsiTreeNode
{
public:
	inline double getProb(){return prob_;}
	inline void setProb(double p){prob_=p;}
	inline SmiScnOsiTreeNode * getParent(){ return parent_;}
	inline void setParent(SmiScnOsiTreeNode * g) {parent_=g;}
	inline void setColOffset(int c) {coffset_ = c;}
	inline int  getColOffset() {return coffset_;}	
	inline void setRowOffset(int r) {roffset_ = r;}
	inline int  getRowOffset() {return roffset_;}
	inline double addProb(double prob){ return prob_+=prob;}
	inline double scaleProb(double prob){ return prob_/=prob;}
	inline SmiScnOsiNode *getNode() {return node_;}
	SmiScnOsiTreeNode(SmiScnOsiNode *&node)
	{node_=node;prob_=0;parent_=NULL;}
	~SmiScnOsiTreeNode(){}
private:
	SmiScnOsiNode *node_;
	SmiScnOsiTreeNode *parent_;  // for convenience in addNode
	double prob_;
	int coffset_;
	int roffset_;
};

class SmiScnOsiModel
{

public:

	// set core model from Osi data
	SmiCoreIndex setCore(OsiSolverInterface *osi, int nstage, 
				SmiStageIndex *cstage, SmiStageIndex *rstage);

	// generate scenario 
	SmiScenarioIndex genScenarioReplaceCoreValues(SmiCoreIndex core, 
				CoinPackedMatrix *matrix,
				CoinPackedVector *dclo, CoinPackedVector *dcup,
				CoinPackedVector *dobj,
				CoinPackedVector *drlo, CoinPackedVector *drup,
				SmiStageIndex branch, SmiScenarioIndex anc, double prob);

		
	// get scenario problem data
	SmiScenarioIndex getNumScenarios(){ return scen_;}
	double getScenarioProb(SmiScenarioIndex ns);

	// getXXX by scenario
	double getObjectiveValue(SmiScenarioIndex ns);
	const double *getColSolution(SmiScenarioIndex ns);
	const double *getRowSolution(SmiScenarioIndex ns);

	// methods implented by OSLSE
	void writeSMPS(const char *c,const char *t, const char* s);
	void readSMPS(const char *c,const char *t, const char* s);

	// OsiSolverInterface
	void setOsiSolverHandle(OsiSolverInterface *osi)
	{
		osiStoch_ = osi->clone(false);
//		osiStoch_->reset();
	}
	OsiSolverInterface * getOsiSolverInterface();

	// load det equiv model into osi and return handle
	OsiSolverInterface * loadOsiSolverData();

	// constructor 
	SmiScnOsiModel(): 
		scen_(-1),solve_synch_(false),totalProb_(0)
	{ }

	// destructor
	~SmiScnOsiModel();

	void addNode(SmiScnOsiTreeNode *node);

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
	std::vector<SmiScnOsiCoreModel *> core_vec_;
	SmiScenarioTree<SmiScnOsiTreeNode *> smiTree_;
};

// function object for addnode loop
class SmiScnOsiModelAddNode{
public:
	void operator() (SmiScnOsiTreeNode *node)
	{
		s_->addNode(node);
	}

	SmiScnOsiModelAddNode(SmiScnOsiModel *s) { s_ = s;}
private:
	SmiScnOsiModel *s_;


};


#endif //#define SmiScnModel_HPP
