// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef SmiScnData_HPP
#define SmiScnData_HPP

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "OsiSolverInterface.hpp"
#include "CoinPackedVector.hpp"

#include <map>
#include <vector>

using namespace std;


typedef int SmiCoreIndex;
typedef int SmiScenarioIndex;
typedef int SmiStageIndex;

class SmiCoreData;

class SmiNodeData
{
public:
	typedef map<int,CoinPackedVector *> SmiRowMap;
	void setCoreNode();
	inline  CoinPackedVector * getRowLower(){return drlo_;}
	inline  CoinPackedVector * getRowUpper(){return drup_;}
	inline  CoinPackedVector * getColLower(){return dclo_;}
	inline  CoinPackedVector * getColUpper(){return dcup_;}
	inline  CoinPackedVector * getObjCoefficients(){return dobj_;}

	CoinPackedVector * getRow(int i) { 
		SmiRowMap::iterator r=rowMap.find(i); 
		if (r!=rowMap.end()) return r->second; 
		else return NULL;}
	inline SmiCoreData * getCore() { return core_;}
	inline int getStage() { return stg_;}
	inline  int getNumElements(){ return nels_; }

	void copyRowLower(double * drlo);
	void copyRowUpper(double * drup);
	void copyColLower(double * dclo);
	void copyColUpper(double * dcup);
	void copyObjCoefficients(double * dobj);


	SmiNodeData(SmiStageIndex stg, SmiCoreData *core,
				 const CoinPackedMatrix *const matrix,
				 CoinPackedVector *dclo, 
				 CoinPackedVector *dcup,
				 CoinPackedVector *dobj,
				 CoinPackedVector *drlo, 
				 CoinPackedVector *drup);
	~SmiNodeData();
private:
	SmiStageIndex stg_;
	SmiRowMap rowMap;
	int nels_;
	CoinPackedVector *drlo_; 
	CoinPackedVector *drup_;
	CoinPackedVector *dobj_;
	CoinPackedVector *dclo_; 
	CoinPackedVector *dcup_;
	double *cdrlo_; 
	double *cdrup_;
	double *cdobj_;
	double *cdclo_; 
	double *cdcup_;
	SmiCoreData *core_;
	bool isCoreNode_;
};


class SmiCoreData
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

	inline SmiNodeData * getNode(SmiStageIndex t){return nodes_[t];}
	SmiCoreData(OsiSolverInterface *osi, int nstag, int *cstag, int *rstag);
	~SmiCoreData();
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
	vector<SmiNodeData*> nodes_;
};

#endif //#define SmiScnData_HPP
