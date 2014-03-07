// Copyright (C) 2014, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef SmiQuadratic_HPP
#define SmiQuadratic_HPP

#include "CoinPragma.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinPackedVector.hpp"

#include <map>
#include <vector>

using namespace std;


//SmiQuadratic.hpp
//
//Classes to implement quadratic SP

class SmiQuadraticData {
protected:
	int _n;
	int *_starts;
	int *_indx;
	double *_els;
	int _coff;
	bool _hasData;
public:
	SmiQuadraticData(int n, int *starts, int *indx, double *els,int coff=0): 
	  _n(n),_starts(starts),_indx(indx),_els(els),_coff(coff){
		  if (n==0 || starts==NULL || indx==NULL || els==NULL)
			  _hasData=false;
		  else
			  _hasData=true;
	  }

	SmiQuadraticData():_n(0),_starts(NULL),_indx(NULL),_els(NULL),_coff(0),_hasData(false){}

	bool hasData() {return _hasData;}
	void setHasData(bool hasData){_hasData=hasData;}

	int *getQDstarts() {return _starts;}
	void setQDstarts(int *starts) {_starts=starts;}
	int *getQDindx()   {return _indx;}
	void setQDindx(int *indx) {_indx=indx;}
	double *getQDels() {return _els;}
	void setQDels( double *els) {_els=els;}
	int getQDncols()   {return _n;}
	void setQDncols(int ncols) {_n=ncols;}
	int getNumEls()		{return _starts[_n] - _starts[0];}

};

class SmiQuadraticDataDC: public SmiQuadraticData{
public:
	SmiQuadraticDataDC(int ncols, int nels){
		this->setQDstarts(new int[ncols+1]);
		memset(this->getQDstarts(),0,(ncols+1)*sizeof(int));
		this->setQDindx(new int[nels]);
		this->setQDels(new double[nels]);
		this->setQDncols(ncols);
	}
	~SmiQuadraticDataDC() {
		delete[] _starts; delete[] _indx; delete[] _els;
	}
};



#endif