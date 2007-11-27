// SmiCoreCombineRule.cpp: implementation of the SmiCoreCombineRule class.
//
//////////////////////////////////////////////////////////////////////

#include "SmiCoreCombineRule.hpp"
#include "CoinPackedVector.hpp"
#include "CoinHelperFunctions.hpp"

//////////////////////////////////////////////////////////////////////
// SmiCoreCombineReplace
//////////////////////////////////////////////////////////////////////
SmiCoreCombineReplace *SmiCoreCombineReplace::_instance = 0;

SmiCoreCombineReplace* SmiCoreCombineReplace::Instance()
{
	if (_instance==0) {
		_instance = new SmiCoreCombineReplace;
	}
	return _instance;
}

CoinPackedVector *  SmiCoreCombineReplace::Process(CoinPackedVector *cr,CoinPackedVector *nr, char *type)
{	
	
	
	CoinPackedVector *newrow=NULL;

	if (!cr && nr)
	{
		newrow = new CoinPackedVector(*nr);
	}

	if (cr && !nr)
	{
		newrow = new CoinPackedVector(*cr);
	}

	
	if (cr && nr)
	{
		// merge using denseVector

		// get max entries
		int maxentries = CoinMax(cr->getMaxIndex(),nr->getMaxIndex());
		
		double* dense = cr->denseVector(maxentries+1);
		double* elt_nr = nr->getElements();
		int* ind_nr = nr->getIndices();
		
		int j;

			for (j=0; j<nr->getNumElements(); ++j)
			{
				dense[ind_nr[j]] = elt_nr[j];
			}

		// generate new packed vector
		newrow = new CoinPackedVector();

		for (int i=0; i<maxentries+1; ++i)
		{
			if (dense[i])
				newrow->insert(i,dense[i]);
		}

		delete [] dense;
	}

	return newrow;

	
}

void SmiCoreCombineReplace::Process(double *d, int o, const CoinPackedVector &cpv, char *type)
{
	const double *cd = cpv.getElements();
	const int *ci = cpv.getIndices();
	
	for (int j=0; j < cpv.getNumElements(); ++j) 
		d[ci[j]-o] = cd[j];

}
int SmiCoreCombineReplace::Process(vector<double> *dr,CoinPackedVector *cpv,double *dels,int *indx)
{

	int numels=0;
	double *cpv_els=cpv->getElements();
	int *cpv_ind=cpv->getIndices();

	for (int j=0; j<dr->size(); ++j)
	{
		dels[numels] = (*dr)[j];
		if (*cpv_ind == j)
		{
			dels[numels] = *cpv_els;
			cpv_ind++;
			cpv_els++;
		}
		if (dels[numels])
		{
			indx[numels]=j;
			numels++;
		}
	}
	return numels;
}
//////////////////////////////////////////////////////////////////////
// SmiCoreCombineReplace
//////////////////////////////////////////////////////////////////////

SmiCoreCombineAdd *SmiCoreCombineAdd::_instance = 0;

SmiCoreCombineAdd* SmiCoreCombineAdd::Instance()
{
	if (_instance==0) {
		_instance = new SmiCoreCombineAdd;
	}
	return _instance;
}

CoinPackedVector *SmiCoreCombineAdd::Process(CoinPackedVector *cr,CoinPackedVector *nr, char *type)
{
	CoinPackedVector *newrow=NULL;

	if (!cr && nr)
	{
		newrow = new CoinPackedVector(*nr);
	}

	if (cr && !nr)
	{
		newrow = new CoinPackedVector(*cr);
	}

	
	if (cr && nr)
	{
		//newrow = new CoinPackedVector(*cr + *nr);
			// merge using denseVector

		// get max entries
		int maxentries = CoinMax(cr->getMaxIndex(),nr->getMaxIndex());
		
		double* dense = cr->denseVector(maxentries+1);
		double* elt_nr = nr->getElements();
		int* ind_nr = nr->getIndices();
		
		int j;

			for (j=0; j<nr->getNumElements(); ++j)
			{
				dense[ind_nr[j]] += elt_nr[j];
			}

		// generate new packed vector
		newrow = new CoinPackedVector();

		for (int i=0; i<maxentries+1; ++i)
		{
			if (dense[i])
				newrow->insert(i,dense[i]);
		}

		delete [] dense;
	}

	return newrow;
}

void SmiCoreCombineAdd::Process(double *d1, int o1, const CoinPackedVector &cpv2, char *type)
{
	const double *cd = cpv2.getElements();
	const int *ci = cpv2.getIndices();
	
	for (int j=0; j < cpv2.getNumElements(); ++j) 
		d1[ci[j]-o1] += cd[j];

}
int SmiCoreCombineAdd::Process(vector<double> *dr,CoinPackedVector *cpv,double *dels,int *indx)
{

	int numels=0;
	double *cpv_els=cpv->getElements();
	int *cpv_ind=cpv->getIndices();

	for (int j=0; j<dr->size(); ++j)
	{
		dels[numels] = (*dr)[j];
		if (*cpv_ind == j)
		{
			dels[numels] += *cpv_els;
			cpv_ind++;
			cpv_els++;
		}
		if (dels[numels])
		{
			indx[numels]=j;
			numels++;
		}
	}
	return numels;
}

