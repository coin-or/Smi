#include "SmiScnOsiModel.hpp"
#include "CoinHelperFunctions.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"

#include <vector>
#include <assert.h>

using namespace std;

// constructor from LP data
// TODO: allow for special node data like integer variables not in core, etc
SmiScnOsiNode::SmiScnOsiNode(SmiStageIndex stg, SmiScnOsiCoreModel *core,
				 const CoinPackedMatrix * const matrix,
				 CoinPackedVector *dclo, 
				 CoinPackedVector *dcup,
				 CoinPackedVector *dobj,
				 CoinPackedVector *drlo, 
				 CoinPackedVector *drup )
{

	core_ = core;
	stg_ = stg;
	nels_ = matrix->getNumElements();
	int i;
	int nrow = core->getNumRows(stg_);
	int ncol = core->getNumCols(stg_);


		
	if (matrix)
	{
		// should already be done but no harm checking
		assert(!matrix->isColOrdered());
		
		// pick up and store rows belonging to node's stage
		// TODO: is this a fast way to do this?
		for (i=core->getRowStart(stg); i<core->getRowStart(stg+1) ; i++)
		{
			int irow = core->getRowExternalIndex(i);
			if (irow < matrix->getNumRows())
			{
				const CoinShallowPackedVector row = matrix->getVector(irow);
				if (row.getNumElements())
				{
					CoinPackedVector *stored = new CoinPackedVector(
						row.getNumElements(), row.getIndices(), row.getElements(),
						false);
					int *indx = stored->getIndices();
					// revise indices
					for(int j=0;j<stored->getNumElements();j++)
					{
						int t= core->getColStage(indx[j]);
						indx[j] = core->getColInternalIndex(indx[j]);
						// TODO: message about row stage incompatible with col stage
						assert(!( t > stg));
					}
					//TODO: this is nice for the addNode function but is it too expensive?
					stored->sortIncrIndex();

					//TODO: is map a good container?
					this->rowMap.insert(SmiRowMap::value_type(i,stored));
				}
			}
		}
	}
	
	const int *ind;
	const double *elt;

	if (dclo)
	{
		this->dclo_ = new CoinPackedVector(false);
		dclo_->reserve(ncol);
		ind = dclo->getIndices();
		elt = dclo->getElements();
		// TODO: is this a fast way to do this?
		for (i=0; i<dclo->getNumElements(); i++)
		{
			int icol = ind[i];
			if ( core->getColStage(icol) == stg)
				this->dclo_->insert(core->getColInternalIndex(icol),elt[i]);
		}
	}
	else
		dclo_ = NULL;
	
	if (dcup)
	{
		this->dcup_ = new CoinPackedVector(false);
		dcup_->reserve(ncol);
		ind = dcup->getIndices();
		elt = dcup->getElements();
		for (i=0; i<dcup->getNumElements(); i++)
		{
			int icol = ind[i];
			if ( core->getColStage(icol) == stg)
				this->dcup_->insert(core->getColInternalIndex(icol),elt[i]);
		}
	}	
	else
		dcup_ = NULL;
	
	if (dobj)
	{
		this->dobj_ = new CoinPackedVector(false);
		dobj_->reserve(ncol);
		ind = dobj->getIndices();
		elt = dobj->getElements();
		for (i=0; i<dobj->getNumElements(); i++)
		{
			int icol = ind[i];
			if ( core->getColStage(icol) == stg)
				this->dobj_->insert(core->getColInternalIndex(icol),elt[i]);
		}
	}	
	else
		dobj_ = NULL;
	
	if (drlo)
	{
		this->drlo_ = new CoinPackedVector(false);
		drlo_->reserve(nrow);
		ind = drlo->getIndices();
		elt = drlo->getElements();
		for (i=0; i<drlo->getNumElements(); i++)
		{
			int icol = ind[i];
			if ( core->getRowStage(icol) == stg)
				this->drlo_->insert(core->getRowInternalIndex(icol),elt[i]);
		}
	}
	else
		drlo_ = NULL;

				
	if (drup)
	{
		this->drup_ = new CoinPackedVector(false);
		drlo_->reserve(nrow);
		ind = drup->getIndices();
		elt = drup->getElements();
		for (i=0; i<drup->getNumElements(); i++)
		{
			int icol = ind[i];
			if ( core->getRowStage(icol) == stg)
				this->drup_->insert(core->getRowInternalIndex(icol),elt[i]);
		}
	}
	else
		drup_ = NULL;
}

SmiScnOsiNode::~SmiScnOsiNode()
{

	if (dclo_)
		delete dclo_;
	if (dcup_)
		delete dcup_;
	if (dobj_)
		delete dobj_;
	if (drlo_)
		delete drlo_;
	if (drup_)
		delete drup_;

}

