#include "SmiScnData.hpp"
#include "CoinHelperFunctions.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"

#include <vector>

using namespace std;

SmiCoreData::SmiCoreData(OsiSolverInterface *osi,int nstag,int *cstag,int *rstag)
{
	int i;
	nrow_ = osi->getNumRows();
	ncol_ = osi->getNumCols();

	// store number stages
	nstag_ = nstag;

	nColInStage_ = new int[nstag_+1];
	nRowInStage_ = new int[nstag_+1];


	colStage_ = new int[ncol_];
	colEx2In_ = new int[ncol_];
	colIn2Ex_ = new int[ncol_];

	rowStage_ = new int[nrow_];
	rowEx2In_ = new int[nrow_];
	rowIn2Ex_ = new int[nrow_];
	
	// store stage maps
	CoinDisjointCopyN(cstag,ncol_,colStage_);
	CoinDisjointCopyN(rstag,nrow_,rowStage_);

	// zero stage counts
	for (i=0;i<nstag_+1;i++)
	{
		nColInStage_[i] = 0;
		nRowInStage_[i] = 0;
	}
		
	// array to point to start of new stage
	stageRowPtr_ = new int[nstag_+1];

	// count rows in each stage
	for (i=0;i<nrow_;i++)
		nRowInStage_[rowStage_[i]]++;

	// set stage pointers
	stageRowPtr_[0] = 0;
	for (i=0;i<nstag_;i++)
		stageRowPtr_[i+1] = stageRowPtr_[i] + nRowInStage_[i];

	// place index into next open position in its stage
	for (i=0;i<nrow_;i++)
	{
		rowEx2In_[i] = stageRowPtr_[rowStage_[i]];
		rowIn2Ex_[rowEx2In_[i]] = i;
		stageRowPtr_[rowStage_[i]]++;
	}

	// reset stage pointers
	stageRowPtr_[0] = 0;
	for (i=0;i<nstag_;i++)
		stageRowPtr_[i+1] = stageRowPtr_[i] + nRowInStage_[i];

			
	// array to point to start of new stage
	stageColPtr_ = new int[nstag_+1];

	// count cols in each stage
	for (i=0;i<ncol_;i++)
		nColInStage_[colStage_[i]]++;

	// set stage pointers
	stageColPtr_[0] = 0;
	for (i=0;i<nstag_;i++)
		stageColPtr_[i+1] = stageColPtr_[i] + nColInStage_[i];

	// place index into next open position in its stage
	for (i=0;i<ncol_;i++)
	{
		colEx2In_[i] = stageColPtr_[colStage_[i]];
		colIn2Ex_[colEx2In_[i]] = i;
		stageColPtr_[colStage_[i]]++;
	}

	// reset stage pointers
	stageColPtr_[0] = 0;
	for (i=0;i<nstag_;i++)
		stageColPtr_[i+1] = stageColPtr_[i] + nColInStage_[i];


	// make nodes

	this->nodes_.reserve(nstag_);

		// TODO: specialize this interface for core nodes
	CoinPackedVector *drlo = new CoinPackedVector(nrow_,osi->getRowLower());
	CoinPackedVector *drup = new CoinPackedVector(nrow_,osi->getRowUpper());
	CoinPackedVector *dclo = new CoinPackedVector(ncol_,osi->getColLower()); 
	CoinPackedVector *dcup = new CoinPackedVector(ncol_,osi->getColUpper());
	CoinPackedVector *dobj = new CoinPackedVector(ncol_,osi->getObjCoefficients());

	CoinPackedMatrix *matrix = new CoinPackedMatrix(*osi->getMatrixByRow());
	matrix->eliminateDuplicates(0.0);


	for (i=0;i<nstag_;i++)
	{
	
		SmiNodeData *node = new SmiNodeData(i,this,
			matrix,dclo,dcup,dobj,drlo,drup);

		node->setCoreNode();

		nodes_.push_back(node);

	}

	
		delete matrix;
		delete drlo;
		delete drup;
		delete dclo;
		delete dcup;
		delete dobj;


}

SmiCoreData::~SmiCoreData()
{
	delete nColInStage_;
	delete nRowInStage_;
	delete stageColPtr_;
	delete stageRowPtr_;
	delete colStage_;
	delete rowStage_;
	delete colEx2In_;
	delete rowEx2In_;
	delete colIn2Ex_;
	delete rowIn2Ex_;
}

void 
SmiNodeData::setCoreNode()
{

	isCoreNode_=true;

	cdclo_ = dclo_->denseVector(core_->getNumCols())+core_->getColStart(stg_);
	cdcup_ = dcup_->denseVector(core_->getNumCols())+core_->getColStart(stg_);
	cdobj_ = dobj_->denseVector(core_->getNumCols())+core_->getColStart(stg_);
	cdrlo_ = drlo_->denseVector(core_->getNumRows())+core_->getRowStart(stg_);
	cdrup_ = drup_->denseVector(core_->getNumRows())+core_->getRowStart(stg_);

}

#include <vector>
#include <assert.h>

using namespace std;

// constructor from LP data
// TODO: allow for special node data like integer variables not in core, etc
SmiNodeData::SmiNodeData(SmiStageIndex stg, SmiCoreData *core,
				 const CoinPackedMatrix * const matrix,
				 CoinPackedVector *dclo, 
				 CoinPackedVector *dcup,
				 CoinPackedVector *dobj,
				 CoinPackedVector *drlo, 
				 CoinPackedVector *drup )
{
	// initialize specialized core node info
	isCoreNode_ = false;
	cdrlo_= NULL; 
	cdrup_= NULL;
	cdobj_= NULL;
	cdclo_= NULL; 
	cdcup_= NULL;

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

void NodeCopyToDouble(double *d, CoinPackedVector *cpv, int o)
{	
	double *cd = cpv->getElements();
	int *ci = cpv->getIndices();
	int j=0;
	while(j < cpv->getNumElements())
		d[ci[j]-o] = cd[j++];
}

void SmiNodeData::copyRowLower(double * drlo)
{
	if (isCoreNode_)
	{
		CoinDisjointCopyN(cdrlo_,core_->getNumRows(stg_),drlo);
	}
	else
	{
		SmiNodeData *cnode = core_->getNode(stg_);
		cnode->copyRowLower(drlo);
		if (drlo_)
			NodeCopyToDouble(drlo,drlo_,core_->getRowStart(stg_));
	}
}

void SmiNodeData::copyRowUpper(double * drup){
	if (isCoreNode_)
	{
		CoinDisjointCopyN(cdrup_,core_->getNumRows(stg_),drup);
	}
	else
	{
		SmiNodeData *cnode = core_->getNode(stg_);
		cnode->copyRowUpper(drup);
		if (drup_)
			NodeCopyToDouble(drup,drup_,core_->getRowStart(stg_));
	}
}

void SmiNodeData::copyColLower(double * dclo){
	if (isCoreNode_)
	{
		CoinDisjointCopyN(cdclo_,core_->getNumCols(stg_),dclo);
	}
	else
	{
		SmiNodeData *cnode = core_->getNode(stg_);
		cnode->copyColLower(dclo);
		if (dclo_)
			NodeCopyToDouble(dclo,dclo_,core_->getColStart(stg_));
	}
}

void SmiNodeData::copyColUpper(double * dcup){
	if (isCoreNode_)
	{
		CoinDisjointCopyN(cdcup_,core_->getNumCols(stg_),dcup);
	}
	else
	{
		SmiNodeData *cnode = core_->getNode(stg_);
		cnode->copyColUpper(dcup);
		if (dcup_)
			NodeCopyToDouble(dcup,dcup_,core_->getColStart(stg_));
	}
}

void SmiNodeData::copyObjCoefficients(double * dobj){
	if (isCoreNode_)
	{
		CoinDisjointCopyN(cdobj_,core_->getNumCols(stg_),dobj);
	}
	else
	{
		SmiNodeData *cnode = core_->getNode(stg_);
		cnode->copyObjCoefficients(dobj);
		if (dobj_)
			NodeCopyToDouble(dobj,dobj_,core_->getColStart(stg_));
	}
}


SmiNodeData::~SmiNodeData()
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

	if (cdclo_)
	{
		cdclo_ -= core_->getColStart(stg_);
		delete cdclo_;
	}
	if (cdcup_)
	{
		cdcup_ -= core_->getColStart(stg_);
		delete cdcup_;
	}
	if (cdobj_)
	{
		cdobj_ -= core_->getColStart(stg_);
		delete cdobj_;
	}
	if (cdrlo_)
	{
		cdrlo_ -= core_->getRowStart(stg_);
		delete cdrlo_;
	}
	if (cdrup_)
	{
		cdrup_ -= core_->getRowStart(stg_);
		delete cdrup_;
	}

}

