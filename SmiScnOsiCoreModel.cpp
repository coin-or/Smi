#include "SmiScnOsiModel.hpp"
#include "CoinHelperFunctions.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"

#include <vector>

using namespace std;

SmiScnOsiCoreModel::SmiScnOsiCoreModel(OsiSolverInterface *osi,int nstag,int *cstag,int *rstag)
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
	for (i=0;i<nstag_;i++)
	{
		// TODO: specialize this interface for core nodes
		CoinPackedVector *drlo = new CoinPackedVector(nrow_,osi->getRowLower());
		CoinPackedVector *drup = new CoinPackedVector(nrow_,osi->getRowUpper());
		CoinPackedVector *dclo = new CoinPackedVector(ncol_,osi->getColLower()); 
		CoinPackedVector *dcup = new CoinPackedVector(ncol_,osi->getColUpper());
		CoinPackedVector *dobj = new CoinPackedVector(ncol_,osi->getObjCoefficients());

		CoinPackedMatrix *matrix = new CoinPackedMatrix(*osi->getMatrixByRow());
		matrix->eliminateDuplicates(0.0);

		SmiScnOsiNode *node = new SmiScnOsiNode(i,this,
			matrix,dclo,dcup,dobj,drlo,drup);

		nodes_.push_back(node);

		delete matrix;
		delete drlo;
		delete drup;
		delete dclo;
		delete dcup;
		delete dobj;
	}

}


