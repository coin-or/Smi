#include "SmiScnModel.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinHelperFunctions.hpp"
#include <assert.h>
#include <algorithm>

using namespace std;


int 
SmiScnNode::getCoreColIndex(int i)
{
	SmiCoreData *core = node_->getCore();
	return core->getColExternalIndex(i-coffset_+core->getColStart(node_->getStage()));
}

int 
SmiScnNode::getCoreRowIndex(int i){
	SmiCoreData *core = node_->getCore();
	return core->getRowExternalIndex(i-roffset_+core->getRowStart(node_->getStage()));
}

SmiScnModel::~SmiScnModel()
{
	delete osiStoch_;
	
	for(unsigned int i=0; i<core_vec_.size() ; i++)
	{
		delete core_vec_[i];
	}

	if (drlo_)
		delete drlo_;
	
	if (drup_)
		delete drup_;
	
	if (dclo_)
		delete dclo_;
	
	if (dcup_)
		delete dcup_;

	if (dobj_)
		delete dobj_;

	if (matrix_)
		delete matrix_;
	

}

SmiScenarioIndex 
SmiScnModel::genScenarioReplaceCoreValues(SmiCoreIndex ic, 
				CoinPackedMatrix *matrix,
				CoinPackedVector *v_dclo, CoinPackedVector *v_dcup,
				CoinPackedVector *v_dobj,
				CoinPackedVector *v_drlo, CoinPackedVector *v_drup,				
				SmiStageIndex branch, SmiScenarioIndex anc, double prob)
{
	
	SmiCoreData *core = core_vec_[ic];
	vector<SmiScnNode *> node_vec;

	node_vec.reserve(core->getNumStages());

	// first scenario
	if (scen_==-1)
	{
		// TODO: warnings if branch and anc are not 0
		anc = 0;
		branch = 0;

		// generate root node
		SmiNodeData *node = core->getNode(0);
		SmiScnNode *tnode = new SmiScnNode(node);
		node_vec.push_back(tnode);
		this->ncol_ = core->getNumCols(0);
		this->nrow_ = core->getNumRows(0);
		this->nels_ = node->getNumElements();

	}
	else
	{
		// TODO: throw error if branch too large
		assert(branch<core->getNumStages());
	}

	// TODO: what to do about duplicate matrix entries?
	// ...can't do the following because eliminates zero entries...
	// matrix->eliminateDuplicates(0.0);

	for (int t=branch+1; t<core->getNumStages(); t++)
	{
		// generate new data node
		SmiNodeData *node = new SmiNodeData(t,core,matrix,
			v_dclo,v_dcup,v_dobj,v_drlo,v_drup);
		// generate new tree node
		SmiScnNode *tnode = new SmiScnNode(node);
		node_vec.push_back(tnode);
		
		this->ncol_ += core->getNumCols(t);
		this->nrow_ += core->getNumRows(t);
		this->nels_ += node->getNumElements();
	}

	scen_ = smiTree_.addPathtoLeaf(anc,branch,node_vec);

	// add probability to all scenario nodes
	SmiTreeNode<SmiScnNode *> *child = smiTree_.getLeaf(scen_);
	SmiTreeNode<SmiScnNode *> *parent = child->getParent();
	SmiTreeNode<SmiScnNode *> *root = smiTree_.getRoot();

	while (child != root)
	{
		SmiScnNode *tnode = child->getDataPtr();
		tnode->addProb(prob);
		tnode->setParent(parent->getDataPtr());
		child = parent;
		parent = child->getParent();
	}
	root->getDataPtr()->addProb(prob);

	this->totalProb_+=prob;

	return scen_;

}


SmiCoreIndex 
SmiScnModel::setCore(OsiSolverInterface *osi, int nstage, 
		SmiStageIndex *cstage, SmiStageIndex *rstage)
		
{

	SmiCoreData *core = new SmiCoreData(osi,nstage,cstage,rstage);

	core_vec_.push_back(core);
	return core_vec_.size()-1;
}

 
OsiSolverInterface *
SmiScnModel::loadOsiSolverData()
{
	osiStoch_->reset();

	// initialize arrays
	this->dclo_ = new double[this->ncol_];
	this->dcup_ = new double[this->ncol_];
	this->dobj_ = new double[this->ncol_];
	this->drlo_ = new double[this->nrow_];
	this->drup_ = new double[this->nrow_];
	// initialize row-ordered matrix with no extragaps or extramajors
	CoinPackedMatrix *matrix = new CoinPackedMatrix(false,0,0);
	matrix->reserve(nrow_,4*nels_);
	this->matrix_=matrix;

	ncol_=0;
	nrow_=0;
	
	// loop to addNodes
	for_each(smiTree_.treeBegin(),smiTree_.treeEnd(),SmiScnModelAddNode(this));

	// pass data to osiStoch
	osiStoch_->loadProblem(CoinPackedMatrix(*matrix_),dclo_,dcup_,dobj_,drlo_,drup_);

	return osiStoch_;
}


CoinPackedVector * ReplaceWithSecond(CoinPackedVector *cr, CoinPackedVector *nr)
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
		
		// replace entries
		for (int j=0; j<nr->getNumElements(); ++j)
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

void 
SmiScnModel::addNode(SmiScnNode *tnode)
{

	SmiNodeData *node = tnode->getNode();

	// set offsets for current node
	tnode->setColOffset(ncol_);
	tnode->setRowOffset(nrow_);

	// OsiSolverInterface *osi = this->osiStoch_;
	SmiCoreData *core = node->getCore();

	// get stage and associated core node
	int stg = node->getStage();
	SmiNodeData *cnode = core->getNode(stg);

	node->copyColLower(dclo_+ncol_);
	node->copyColUpper(dcup_+ncol_);
	node->copyObjCoefficients(dobj_+ncol_);
	node->copyRowLower(drlo_+nrow_);
	node->copyRowUpper(drup_+nrow_);
	
	// multiply obj coeffs by node probability and normalize
	double prob = tnode->getProb()/this->totalProb_;
	tnode->setModelProb(prob);

	for(int j=ncol_; j<ncol_+core->getNumCols(stg); ++j)
		dobj_[j] *= prob;
	
	
	// add rows to matrix
	for (int i=core->getRowStart(stg); i<core->getRowStart(stg+1) ; i++)
	{
		// get node rows
		CoinPackedVector *cr = cnode->getRow(i);
		if (stg)
		{
			CoinPackedVector *newrow = ReplaceWithSecond(cr,node->getRow(i));

			// TODO: this is probably a throwable error
			if (!newrow)
				continue;
			
			
			// coefficients of new row
			int *indx = newrow->getIndices();
			
			// stage starts
			int t=stg;
			int jlo=core->getColStart(stg);
			
			// net offset to be added to indices
			int coff = ncol_-jlo;
			
			if(coff)
			{
				// parent node
				SmiScnNode *pnode=tnode;
				
				// main loop iterates backwards through indices
				for (int j=newrow->getNumElements()-1; j>-1;--j)
				{
					// get new offset from parent node when crossing stage bndy
					while (indx[j]<jlo)
					{
						jlo = core->getColStart(--t);
						pnode=pnode->getParent();
						coff = pnode->getColStart()-jlo;
					}
					
					// add offset to index
					indx[j]+=coff;
				}
			}
			matrix_->appendRow(*newrow);
		}
		else
			matrix_->appendRow(*cr);

	}

	// update row, col counts
	ncol_ += core->getNumCols(stg);
	nrow_ += core->getNumRows(stg);

	// for debug sanity
	int mnrow = matrix_->getNumRows();
	int mncol = matrix_->getNumCols();
	assert(mnrow == nrow_);
	assert(mncol == ncol_);
}

OsiSolverInterface *
SmiScnModel::getOsiSolverInterface()
{
	return osiStoch_;
}

const double *
SmiScnModel::getColSolution(int ns)
{
	return NULL;
}

const double *
SmiScnModel::getRowSolution(int ns)
{
	return NULL;
}

