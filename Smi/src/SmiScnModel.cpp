#include "SmiDiscreteDistribution.hpp"
#include "SmiScenarioTree.hpp"
#include "SmiScnModel.hpp"
#include "SmiSmpsIO.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinError.hpp"
#include "CoinPackedVector.hpp"
#include <assert.h>
#include <algorithm>

using namespace std;


int
SmiScnNode::getCoreColIndex(int i)
{
	SmiCoreData *core = node_->getCore();
	int coffset = this->getColStart();
	return core->getColExternalIndex(i-coffset+core->getColStart(node_->getStage()));
}

int
SmiScnNode::getCoreRowIndex(int i){
	SmiCoreData *core = node_->getCore();
	int roffset = this->getRowStart();
	return core->getRowExternalIndex(i-roffset+core->getRowStart(node_->getStage()));
}

SmiScnModel::~SmiScnModel()
{
	if (osiStoch_)
		delete osiStoch_;

	if (roffset_)
		delete roffset_;

	if (coffset_)
		delete coffset_;

	if (core_)
		delete core_;

	if (drlo_)
		delete [] drlo_;

	if (drup_)
		delete [] drup_;

	if (dclo_)
		delete [] dclo_;

	if (dcup_)
		delete [] dcup_;

	if (dobj_)
		delete [] dobj_;

	if (matrix_)
		delete matrix_;


}

SmiScenarioIndex
SmiScnModel::generateScenario(SmiCoreData *core,
				CoinPackedMatrix *matrix,
				CoinPackedVector *v_dclo, CoinPackedVector *v_dcup,
				CoinPackedVector *v_dobj,
				CoinPackedVector *v_drlo, CoinPackedVector *v_drup,
				SmiStageIndex branch, SmiScenarioIndex anc, double prob,
				SmiCoreCombineRule *r)
{

	// this coding takes branch to be the node that the scenario branches *from*
	--branch;

	vector<SmiScnNode *> node_vec;

	node_vec.reserve(core->getNumStages());

	// first scenario
	if (this->getNumScenarios()==0)
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
		this->nels_ = node->getNumMatrixElements();

	}
	else
	{
		// TODO: throw error if branch too large
		assert(branch<core->getNumStages());
	}

	// TODO: what to do about duplicate matrix entries?
	// ...can't do the following because eliminates zero entries...
	// matrix->eliminateDuplicates(0.0);

	int t;
	for (t=branch+1; t<core->getNumStages(); t++)
	{
		// generate new data node
		SmiNodeData *node = new SmiNodeData(t,core,matrix,
			v_dclo,v_dcup,v_dobj,v_drlo,v_drup);
		node->setCoreCombineRule(r);

		// generate new tree node
		SmiScnNode *tnode = new SmiScnNode(node);
	
		node_vec.push_back(tnode);

		this->ncol_ += core->getNumCols(t);
		this->nrow_ += core->getNumRows(t);
		this->nels_ += core->getNode(t)->getNumMatrixElements() + node->getNumMatrixElements();
	}

	int scen = smiTree_.addPathtoLeaf(anc,branch,node_vec);

	// add probability to all scenario nodes in path
	SmiTreeNode<SmiScnNode *> *child = smiTree_.getLeaf(scen);
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

	return scen;

}

SmiScenarioIndex
SmiScnModel::generateScenario(SmiCoreData *core,
				CoinPackedMatrix *matrix,
				CoinPackedVector *v_dclo, CoinPackedVector *v_dcup,
				CoinPackedVector *v_dobj,
				CoinPackedVector *v_drlo, CoinPackedVector *v_drup,
				vector<int> labels, double prob,
				SmiCoreCombineRule *r)
{

	// this code assumes that full path data (incl root node data)
	// is passed in.

	vector<SmiScnNode *> node_vec;

	node_vec.reserve(core->getNumStages());


	// TODO: what to do about duplicate matrix entries?
	// ...can't do the following because eliminates zero entries...
	// matrix->eliminateDuplicates(0.0);

	int t;
	for (t=0; t<core->getNumStages(); t++)
	{
		// generate new data node
		SmiNodeData *node = new SmiNodeData(t,core,matrix,
			v_dclo,v_dcup,v_dobj,v_drlo,v_drup);

		node->setCoreCombineRule(r);
		// generate new tree node
		SmiScnNode *tnode = new SmiScnNode(node);
		node_vec.push_back(tnode);

		this->ncol_ += core->getNumCols(t);
		this->nrow_ += core->getNumRows(t);
		this->nels_ += core->getNode(t)->getNumMatrixElements() + node->getNumMatrixElements();
	}

	SmiTreeNode<SmiScnNode *> *node = smiTree_.find(labels);
	int scen = smiTree_.addPathtoLeaf(node->scenario(),node->depth(),node_vec);
	smiTree_.setChildLabels(node,labels);

	// add probability to all scenario nodes in path
	SmiTreeNode<SmiScnNode *> *child = smiTree_.getLeaf(scen);
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

	return scen;

}





OsiSolverInterface *
SmiScnModel::loadOsiSolverData()
{
	this->setupOsiSolverData();
	this->gutsofloadOsiSolverData();
	return osiStoch_;
}

void
SmiScnModel::gutsofloadOsiSolverData()
{
	// loop to addNodes
	for_each(smiTree_.treeBegin(),smiTree_.treeEnd(),SmiScnModelAddNode(this));

	matrix_ = new CoinPackedMatrix(false,0,0);
	int *len=NULL;
	matrix_->assignMatrix(false,ncol_,nrow_,nels_,
		dels_,indx_,rstrt_,len);

	// pass data to osiStoch
	// *THINK* why do I pass the matrix this way?
	osiStoch_->loadProblem(CoinPackedMatrix(*matrix_),dclo_,dcup_,dobj_,drlo_,drup_);
}
void SmiScnModel::setupOsiSolverData()
{
	osiStoch_->reset();

	// initialize arrays
	this->dclo_ = new double[this->ncol_];
	this->dcup_ = new double[this->ncol_];
	this->dobj_ = new double[this->ncol_];
	this->drlo_ = new double[this->nrow_];
	this->drup_ = new double[this->nrow_];

	this->roffset_ = new int[this->smiTree_.getNumNodes()];
	this->coffset_ = new int[this->smiTree_.getNumNodes()];

	// initialize row-ordered matrix arrays
	this->dels_ = new double[this->nels_];
	this->indx_ = new int[this->nels_];
	this->rstrt_ = new int[this->nrow_+1];
	this->rstrt_[0] = 0;
	this->nels_max = nels_;

	ncol_=0;
	nrow_=0;
	nels_=0;
	nnodes_=0;
}


void
SmiScnModel::addNode(SmiScnNode *tnode)
{

	SmiNodeData *node = tnode->getNode();

	tnode->setNodeId(nnodes_);
	tnode->setModel(this);
	nnodes_++;

	// set offsets for current node
	this->coffset_[tnode->getNodeId()] = ncol_;
	this->roffset_[tnode->getNodeId()] = nrow_;

	// OsiSolverInterface *osi = this->osiStoch_;
	SmiCoreData *core = node->getCore();

	// get stage and associated core node
	int stg = node->getStage();
	SmiNodeData *cnode = core->getNode(stg);

	// pretty sure this is an error?
	core->copyRowLower(drlo_+nrow_,stg);
	core->copyRowUpper(drup_+nrow_,stg);
	core->copyColLower(dclo_+ncol_,stg);
	core->copyColUpper(dcup_+ncol_,stg);
	core->copyObjective(dobj_+ncol_,stg);


	node->copyColLower(dclo_+ncol_);
	node->copyColUpper(dcup_+ncol_);
	node->copyObjective(dobj_+ncol_);
	node->copyRowLower(drlo_+nrow_);
	node->copyRowUpper(drup_+nrow_);

	// multiply obj coeffs by node probability and normalize
	double prob = tnode->getProb()/this->totalProb_;
	tnode->setModelProb(prob);

	for(int j=ncol_; j<ncol_+core->getNumCols(stg); ++j)
		dobj_[j] *= prob;


	vector<int> stochColStart(stg+1);
	SmiScnNode *pnode=tnode;

	stochColStart[stg]=ncol_;
	for (int t=stg-1; t>0; t--)
	{
		pnode=pnode->getParent();
		stochColStart[t] = pnode->getColStart();
	}

	// row counter
	int rowCount=nrow_;

	// add rows to matrix
	for (int i=core->getRowStart(stg); i<core->getRowStart(stg+1) ; i++)
	{
		if (stg)
		{
			// build row explicitly into sparse arrays
			int rowStart=this->rstrt_[rowCount];
			int rowNumEls=0;
			if (node->getRowLength(i))
			{
				vector<double> *denseCoreRow = cnode->getDenseRow(i);
				rowNumEls=node->combineWithDenseCoreRow(denseCoreRow,node->getRowLength(i),node->getRowIndices(i),node->getRowElements(i),dels_+rowStart,indx_+rowStart);
			}
			else
			{
				const double *cels=cnode->getRowElements(i);
				const int *cind=cnode->getRowIndices(i);
				const int len=cnode->getRowLength(i);
				memcpy(dels_+rowStart,cels,sizeof(double)*len);
				memcpy(indx_+rowStart,cind,sizeof(int)*len);
				rowNumEls=len;
			}


			rowCount++;
			nels_+=rowNumEls;
			this->rstrt_[rowCount] = nels_;

			// coefficients of new row
			int *indx = indx_+rowStart;

			// stage starts
			int t=stg;
			int jlo=core->getColStart(stg);

			// net offset to be added to indices
			int coff= stochColStart[stg]-jlo;

			if(coff)
			{
				// TODO -- decide if better to sort COL indices in core
#if 1

				// main loop iterates backwards through indices
				for (int j=rowNumEls-1; j>-1;--j)
				{
					// get new offset from parent node when crossing stage bndy
					while (indx[j]<jlo)
					{
						jlo = core->getColStart(--t);
						coff = stochColStart[t] - jlo;
					}

					// add offset to index
					indx[j]+=coff;
				}
#else
				for (int j=0; j<rowNumEls; ++j)
				{
					t = core->getColStage(indx[j]);
					indx[j] += stochColStart[t] - core->getColStart(t);
				}
#endif
			}
		}
		else
		{
			// build row explicitly into sparse arrays
			const double *els = cnode->getRowElements(i);
			const int *ind = cnode->getRowIndices(i);
			const int len = cnode->getRowLength(i);

			int rowStart=this->rstrt_[rowCount];

			memcpy(dels_+rowStart,els,sizeof(double)*len);
			memcpy(indx_+rowStart,ind,sizeof(int)*len);

			rowCount++;
			nels_+=len;
			this->rstrt_[rowCount] = nels_;

		}


	}

	// update row, col counts
	ncol_ += core->getNumCols(stg);
	nrow_ += core->getNumRows(stg);

	// sanity check
	assert(! ( this->nels_ > this->nels_max ) );

}

OsiSolverInterface *
SmiScnModel::getOsiSolverInterface()
{
	return osiStoch_;
}


double *
SmiScnModel::getColSolution(int ns, int *length)
{
	//CoinPackedVector *soln=new CoinPackedVector();
	const double * osiSoln = this->getOsiSolverInterface()->getColSolution();
	int numcols=0;

	assert( ns < this->getNumScenarios() );

	// start with leaf node
	SmiScnNode *node = this->getLeafNode(ns);
	while (node != NULL){
		// accumulate number of rows along scenario
		numcols+=node->getNumCols();

		// get parent of node
		node = node->getParent();
	}

	// malloc vector
	double *dsoln = (double *)calloc(numcols,sizeof(double));

	// start with leaf node
	node = this->getLeafNode(ns);
	while (node != NULL){
		// copy entries
		// getColStart returns the starting index of node in OSI model
		for(int j=node->getColStart(); j<node->getColStart()+node->getNumCols(); ++j){
				// getCoreRowIndex returns the corresponding Core index
				// in the original (user's) ordering
				dsoln[node->getCoreColIndex(j)] = osiSoln[j];
		}
		// get parent of node
		node = node->getParent();
	}
	*length=numcols;
	return dsoln;
}

double *
SmiScnModel::getRowSolution(int ns, int *length)
{
	//CoinPackedVector *soln=new CoinPackedVector();
	const double * osiSoln = this->getOsiSolverInterface()->getRowActivity();
	int numrows=0;

	assert( ns < this->getNumScenarios() );

	// start with leaf node
	SmiScnNode *node = this->getLeafNode(ns);
	while (node != NULL){
		// accumulate number of rows along scenario
		numrows+=node->getNumRows();

		// get parent of node
		node = node->getParent();
	}

	// malloc vector
	double *dsoln = (double *)calloc(numrows,sizeof(double));

	// start with leaf node
	node = this->getLeafNode(ns);
	while (node != NULL){
		// copy entries
		// getRowStart returns the starting index of node in OSI model
		for(int j=node->getRowStart(); j<node->getRowStart()+node->getNumRows(); ++j){
				// getCoreRowIndex returns the corresponding Core index
				// in the original (user's) ordering
				dsoln[node->getCoreRowIndex(j)] = osiSoln[j];
		}
		// get parent of node
		node = node->getParent();
	}
	*length=numrows;
	return dsoln;
}

int
SmiScnModel::readSmps(const char *c, SmiCoreCombineRule *r)
{
	int i;
	SmiSmpsIO *smiSmpsIO=NULL;
	string fname(c);

	const char* core_ext[] = {"cor","core"};
	for (i = sizeof(core_ext)/sizeof(const char*) - 1; i >= 0; --i) {
		string ext(core_ext[i]);
		string fullname=fname+"."+ext;
		if (fileCoinReadable(fullname))
			break;
	}
	if (i == -1)
	{
		cerr << "SmiScnModel::readSmps() - No file "<< c <<" with extensions .core or .cor were found." << endl;
		return -1;
	}

	smiSmpsIO = new SmiSmpsIO();

	if (r != NULL)
		smiSmpsIO->setCoreCombineRule(r);

	if (smiSmpsIO->readMps(c,core_ext[i]) == -1)
	{
		delete smiSmpsIO;
		return -1;
	}

	SmiCoreData *smiCore = NULL;
	const char* time_ext[] = {"tim", "time"};

	for (i = sizeof(time_ext)/sizeof(const char*) - 1; i >= 0; --i) {
		string ext(time_ext[i]);
		string fullname=fname+"."+ext;
		if (fileCoinReadable(fullname))
			break;
	}
	if (i == -1)
	{
		cerr << "SmiScnModel::readSmps() - No file "<< c <<" with extensions .time or .tim were found." << endl;
		delete smiSmpsIO;
		return -1;
	}

	smiCore = smiSmpsIO->readTimeFile(this,c,time_ext[i]);
	if (!smiCore)
	{
		delete smiSmpsIO;
		return -1;
	}

	const char* stoch_ext[] = {"sto", "stoc","stoch"};
	for (i = sizeof(stoch_ext)/sizeof(const char*) - 1; i >= 0; --i) {
		string ext(stoch_ext[i]);
		string fullname=fname+"."+ext;
		if (fileCoinReadable(fullname))
			break;
	}
	if (i == -1)
	{
		cerr << "SmiScnModel::readSmps() - No file "<< c <<" with extensions .stoch, .stoc, or .sto were found." << endl;
		delete smiSmpsIO;
		return -1;
	}
	if (smiSmpsIO->readStochFile(this,smiCore,c,stoch_ext[i]) == -1)
	{
		delete smiSmpsIO;
		return -1;
	}

	delete smiSmpsIO;
	return 0;
}

void replaceFirstWithSecond(CoinPackedVector &dfirst, const CoinPackedVector &dsecond)
{
	double *delt1 = dfirst.getElements();
	const double *delt2 = dsecond.getElements();
	const int *indx2 = dsecond.getIndices();
	for(int j=0;j<dsecond.getNumElements();++j)
				delt1[dfirst.findIndex(indx2[j])] = delt2[j];
}

void
SmiScnModel::processDiscreteDistributionIntoScenarios(SmiDiscreteDistribution *smiDD, bool test)

{
	SmiCoreData *core=smiDD->getCore();

	int nindp = smiDD->getNumRV();
	int nstages = 1;

	if (test)
	{
		nstages = 3;
		assert(nindp==4);
	}
	else
	{
		nstages = core->getNumStages();
		assert(nindp > 0);
	}


	int ns=1;
	double dp=1.0;

	CoinPackedMatrix matrix ;
	CoinPackedVector cpv_dclo ;
	CoinPackedVector cpv_dcup ;
	CoinPackedVector cpv_dobj ;
	CoinPackedVector cpv_drlo ;
	CoinPackedVector cpv_drup ;

	cpv_dclo.setTestForDuplicateIndex(true);
	cpv_dcup.setTestForDuplicateIndex(true);
	cpv_dobj.setTestForDuplicateIndex(true);
	cpv_drlo.setTestForDuplicateIndex(true);
	cpv_drup.setTestForDuplicateIndex(true);

	// initialize data for first scenario
	vector<int> indx(nindp);
	vector<int> nsamp(nindp);
	vector<int> label(nstages);
	vector<int>::iterator iLabel;

	for (iLabel=label.begin(); iLabel<label.end(); ++iLabel)
		*iLabel=0;

	int jj;
	for (jj=0;jj<nindp;jj++) {

		SmiDiscreteRV *smiRV = smiDD->getDiscreteRV(jj);

		indx[jj] = 0;
		nsamp[jj] = smiRV->getNumEvents();

		assert( COIN_INT_MAX / ns > nsamp[jj] );
		ns *= nsamp[jj];

		dp *= smiRV->getEventProb(indx[jj]);

		if (test)
		{
			double p;
			p=0.5*nsamp[jj]*(nsamp[jj]+1);
			assert(smiRV->getEventProb(indx[jj])==(indx[jj]+1)/p);
		}


		cpv_dclo.append(smiRV->getEventColLower(indx[jj]));
		cpv_dcup.append(smiRV->getEventColUpper(indx[jj]));
		cpv_dobj.append(smiRV->getEventObjective(indx[jj]));
		cpv_drlo.append(smiRV->getEventRowLower(indx[jj]));
		cpv_drup.append(smiRV->getEventRowUpper(indx[jj]));

		//TODO test smiModel code
		CoinPackedMatrix m = smiRV->getEventMatrix(indx[jj]);
		if (m.getNumElements()) assert(!m.isColOrdered());

		if (matrix.getNumElements())
		{
			for (int i=0; i<m.getNumRows(); ++i)
			{
				CoinPackedVector row=m.getVector(i);
				CoinPackedVector rrow=matrix.getVector(i);
				for (int j=m.getVectorFirst(i); j<m.getVectorLast(j); ++j)
				{
					matrix.modifyCoefficient(i,j,row[j],true);
				}
			}
		}
		else
			matrix = m;

    }


	// first scenario
	int anc = 0;
	int branch = 1;
	int	is = 0;

	if (!test)
		is=this->generateScenario(core,&matrix,&cpv_dclo,&cpv_dcup,&cpv_dobj,
						&cpv_drlo,&cpv_drup,branch,anc,dp);
	else
	{
		assert(matrix.getNumElements()==4);
		assert(cpv_dclo.getNumElements()==4);
		for (int j=0;j<nindp;j++)
		{
			assert(cpv_dclo.getIndices()[j]==j);
			assert(cpv_drlo.getIndices()[j]==indx[j]);
			assert(matrix.getCoefficient(indx[j],j)==(double)(j*indx[j]));
		}
	}



	SmiTreeNode<SmiScnNode *> *root = this->smiTree_.getRoot();
	this->smiTree_.setChildLabels(root,label);

	/* sample space increment initialized to 1 */
    int *incr = (int *) malloc( nindp*sizeof(int) );
    for (jj=0;jj<nindp;jj++) incr[jj] = 1;

	/***** ...main loop to generate scenarios from discrete random variables
	For each scenario index ii:
	If the sample size nsamp[jj] divides the scenario index ii,
	reverse the increment direction incr[jj]
	and increase the random variable index jj by 1.
	Increment the jj'th random variable by incr[jj]
	and generate new sample data.
    ***** */

    for (int iss=1;iss<ns;iss++) {
		int iii=iss; jj=0;
		while ( !(iii%nsamp[jj]) ) {
			iii /= nsamp[jj];
			incr[jj] = -incr[jj];
			jj++;
		}

		SmiDiscreteRV *smiRV = smiDD->getDiscreteRV(jj);

		dp /= smiRV->getEventProb(indx[jj]);
		indx[jj] += incr[jj];
		dp *= smiRV->getEventProb(indx[jj]);

		if (test)
		{
			double p;
			p=0.5*nsamp[jj]*(nsamp[jj]+1);
			assert(smiRV->getEventProb(indx[jj])==(indx[jj]+1)/p);
		}

		for (iLabel=label.begin(); iLabel<label.end(); ++iLabel)
			*iLabel=0;

		for(int jjj=0; jjj<nindp; jjj++)
		{
			SmiDiscreteRV *s = smiDD->getDiscreteRV(jjj);

			label[s->getStage()] *= s->getNumEvents();
			label[s->getStage()] += indx[jjj];
		}


		// set data
		//TODO -- should we declare NULL entries to have 0 entries?
		//this would eliminate these tests
		replaceFirstWithSecond(cpv_dclo,smiRV->getEventColLower(indx[jj]));
		replaceFirstWithSecond(cpv_dcup,smiRV->getEventColUpper(indx[jj]));
		replaceFirstWithSecond(cpv_dobj,smiRV->getEventObjective(indx[jj]));
		replaceFirstWithSecond(cpv_drlo,smiRV->getEventRowLower(indx[jj]));
		replaceFirstWithSecond(cpv_drup,smiRV->getEventRowUpper(indx[jj]));

		//TODO test this code
		CoinPackedMatrix m = smiRV->getEventMatrix(indx[jj]);
		if (m.getNumElements()) assert(!m.isColOrdered());
		if (matrix.getNumElements())
		{
			for (int i=0; i<m.getNumRows(); ++i)
			{
				CoinPackedVector row=m.getVector(i);
				CoinPackedVector rrow=matrix.getVector(i);
				for (int j=m.getVectorFirst(i); j<m.getVectorLast(j); ++j)
				{
					matrix.modifyCoefficient(i,j,row[j],true);
				}
			}
		}
		else
			matrix = m;

		// find ancestor node
		SmiTreeNode<SmiScnNode *> *tnode = this->smiTree_.find(label);
		anc = tnode->scenario();
		branch = tnode->depth();
		if (!test)
		{
			is = this->generateScenario(core,&matrix,&cpv_dclo,&cpv_dcup,&cpv_dobj,&cpv_drlo,&cpv_drup,branch,anc,dp);
		}
		else
		{
			assert(matrix.getNumElements()==4);
			assert(cpv_dclo.getNumElements()==4);
			for (int j=0;j<nindp;j++)
			{
				assert(cpv_dclo.getIndices()[j]==j);
				assert(cpv_drlo.getIndices()[j]==indx[j]);
				assert(matrix.getCoefficient(j,indx[j])==(double)(j*indx[j]));
			}
		}

		this->smiTree_.setChildLabels(tnode,label);

	}

	free (incr);
}

double SmiScnModel::getObjectiveValue(SmiScenarioIndex ns)
{
	const double *dsoln = this->getOsiSolverInterface()->getColSolution();
	const double *dobj  = this->getOsiSolverInterface()->getObjCoefficients();

	/* calculate the scenario objective value */
	double scenSum = 0.0;

	// start with leaf node
	SmiScnNode *node = this->getLeafNode(ns);

	while (node != NULL)
	{
		double nodeSum = 0.0;
		double nodeProb = node->getModelProb();

		assert(nodeProb>0);

		// getColStart returns the starting index of node in OSI model
		for(int j=node->getColStart(); j<node->getColStart()+node->getNumCols(); ++j)
		{
			nodeSum += dobj[j]*dsoln[j];
		}
		nodeSum /= nodeProb;
		scenSum += nodeSum;

		// get parent of node
		node = node->getParent();
	}

	return scenSum;
}





