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

int SmiScnNode::getCoreColIndex(int i)
{
    SmiCoreData *core = node_->getCore();
    return core->getColExternalIndex(i-coffset_+core->getColStart(node_->getStage()));
}
//int
//SmiScnNode::getSolverColIndex(int i)
//{
//	SmiCoreData *core = node_->getCore();
//	return core->getColExternalIndex(i-coffset_+core->getColStart(node_->getStage()));
//}

int SmiScnNode::getCoreRowIndex(int i){
    SmiCoreData *core = node_->getCore();
    return core->getRowExternalIndex(i-roffset_+core->getRowStart(node_->getStage()));
}
//int
//SmiScnNode::getSolverRowIndex(int i)
//{
//	SmiCoreData *core = node_->getCore();
//	return core->getRowInternalIndex(i-coffset_+core->getColStart(node_->getStage()));
//}

SmiScnModel::~SmiScnModel()
{
    // loop to deleteNodes
    for_each(smiTree_.treeBegin(),smiTree_.treeEnd(),SmiScnModelDeleteNode(this));

    if (osiStoch_)
        delete osiStoch_;

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

    if (maxNelsPerScenInStage)
        delete [] maxNelsPerScenInStage;

	if (this->columnNode)
		free(columnNode);
		//delete [] columnNode;

	if (this->rowNode)
		free(rowNode);
		//delete [] rowNode;

}

//Generates Tree Nodes with Data for given scenario
// branch: Stage where the new scenario differs from previous scenario
// anc: the ancestor of this scenario, this is the one it branches from
// prob: unconditional probability of this scenario (leaf node probability)
SmiScenarioIndex SmiScnModel::generateScenario(
                              CoinPackedMatrix *matrix,
                              CoinPackedVector *v_dclo, CoinPackedVector *v_dcup,
                              CoinPackedVector *v_dobj,
                              CoinPackedVector *v_drlo, CoinPackedVector *v_drup,
                              SmiStageIndex branch, SmiScenarioIndex anc, double prob,
                              SmiCoreCombineRule *r)
{
	return generateScenario(this->core_,matrix,v_dclo,v_dcup,v_dobj,v_drlo,v_drup,branch,anc,prob,r);
}
SmiScenarioIndex SmiScnModel::generateScenario(SmiCoreData *core,
                              CoinPackedMatrix *matrix,
                              CoinPackedVector *v_dclo, CoinPackedVector *v_dcup,
                              CoinPackedVector *v_dobj,
                              CoinPackedVector *v_drlo, CoinPackedVector *v_drup,
                              SmiStageIndex branch, SmiScenarioIndex anc, double prob,
                              SmiCoreCombineRule *r)
{

    // this coding takes branch to be the node that the scenario branches *from*
    --branch;

	// if core is null we use the SmiScnModel::core_
	if (!core)
		core = core_;

	assert(core != 0);

    vector<SmiScnNode *> node_vec;

    node_vec.reserve(core->getNumStages());

    // If we add the first scenario..
    if (this->getNumScenarios()==0)
    {

        assert(anc == 0);
        assert(branch == 0);
        anc = 0;
        branch = 0;
        //Add core node pointer to SmiScnModel for proper cleanup. Kind of unnice to do it in this way.
        //if (core_ == 0)
        //    core_ = core;


        // generate root node
        SmiNodeData *node = core->getNode(0);
        SmiScnNode *tnode = new SmiScnNode(node);
        tnode->setScenarioIndex(0);
        node_vec.push_back(tnode);
        this->ncol_ = core->getNumCols(0);
        this->nrow_ = core->getNumRows(0);
        this->nels_ = node->getNumMatrixElements();

		if (node->hasQdata())
		{
			this->nqels_=node->getQdata()->getNumEls();
		}
		else
		{
			this->nqels_=0;
		}

        // prepare array for counting max nels per stage and scenario
        if (this->maxNelsPerScenInStage)
            delete[] maxNelsPerScenInStage;
        this->maxNelsPerScenInStage = new int[core->getNumStages()];
        std::fill_n(maxNelsPerScenInStage, core->getNumStages(), 0); //Initialize before using it. 
        this->maxNelsPerScenInStage[0] = this->nels_;
    }
    else
    {
        // TODO: throw error if branch too large
        assert(branch<core->getNumStages());
    }

    // TODO: what to do about duplicate matrix entries?
    // ...can't do the following because eliminates zero entries...
    // matrix->eliminateDuplicates(0.0);

    //Christian: This works as follows
    // First: All nodes for the scenario from the branching stage to the leaf get created
    // Second: Connect the nodes with the present tree via the branching node, generates scenario index
    // Third: Add correct probabilities and set correct scenario index

    int t;
    for (t=branch+1; t<core->getNumStages(); t++) //Christian: Generate SmiScnNodes until the leaf nodes are reached
    {
        // generate new data node for given stage
        SmiNodeData *node = new SmiNodeData(t,core,matrix,
            v_dclo,v_dcup,v_dobj,v_drlo,v_drup);
        node->setCoreCombineRule(r);
        // generate new tree node
        SmiScnNode *tnode = new SmiScnNode(node);
        node_vec.push_back(tnode);

        this->ncol_ += core->getNumCols(t); //TODO: Why is this done for every stage? What if a new scenario branches from e.g. root node? 
        this->nrow_ += core->getNumRows(t); //Reson seems to be intersection between Tree Generation and Preparations for Det. Eq. Generation in SmiScnModel.. 
        this->nels_ += core->getNode(t)->getNumMatrixElements() + node->getNumMatrixElements(); //Maybe we have to split that..
        this->maxNelsPerScenInStage[t] = max(this->maxNelsPerScenInStage[t], core->getNode(t)->getNumMatrixElements() + node->getNumMatrixElements());

		//Qdata comes only from core nodes
		if (core->getNode(t)->hasQdata())
		{
			this->nqels_+=core->getNode(t)->getQdata()->getNumEls();
		}
    }

    //Christian: What is this method doing? Connects the newly created nodes above
    SmiScenarioIndex scen = smiTree_.addPathtoLeaf(anc,branch,node_vec);

    // add probability to all scenario nodes in path: TODO: Outsource to own method
    SmiTreeNode<SmiScnNode *> *child = smiTree_.getLeaf(scen);
    SmiTreeNode<SmiScnNode *> *parent = child->getParent();
    SmiTreeNode<SmiScnNode *> *root = smiTree_.getRoot();

    while (child != root)
    {
        SmiScnNode *tnode = child->getDataPtr();
        tnode->addProb(prob);
        tnode->setParent(parent->getDataPtr());
        if (tnode->getScenarioIndex()==-1)
            tnode->setScenarioIndex(scen);
        child = parent;
        parent = child->getParent();
    }
    root->getDataPtr()->addProb(prob);

    this->totalProb_+=prob;

    return scen;
}

SmiScenarioIndex SmiScnModel::generateScenario(
                              CoinPackedMatrix *matrix,
                              CoinPackedVector *v_dclo, CoinPackedVector *v_dcup,
                              CoinPackedVector *v_dobj,
                              CoinPackedVector *v_drlo, CoinPackedVector *v_drup,
                              vector<int> labels, double prob,
                              SmiCoreCombineRule *r)
{
	return generateScenario(this->core_,matrix,v_dclo,v_dcup,v_dobj,v_drlo,v_drup,labels,prob,r);
}

SmiScenarioIndex SmiScnModel::generateScenario(SmiCoreData *core,
                              CoinPackedMatrix *matrix,
                              CoinPackedVector *v_dclo, CoinPackedVector *v_dcup,
                              CoinPackedVector *v_dobj,
                              CoinPackedVector *v_drlo, CoinPackedVector *v_drup,
                              vector<int> labels, double prob,
                              SmiCoreCombineRule *r)
{
	
	SmiTreeNode<SmiScnNode *> *node = smiTree_.find(labels);
	SmiScenarioIndex s=generateScenario(core,matrix,v_dclo,v_dcup,v_dobj,v_drlo,v_drup,node->depth(),node->scenario(),prob,r);
	smiTree_.setChildLabels(node,labels);
	return s;
}

SmiScenarioIndex SmiScnModel::generateScenarioFromCore(SmiCoreData *core, double prob,
        SmiCoreCombineRule *r)
{
	return this->generateScenario(core,
		NULL,NULL,NULL,NULL,NULL,NULL,1,0,prob,r);
}
SmiScenarioIndex SmiScnModel::generateScenarioFromCore(double prob,
        SmiCoreCombineRule *r)
{
	return this->generateScenario(NULL,NULL,NULL,NULL,NULL,NULL,1,0,prob,r);
}

// Main public interfaces for getting WS, EV, and EEV values (no coefficients whatsoever).
double SmiScnModel::getWSValue(OsiSolverInterface *osiSolver, double objSense ){
    //Test if wsValues is empty, if this is the case, run solveWS
    std::vector< std::pair<double,double> > tempVector = solveWS(osiSolver, objSense);
    assert( tempVector.size() == smiTree_.getNumScenarios()); 
    double result = 0;
    for (unsigned int i = 0; i < tempVector.size(); i++){
        result += tempVector[i].first*tempVector[i].second;
    }
    return result;
}

double SmiScnModel::getEVValue(OsiSolverInterface *osiSolver, double objSense ){
    //Test if wsValues is empty, if this is the case, run solveWS
    std::pair<double, double*> results(solveEV(osiSolver,objSense));
    delete[] results.second;
    return results.first;
}

double SmiScnModel::getEEVValue(OsiSolverInterface* osiSolver, double objSense){
    return solveEEV(osiSolver,objSense);
}

std::pair<double,double*> SmiScnModel::solveEV(OsiSolverInterface *osiSolver, double objSense )
{
    // We create a new osiSolver Handle at this point. 
    OsiSolverInterface* tempPtr = osiStoch_;
    osiStoch_ = osiSolver->clone(false);

    //Clean up previous stuff.
    delete[] dclo_;
    delete[] dcup_;
    delete[] dobj_;
    delete[] drlo_;
    delete[] drup_;
    delete matrix_;

    // initialize arrays
    this->dclo_ = new double[this->core_->getNumCols()];
    this->dcup_ = new double[this->core_->getNumCols()];
    this->dobj_ = new double[this->core_->getNumCols()];
    this->drlo_ = new double[this->core_->getNumRows()];
    this->drup_ = new double[this->core_->getNumRows()];

    int tempNels = 0; //= nels_ ? nels_: this->core_->getNumCols()*this->core_->getNumRows(); //nels_ got set by det.eq. generations. If this is not done, we have a problem.
    for (int i = 0; i < this->core_->getNumStages(); i++) //Store this value in Core..
        tempNels += this->maxNelsPerScenInStage[i];

    // We have an empty model (zero cols,rows,els)
    ncol_=0;
    nrow_=0;
    nels_=0;
    // initialize row-ordered matrix arrays
    this->dels_ = new double[tempNels]; //TODO: Remove the +1 as it should work without it.. something is wrong there during addition of the elements..
    this->indx_ = new int[tempNels];
    this->rstrt_ = new int[this->core_->getNumRows()+1];
    this->rstrt_[0] = 0;

    // Call addNode on all these nodes
    for (int t = 0; t < this->core_->getNumStages(); t++) {
        addNode(this->core_->getNode(t));
    }

    matrix_ = new CoinPackedMatrix(false,0,0);
    int *len=NULL;
    // Assign values from current arrays (which get nulled thereafter)
    matrix_->assignMatrix(false,ncol_,nrow_,nels_,
        dels_,indx_,rstrt_,len);
    // pass data to osiStoch
    osiStoch_->loadProblem(*matrix_,dclo_,dcup_,dobj_,drlo_,drup_); //This works only for same-dimensional subproblems.

    // load integer values in solver
    for (unsigned int i = 0; i < intIndices.size(); i++) {
        osiStoch_->setInteger(intIndices[i]);
    }
    //Set objSense
    osiStoch_->setObjSense(objSense);
    osiStoch_->initialSolve(); //Solve this problem. We need objSense..

    //Delete new-ed objects
    delete matrix_;
    delete [] this->dels_;
    delete [] this->indx_;
    delete [] this->rstrt_;
    matrix_ = 0; dels_ = 0; indx_ = 0; rstrt_ = 0;

    //Print out Matrix.. if in Debug Mode
    //TODO: Define proper Debug Mode..
    //TODO: Printout obj-functions
#ifdef NDEBUG
    for (int mi = 0; mi < osiStoch_->getMatrixByRow()->getMinorDim();mi++){
        std::cout << osiStoch_->getObjCoefficients()[mi] << "*x" << mi << " + ";
    }
    std::cout << std::endl;
    for (int mi = 0; mi < osiStoch_->getMatrixByRow()->getMajorDim(); mi++){
        printf("\n%g <= ", osiStoch_->getRowLower()[mi]);
        for (int mj = 0; mj < osiStoch_->getMatrixByRow()->getMinorDim(); mj++){
            printf("%g ",osiStoch_->getMatrixByRow()->getCoefficient(mi,mj));
        }
        printf("<= %g",osiStoch_->getRowUpper()[mi]);
    }
    printf("\n Matrix printed \n");
    for (int mi = 0; mi < osiStoch_->getNumCols(); mi++) {
        printf("\n%g <= %d <= %g",osiStoch_->getColLower()[mi],mi,osiStoch_->getColUpper()[mi]);
    }
#endif
    
    double * colSolution = new double[core_->getColStart(1)];
    memcpy(colSolution, osiStoch_->getColSolution(), (core_->getColStart(1) * sizeof(double)) );
    double objValue = osiStoch_->getObjValue();
    delete osiStoch_;
    osiStoch_ = tempPtr;
    return make_pair(objValue,colSolution );
}

double SmiScnModel::solveEEV(OsiSolverInterface *osiSolver, double objSense)
{
    OsiSolverInterface* tempPtr = osiStoch_;
    osiStoch_ = osiSolver->clone(false);
    std::pair<double, double* > evResult(solveEV(osiSolver, objSense)); //We create a new osiSolver object in here (osiStoch_)
    

    // save the column solution for stage-1-columns
    int numStage1Cols = core_->getColStart(1);
    double * colSolution = evResult.second;


#pragma region EEVP
    double eev = 0; // the value of the EEVP

    //Clean up previous stuff.
    delete[] dclo_;
    delete[] dcup_;
    delete[] dobj_;
    delete[] drlo_;
    delete[] drup_;
    delete matrix_;

    // initialize arrays
    this->dclo_ = new double[this->core_->getNumCols()];
    this->dcup_ = new double[this->core_->getNumCols()];
    this->dobj_ = new double[this->core_->getNumCols()];
    this->drlo_ = new double[this->core_->getNumRows()];
    this->drup_ = new double[this->core_->getNumRows()];

    int tempNels = 0; //= nels_ ? nels_: this->core_->getNumCols()*this->core_->getNumRows(); //nels_ got set by det.eq. generations. If this is not done, we have a problem.
    for (int i = 0; i < this->core_->getNumStages(); i++)
        tempNels += this->maxNelsPerScenInStage[i];

    // loop over all scenarios and solve each of it individually TODO: This can get easily parallelized.
    for( int i = 0; i < this->smiTree_.getNumScenarios(); i++) {
        // We have an empty model (zero cols,rows,els)
        ncol_=0;
        nrow_=0;
        nels_=0;
        osiStoch_->reset();
        // initialize row-ordered matrix arrays
        this->dels_ = new double[tempNels+1]; //TODO: Remove the +1 as it should work without it.. something is wrong there during addition of the elements..
        this->indx_ = new int[tempNels];
        this->rstrt_ = new int[this->core_->getNumRows()+1];
        this->rstrt_[0] = 0;
        // Get all nodes for this scenario
        std::vector<SmiScnNode*>& nodes = smiTree_.getScenario(i);
        // Call addNode on all these nodes
        for (unsigned int i = 0; i < nodes.size(); i++) {
            addNode(nodes[i],true);
        }

        matrix_ = new CoinPackedMatrix(false,0,0);
        int *len=NULL;
        // Assign values from current arrays (which get nulled thereafter)
        matrix_->assignMatrix(false,ncol_,nrow_,nels_,
            dels_,indx_,rstrt_,len);
        // pass data to osiStoch
        osiStoch_->loadProblem(*matrix_,dclo_,dcup_,dobj_,drlo_,drup_); //This works only for same-dimensional subproblems.

        // load integer values in solver
        for (unsigned int i = 0; i < intIndices.size(); i++) {
            osiStoch_->setInteger(intIndices[i]);
        }

        // fix stage-1-columns to the values of the evp
        for (int col = 0; col < numStage1Cols; col++) {
            if (this->core_->getColStage(col) == 0) {
                osiStoch_->setColBounds(col, colSolution[col], colSolution[col]);
            }
        }

        //Set objSense
        osiStoch_->setObjSense(objSense);
        osiStoch_->initialSolve(); //Solve this problem. We need objSense..

        // sum up the solution values (multiplied with probabilities)
        eev += osiStoch_->getObjValue() * nodes[nodes.size()-1]->getProb();

        //Delete new-ed objects
        delete matrix_;
        delete [] this->dels_;
        delete [] this->indx_;
        delete [] this->rstrt_;
        matrix_ = 0; dels_ = 0; indx_ = 0; rstrt_ = 0;

        //Print out Matrix.. if in Debug Mode
        //TODO: Define proper Debug Mode..
        //TODO: Printout obj-functions
#ifdef NDEBUG
        for (int mi = 0; mi < osiStoch_->getMatrixByRow()->getMinorDim();mi++){
            std::cout << osiStoch_->getObjCoefficients()[mi] << "*x" << mi << " + ";
        }
        std::cout << std::endl;
        for (int mi = 0; mi < osiStoch_->getMatrixByRow()->getMajorDim(); mi++){
            printf("\n%g <= ", osiStoch_->getRowLower()[mi]);
            for (int mj = 0; mj < osiStoch_->getMatrixByRow()->getMinorDim(); mj++){
                printf("%g ",osiStoch_->getMatrixByRow()->getCoefficient(mi,mj));
            }
            printf("<= %g",osiStoch_->getRowUpper()[mi]);
        }
        printf("\n Matrix printed \n");
        for (int mi = 0; mi < osiStoch_->getNumCols(); mi++) {
            printf("\n%g <= %d <= %g",osiStoch_->getColLower()[mi],mi,osiStoch_->getColUpper()[mi]);

        }
        printf("\n");
#endif
    }
#pragma endregion EEVP;

    delete [] colSolution;
    delete osiStoch_;
    osiStoch_ = tempPtr;
    return eev;
}

std::vector< std::pair<double,double> > SmiScnModel::solveWS(OsiSolverInterface *osiSolver, double objSense) 
{
    OsiSolverInterface* tempPtr = osiStoch_;
    osiStoch_ = osiSolver->clone(false); // We copy the existing solverHandle, so we have to delete it at the end?

    //Clean up previous stuff.
    delete[] dclo_;
    delete[] dcup_;
    delete[] dobj_;
    delete[] drlo_;
    delete[] drup_;
    delete matrix_;

    // initialize arrays
    this->dclo_ = new double[this->core_->getNumCols()];
    this->dcup_ = new double[this->core_->getNumCols()];
    this->dobj_ = new double[this->core_->getNumCols()];
    this->drlo_ = new double[this->core_->getNumRows()];
    this->drup_ = new double[this->core_->getNumRows()];

    int tempNels = 0; //= nels_ ? nels_: this->core_->getNumCols()*this->core_->getNumRows(); //nels_ got set by det.eq. generations. If this is not done, we have a problem.
    for (int i = 0; i < this->core_->getNumStages(); i++)
        tempNels += this->maxNelsPerScenInStage[i];

    std::vector<std::pair<double,double> > solutionValues;
    solutionValues.reserve(this->smiTree_.getNumScenarios());


    // loop over all scenarios and solve each of it individually
    for( int i = 0; i < this->smiTree_.getNumScenarios(); i++) {
        // We have an empty model (zero cols,rows,els)
        ncol_=0;
        nrow_=0;
        nels_=0;
        //TODO: Create a new SolverObject by using the type of osiStoch. Need to load osiCoreData in the newly created osiStoch..
        osiStoch_->reset();
        // initialize row-ordered matrix arrays
        this->dels_ = new double[tempNels]; 
        this->indx_ = new int[tempNels];
        this->rstrt_ = new int[this->core_->getNumRows()+1];
        this->rstrt_[0] = 0;
        // Get all nodes for this scenario
        std::vector<SmiScnNode*>& nodes = smiTree_.getScenario(i);
        // Call addNode on all these nodes
        for (unsigned int i = 0; i < nodes.size(); i++) {
            addNode(nodes[i],true);
        }

        matrix_ = new CoinPackedMatrix(false,0,0);
        int *len=NULL;
        // Assign values from current arrays (which get nulled thereafter)
        matrix_->assignMatrix(false,ncol_,nrow_,nels_,
            dels_,indx_,rstrt_,len);
        // pass data to osiStoch
        osiStoch_->loadProblem(*matrix_,dclo_,dcup_,dobj_,drlo_,drup_); //This works only for same-dimensional subproblems.

        // load integer values in solver
        for (unsigned int i = 0; i < intIndices.size(); i++) {
            osiStoch_->setInteger(intIndices[i]);
        }
        //Set objSense
        osiStoch_->setObjSense(objSense);
        osiStoch_->initialSolve(); //Solve this problem. We need objSense..
        solutionValues.push_back(make_pair<double,double>(osiStoch_->getObjValue(),nodes.back()->getProb()));

        //Delete new-ed objects
        delete matrix_;
        delete [] this->dels_;
        delete [] this->indx_;
        delete [] this->rstrt_;
        matrix_ = 0; dels_ = 0; indx_ = 0; rstrt_ = 0;

        //Print out Matrix.. if in Debug Mode
        //TODO: Define proper Debug Mode..
        //TODO: Printout obj-functions
#ifdef NDEBUG
        for (int mi = 0; mi < osiStoch_->getMatrixByRow()->getMinorDim();mi++){
            std::cout << osiStoch_->getObjCoefficients()[mi] << "*x" << mi << " + ";
        }
        std::cout << std::endl;
        for (int mi = 0; mi < osiStoch_->getMatrixByRow()->getMajorDim(); mi++){
            printf("\n%g <= ", osiStoch_->getRowLower()[mi]);
            for (int mj = 0; mj < osiStoch_->getMatrixByRow()->getMinorDim(); mj++){
                printf("%g ",osiStoch_->getMatrixByRow()->getCoefficient(mi,mj));
            }
            printf("<= %g",osiStoch_->getRowUpper()[mi]);
        }
        printf("\n Matrix printed \n");
        for (int mi = 0; mi < osiStoch_->getNumCols(); mi++) {
            printf("\n%g <= %d <= %g",osiStoch_->getColLower()[mi],mi,osiStoch_->getColUpper()[mi]);
        }
        printf("\n");
#endif
    }
    delete osiStoch_;
    osiStoch_ = tempPtr;
    return solutionValues;
}

#if COIN_HAS_CLP == 1
ClpModel * SmiScnModel::loadQuadraticSolverData()
{
	generateSolverArrays();
	clp_->loadProblem(*matrix_,dclo_,dcup_,dobj_,drlo_,drup_);
	clp_->loadQuadraticObjective(matrix_->getNumCols(),this->qstart_,this->qindx_,this->qdels_);
	// load integer values in solver
    for (unsigned int i = 0; i < intIndices.size(); i++) {
        clp_->setInteger(intIndices[i]);
    }
	return clp_;
}
#endif

//Christian: Generates Deterministic Equivalent (Big Matrix)
OsiSolverInterface *
SmiScnModel::loadOsiSolverData()
{
    osiStoch_->reset();
	generateSolverArrays();
	    // pass data to osiStoch
    osiStoch_->loadProblem(*matrix_,dclo_,dcup_,dobj_,drlo_,drup_);

    // load integer values in solver
    for (unsigned int i = 0; i < intIndices.size(); i++) {
        osiStoch_->setInteger(intIndices[i]);
    }
	return osiStoch_;
}

void SmiScnModel::generateSolverArrays()
{

    delete[] dclo_;
    delete[] dcup_;
    delete[] dobj_;
    delete[] drlo_;
    delete[] drup_;
    delete matrix_;

    // initialize arrays
    this->dclo_ = new double[this->ncol_];
    this->dcup_ = new double[this->ncol_];
    this->dobj_ = new double[this->ncol_];
    this->drlo_ = new double[this->nrow_];
    this->drup_ = new double[this->nrow_];

    // initialize row-ordered matrix arrays
    this->dels_ = new double[this->nels_+1]; //TODO: Remove the +1 as it should work without it.. something is wrong there during addition of the elements..
    this->indx_ = new int[this->nels_];
    this->rstrt_ = new int[this->nrow_+1];
    this->rstrt_[0] = 0;
    this->nels_max = nels_;

	int nqels_max;
	if (this->nqels_)
	{
		this->qdels_= new double[this->nqels_];
		this->qindx_= new int[this->nqels_];
		this->qstart_= new int[this->ncol_+1];
		nqels_max=this->nqels_;
		this->nqels_=0;
		memset(this->qstart_,0,(this->ncol_+1)*sizeof(int));
		
	}

    ncol_=0;
    nrow_=0;
    nels_=0;

    // loop to addNodes
    for_each(smiTree_.treeBegin(),smiTree_.treeEnd(),SmiScnModelAddNode(this));

    //What happens if no ScenarioTree is present, but a core model is?
    matrix_ = new CoinPackedMatrix(false,0,0);
    int *len=NULL;
    matrix_->assignMatrix(false,ncol_,nrow_,nels_,
        dels_,indx_,rstrt_,len);



    //Print out Matrix.. if in Debug Mode
    //TODO: Define proper Debug Mode..
    //TODO: Printout obj-functions
#ifndef NDEBUG 
    //   for (int mi = 0; mi < osiStoch_->getMatrixByRow()->getMinorDim();mi++){
    //       std::cout << osiStoch_->getObjCoefficients()[mi] << "*x" << mi << " + ";
    //   }
    //   std::cout << std::endl;
    //for (int mi = 0; mi < osiStoch_->getMatrixByRow()->getMajorDim(); mi++){
    //	printf("\n%g <= ", osiStoch_->getRowLower()[mi]);
    //	for (int mj = 0; mj < osiStoch_->getMatrixByRow()->getMinorDim(); mj++){
    //		printf("%g ",osiStoch_->getMatrixByRow()->getCoefficient(mi,mj));
    //	}
    //	printf("<= %g",osiStoch_->getRowUpper()[mi]);
    //       //printf(" %g",this->core_->getRowStage(mi)); //We have to transform indices back to core model..
    //}
    //printf("\n Matrix printed \n");
#endif
}

//Christian: Generate Submodel Specified by Stage and Scenario Number
// This Model contains only constraints for the specified stage and scenario numbers..
//TODO: Not yet implemented
OsiSolverInterface *
SmiScnModel::loadOsiSolverDataForSubproblem(int stage, int scenStart)
{

    //TODO Declare new SolverInterface?
    osiStoch_->reset();

    delete[] dclo_;
    delete[] dcup_;
    delete[] dobj_;
    delete[] drlo_;
    delete[] drup_;
    delete matrix_;

    // initialize arrays with correct dimensions
    //this->core_->nColInStage_
    this->dclo_ = new double[this->ncol_];
    this->dcup_ = new double[this->ncol_];
    this->dobj_ = new double[this->ncol_];
    this->drlo_ = new double[this->nrow_];
    this->drup_ = new double[this->nrow_];

    // initialize row-ordered matrix arrays
    this->dels_ = new double[this->nels_];
    this->indx_ = new int[this->nels_];
    this->rstrt_ = new int[this->nrow_+1];
    this->rstrt_[0] = 0;
    this->nels_max = nels_;

    ncol_=0;
    nrow_=0;
    nels_=0;

    // loop to addNodes
    for_each(smiTree_.treeBegin(),smiTree_.treeEnd(),SmiScnModelAddNode(this));

    matrix_ = new CoinPackedMatrix(false,0,0);
    int *len=NULL;
    matrix_->assignMatrix(false,ncol_,nrow_,nels_,
        dels_,indx_,rstrt_,len);

    // pass data to osiStoch
    osiStoch_->loadProblem(CoinPackedMatrix(*matrix_),dclo_,dcup_,dobj_,drlo_,drup_);

    return osiStoch_;
    ////	SmiScnNode tnode;
    //	SmiNodeData *node = tnode->getNode();
    //
    //	// set offsets for current node
    //	tnode->setColOffset(ncol_);
    //	tnode->setRowOffset(nrow_);
    //
    //	if (tnode->isVirtualNode())
    //	{
    //		// update column count
    //		ncol_ += tnode->getNumCols();
    //		return;
    //	}
    //
    //	// OsiSolverInterface *osi = this->osiStoch_;
    //	SmiCoreData *core = node->getCore();
    //
    //	// get stage and associated core node
    //	int stg = node->getStage();
    //	SmiNodeData *cnode = core->getNode(stg);
    //
    //	// pretty sure this is an error? //Christian: At a first glance, looks like it..
    //	core->copyRowLower(drlo_+nrow_,stg);
    //	core->copyRowUpper(drup_+nrow_,stg);
    //	core->copyColLower(dclo_+ncol_,stg);
    //	core->copyColUpper(dcup_+ncol_,stg);
    //	core->copyObjective(dobj_+ncol_,stg);
    //
    //
    //	node->copyColLower(dclo_+ncol_);
    //	node->copyColUpper(dcup_+ncol_);
    //	node->copyObjective(dobj_+ncol_);
    //	node->copyRowLower(drlo_+nrow_);
    //	node->copyRowUpper(drup_+nrow_);
    //
    //	// multiply obj coeffs by node probability and normalize
    //	double prob = tnode->getProb()/this->totalProb_;
    //	tnode->setModelProb(prob);
    //
    //	for(int j=ncol_; j<ncol_+core->getNumCols(stg); ++j)
    //		dobj_[j] *= prob;
    //
    //
    //	vector<int> stochColStart(stg+1);
    //	SmiScnNode *pnode=tnode;
    //
    //	stochColStart[stg]=ncol_;
    //	for (int t=stg-1; t>0; t--)
    //	{
    //		pnode=pnode->getParent();
    //		stochColStart[t] = pnode->getColStart();
    //	}
    //
    //	// row counter
    //	int rowCount=nrow_;
    //
    //	// add rows to matrix
    //	for (int i=core->getRowStart(stg); i<core->getRowStart(stg+1) ; i++)
    //	{
    //		if (stg)
    //		{
    //			// build row explicitly into sparse arrays
    //			int rowStart=this->rstrt_[rowCount];
    //			int rowNumEls=0;
    //			if (node->getRowLength(i))
    //			{
    //				double *denseCoreRow = cnode->getDenseRow(i);
    //				rowNumEls=node->combineWithDenseCoreRow(denseCoreRow,node->getRowLength(i),node->getRowIndices(i),node->getRowElements(i),dels_+rowStart,indx_+rowStart);
    //			}
    //			else
    //			{
    //				const double *cels=cnode->getRowElements(i);
    //				const int *cind=cnode->getRowIndices(i);
    //				const int len=cnode->getRowLength(i);
    //				memcpy(dels_+rowStart,cels,sizeof(double)*len);
    //				memcpy(indx_+rowStart,cind,sizeof(int)*len);
    //				rowNumEls=len;
    //			}
    //
    //
    //			rowCount++;
    //			nels_+=rowNumEls;
    //			this->rstrt_[rowCount] = nels_;
    //
    //			// coefficients of new row
    //			int *indx = indx_+rowStart;
    //
    //			// stage starts
    //			int t=stg;
    //			int jlo=core->getColStart(stg);
    //
    //			// net offset to be added to indices
    //			int coff= stochColStart[stg]-jlo;
    //
    //			if(coff)
    //			{
    //				// TODO -- decide if better to sort COL indices in core
    //#if 1
    //
    //				// main loop iterates backwards through indices
    //				for (int j=rowNumEls-1; j>-1;--j)
    //				{
    //					// get new offset from parent node when crossing stage bndy
    //					while (indx[j]<jlo)
    //					{
    //						jlo = core->getColStart(--t);
    //						coff = stochColStart[t] - jlo;
    //					}
    //
    //					// add offset to index
    //					indx[j]+=coff;
    //				}
    //#else
    //				for (int j=0; j<rowNumEls; ++j)
    //				{
    //					t = core->getColStage(indx[j]);
    //					indx[j] += stochColStart[t] - core->getColStart(t);
    //				}
    //#endif
    //			}
    //		}
    //		else
    //		{
    //			// build row explicitly into sparse arrays
    //			const double *els = cnode->getRowElements(i);
    //			const int *ind = cnode->getRowIndices(i);
    //			const int len = cnode->getRowLength(i);
    //
    //			int rowStart=this->rstrt_[rowCount];
    //
    //			memcpy(dels_+rowStart,els,sizeof(double)*len);
    //			memcpy(indx_+rowStart,ind,sizeof(int)*len);
    //
    //			rowCount++;
    //			nels_+=len;
    //			this->rstrt_[rowCount] = nels_;
    //
    //		}
}


void
SmiScnModel::deleteNode(SmiScnNode *tnode)
{
    //cout << "Deleting node from scenario " << tnode->getScenarioIndex() << " Stage " <<  tnode->getNode()->getStage() << endl;
    delete tnode;
}

void
SmiScnModel::addNode(SmiScnNode *tnode,bool notDetEq /* = false */)
{

    SmiNodeData *node = tnode->getNode();
		
    // OsiSolverInterface *osi = this->osiStoch_;
    SmiCoreData *core = node->getCore();

	// get stage and associated core node 
    int stg = node->getStage();
    SmiNodeData *cnode = core->getNode(stg);

	this->numNodes++;
	node->setNodeIndex(this->numNodes);
    if (this->columnNode==NULL) {
	    this->columnNode = (int *)malloc(sizeof(int)*(ncol_+core->getNumCols(stg)));
    }
    else {
	    this->columnNode = (int *)realloc(this->columnNode,sizeof(int)*(ncol_+core->getNumCols(stg)));
    }
    if (this->rowNode==NULL) {
	    this->rowNode = (int *)malloc(sizeof(int)*(nrow_+core->getNumRows(stg)));
    }
    else {
	    this->rowNode = (int *)realloc(this->rowNode,sizeof(int)*(nrow_+core->getNumRows(stg)));
    }

	
    for(int j=ncol_; j<ncol_+core->getNumCols(stg); ++j)
        columnNode[j] = node->getNodeIndex();
		
    for(int j=nrow_; j<nrow_+core->getNumRows(stg); ++j)
        rowNode[j] = node->getNodeIndex();


    // set offsets for current node
    tnode->setColOffset(ncol_);
    tnode->setRowOffset(nrow_);
	

    if (tnode->isVirtualNode())
    {
        // update column count
        ncol_ += tnode->getNumCols();
        return;
    }




    //TODO: Generate integer information from information given in core node (offset manipulation)
    // Get index list of integers for current stage
    // Iterate through list 
    //    add offset
    //    add to overall integer list
    //cout << "Stage " << stg << endl;
    vector<int> intCols = core->getIntCols(stg);
    for (unsigned int i = 0; i < intCols.size(); i++) {
        //cout << intCols[i] << " + " << ncol_ << endl; // Prints out stage information for given column
        this->addIntIndice(intCols[i] + ncol_);
    }

    // pretty sure this is an error? //Christian: At a first glance, looks like it..
    //Christian: Copy and Combine with Core-Arrays to create Arrays for current scenario and stage
    //core->copyRowLower(drlo_+nrow_,stg);
    //core->copyRowUpper(drup_+nrow_,stg);
    //core->copyColLower(dclo_+ncol_,stg);
    //core->copyColUpper(dcup_+ncol_,stg);
    //core->copyObjective(dobj_+ncol_,stg);
    //Copy Core Data and Replace Stochastic Data to obtain a copy for the det. eq.
    node->copyColLower(dclo_+ncol_);
    node->copyColUpper(dcup_+ncol_);
    node->copyObjective(dobj_+ncol_);
    node->copyRowLower(drlo_+nrow_);
    node->copyRowUpper(drup_+nrow_);


    // multiply obj coeffs by node probability and normalize
    double prob = !notDetEq ? tnode->getProb()/this->totalProb_ : 1; //Christian: TODO: totalProb_ = Summer aller absoluten W'keiten (am Ende muesste die bei 1 sein..?!)
    tnode->setModelProb(prob);

    for(int j=ncol_; j<ncol_+core->getNumCols(stg); ++j)
        dobj_[j] *= prob;

	if (cnode->hasQdata())
	{
		//sanity check
		assert(this->qstart_[ncol_]==nqels_);

		// column counter
		int nqcol = ncol_;

		//get quadratic data for core node
		SmiQuadraticData *qdata = cnode->getQdata();

		
		//calculate the core node indices offset for this stage
		int coff = ncol_-core->getColStart(stg);

		//loop over core columns in the current stage
		for (int j=core->getColStart(stg);j<core->getColStart(stg+1);++j)
		{
			//loop over the elements in the quadratic data array for core column j
			for (int iels=qdata->getQDstarts()[j];iels<qdata->getQDstarts()[j+1];++iels)
			{
				// copy the element multiplied by the probability
				this->qdels_[nqels_]=qdata->getQDels()[iels]*prob;

				// copy the index with offset
				this->qindx_[nqels_]=qdata->getQDindx()[iels] + coff;

				// increment the qelement pointer
				++nqels_;
			}

			//increment the qcolumn pointer and do a column counting sanity check
			++nqcol;
			assert(nqcol == coff+j+1);

			//record the starting element for the next column 
			this->qstart_[nqcol]=nqels_;
		}
	}
		
    //Christian: Get Column Start Values for every stage less than current stage. This depends on the ordering of nodes in the nodes array..
    vector<int> stochColStart(stg+1);
    SmiScnNode *pnode=tnode;

    stochColStart[stg]=ncol_;
    for (int t=stg-1; t>0; t--)
    {
        pnode=pnode->getParent();
        stochColStart[t] = pnode->getColStart();
    }

    // row counter : initialized with number of rows counted so far
    int rowCount=nrow_;

    // add rows to det. eq. matrix for current stage
    for (int i=core->getRowStart(stg); i<core->getRowStart(stg+1) ; i++)
    {

        // build row explicitly into sparse arrays
        int rowStart=this->rstrt_[rowCount];
        int rowNumEls=0;
        //Christian: In my opinion the if/else is the same for all rows in a given stage
        //Christian: If row has stochastic coefficients (that implies that we are not in stage one which means 0 here)
        if (stg && node->getRowLength(i))
        {	//Christian: If I understood it correctly, a dense row is inserted (a row with values for all collumns)
            double *denseCoreRow = cnode->getDenseRow(i); //Christian: Why is the dense core row needed?? I don't get it..
            //I think the DenseCoreRow is needed because it simplifies things a little bit on the cost of performance
            //TODO: Change methods to use CompressedRowStorage for the core row also.. (Performance)
            //Christian: Returned row is a row that contains only non-zero elements, so it is not dense anymore
            rowNumEls=node->combineWithDenseCoreRow(denseCoreRow,node->getRowLength(i),node->getRowIndices(i),node->getRowElements(i),dels_+rowStart,indx_+rowStart);
        }
        //Christian: If row does not have stochastic coefficients, copy values from core node for current stage (no combination needed..)
        else
        {
            const double *cels=cnode->getRowElements(i);
            const int *cind=cnode->getRowIndices(i);
            const int len=cnode->getRowLength(i);
            memcpy(dels_+rowStart,cels,sizeof(double)*len);
            memcpy(indx_+rowStart,cind,sizeof(int)*len);
            rowNumEls=len;
        }

        //preparation for the next row
        rowCount++;
        nels_+=rowNumEls;
        this->rstrt_[rowCount] = nels_;
        //done with preparations for the next row

        //Coefficients of the newly added row needs to get adjusted, as they differ from the core model
        //We introduce the new variables for the det. eq. via the indx. array
        //This can only be done if we are not in the first stage aka stage 0
        if (stg) {

            int *indx = indx_+rowStart;

            // We start with the current stage
            int t=stg;
            int jlo=core->getColStart(stg);

            // net offset to be added to indices for the current stage
            int coff= stochColStart[stg]-jlo;
            if(coff)
            {
                // TODO -- decide if better to sort COL indices in core
#if 1

                // main loop iterates backwards through indices
                for (int j=rowNumEls-1; j>-1;--j)
                {
                    // get new offset from parent node when crossing stage boundary
                    // This happens, if the index of the coefficient at hand is less than the index of the first coefficient belonging to the current stage.
                    // In this case the index at hand is actually from an variable of an lower stage. Repeat procedure, until index is at the correct stage.
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

    }
    // update row, col counts
    ncol_ += core->getNumCols(stg);
    nrow_ += core->getNumRows(stg);

    // sanity check
    //assert(! ( this->nels_ > this->nels_max ) );

}

void
SmiScnModel::addNode(SmiNodeData *node)
{

    SmiCoreData *core = node->getCore();
    int stg = node->getStage();

    //Generate integer information from information given in core node (offset manipulation)
    vector<int> intCols = core->getIntCols(stg);
    for (unsigned int i = 0; i < intCols.size(); i++) {
        //cout << intCols[i] << " + " << ncol_ << endl; // Prints out stage information for given column
        this->addIntIndice(intCols[i] + ncol_);
    }


    core->copyRowLower(drlo_+nrow_,stg);
    core->copyRowUpper(drup_+nrow_,stg);
    core->copyColLower(dclo_+ncol_,stg);
    core->copyColUpper(dcup_+ncol_,stg);
    core->copyObjective(dobj_+ncol_,stg);


    // row counter : initialized with number of rows counted so far
    int rowCount=nrow_;

    // add rows to det. eq. matrix for current stage
    for (int i=core->getRowStart(stg); i<core->getRowStart(stg+1) ; i++)
    {

        // build row explicitly into sparse arrays
        int rowStart=this->rstrt_[rowCount];
        int rowNumEls=0;

        const double *cels=node->getRowElements(i);
        const int *cind=node->getRowIndices(i);
        const int len=node->getRowLength(i);
        memcpy(dels_+rowStart,cels,sizeof(double)*len);
        memcpy(indx_+rowStart,cind,sizeof(int)*len);
        rowNumEls=len;

        //preparation for the next row
        rowCount++;
        nels_+=rowNumEls;
        this->rstrt_[rowCount] = nels_;
        //done with preparations for the next row

    }
    // update row, col counts
    ncol_ += core->getNumCols(stg);
    nrow_ += core->getNumRows(stg);

    // sanity check
    //assert(! ( this->nels_ > this->nels_max ) );

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
	return this->getColValue(osiSoln,ns,length);
}
double * SmiScnModel::getColValue(const double *osiSoln, int ns, int *length)
{
    int numcols=0;

    assert( ns < this->getNumScenarios() );

    // start with leaf node
    SmiScnNode *node = this->getLeafNode(ns);
    while (node != NULL){
        // accumulate number of collumns along scenario
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
            // getCoreColIndex returns the corresponding Core index
            // in the original (user's) ordering
            dsoln[node->getCoreColIndex(j)] = osiSoln[j];
        }
        // get parent of node
        node = node->getParent();
    }
    *length=numcols;
    return dsoln;
}

double 
SmiScnModel::getColSolution(int ns, int stage, int colIndex)
{
    //CoinPackedVector *soln=new CoinPackedVector();
    const double * osiSoln = this->getOsiSolverInterface()->getColSolution();
	return this->getColValue(osiSoln,ns,stage,colIndex);
}
double SmiScnModel::getColValue(const double *osiSoln,int ns, int stage, int colIndex)
{

    assert( ns < this->getNumScenarios() );
    assert( stage < this->core_->getNumStages() );

    SmiScnNode *node = this->getLeafNode(ns);
    int curStage = this->core_->getNumStages() - stage-1; // Number of stages
    while (node != NULL && curStage){ // Go to the correct node. If Stage is reduced to zero we are at the correct node.
        // get parent of node
        node = node->getParent();
        --curStage;
    }

    assert( colIndex < core_->getColStart(stage)+node->getNumCols() && colIndex >= core_->getColStart(stage));
    return osiSoln[node->getColStart()+colIndex-core_->getColStart(stage)];// Return element in matrix that is at nodeColStart + colIndex for stage. Where to get colIndex for stage?
}

double *
SmiScnModel::getRowSolution(int ns, int *length)
{
    bool isDual=false;
    const double * osiSoln = this->getOsiSolverInterface()->getRowActivity();
	return getRowValue(osiSoln,ns,length, isDual);
}
double *SmiScnModel::getRowDuals(int ns, int *length){
	bool isDual=true;
	const double *osiDuals = this->getOsiSolverInterface()->getRowPrice();
	return getRowValue(osiDuals,ns,length,isDual);
}
double *SmiScnModel::getRowValue(const double *osiSoln, int ns, int *length, bool isDual)
{
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
			int jj= node->getCoreRowIndex(j);
			// if the Duals are called for then need to extract the probability
			// because the linear form for the duality principle incorporates probability explicitly
            double dprob=node->getProb();
			if (isDual && dprob)
				dsoln[jj] = osiSoln[j]/dprob;
			else
				dsoln[jj] = osiSoln[j];
        }
        // get parent of node
        node = node->getParent();
    }
    *length=numrows;
    return dsoln;
}

double 
SmiScnModel::getRowSolution(int ns, int stage, int rowIndex)
{
    //CoinPackedVector *soln=new CoinPackedVector();
    const double * osiSoln = this->getOsiSolverInterface()->getRowActivity();
	return this->getRowValue(osiSoln,ns,stage,rowIndex,false);
}
double SmiScnModel::getRowDuals(int ns, int stage, int rowIndex)
{
    //CoinPackedVector *soln=new CoinPackedVector();
    const double * osiSoln = this->getOsiSolverInterface()->getRowPrice();
	return this->getRowValue(osiSoln,ns,stage,rowIndex, true);
}
double SmiScnModel::getRowValue(const double *osiSoln, int ns, int stage, int rowIndex, bool isDual)
{
    assert( ns < this->getNumScenarios() );
    assert( stage < this->core_->getNumStages() );
    // start with leaf node
    SmiScnNode *node = this->getLeafNode(ns);
    int curStage = this->core_->getNumStages() - stage-1; // Number of stages
    while (node != NULL && curStage){ // Go to the correct node. If Stage is reduced to zero we are at the correct node.
        // get parent of node
        node = node->getParent();
        --curStage;
    }

    assert( rowIndex < core_->getRowStart(stage)+node->getNumRows() && rowIndex >= core_->getRowStart(stage));
	double dvalue = osiSoln[node->getRowStart()+rowIndex-core_->getRowStart(stage)];
	if (isDual && node->getProb())
		dvalue /= node->getProb();

    return dvalue;
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
    //if (osiStoch_) 
    //    smiSmpsIO->setInfinity(osiStoch_->getInfinity()); // Set solver infinity..

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

    core_ = smiCore;
    return 0;
}

int SmiScnModel::writeSmps(const char *name, bool winFileExtensions, bool strictFormat) {
    if (this != NULL) {
        SmiSmpsIO * smpsIO = new SmiSmpsIO(this->getCore(), this->getSmiTree());
        if (osiStoch_)
            smpsIO->setSolverInfinity(osiStoch_->getInfinity());
        else if (core_)
            smpsIO->setSolverInfinity(core_->getInfinity() );
        smpsIO->writeSmps(name, winFileExtensions, strictFormat);
        delete smpsIO;
        return 0;
    } else {
        std::cerr << "SMPS files for " << name << " cant be written - no stochastic model.";
        return -1;
    }
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
        nsamp[jj] = static_cast<int>(smiRV->getNumEvents());

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

        if (matrix.getNumElements()) //Christian: What is this about?! This can not work, as matrix is not yet initialized..
        {                            //Alan: matrix is initialized with CoinPackedMatrix(), so getNumElements returns 0.
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

            label[s->getStage()] *= static_cast<int>(s->getNumEvents());
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
        branch = tnode->depth()+1;
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

int
SmiScnModel::addNodeToSubmodel(SmiScnNode * smiScnNode)
{
    // stage
    int stg=smiScnNode->getStage();

    // vector of nodes
    SmiScnNode ** vecNode = (SmiScnNode **)malloc((stg+1)*sizeof(SmiScnNode*));

    // vector of labels -- the node's scenario index
    int * label = (int *) malloc((stg+1)*sizeof(int));

    // temporaries
    SmiScnNode *s=smiScnNode;

    // place incoming label into vectors offset by stage location
    vecNode[stg] = smiScnNode;
    label[stg]   = smiScnNode->getScenarioIndex();

    // fill in vector up to root
    int t=stg;
    while ((s = s->getParent()) != NULL )
    {
        vecNode[--t] = s;
        label[t]   = s->getScenarioIndex();
    }

    //assert(t==-1);

    // Search the tree for the labels
    SmiTreeNode<SmiScnNode *> *tnode = this->smiTree_.find(label,stg+1);

    // Set starting stage for new vector of SmiScnNodes
    int start_stg=0;
    if (tnode)
        start_stg=tnode->getDataPtr()->getStage()+1;

    if (start_stg<stg+1)
    {
        vector<SmiScnNode *> nwVecNode;
        SmiScnNode *nwNode;
        for (t=start_stg; t<stg+1; t++)
        {
            nwNode = new SmiScnNode(vecNode[t]->getNode());
            nwNode->setProb(vecNode[t]->getProb());
            nwVecNode.push_back(nwNode);
            nwNode->setIncludeOff();
            this->ncol_ += nwNode->getNumCols();
            this->nrow_ += nwNode->getNumRows();
            if (t)
                this->nels_ += nwNode->getNode()->getNumMatrixElements() + nwNode->getNode()->getCore()->getNode(t)->getNumMatrixElements();
            else
                this->nels_ += smiScnNode->getNode()->getNumMatrixElements();


        }
        this->smiTree_.addNodesToTree(tnode,0,nwVecNode,0);

        // Last node in path is the copy of the incoming node -- set Include on.
        nwNode->setIncludeOn();
    }

    free(vecNode);
    free(label);

    return (stg+1 - start_stg);
}





