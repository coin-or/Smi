#include "SmiModel.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinHelperFunctions.hpp"
#include <assert.h>


void
SmiModel::gutsOfDelete()
{
	ekks_deleteStoch(ekkStoch_);
	ekks_endContext(env_);
	delete osiStoch_;
	
	for(unsigned int i=0; i<core_vec_.size() ; i++)
	{
		
		if (core_vec_[i]->rowIndex)
			free(core_vec_[i]->rowIndex);
		if (core_vec_[i]->colIndex)
			free(core_vec_[i]->colIndex);
		if (core_vec_[i]->colSol)
			free(core_vec_[i]->colSol);
		if (core_vec_[i]->rowSol)
			free(core_vec_[i]->rowSol);
		free(core_vec_[i]);
	}
	

}

SmiScenarioIndex 
SmiModel::genScenarioReplaceCoreValues(SmiCoreIndex ic, 
				CoinPackedMatrix *matrix,
				CoinPackedVector *v_drlo, CoinPackedVector *v_drup,
				CoinPackedVector *v_dobj,
				CoinPackedVector *v_dclo, CoinPackedVector *v_dcup,
				SmiStageIndex branch, SmiScenarioIndex anc, double prob)
{
	if (ic != 0)
	{
		//message multiple core models_ not supported in OSLSE
		abort();
	}
	SmiCoreModel *core = core_vec_[ic];
	
	int ncol = core->ncol;
	int nrow = core->nrow;
	int *mrow=NULL,*mcol=NULL;
	double *dels=NULL;
	int nels = 0;

	double *drlo = (double *)malloc(nrow*sizeof(double));
	double *drup = (double *)malloc(nrow*sizeof(double));
	double *dclo = (double *)malloc(ncol*sizeof(double));
	double *dcup = (double *)malloc(ncol*sizeof(double));
	double *dobj = (double *)malloc(ncol*sizeof(double));

	memset(drlo,0,nrow*sizeof(double));
	memset(drup,0,nrow*sizeof(double));
	memset(dclo,0,ncol*sizeof(double));
	memset(dcup,0,ncol*sizeof(double));
	memset(dobj,0,ncol*sizeof(double));



	if (matrix)
	{
		nels = matrix->getNumElements();
		mrow = (int*)malloc(nels*sizeof(int));
		dels = (double*)malloc(nels*sizeof(double));
		mcol = (int *)malloc(nels*sizeof(int));
		
		const CoinPackedMatrix *core_matrix = core->osi->getMatrixByCol();

		// assume minrow_ equals one or zero and that this value is the same for columns
		int iels = 0,j,i;
		switch(minrow_)
		{
		case 0:
			
			for (j=0; j<matrix->getNumCols(); j++)
			{
				CoinShallowPackedVector core_col = core_matrix->getVector(j);
				for (i=matrix->getVectorFirst(j); i<matrix->getVectorLast(j); i++)
				{
					int rowindx = matrix->getIndices()[i];
					mrow[iels] = rowindx + 1;
					mcol[iels] = j+1;
					dels[iels] = matrix->getElements()[i] - core_col[rowindx];
					iels++;
				}
			}
			break;
		case 1:
			
			for (j=0; j<matrix->getNumCols(); j++)
			{
				CoinShallowPackedVector core_col = core_matrix->getVector(j);
				for (i=matrix->getVectorFirst(j); i<matrix->getVectorLast(j); i++)
				{
					int rowindx = matrix->getIndices()[i];
					mrow[iels] = rowindx ;
					mcol[iels] = j;
					dels[iels] = matrix->getElements()[i] - core_col[rowindx];
					iels++;
				}
			}
			
			break;
		default:
			break;
		}
		
		assert(iels==nels);
	}

	int i;
	if (v_drlo)
	{
		for(i=0;i<v_drlo->getNumElements();i++)
		{
			int j=v_drlo->getIndices()[i];
			drlo[j] = v_drlo->getElements()[i] - core->osi->getRowLower()[j];
		}
	}
	
	if (v_drup)
	{
		for(i=0;i<v_drup->getNumElements();i++)
		{
			int j=v_drup->getIndices()[i];
			drup[j] = v_drup->getElements()[i] - core->osi->getRowUpper()[j];
		}
	}
	if (v_dclo)
	{
		for(i=0;i<v_dclo->getNumElements();i++)
		{
			int j=v_dclo->getIndices()[i];
			dclo[j] = v_dclo->getElements()[i] - core->osi->getColLower()[j];
		}
	}
	
	if (v_dcup)
	{
		for(i=0;i<v_dcup->getNumElements();i++)
		{
			int j=v_dcup->getIndices()[i];
			dcup[j] = v_dcup->getElements()[i] - core->osi->getColUpper()[j];
		}
	}
	
	if (v_dobj)
	{
		for(i=0;i<v_dobj->getNumElements();i++)
		{
			int j=v_dobj->getIndices()[i];
			dobj[j] = v_dobj->getElements()[i] - core->osi->getObjCoefficients()[j];
		}
	}

	ekks_messagePrintOn(ekkStoch_,404);

	scen_++;
	ekks_addScenario(ekkStoch_,scen_,anc,branch,prob,1,nrow,ncol,nels,
		dobj,drlo,drup,dclo,dcup,mrow,mcol,dels,ADD_TO_CORE_VALUES);
	
	
	ekks_messagePrintOff(ekkStoch_,404);

	if(mcol) delete(mcol);
	if(mrow) delete(mrow);
	if(dels) delete(dels);
	delete(drlo);
	delete(drup);
	delete(dobj);
	delete(dcup);
	delete(dclo);
	
	return scen_;
}


SmiCoreIndex 
SmiModel::setCore(OsiSolverInterface *osi, int nstage, 
		SmiStageIndex *cstage, SmiStageIndex *rstage)
		
{
	if (core_vec_.size() > 0)
	{
		// message: multiple core models_ not supported by OSLSE
		abort();
	}

	// process matrix
	const CoinPackedMatrix *matrix = osi->getMatrixByCol();
	int nels = matrix->getNumElements();
	int *mrow= (int*)malloc(nels*sizeof(int));
	double *dels = (double*)malloc(nels*sizeof(double));
	int *mcol = (int *)malloc(nels*sizeof(int));

	int ncol = matrix->getNumCols();
	int nrow = matrix->getNumRows();
	
	// find min row index

	int j;
	minrow_ = matrix->getVector(0).getMinIndex();
	for (j=1; j<ncol; j++)
		minrow_ = min(matrix->getVector(j).getMinIndex(),minrow_);
	
	// assume minrow_ equals one or zero and that this value is the same for columns
	int iels = 0;
	switch(minrow_)
	{
	case 0:
		
		for (j=0; j<ncol; j++)
		{
			for (int i=matrix->getVectorFirst(j); i<matrix->getVectorLast(j); i++)
			{
				mrow[iels] = matrix->getIndices()[i] + 1;
				mcol[iels] = j+1;
				dels[iels] = matrix->getElements()[i];
				iels++;
			}
		}
		break;
	case 1:
		
		for (j=0; j<ncol; j++)
		{
			for (int i=matrix->getVectorFirst(j); i<matrix->getVectorLast(j); i++)
			{
				mrow[iels] = matrix->getIndices()[i] ;
				mcol[iels] = j;
				dels[iels] = matrix->getElements()[i];
				iels++;
			}
		}
		ncol--;
		nrow--;
		break;
	default:
		break;
	}
	
	assert(iels==nels);

	
	// allocate and copy core structures.
	SmiCoreModel *core = (SmiCoreModel *)malloc(sizeof(SmiCoreModel));

	OsiSolverInterface *model = core->osi = osi->clone(true);
	core->nstag = nstage;
	core->cstag = (int *)malloc(ncol*sizeof(int));
	CoinDisjointCopyN(cstage,ncol,core->cstag);
	core->rstag = (int *)malloc(nrow*sizeof(int));
	CoinDisjointCopyN(rstage,nrow,core->rstag);
	core->ncol = ncol;
	core->nrow = nrow;

	core->rowIndex = NULL;
	core->colIndex = NULL;
	core->rowSol = NULL;
	core->colSol = NULL;


	core_vec_.push_back(core);
	
	double *dobj = (double *)model->getObjCoefficients();
	double *drlo = (double *)model->getRowLower();
	double *drup = (double *)model->getRowUpper();
	double *dclo = (double *)model->getColLower();
	double *dcup = (double *)model->getColUpper();
	
	int i=0;
	
	// cleanup infinities
	double infty=model->getInfinity();
	for (i=0;i<nrow;i++)
	{
		
		if (drlo[i] == - infty)
			drlo[i] = -1.0e31;
		if (drup[i] == infty)
			drup[i] = 1.0e31;
	}
	
	int *intnums = (int *)malloc(ncol*sizeof(int));
	int *intType = (int *)malloc(ncol*sizeof(int));
	int numints=0;
	for (i=0;i<ncol;i++)
	{
		
		if (dclo[i] == - infty)
			dclo[i] = -1.0e31;
		if (dcup[i] == infty)
			dcup[i] = 1.0e31;
		if (model->isInteger(i))
		{
			intnums[numints] = i+1;
			intType[numints] = 1;
			numints++;
		}
		
	}
	
	core->numints = numints;
	
	ekks_createCore(ekkStoch_,nstage,rstage,cstage,nrow,ncol,
		nels,dobj,drlo,drup,dclo,dcup,mrow,mcol,dels,NULL);
	
	if (numints)
		ekks_setIntegersAtCore(ekkStoch_, numints, intnums, intType);
	
	
	delete(intnums);
	delete(intType);

	delete(mcol);
	delete(mrow);
	delete(dels);

	return core_vec_.size()-1;
}


void
SmiModel::gutsOfConstructor(int numScen)
{

  // initialize OSLSE
  env_ = ekks_initializeContext();
  ekk_messagesPrintOff(ekk_baseModel(env_),1,9000);
  ekkStoch_=ekks_newStoch(env_,"SmiModel",numScen);

}
  
OsiSolverInterface *
SmiModel::loadOsiSolverData()
  {
	  EKKStoch *stoch = ekkStoch_;
	  ekks_describeFullModel(stoch,1); 

	  SmiCoreModel *core = core_vec_[0];
	  if (core->numints)
		  ekks_markIntegers(stoch);
	  
	  EKKModel *ekkm=ekkse_getCurrentModel(stoch);

	  ekk_mergeBlocks(ekkm,1);

	  OsiSolverInterface *osi =  osiStoch_;

	  int ncol = ekk_getInumcols(ekkm);
	  int nrow = ekk_getInumrows(ekkm);
	  const int *mcol = ekk_blockColumn(ekkm,0);
	  const int *mrow_ = ekk_blockRow(ekkm,0);
	  const double *dels_ = ekk_blockElement(ekkm,0);

	  osi->loadProblem(ncol,nrow,mcol,mrow_,dels_,
		  ekk_collower(ekkm),ekk_colupper(ekkm),ekk_objective(ekkm),
		  ekk_rowlower(ekkm),ekk_rowupper(ekkm));

	  osi->setInteger(ekk_listOfIntegers(ekkm),ekk_getInumints(ekkm));

	  solve_synch_ = false;

	  return osi;
	 
  }

void
SmiModel::decompSolve()
{
	// solve with L-shaped
	EKKStoch *stoch = ekkStoch_;
	
	
	
	ekks_messagePrintOn(stoch,420);
	ekks_messagePrintOn(stoch,1);
	ekks_messagePrintOn(stoch,445);
	ekks_messagePrintOn(stoch,446);
	ekks_messagePrintOn(stoch,447);
	ekks_messagePrintOn(stoch,425);
	ekks_messagePrintOn(stoch,452);
	ekks_messagePrintOn(stoch,453);
	
	
	ekks_describeFullModel(stoch,0);                      /* Describe Full model and populate with data */
	
	ekks_setScaleOn(stoch);
	ekks_setPresolve(stoch,3);
	ekks_setCrash(stoch,3);
	ekks_setFinalSimplexSolverOff(stoch);
	ekks_setMaxsubmodels(stoch,200);
	ekks_setMinNumRows(stoch,500);
	
	ekkse_nestedLSolve(stoch,CUTBYROWSIZE,0);                /* Solve using nestedLSolve() -- (parm2=1)CUTBYROWSIZE, (parm3=0) No Hotstart */
	
	EKKModel *ekkm=ekkse_getCurrentModel(ekkStoch_);
	
	double *ekksol = ekk_getColsol(ekkm);
	
	osiStoch_->setColSolution(ekksol);
	ekk_free(ekksol);
	
	solve_synch_ = true;

	ekks_messagePrintOff(stoch,420);
	ekks_messagePrintOff(stoch,1);
	ekks_messagePrintOff(stoch,445);
	ekks_messagePrintOff(stoch,446);
	ekks_messagePrintOff(stoch,447);
	ekks_messagePrintOff(stoch,425);
	ekks_messagePrintOff(stoch,452);
	ekks_messagePrintOff(stoch,453);

	
}


OsiSolverInterface *
SmiModel::getOsiSolverInterface()
{
	return osiStoch_;
}

const double *
SmiModel::getColSolution(int ns)
{
	SmiCoreModel *core = core_vec_[0];
	
	if (core->colIndex == NULL)
	{
		// malloc arrays for getScenarioSolution
		core->colIndex = (int *)malloc(core->ncol*sizeof(int));
		core->colSol = (double *)malloc(core->ncol*sizeof(double));
	}
	
	if (!solve_synch_)
	{
		//kludge:
		// update solution vector in ekkStoch_ after solve
		
		EKKModel *ekkm=ekkse_getCurrentModel(ekkStoch_);
		double *ekksol = ekk_getColsol(ekkm);
		memcpy(ekksol,osiStoch_->getColSolution(),sizeof(double)*osiStoch_->getNumCols());
		ekk_setColsol(ekkm,ekksol);
		ekk_free(ekksol);
		
		solve_synch_ = true;
	}

	ekks_getScenarioSolution(ekkStoch_,ns,0,core->colSol,core->colIndex);
	return core->colSol;
}

const double *
SmiModel::getRowSolution(int ns)
{
	SmiCoreModel *core = core_vec_[0];

	if (core->rowIndex == NULL)
	{
		// malloc arrays for getScenarioSolution
		core->rowIndex = (int *)malloc(core->nrow*sizeof(int));
		core->rowSol = (double *)malloc(core->nrow*sizeof(double));
	}

	if (!solve_synch_)
	{
		//kludge:
		// update solution vector in ekkStoch_ after solve
		
		EKKModel *ekkm=ekkse_getCurrentModel(ekkStoch_);
		double *ekksol = ekk_getColsol(ekkm);
		memcpy(ekksol,osiStoch_->getColSolution(),sizeof(double)*osiStoch_->getNumCols());
		ekk_setColsol(ekkm,ekksol);
		ekk_free(ekksol);
		
		solve_synch_ = true;
	}

	ekks_getScenarioSolution(ekkStoch_,ns,1,core->rowSol,core->rowIndex);

	return core->rowSol;
}


void 
SmiModel::writeSMPS(const char *c,const char *t, const char* s)
{	  ekks_outMatrixSMPS(ekkStoch_,SCENARIOS,REPLACE_CORE_VALUES,c,t,s);
}

void 
SmiModel::readSMPS(const char *c,const char *t, const char* s)
{	  ekks_readSMPSData(ekkStoch_,c,t,s);
}
