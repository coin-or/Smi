#include "SmiScnOsiModel.hpp"
#include "OsiClpSolverInterface.hpp"
//#include "scn_c_api.h"

#define INF 1.0e31
#define COLUMNS 0
#define ROWS    1


int main()
{
	
    /* Model dimensions */
    int ncol=27, nrow=9, nels=44;
	
	/* Sparse matrix data...organized by row */
    int *mrow,cmrow[]={ 0, 0, 0, 0, 0,
		1, 1, 1, 1,
		2, 2, 2,
		3, 3, 3, 3, 3,
		4, 4, 4, 4,
		5, 5, 5, 5, 5, 5,
		6, 6, 6, 6, 6,
		7, 7, 7, 7, 7, 7,
		8, 8, 8, 8, 8, 8 };
	  int *mcol,cmcol[]={ 0, 1, 2, 3, 4,
		5, 6, 7, 8,
		9,10, 11, 
		12, 13, 14, 15, 16, 
		0,        12, 17, 18,
		1, 5, 9,  13, 19, 20,
		2, 6,     14, 21, 22,
		3, 7, 10, 15, 23, 24,
		4, 8, 11, 16, 25, 26 };

    double dels[] = { 1.0, 1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0,
		1.0, 1.0, 1.0, 1.0, 1.0,
		16.0,              9.0, -1.0, 1.0,
		15.0, 10.0,  5.0, 11.0, -1.0, 1.0,
		28.0, 14.0,       22.0, -1.0, 1.0,
		23.0, 15.0,  7.0, 17.0, -1.0, 1.0,
		81.0, 57.0, 29.0, 55.0, -1.0, 1.0 };
	
    /* Objective */
    double *dobj,cdobj[]={ 18.0, 21.0, 18.0, 16.0, 10.0, 15.0, 16.0, 14.0, 9.0,
		10.0,  9.0,  6.0, 17.0, 16.0, 17.0, 15.0, 10.0, 0.0,
		13.0,  0.0, 13.0,  0.0,  7.0,  0.0,  7.0,  0.0, 1.0 };
	
    /* Column bounds */
    double *dclo,cdclo[]={ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
		0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
		0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 };
    double *dcup,cdcup[]={ INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,
		INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,
		INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF };
	
    /* Row bounds */
    double *drlo,cdrlo[]={ -INF, -INF, -INF, -INF,  0.0, 4.0, 0.0, 8.0, 10.0 };
    double *drup,cdrup[]={ 10.0, 19.0, 25.0, 15.0,  0.0, 7.0, 0.0, 8.0, 90.0 };
	
    /* Stages */
    int nstg=2;
    int n_first_stg_rows=4;
	int n_second_stg_rows=4;
    int n_first_stg_cols=17;
	int n_second_stg_cols=8;
    int *rstg,crstg[]={ 0,0,0,0,1,1,1,1,2 };
    int *cstg,ccstg[]={ 0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,1,
		1,1,1,1,1,1,1,2,2 };
	
    /* Stochastic data */
    int nindp=5;
    int nsamp[]={ 5, 2, 5, 5, 3 };
    double demand[]={ 200, 220, 250, 270, 300,
		50, 150,
		140, 160, 180, 200, 220,
		10, 50, 80, 100, 340,
		580, 600, 620 };
    double dprobs[]={ 0.2, 0.05, 0.35, 0.2, 0.2,
		0.3, 0.7,
		0.1, 0.2, 0.4, 0.2, 0.1,
		0.2, 0.2, 0.3, 0.2, 0.1,
		0.1, 0.8, 0.1 };
	
	/* scramble */

	int irow[]={ 1,2,7,8,0,3,4,5,6};
    int icol[]={ 9,2,3,4,5,6,7,8,1,
		19,21,23,25,0,26,24,22,20,
		10,11,12,13,14,15,16,17,18 };

    /* local variables */
    int ns=1,ii,iii,jj,*indx,*incr;
    double dp=1.0,dd;

    for (ii=0;ii<nindp;ii++) ns *= nsamp[ii];     /* Compute number of scenarios */

	// debug small sample
	// ns = 3;
	
	// initialize SmiModel
	SmiScnOsiModel *smiModel = new SmiScnOsiModel();

//	OsiClpSolverInterface osiClp;
//	OsiSolverInterface *osi = osiClp.clone(false);
	OsiClpSolverInterface *osiClp1 = new OsiClpSolverInterface();
	OsiSolverInterface *osi1 = osiClp1->clone(false);

	smiModel->setOsiSolverHandle(osi1);
	
	/* scramble LP entries */
	mrow = (int*)malloc(nels*sizeof(int));
	mcol = (int*)malloc(nels*sizeof(int));
	for (ii=0;ii<nels;ii++)
	{
		mcol[ii] = icol[cmcol[ii]];
		mrow[ii] = irow[cmrow[ii]];
	}
	
	drlo = (double *)malloc(nrow*sizeof(double));
	drup = (double *)malloc(nrow*sizeof(double));
	rstg = (int *)malloc(nrow*sizeof(int));
	for (ii=0;ii<nrow;ii++)
	{
		drlo[irow[ii]] = cdrlo[ii];
		drup[irow[ii]] = cdrup[ii];
		rstg[irow[ii]] = crstg[ii];
	}
	
	dclo = (double *)malloc(ncol*sizeof(double));
	dcup = (double *)malloc(ncol*sizeof(double));
	dobj = (double *)malloc(ncol*sizeof(double));
	cstg = (int *)malloc(ncol*sizeof(int));
	for (ii=0;ii<ncol;ii++)
	{
		dclo[icol[ii]] = cdclo[ii];
		dcup[icol[ii]] = cdcup[ii];
		dobj[icol[ii]] = cdobj[ii];
		cstg[icol[ii]] = ccstg[ii];
	}

	// this to test the matrix update stanza in genScenario
	CoinPackedMatrix *origmat = new CoinPackedMatrix(false,mrow,mcol,dels,nels);
	int corenels = nels - 4;
	
	
	// set core model using Osi interface
	OsiClpSolverInterface ocsi;
	ocsi.loadProblem(CoinPackedMatrix( 1,mrow,mcol,dels,corenels),dclo,dcup,dobj,drlo,drup);
	
	OsiSolverInterface *ohoh= ocsi.clone();
	SmiCoreIndex ic = smiModel->setCore(ohoh,3,cstg,rstg);

	// test Core Model
	SmiScnOsiCoreModel *osiCore = new SmiScnOsiCoreModel(ohoh,3,cstg,rstg);

	assert(osiCore->getNumCols(0) == n_first_stg_cols);
	assert(osiCore->getNumCols(1) == n_second_stg_cols);
	assert(osiCore->getNumCols(2) == ncol - n_first_stg_cols - n_second_stg_cols);

	assert(osiCore->getNumRows(0) == n_first_stg_rows);
	assert(osiCore->getNumRows(1) == n_second_stg_rows);
	assert(osiCore->getNumRows(2) == nrow - n_first_stg_rows - n_second_stg_rows);

	assert(osiCore->getColStart(0) == 0);
	assert(osiCore->getColStart(1) == n_first_stg_cols );
	assert(osiCore->getColStart(2) == n_first_stg_cols + n_second_stg_cols);
	assert(osiCore->getColStart(3) == ncol);
	
	assert(osiCore->getRowStart(0) == 0);
	assert(osiCore->getRowStart(1) == n_first_stg_rows );
	assert(osiCore->getRowStart(2) == n_first_stg_rows + n_second_stg_rows);
	assert(osiCore->getRowStart(3) == nrow);

	for (ii = 0; ii < n_first_stg_cols ; ii++)
		assert(cstg[osiCore->getColExternalIndex(ii)] == 0);
	for (ii = n_first_stg_cols; ii < n_first_stg_cols + n_second_stg_cols ; ii++)
		assert(cstg[osiCore->getColExternalIndex(ii)] == 1);
	for (ii = n_first_stg_cols + n_second_stg_cols; ii < ncol ; ii++)
		assert(cstg[osiCore->getColExternalIndex(ii)] == 2);

	for (ii = 0; ii < n_first_stg_rows ; ii++)
		assert(rstg[osiCore->getRowExternalIndex(ii)] == 0);
	for (ii = n_first_stg_rows; ii < n_first_stg_rows + n_second_stg_rows ; ii++)
		assert(rstg[osiCore->getRowExternalIndex(ii)] == 1);
	for (ii = n_first_stg_rows + n_second_stg_rows; ii < nrow ; ii++)
		assert(rstg[osiCore->getRowExternalIndex(ii)] == 2);

	const CoinPackedMatrix *origCore = ohoh->getMatrixByRow();
	for (int t=0;t<3;t++)
	{
		CoinPackedVector *cpvdrlo = osiCore->getRowLower(t);
		CoinPackedVector *cpvdrup = osiCore->getRowUpper(t);
		CoinPackedVector *cpvdclo = osiCore->getColLower(t);
		CoinPackedVector *cpvdcup = osiCore->getColUpper(t);
		CoinPackedVector *cpvdobj = osiCore->getObjCoefficients(t);

		double elt1,elt2;
		int ic;
		for(ii=osiCore->getColStart(t);ii<osiCore->getColStart(t+1);ii++)
		{
			ic = osiCore->getColExternalIndex(ii);
			elt1 = (*cpvdclo)[ii];
			elt2 = dclo[ic];
			assert(elt1==elt2);
			elt1 = (*cpvdcup)[ii];
			elt2 = dcup[ic];
			assert(elt1==elt2);
			elt1 = (*cpvdobj)[ii];
			elt2 = dobj[ic];
			assert(elt1==elt2);
		}
		for(ii=osiCore->getRowStart(t);ii<osiCore->getRowStart(t+1);ii++)
		{
			assert((*cpvdrlo)[ii]==drlo[osiCore->getRowExternalIndex(ii)]);
			assert((*cpvdrup)[ii]==drup[osiCore->getRowExternalIndex(ii)]);

			CoinPackedVector *row1 = osiCore->getMatrixRow(t,ii);
			const CoinPackedVector row2 = origCore->getVector(osiCore->getRowExternalIndex(ii));
			assert(row1->getNumElements() == row2.getNumElements());
			int *indx = row1->getIndices();
			double *els = row1->getElements();
			for (int j=0; j<row1->getNumElements(); j++)
			{
				elt1 = els[j];
				ic = osiCore->getColExternalIndex(indx[j]);
				elt2 = row2[ic];
				assert(elt1==elt2);
			}
		}		
	}

	printf(" *** Successfully tested SmiScnCoreModel.\n");
	
	
	
	// Coin structures for scenario updates to right hand sides
	CoinPackedVector cpv_rlo;
	CoinPackedVector cpv_rup;

	// Coin structure for scenario "updates" to core matrix
	// ..row-ordered
	CoinPackedMatrix *cpm_mat = new CoinPackedMatrix(false,mrow+corenels,mcol+corenels,dels+corenels,nels-corenels);
		
    // initialize right hand side data for first scenario
    indx = (int *) malloc( (1+nindp)*sizeof(int) );
    memset( indx,0,(1+nindp)*sizeof(int));
    for (jj=0;jj<nindp;jj++) {
		indx[jj+1] += indx[jj] + nsamp[jj];
		dp *= dprobs[ indx[jj] ];

		drlo[irow[n_first_stg_rows + jj]] = demand[ indx[jj] ];
		drup[irow[n_first_stg_rows + jj]] = demand[ indx[jj] ];
		
		cpv_rlo.insert(irow[n_first_stg_rows + jj],demand[ indx[jj] ]);
		cpv_rup.insert(irow[n_first_stg_rows + jj],demand[ indx[jj] ]);
    }
	
	// first scenario
	int	is = smiModel->genScenarioReplaceCoreValues(ic,cpm_mat,NULL,NULL,NULL,
									&cpv_rlo,&cpv_rup,0,0,dp);	
	

	// test first scenario

	// load problem data into OsiSolver
	smiModel->loadOsiSolverData();
	// get Osi pointer
	OsiSolverInterface *smiOsi1 = smiModel->getOsiSolverInterface();

	int nStochCol = smiOsi1->getNumCols();
	int nStochRow = smiOsi1->getNumRows();
	double totalProb = dp;

	// get arrays	
	const double *stochdrlo = smiOsi1->getRowLower();
	const double *stochdrup = smiOsi1->getRowUpper();
	const double *stochdclo = smiOsi1->getColLower();
	const double *stochdcup = smiOsi1->getColUpper();
	const double *stochdobj = smiOsi1->getObjCoefficients();

	// get matrix
	const CoinPackedMatrix *stochmat = smiOsi1->getMatrixByRow();
	for (t=0;t<3;t++)
	{
		double elt1,elt2;
		int ic;
		for(ii=osiCore->getColStart(t);ii<osiCore->getColStart(t+1);ii++)
		{
			ic = osiCore->getColExternalIndex(ii);
			elt1 = stochdclo[ii];
			elt2 = dclo[ic];
			assert(elt1==elt2);
			elt1 = stochdcup[ii];
			elt2 = dcup[ic];
			assert(elt1==elt2);
			elt1 = stochdobj[ii];
			elt2 = dobj[ic];
			assert(elt1==elt2);
		}
		int ir;
		for(ii=osiCore->getRowStart(t);ii<osiCore->getRowStart(t+1);ii++)
		{
		
			ir = osiCore->getRowExternalIndex(ii);

			assert(stochdrlo[ii]==drlo[ir]);
			assert(stochdrup[ii]==drup[ir]);

			const CoinPackedVector row1 = stochmat->getVector(ii);
			const CoinPackedVector row2 = origmat->getVector(ir);
			assert(row1.getNumElements() == row2.getNumElements());
			const int *indx = row1.getIndices();
			const double *els = row1.getElements();
			for (int j=0; j<row1.getNumElements(); j++)
			{
				elt1 = els[j];
				ic = osiCore->getColExternalIndex(indx[j]);
				elt2 = row2[ic];
				assert(elt1==elt2);
			}
		}		
	}

	printf(" *** Successfully tested scenario %d.\n",is);

	
	/***** ...main loop to generate scenarios from discrete random variables
		For each scenario index ii:
        If the sample size nsamp[jj] divides the scenario index ii,
		reverse the increment direction incr[jj]
		and increase the random variable index jj by 1.
        Increment the jj'th random variable by incr[jj]
		and generate new sample data.
    ***** */
	
    /* sample space increment initialized to 1 */
    incr = (int *) malloc( nindp*sizeof(int) );
    for (jj=0;jj<nindp;jj++) incr[jj] = 1;
	
    for (int iss=1;iss<ns;iss++) {
		iii=iss; jj=0;
		while ( !(iii%nsamp[jj]) ) {
			iii /= nsamp[jj];
			incr[jj] = -incr[jj];
			jj++;
		}
		dp /= dprobs[ indx[jj] ];
		indx[jj] += incr[jj];
		dp *= dprobs[ indx[jj] ];

		// set data
		drlo[irow[n_first_stg_rows + jj]] = demand[ indx[jj] ];
		drup[irow[n_first_stg_rows + jj]] = demand[ indx[jj] ];

		cpv_rlo.setElement(cpv_rlo.findIndex(irow[n_first_stg_rows + jj]),demand[ indx[jj] ]);
		cpv_rup.setElement(cpv_rup.findIndex(irow[n_first_stg_rows + jj]),demand[ indx[jj] ]);
		
		// genScenario
		is = smiModel->genScenarioReplaceCoreValues(ic,cpm_mat,NULL,NULL,NULL,
			&cpv_rlo,&cpv_rup,0,0,dp);	
		
		// test scenario
		
		// load problem data into OsiSolver
		smiModel->loadOsiSolverData();
		// get Osi pointer
		OsiSolverInterface *smiOsi1 = smiModel->getOsiSolverInterface();

		totalProb += dp;
		
		// get arrays	
		const double *stochdrlo = smiOsi1->getRowLower()+nStochRow;
		const double *stochdrup = smiOsi1->getRowUpper()+nStochRow;
		const double *stochdclo = smiOsi1->getColLower()+nStochCol;
		const double *stochdcup = smiOsi1->getColUpper()+nStochCol;
		const double *stochdobj = smiOsi1->getObjCoefficients()+nStochCol;
		
		// get matrix
		const CoinPackedMatrix *stochmat = smiOsi1->getMatrixByRow();
		for (t=1;t<3;t++)
		{
			double elt1,elt2;
			int ic;
			int colOff = osiCore->getColStart(1);
			for(ii=osiCore->getColStart(t);ii<osiCore->getColStart(t+1);ii++)
			{
				ic = osiCore->getColExternalIndex(ii);
				elt1 = stochdclo[ii-colOff];
				elt2 = dclo[ic];
				assert(elt1==elt2);
				elt1 = stochdcup[ii-colOff];
				elt2 = dcup[ic];
				assert(elt1==elt2);
				elt1 = stochdobj[ii-colOff];
				elt2 = dobj[ic];
				assert(fabs(elt1 - (elt2*dp/totalProb)) < 1.0e-8);
			}
			int ir;
			int rowOff = osiCore->getRowStart(1);
			for(ii=osiCore->getRowStart(t);ii<osiCore->getRowStart(t+1);ii++)
			{
				
				ir = osiCore->getRowExternalIndex(ii);
				
				assert(stochdrlo[ii-rowOff]==drlo[ir]);
				assert(stochdrup[ii-rowOff]==drup[ir]);
				
				const CoinPackedVector row1 = stochmat->getVector(ii);
				const CoinPackedVector row2 = origmat->getVector(ir);
				assert(row1.getNumElements() == row2.getNumElements());
				const int *indx = row1.getIndices();
				const double *els = row1.getElements();
				for (int j=0; j<row1.getNumElements(); j++)
				{
					elt1 = els[j];
					ic = osiCore->getColExternalIndex(indx[j]);
					elt2 = row2[ic];
					assert(elt1==elt2);
				}
			}		
		}
				
		nStochCol = smiOsi1->getNumCols();
		nStochRow = smiOsi1->getNumRows();

		printf(" *** Successfully tested scenario %d.\n",is);
		
	}
	
	// solve with decomp solver
//	smiModel->decompSolve();

	// load problem data into OsiSolver
	smiModel->loadOsiSolverData();
	// get Osi pointer
	OsiSolverInterface *smiOsi = smiModel->getOsiSolverInterface();
	// solve using Osi Solver
	smiOsi->initialSolve();
	// test optimal value
    assert(fabs(smiOsi->getObjValue()-1566.042)<0.01);
	return 0;
}
