#include "SmiModel.hpp"
#include "OsiClpSolverInterface.hpp"
#include "scn_c_api.h"

#define INF 1.0e31
#define COLUMNS 0
#define ROWS    1

#if 0
int main ()
{
  
    /* Model dimensions */
    int ncol=27, nrow=9, nels=44;

    /* Sparse matrix data...organized by row */
    int *mrow,cmrow[]={ 1, 1, 1, 1, 1,
            2, 2, 2, 2,
            3, 3, 3,
            4, 4, 4, 4, 4,
            5, 5, 5, 5,
            6, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7,
            8, 8, 8, 8, 8, 8,
            9, 9, 9, 9, 9, 9 };
    int *mcol,cmcol[]={ 1, 2, 3, 4, 5,
            6, 7, 8, 9,
            10, 11, 12,
            13, 14, 15, 16, 17,
            1,        13, 18, 19,
            2, 6, 10, 14, 20, 21,
            3, 7,     15, 22, 23,
            4, 8, 11, 16, 24, 25,
            5, 9, 12, 17, 26, 27 };
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
    double *drlo,cdrlo[]={ -INF, -INF, -INF, -INF,  1.0, 2.0, 4.0, 6.0, 80.0 };
    double *drup,cdrup[]={ 10.0, 19.0, 25.0, 15.0,  1.0, 3.0, 5.0, 7.0, 90.0 };

    /* Stages */
    int nstg=2;
    int n_first_stg_rows=4;
    int n_first_stg_cols=17;
    int *rstg,crstg[]={ 1,1,1,1,2,2,2,2,2 };
    int *cstg,ccstg[]={ 1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,2,
            2,2,2,2,2,2,2,2,2 };

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

    int irow[]={ 9,8,7,6,5,4,3,2,1 };

    int icol[]={ 9,2,3,4,5,6,7,8,1,

                                 19,21,23,25,27,26,24,22,20,

                                 10,11,12,13,14,15,16,17,18 };
#if 0
        /* scramble */
    int irow[]={ 5,4,3,2,9,8,7,6,1 };
    int icol[]={ 9,2,3,4,5,6,7,8,1,23,25,13,14,15,16,
				27,26,19,21,24,22,20,10,11,12,17,18 };
#endif

    /* local variables */
    int ns=1,ii,iii,jj,*indx,*incr;
    double dd,dp=1.0;
	
	

	// compute number of scenarios
    for (ii=0;ii<nindp;ii++) ns *= nsamp[ii]; 
	
	// initialize SmiModel and set solver
	SmiModel *stoch = new SmiModel(ns);
	OsiClpSolverInterface *osi = new OsiClpSolverInterface();
	stoch->setOsiSolverHandle(osi);

	/* scramble indices to test re-ordering */
	mrow = (int*)malloc(nels*sizeof(int));
	mcol = (int*)malloc(nels*sizeof(int));
	for (ii=0;ii<nels;ii++)
	{
		mcol[ii] = icol[cmcol[ii]-1];
		mrow[ii] = irow[cmrow[ii]-1];
	}
	
	drlo = (double *)malloc(nrow*sizeof(double));
	drup = (double *)malloc(nrow*sizeof(double));
	rstg = (int *)malloc(nrow*sizeof(int));
	for (ii=0;ii<nrow;ii++)
	{
		drlo[irow[ii]-1] = cdrlo[ii];
		drup[irow[ii]-1] = cdrup[ii];
		rstg[irow[ii]-1] = crstg[ii];
	}
	
	dclo = (double *)malloc(ncol*sizeof(double));
	dcup = (double *)malloc(ncol*sizeof(double));
	dobj = (double *)malloc(ncol*sizeof(double));
	cstg = (int *)malloc(ncol*sizeof(int));
	
	for (ii=0;ii<ncol;ii++)
	{
		dclo[icol[ii]-1] = cdclo[ii];
		dcup[icol[ii]-1] = cdcup[ii];
		dobj[icol[ii]-1] = cdobj[ii];
		cstg[icol[ii]-1] = ccstg[ii];
	}
	
	// set core problem
	OsiClpSolverInterface *ocsi = new OsiClpSolverInterface();
	ocsi->loadProblem(CoinPackedMatrix( 1,mrow,mcol,dels,nels),dclo,dcup,dobj,drlo,drup);
	
	SmiCoreIndex ic = stoch->setCore(ocsi,2,cstg,rstg);

    /* ...initialization -- indx points to first sample of each rv */
    indx = (int *) malloc( (1+nindp)*sizeof(int) );
    memset( indx,0,(1+nindp)*sizeof(int));
	
	CoinPackedVector cpv_rlo;
	CoinPackedVector cpv_rup;

    for (jj=0;jj<nindp;jj++) {
		indx[jj+1] += indx[jj] + nsamp[jj];
		dp *= dprobs[ indx[jj] ];
		
		cpv_rlo.insert(irow[n_first_stg_rows + jj]-1,demand[ indx[jj] ]);
		cpv_rup.insert(irow[n_first_stg_rows + jj]-1,demand[ indx[jj] ]);
		
    }


	// add first scenario
	int	is = smiModel->genScenarioReplaceCoreValues(ic,cpm_mat,&cpv_rlo,&cpv_rup,
													NULL,NULL,NULL,2,0,dp);	

	

    /***** ...main loop
	The logic of this loop is as follows:
	While the sample size nsamp[jj] divides the scenario index ii,
	reverse the increment direction incr[jj],
	and increase the random variable index jj by 1.
	Increment the jj'th random variable by incr[jj]
	and call ekks_addScenario with the new sample.
    ***** */

	
	
    /* ...set sample space increment for each rv to 1 */
    incr = (int *) malloc( nindp*sizeof(int) );
    for (jj=0;jj<nindp;jj++) incr[jj] = 1;
	
    for (ii=1;ii<ns;ii++) {
		iii=ii; jj=0;
		while ( !(iii%nsamp[jj]) ) {
			iii /= nsamp[jj];
			incr[jj] = -incr[jj];
			jj++;
		}
		dd = dprobs[ indx[jj] ];
		indx[jj] += incr[jj];
		dp = dp* dprobs[ indx[jj] ]/dd;
		
		drlo[ irow[n_first_stg_rows + jj] - 1 ] = demand[ indx[jj] ];
		drup[ irow[n_first_stg_rows + jj] - 1 ] = demand[ indx[jj] ];
		
	 	cpv_rlo.setElement(cpv_rlo.findIndex(irow[n_first_stg_rows + jj]-1),demand[ indx[jj] ]);
		cpv_rup.setElement(cpv_rup.findIndex(irow[n_first_stg_rows + jj]-1),demand[ indx[jj] ]);

		is = smiModel->genScenarioReplaceCoreValues(ic,cpm_mat,&cpv_rlo,&cpv_rup,
													NULL,NULL,NULL,2,0,dp);	

    }
	
	assert (is==ns);
	
	stoch->load();
	stoch->initialSolve();
	
	assert(fabs(stoch->getOsiSolverInterface().getObjValue()-1566.042)<0.01);

	return 0;

#ifdef MSH_MALLOC_DEBUG
MshMemStats(1,0);
#endif


}
#endif

int main()
{
	
    /* Model dimensions */
    int ncol=27, nrow=9, nels=44;
	
    /* Sparse matrix data...organized by row */
    int *mrow,cmrow[]={ 1, 1, 1, 1, 1,
		2, 2, 2, 2,
		3, 3, 3,
		4, 4, 4, 4, 4,
		5, 5, 5, 5,
		6, 6, 6, 6, 6, 6,
		7, 7, 7, 7, 7,
		8, 8, 8, 8, 8, 8,
		9, 9, 9, 9, 9, 9 };
    int *mcol,cmcol[]={ 1, 2, 3, 4, 5,
		6, 7, 8, 9,
		10, 11, 12,
		13, 14, 15, 16, 17,
		1,        13, 18, 19,
		2, 6, 10, 14, 20, 21,
		3, 7,     15, 22, 23,
		4, 8, 11, 16, 24, 25,
		5, 9, 12, 17, 26, 27 };
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
    int n_first_stg_cols=17;
    int *rstg,crstg[]={ 1,1,1,1,2,2,2,2,2 };
    int *cstg,ccstg[]={ 1,1,1,1,1,1,1,1,1,
		1,1,1,1,1,1,1,1,2,
		2,2,2,2,2,2,2,2,2 };
	
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
	//   int irow[]={ 9,8,7,6,5,4,3,2,1 };
	int irow[]={ 1,2,7,8,9,3,4,5,6};
    int icol[]={ 9,2,3,4,5,6,7,8,1,
		19,21,23,25,27,26,24,22,20,
		10,11,12,13,14,15,16,17,18 };

    /* local variables */
    int ns=1,ii,iii,jj,*indx,*incr;
    double dp=1.0;

    for (ii=0;ii<nindp;ii++) ns *= nsamp[ii];     /* Compute number of scenarios */

	
	// initialize SmiModel
	SmiModel *smiModel = new SmiModel(ns);
	OsiClpSolverInterface *osi = new OsiClpSolverInterface();
	smiModel->setOsiSolverHandle(osi);
	
	/* scramble LP entries */
	mrow = (int*)malloc(nels*sizeof(int));
	mcol = (int*)malloc(nels*sizeof(int));
	for (ii=0;ii<nels;ii++)
	{
		mcol[ii] = icol[cmcol[ii]-1];
		mrow[ii] = irow[cmrow[ii]-1];
	}
	
	drlo = (double *)malloc(nrow*sizeof(double));
	drup = (double *)malloc(nrow*sizeof(double));
	rstg = (int *)malloc(nrow*sizeof(int));
	for (ii=0;ii<nrow;ii++)
	{
		drlo[irow[ii]-1] = cdrlo[ii];
		drup[irow[ii]-1] = cdrup[ii];
		rstg[irow[ii]-1] = crstg[ii];
	}
	
	dclo = (double *)malloc(ncol*sizeof(double));
	dcup = (double *)malloc(ncol*sizeof(double));
	dobj = (double *)malloc(ncol*sizeof(double));
	cstg = (int *)malloc(ncol*sizeof(int));
	for (ii=0;ii<ncol;ii++)
	{
		dclo[icol[ii]-1] = cdclo[ii];
		dcup[icol[ii]-1] = cdcup[ii];
		dobj[icol[ii]-1] = cdobj[ii];
		cstg[icol[ii]-1] = ccstg[ii];
	}

	// this to test the matrix update stanza in genScenario
	int corenels = nels - 4;
	
	// set core model using Osi interface
	OsiClpSolverInterface *ocsi = new OsiClpSolverInterface();
	ocsi->loadProblem(CoinPackedMatrix( 1,mrow,mcol,dels,corenels),dclo,dcup,dobj,drlo,drup);
	
	SmiCoreIndex ic = smiModel->setCore(ocsi,2,cstg,rstg);
	
	
	// Coin structures for scenario updates to right hand sides
	CoinPackedVector cpv_rlo;
	CoinPackedVector cpv_rup;

	// Coin structure for scenario "updates" to core matrix
	CoinPackedMatrix *cpm_mat = new CoinPackedMatrix(1,mrow+corenels,mcol+corenels,dels+corenels,nels-corenels);
		
    // initialize right hand side data for first scenario
    indx = (int *) malloc( (1+nindp)*sizeof(int) );
    memset( indx,0,(1+nindp)*sizeof(int));
    for (jj=0;jj<nindp;jj++) {
		indx[jj+1] += indx[jj] + nsamp[jj];
		dp *= dprobs[ indx[jj] ];
		
		cpv_rlo.insert(irow[n_first_stg_rows + jj]-1,demand[ indx[jj] ]);
		cpv_rup.insert(irow[n_first_stg_rows + jj]-1,demand[ indx[jj] ]);
    }
	
	// first scenario
	int	is = smiModel->genScenarioReplaceCoreValues(ic,cpm_mat,&cpv_rlo,&cpv_rup,
		NULL,NULL,NULL,2,0,dp);	
	
	
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
	
    for (ii=1;ii<ns;ii++) {
		iii=ii; jj=0;
		while ( !(iii%nsamp[jj]) ) {
			iii /= nsamp[jj];
			incr[jj] = -incr[jj];
			jj++;
		}
		dp /= dprobs[ indx[jj] ];
		indx[jj] += incr[jj];
		dp *= dprobs[ indx[jj] ];
		
		// set data
		cpv_rlo.setElement(cpv_rlo.findIndex(irow[n_first_stg_rows + jj]-1),demand[ indx[jj] ]);
		cpv_rup.setElement(cpv_rup.findIndex(irow[n_first_stg_rows + jj]-1),demand[ indx[jj] ]);
		
		// genScenario
		is = smiModel->genScenarioReplaceCoreValues(ic,cpm_mat,&cpv_rlo,&cpv_rup,
			NULL,NULL,NULL,2,0,dp);	
    }

	// solve with decomp solver
	smiModel->decompSolve();

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
