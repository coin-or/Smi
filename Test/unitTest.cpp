// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iostream>

#include "SmiScenarioTree.hpp"
#include "SmiScnModel.hpp"
#include "OsiClpSolverInterface.hpp"

//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif

void
SmiScenarioTreeUnitTest()
{

	int  *i1 = new int(1);
	int  *i2= new int(12);
	int  *i3= new int(13);
	int  *i4= new int(24);
	int  *i5= new int(25);
	int  *i6= new int(36);
	int  *i7= new int(37);
	int  *i8= new int(38);

	vector<int *> path4(3);
	path4[0]= i1; 
	path4[1]= i2; 
	path4[2]= i4; 

	vector<int *> path5(1, i5);
	
	vector<int *> path6(2);
	path6[0]= i3; 
	path6[1]= i6; 

	vector<int *> path7(1, i7);
	
	vector<int *> path8(1, i8);
	
	
	SmiScenarioTree<int*> s;

	int is1 = s.addPathtoLeaf(0,0,path4);
	int is2 = s.addPathtoLeaf(is1,1,path5);
	int is3 = s.addPathtoLeaf(is2,0,path6);
	int is4 = s.addPathtoLeaf(is3,1,path7);
	int is5 = s.addPathtoLeaf(is3,1,path8);

	int  *ii1 = new int(1);
	int  *ii2= new int(12);
	int  *ii3= new int(13);
	int  *ii4= new int(24);
	int  *ii5= new int(25);
	int  *ii6= new int(36);
	int  *ii7= new int(37);
	int  *ii8= new int(38);

	
	vector<int *> ppath4(3);
	vector<int> label4(3);
	ppath4[0]= ii1; label4[0]=*ii1;
	ppath4[1]= ii2; label4[1]=*ii2;
	ppath4[2]= ii4; label4[2]=*ii4;

	vector<int *> ppath5(3);
	ppath5[0] = ppath4[0];
	ppath5[1] = ppath4[1];
	ppath5[2] = ii5;
	vector<int> label5(3);
	label5[0] = label4[0];
	label5[1] = label4[1];
	label5[2] = *ii5;

	vector<int *> ppath6(3);
	ppath6[0] = ppath5[4];
	vector<int> label6(3);
	label6[0] = label5[0];
	ppath6[1]= ii3; label6[1] = *ii3;
	ppath6[2]= ii6; label6[2] = *ii6;

	vector<int *> ppath7(3);
	ppath7[0] = ppath6[0];
	ppath7[1] = ppath6[1];
	ppath7[2] = ii7;

	vector<int> label7(3);
	label7[0] = label6[0];
	label7[1] = label6[1];
	label7[2] = *ii7;

	vector<int *> ppath8(3);
	ppath8[0] = ppath6[0];
	ppath8[1] = ppath6[1];
	ppath8[2] = ii8;
	vector<int> label8(3);
	label8[0] = label6[0];
	label8[1] = label6[1];
	label8[2] = *ii8;

	// tree s1 is identical to s just uses labels to id nodes 
	SmiScenarioTree<int*> s1;
	int is11 = s1.addPathtoLeaf(label4,ppath4);
	int is12 = s1.addPathtoLeaf(label5,ppath5);
	int is13 = s1.addPathtoLeaf(label6,ppath6);
	int is14 = s1.addPathtoLeaf(label7,ppath7);
	int is15 = s1.addPathtoLeaf(label8,ppath8);


	vector<int *> vec1 = s.getScenario(is1);
	assert( vec1[0] == i1 );
	assert( vec1[1] == i2 );
	assert( vec1[2] == i4 );

	vector<int *> vec11 = s1.getScenario(is11);
	assert( vec11[0] == ii1 );
	assert( vec11[1] == ii2 );
	assert( vec11[2] == ii4 );

	vector<int *>::iterator vbeg = s.scenBegin(is1);
	vector<int *>::iterator vend = s.scenEnd(is1);
	assert( *vbeg++ == i1 );
	assert( *vbeg++ == i2 );
	assert( *vbeg++ == i4 );
	assert( vbeg == vend );

	vector<int *>::iterator vbeg1 = s1.scenBegin(is11);
	vector<int *>::iterator vend1 = s1.scenEnd(is11);
	assert( *vbeg1++ == ii1 );
	assert( *vbeg1++ == ii2 );
	assert( *vbeg1++ == ii4 );
	assert( vbeg1 == vend1 );

	
	vector<int *> vec2 = s.getScenario(is2);
	assert( vec2[0] == i1 );
	assert( vec2[1] == i2 );
	assert( vec2[2] == i5 );

		
	vector<int *> vec12 = s1.getScenario(is12);
	assert( vec12[0] == ii1 );
	assert( vec12[1] == ii2 );
	assert( vec12[2] == ii5 );

	vector<int *> vec3 = s.getScenario(is3);
	assert( vec3[0] == i1 );
	assert( vec3[1] == i3 );
	assert( vec3[2] == i6 );

	vector<int *> vec13 = s1.getScenario(is13);
	assert( vec13[0] == ii1 );
	assert( vec13[1] == ii3 );
	assert( vec13[2] == ii6 );

	vector<int *> vec4 = s.getScenario(is4);
	assert( vec4[0] == i1 );
	assert( vec4[1] == i3 );
	assert( vec4[2] == i7 );

	
	vector<int *> vec14 = s1.getScenario(is14);
	assert( vec14[0] == ii1 );
	assert( vec14[1] == ii3 );
	assert( vec14[2] == ii7 );

	vector<int *> vec5 = s.getScenario(is5);
	assert( vec5[0] == i1 );
	assert( vec5[1] == i3 );
	assert( vec5[2] == i8 );

	
	vector<int *> vec15 = s1.getScenario(is15);
	assert( vec15[0] == ii1 );
	assert( vec15[1] == ii3 );
	assert( vec15[2] == ii8 );

	vector<int*>::iterator i;
	i=s.treeBegin();
	assert(*i==i1);
	i++;
	assert(*i==i2);
	i = s.treeEnd();
	i--;
	assert(*i==i8);

	vector<int*>::iterator ii;
	ii=s1.treeBegin();
	assert(*ii==ii1);
	ii++;
	assert(*ii==ii2);
	ii = s1.treeEnd();
	ii--;
	assert(*ii==ii8);

	vbeg = s.scenBegin(is1);
	vend = s.scenEnd(is1);
	assert( *vbeg++ == i1 );
	assert( *vbeg++ == i2 );
	assert( *vbeg++ == i4 );
	assert( vbeg == vend );

	vbeg1 = s1.scenBegin(is11);
	vend1 = s1.scenEnd(is11);
	assert( *vbeg1++ == ii1 );
	assert( *vbeg1++ == ii2 );
	assert( *vbeg1++ == ii4 );
	assert( vbeg1 == vend1 );

}

void
SmiTreeNodeUnitTest()
{

	int  *i1 = new int(1);
	int  *i2= new int(12);
	int  *i3= new int(13);
	int  *i4= new int(24);
	int  *i5= new int(25);
	int  *i6= new int(36);
	int  *i7= new int(37);
	int  *i8= new int(38);

	SmiTreeNode<int *> *n1 = new SmiTreeNode<int *>(i1);

	SmiTreeNode<int *> *n2 = n1->addChild(i2);
	SmiTreeNode<int *> *n3 = n1->addChild(i3);
	SmiTreeNode<int *> *n4 = n2->addChild(i4);
	SmiTreeNode<int *> *n5 = n2->addChild(i5);
	SmiTreeNode<int *> *n6 = n3->addChild(i6);
	SmiTreeNode<int *> *n7 = n3->addChild(i7);
	SmiTreeNode<int *> *n8 = n3->addChild(i8);

	assert( n1->depth() == 0 );
	assert( n2->depth() == 1 );
	assert( n3->depth() == 1 );
	assert( n4->depth() == 2 );
	assert( n5->depth() == 2 );
	assert( n6->depth() == 2 );
	assert( n7->depth() == 2 );
	assert( n8->depth() == 2 );

	// parents point to last children
	assert( n1->getChild() == n3 );
	assert( n3->getChild() == n8 );

	// siblings
	assert( n3->hasSibling() );
	assert( n3->getSibling() == n2 );
	assert( !n2->hasSibling() );

	// sibling and cousin pointers are a linked list
	// for same level nodes
	assert( n8->hasSibling() );
	assert( n8->getSibling() == n7 );
	assert( n7->hasSibling() );
	assert( n7->getSibling() == n6 );
	assert( !n6->hasSibling());

	assert( n6->getParent()->getSibling()  == n2 );
	assert( n2->getChild() == n5 );
	assert( n5->hasSibling() );
	assert( n5->getSibling() == n4 );
	// last element of same level list
	assert( !n4->hasSibling());

	vector<SmiTreeNode<int *> *> vec1 = n1->getChildren();
	assert (vec1[0] == n2 );
	assert (vec1[1] == n3 );

}



 
// Display message on stdout and stderr
void testingMessage( const char * const msg )
{
//  std::cerr <<msg;
  cout <<endl <<"*****************************************"
       <<endl <<msg <<endl;
}


void SmiScnSmpsIOUnitTest()
{
	{
		// test SMPS files app0110R
		SmiScnModel smi;	
		smi.readSmps("../../Mps/Stochastic/app0110R");		
		OsiClpSolverInterface *clp = new OsiClpSolverInterface();
		smi.setOsiSolverHandle(*clp);	
		OsiSolverInterface *osiStoch = smi.loadOsiSolverData();
		
		int nrows = osiStoch->getNumRows();
		int ncols = osiStoch->getNumCols();
		assert(nrows==129);
		assert(ncols==268);
		
		osiStoch->initialSolve();		
		assert(fabs(osiStoch->getObjValue()-44.66666) < 0.0001);
		printf(" *** Successfully tested SMPS interfaces on app0110R.\n");
	}
	
	{
		// test SMPS files from Watson test suite (Cambridge, UK)
		
		SmiScnModel smi;
		smi.readSmps("../../Mps/Stochastic/wat_10_C_32");
		OsiClpSolverInterface *clp = new OsiClpSolverInterface();
		smi.setOsiSolverHandle(*clp);	
		OsiSolverInterface *osiStoch = smi.loadOsiSolverData();
		osiStoch->setHintParam(OsiDoPresolveInInitial,true);
		osiStoch->setHintParam(OsiDoScale,true);
		osiStoch->setHintParam(OsiDoCrash,true);
		osiStoch->initialSolve();		
		assert(fabs(osiStoch->getObjValue()+2622.062) < 0.01);
		printf(" *** Successfully tested SMPS interfaces on wat_10_32_C.\n");
	}
}	
void SmiScnModelUnitTest()
{
	// exhaustive test of direct interfaces
	
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

#define INF 1.0e31

    double *dcup,cdcup[]={ INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,
		INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,
		INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF };
	
    /* Row bounds */
    double *drlo,cdrlo[]={ -INF, -INF, -INF, -INF,  0.0, 4.0, 0.0, 8.0, 10.0 };
    double *drup,cdrup[]={ 10.0, 19.0, 25.0, 15.0,  0.0, 7.0, 0.0, 8.0, 90.0 };
	
    /* Stages */
    //int nstg=2;
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
    double dp=1.0;

    for (ii=0;ii<nindp;ii++) ns *= nsamp[ii];     /* Compute number of scenarios */

	// debug small sample
	// ns = 3;
	
	// initialize SmiModel
	SmiScnModel *smiModel = new SmiScnModel();


	OsiClpSolverInterface *osiClp1 = new OsiClpSolverInterface();

	smiModel->setOsiSolverHandle(*osiClp1);
	
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
	SmiCoreData *osiCore = new SmiCoreData(ohoh,3,cstg,rstg);

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
  int t;
	for ( t=0;t<3;t++)
	{
		CoinPackedVector *cpvdrlo = osiCore->getRowLower(t);
		CoinPackedVector *cpvdrup = osiCore->getRowUpper(t);
		CoinPackedVector *cpvdclo = osiCore->getColLower(t);
		CoinPackedVector *cpvdcup = osiCore->getColUpper(t);
		CoinPackedVector *cpvdobj = osiCore->getObjCoefficients(t);

		double *core_drlo = new double[osiCore->getNumRows(t)];
		osiCore->getNode(t)->copyRowLower(core_drlo);

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
			assert(core_drlo[ii-osiCore->getRowStart(t)] ==drlo[osiCore->getRowExternalIndex(ii)]);

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
	int anc = 0;
	int branch = 1;
	int	is = smiModel->genScenarioReplaceCoreValues(ic,cpm_mat,NULL,NULL,NULL,
									&cpv_rlo,&cpv_rup,branch,anc,dp);	
	

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
	for ( t=0;t<3;t++)
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

	printf(" *** Successfully tested problem with scenario %d.\n",is);

	
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
			&cpv_rlo,&cpv_rup,branch,anc,dp);	
		
		
		if (is < 3)
		{		
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
			int t;
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
			
			printf(" *** Successfully tested problem with scenario %d.\n",is);
		}	
	}
	
	// solve with decomp solver
	
	// load problem data into OsiSolver
	smiModel->loadOsiSolverData();
	// get Osi pointer
	OsiSolverInterface *smiOsi = smiModel->getOsiSolverInterface();
	// set some parameters
	smiOsi->setHintParam(OsiDoPresolveInInitial,true);
	smiOsi->setHintParam(OsiDoScale,true);
	smiOsi->setHintParam(OsiDoCrash,true);
	// solve using Osi Solver
	smiOsi->initialSolve();
	// test optimal value
    assert(fabs(smiOsi->getObjValue()-1566.042)<0.01);

	printf(" *** Successfully tested direct interfaces.\n");
	// test solutions
	const double *dsoln = smiOsi->getColSolution();
	double objSum = 0.0;

	/* The canonical way to traverse the tree:
	   For each scenario, get the leaf node.
	   Then get the parent.  Repeat until parent is NULL.
	   (Only the root node has a NULL parent.)
	 */
	
	for(is=0; is<ns; ++is)
	{
		/* this loop calculates the scenario objective value */
		double scenSum = 0.0;

		// start with leaf node
		SmiScnNode *node = smiModel->getLeafNode(is);

		// leaf node probability is the scenario probability
		double scenprob = node->getModelProb();
	
		while (node != NULL)
		{

			// getColStart returns the starting index of node in OSI model
			for(int j=node->getColStart(); j<node->getColStart()+node->getNumCols(); ++j)
			{
				// getCoreColIndex returns the corresponding Core index
				// in the original (user's) ordering
				scenSum += dobj[node->getCoreColIndex(j)]*dsoln[j];				
			}			
			// get parent of node
			node = node->getParent();
		}
		objSum += scenSum*scenprob;
	}

	assert(fabs(smiOsi->getObjValue()-objSum) < 0.01);

	printf(" *** Completed all tests for SmiScnModel.\n");

}

int main()
{
  

	testingMessage( "Testing SmiTreeNode \n");
	SmiTreeNodeUnitTest();

	testingMessage( "Testing SmiScenarioTree\n" );
	SmiScenarioTreeUnitTest();

	testingMessage( "Testing SmiScnSmpsIO\n" );
	SmiScnSmpsIOUnitTest();

	testingMessage( "Testing SmiScnModel\n" );
	SmiScnModelUnitTest();


	testingMessage( "*** Done! *** \n");
	

  return 0;
}
