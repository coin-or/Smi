// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <string>
using namespace std;

#include <cassert>
#include <iostream>
#include <vector>
#include <map>

#include "CoinHelperFunctions.hpp"

#include "SmiScnModel.hpp"
#include "SmiScnData.hpp"
#include "OsiClpSolverInterface.hpp"



class SmiLShapedModel : public SmiScnModel
{
private:
	SmiLShapedModel *parent;
	typedef vector<SmiLShapedModel *> Children;
	Children children;
	vector<SmiScnNode *> vecNaturalNodes;
	vector<SmiScnNode *> vecVirtualNodes;
	vector<int> vecVirtualNodeOffset;
	vector<int> vecParentNodeOffset;
	int numVirtualCols_;
	int *virt2nat;
public:
	SmiLShapedModel(): parent(NULL),SmiScnModel(){}

	void setParent(SmiScnModel *s) {parent = (SmiLShapedModel*)(s);}
	SmiLShapedModel *getParent() { return parent;}
	void addChild(SmiLShapedModel *s) {children.push_back(s);}
	Children * getChildren() {return &children;}
	void addNaturalNode(SmiScnNode *node) { this->vecNaturalNodes.push_back(node);}
	void addVirtualNode(SmiScnNode *node) { this->vecVirtualNodes.push_back(node); numVirtualCols_+=node->getNumCols();}
	int getNumberNaturalNodes(){ return (int) vecNaturalNodes.size(); }
	int getNumberVirtualNodes(){ return (int) vecVirtualNodes.size(); }

	OsiSolverInterface *loadOsiSolverData()
	{
		this->setupOsiSolverData();
		if (parent!=NULL)
		{
			virt2nat = (int *) malloc (sizeof(int)*this->numVirtualCols_);
			for (int i=0; i<(int)vecVirtualNodes.size(); i++)
			{
				SmiScnNode *node = vecVirtualNodes[i];
	
				vecVirtualNodeOffset.push_back(this->ncol_);
				vecParentNodeOffset.push_back(node->getColStart());
				this->ncol_+=node->getNumCols();

				CoinIota(virt2nat+vecVirtualNodeOffset[i],
					virt2nat+this->ncol_,node->getColStart());

			}
		}
		this->gutsofloadOsiSolverData();
		return this->osiStoch_;
	}


};

int main()
{

		const char* name="../../Data/Stochastic/wat_10_C_32";		
		SmiScnModel smi;	

		// read SMPS model from files
		//	<name>.core, <name>.time, and <name>.stoch
		smi.readSmps(name);		

		// set solver object for SmiScnModel to be Clp
		smi.setOsiSolverHandle(new OsiClpSolverInterface());

/* L-Shaped Method */

		// tree of LShaped Models
		vector<SmiLShapedModel *> vecLSModel;
		// generic LSModel pointer
		SmiLShapedModel *smiLSModel;


		// assign nodes to LShaped Models
		// here we just assign based on cut at first time stage.

		// temp vector to hold nodes
		vector<SmiScnNode *> vecNodes;

		// Add the first scenario to the root model
		smiLSModel = new SmiLShapedModel();
		vecLSModel.push_back(smiLSModel);

		// start from leaf and work back to root
		SmiScnNode *node = smi.getLeafNode(0);
		while (node != NULL)
		{
				vecNodes.push_back(node);
				node=node->getParent();
		}
		// add nodes in reverse order		
		vector<SmiScnNode *>::reverse_iterator vFirst(vecNodes.rbegin());
		vector<SmiScnNode *>::reverse_iterator vLast(vecNodes.rend());

		// add nodes to model, from leaf back towards root
		while (vFirst != vLast)
		{
			node=*vFirst;
			smiLSModel->addNaturalNode(node);
			node->setModel(smiLSModel);
			vFirst++;
		}


		// rest of scenarios
		for (SmiScenarioIndex i=1; i<smi.getNumScenarios(); ++i)
		{
			SmiScnNode *parent=NULL;
			SmiScnNode *node = smi.getLeafNode(i);
			vecNodes.clear();
			while (node->getModel() == NULL) //form vector of unallocated nodes
			{
				vecNodes.push_back(node);
				parent=node->getParent();
				if (parent==smi.getRootNode()) 
					break;
				else
					node=parent;
			}
			if (parent==smi.getRootNode()) //last node added is child of root node
				                           //so create a new submodel
			{
				// create new model
				smiLSModel = new SmiLShapedModel();
				vecLSModel.push_back(smiLSModel);
				// set parent
				smiLSModel->setParent(parent->getModel());
				// add self to children of parent
				SmiLShapedModel *parentModel = (SmiLShapedModel *)parent->getModel();
				assert(parentModel != NULL && "Parent Model bad reference");
				parentModel->addChild(smiLSModel);

			} else				//parent is in a model
			{
				smiLSModel = (SmiLShapedModel *)parent->getModel();
				assert(smiLSModel != NULL && "Parent Model bad reference");
	
			}
			
			// record branch stage from parent
			int branch_stage = node->getNode()->getStage();
						
			// check to see if we have to add virtual nodes
			if (parent->getModel() != node->getModel())
			{
				// push virtual nodes onto the vector
				node=parent;
				while (node != NULL)
				{	
					vecNodes.push_back(node);
					node = node->getParent();
				}
			}

			// add nodes to model, from leaf back towards root
			vFirst=vecNodes.rbegin();
			vLast=vecNodes.rend();
			while (vFirst != vLast)
			{
				node=*vFirst;
				if (node->getNode()->getStage() < branch_stage)
					smiLSModel->addVirtualNode(node);
				else
				{
					smiLSModel->addNaturalNode(node);
					node->setModel(smiLSModel);
				}
				vFirst++;
			}
		}

		cout << "There are "<< vecLSModel.size() << " models" << endl;

		for (int i=0; i< vecLSModel.size(); i++)
		{
			smiLSModel = vecLSModel[i];
			cout << "Model " << i << " has "
				<< smiLSModel->getNumberNaturalNodes()<< " natural nodes, and "
				<< smiLSModel->getNumberVirtualNodes() << " virtual nodes." 
				<< endl;
		}

		for (int i=0; i < vecLSModel.size(); i++)
		{
			smiLSModel = vecLSModel[i];
			smiLSModel->setOsiSolverHandle(new OsiClpSolverInterface());
			smiLSModel->loadOsiSolverData();
		}



		// load solver data
		// 	this step generates the deterministic equivalent 
		//	and returns an OsiSolver object 
		OsiSolverInterface *osiStoch = smi.loadOsiSolverData();

		// set some nice Hints to the OSI solver
		osiStoch->setHintParam(OsiDoPresolveInInitial,true);
		osiStoch->setHintParam(OsiDoScale,true);
		osiStoch->setHintParam(OsiDoCrash,true);

		// solve
		

		// print results
		printf("Solved stochastic program %s\n", name);
		printf("Number of rows: %d\n",osiStoch->getNumRows());
		printf("Number of cols: %d\n",osiStoch->getNumCols());
		printf("Optimal value: %g\n",osiStoch->getObjValue());		

}	
