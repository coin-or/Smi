// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef SmiScenarioTree_H
#define SmiScenarioTree_H

/** Scenario Tree

  This class is used for storing and accessing scenario trees.

*/
#include <list>
#include <vector>
#include <assert.h>
#include "SmiTreeNode.hpp"

using namespace std;

template<class T>
class SmiScenarioTree  {
   friend void SmiScenarioTreeUnitTest();
public:

//---------------------------------------------------------------------------
  /**@name Iterators */
  //@{


	/** begin */
	vector<T>::iterator treeBegin(){ return node_data.begin(); }
	/** end */
	vector<T>::iterator treeEnd(){ return node_data.end();}
	/** whole tree */
	vector<T> &wholeTree(){ return node_data;}

	/** scenario iterators
	TODO: native code for these iterators that does not
	depend on copying. */
	vector<T>::iterator scenBegin(int s){
		getScenario(s);
		return scen_data.begin();
	}

	vector<T>::iterator scenEnd(int s){
		getScenario(s);
		return scen_data.begin()+leaf_[s]->depth()+1;
	}




//---------------------------------------------------------------------------
  /**@name Query members */
  //@{

	/** Get root node.	*/
	SmiTreeNode<T> *getRoot() { return root_; }

	/** Get leaf node. */
	SmiTreeNode<T> *getLeaf(int scn) { return leaf_[scn];}

	/** Get node identified by scenario/stage.	*/
	SmiTreeNode<T> &find(int scenario, int stage)
	{
		assert (scenario < leaf_.size());		
		SmiTreeNode<T> * n = leaf_[scenario];
		assert (stage < n->depth() + 1);
		while (stage < n->depth())
			n = n->getParent();
		return *n;
	}

	/** Get vector of node data for given scenario */
	vector<T> &getScenario(int scenario)
	{
		assert (scenario < leaf_.size());
		SmiTreeNode<T> * n = leaf_[scenario];

//		if ( n->getDataPtr()==scen_data[n->depth()] ) return scen_data;

		int i=n->depth()+1;
		while(i>0)
		{
			scen_data[--i] = n->getDataPtr();
			n = n->getParent();
		}
		return scen_data;
	}


  //@}

//---------------------------------------------------------------------------
  /**@name Tree modification members */
  //@{
    /** Add path to leaf.
	    Responsibility for memory management of SmiTreeNodeData elements
		is assigned to SmiScenarioTree.
		SmiTreeNodeData elements must be created with "new" operator.
	*/
	int addPathtoLeaf(int brscenario, int stage, vector<T> &pathdata)
	{
		SmiTreeNode<T> *parent = NULL;		
		if (leaf_.size())
			parent = &find(brscenario,stage);
		for ( int i = 0; i < pathdata.size(); i++)
		{	
			if (parent)
				parent = parent->addChild(pathdata[i]);
			else
			{
				parent = root_ = new SmiTreeNode<T>(pathdata[0]);
			}
			node_data.push_back(pathdata[i]);
			if (i+stage < scen_data.size())
				scen_data[i+stage] = pathdata[i];
			else
				scen_data.push_back(pathdata[i+stage]);
		}
		if (pathdata.size())
		{
			leaf_.push_back(parent);
			
		}
		return leaf_.size()-1;
	}
  //@}

//--------------------------------------------------------------------------
   /**@name Constructors, destructors and major modifying methods*/
   //@{
   /// Default Constructor creates an empty scenario tree
	SmiScenarioTree<T>(): leaf_(0),root_(NULL) {};

   /// Destructor 
   virtual ~SmiScenarioTree<T>() {delete root_;}
   //@}

private:
	vector<T> node_data;
	vector<T> scen_data;
	vector<SmiTreeNode<T> *> leaf_;
	SmiTreeNode<T> *root_;
};

//#############################################################################
/** A function that tests the methods in the SmiScenarioTree class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void
SmiScenarioTreeUnitTest();

#endif //SmiScenarioTree_H
