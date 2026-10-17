// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "SmiSmpsIO.hpp"
#include <cassert>
#include <cmath>
#include <cfloat>
#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "CoinMpsIO.hpp"
#include "CoinMessage.hpp"
#include "CoinError.hpp"


#if 1
const static char *section[] = {
  "", "NAME", "ENDATA", " ", "PERIODS", "SCENARIOS", "INDEPENDENT"," "
};

const static char *smpsType[] = {
	"SC","BL","  ","DISCRETE","ADD","REPLACE",""
};
#endif


//#############################################################################

SmiCoreData *
SmiSmpsIO::readTimeFile(SmiScnModel *smi, const char *c, const char *ext)
{

        CoinFileInput *input = 0;
        int returnCode = dealWithFileName(c,ext,input);
        if (returnCode<0) {
          return NULL;
        } else if (returnCode>0) {
          // delete cardReader_;
          smpsCardReader_ = new SmiSmpsCardReader ( input, this);
        }
	
	smpsCardReader_->readToNextSection();
	if ( smpsCardReader_->whichSection (  ) == COIN_NAME_SECTION ) {
		iftime = true;
		// check problem name and issue warning if needed
		if (strcmp(problemName_,smpsCardReader_->columnName()))
		{
			printf("Warning: Time file name %s does not match problem file name %s\n",
				smpsCardReader_->columnName(),problemName_);
		}
	} else if ( smpsCardReader_->whichSection (  ) == COIN_UNKNOWN_SECTION ) {
		handler_->message(COIN_MPS_BADFILE1,messages_)<<smpsCardReader_->card()
			<<1
			<<fileName_
			<<CoinMessageEol;
/*
#ifdef COIN_HAS_ZLIB
		if (!smpsCardReader_->filePointer()) 
			handler_->message(COIN_MPS_BADFILE2,messages_)<<CoinMessageEol;
		
#endif
*/
		return NULL;
	} else if ( smpsCardReader_->whichSection (  ) != COIN_EOF_SECTION ) {
		// save name of section
		free(problemName_);
		problemName_=strdup(smpsCardReader_->card());
	} else {
		handler_->message(COIN_MPS_EOF,messages_)<<fileName_
			<<CoinMessageEol;
		return NULL;
	}

	if (iftime)
	{
		int rowEnd, rowStart = 0;
		int colEnd, colStart = 0;
		if (	smpsCardReader_->nextSmpsField() != SMI_TIME_SECTION ) // PERIODS card
			return NULL;
	
		if (	smpsCardReader_->nextSmpsField() == SMI_TIME_SECTION ) // first period
		{
		
			cstag_ = new int[this->getNumCols()];
			rstag_ = new int[this->getNumRows()];
			std::string SunStudioNeedsThis = smpsCardReader_->periodName() ;
			periodMap_.insert(make_pair(SunStudioNeedsThis,nstag_));
		}
		else
			return NULL;

		// don't handle the unordered case yet
		if (smpsCardReader_->whichSmpsType() != SMI_TIME_ORDERED_CORE_TYPE)
			return NULL;

		while( smpsCardReader_->nextSmpsField (  ) == SMI_TIME_SECTION ) 
		{
			if (smpsCardReader_->whichSmpsType() == SMI_TIME_ORDERED_CORE_TYPE)
			{
				
				colEnd = this->columnIndex(smpsCardReader_->columnName());
				rowEnd = this->rowIndex(smpsCardReader_->rowName());
				
				std::string SunStudioNeedsThis = smpsCardReader_->periodName() ;
				periodMap_.insert(make_pair(SunStudioNeedsThis,nstag_+1));
				
				int i;
				for(i=colStart;i<colEnd;++i)
					cstag_[i]=nstag_;
				for(i=rowStart;i<rowEnd;++i)
					rstag_[i]=nstag_;
				++nstag_;
				colStart = colEnd;
				rowStart = rowEnd;
			}
			else
				return NULL;
		}
		
		if ( smpsCardReader_->whichSmpsSection() == SMI_ENDATA_SECTION
			&& smpsCardReader_->whichSmpsType() == SMI_TIME_ORDERED_CORE_TYPE)
		{
			int i;

			for(i=std::max(colStart, 0);i<this->getNumCols();++i)
				cstag_[i]=nstag_;
			for(i=std::max(rowStart, 0);i<this->getNumRows();++i)
				rstag_[i]=nstag_;
		}
		else
			return NULL;

	}

    // Aenderung: IntegerType
    int intLen=0;
    int binLen=0;
    const double* colLower = this->getColLower();
    const double* colUpper = this->getColUpper();
    for(int i=0;i<this->getNumCols();i++ ){
        if(this->isInteger(i)){
            intLen++;
            if(colUpper[i]==1&&colLower[i]==0){
                binLen++;
            }
        }
    }
    int* intInd = new int[intLen];
    int* binInd = new int[binLen];
    int cint = 0;
    int cbin = 0;
    for (int i=0; i<this->getNumCols();i++ ){
        if(this->isInteger(i)){
            intInd[cint]=i;
            cint++;
            if(colUpper[i]==1&&colLower[i]==0){
                binInd[cbin]=i;
                cbin++;
            }
        }
    }
    // Aenderung: IntegerType


	SmiCoreData *smiCoreData = new SmiCoreData(this,nstag_+1,cstag_,rstag_,intInd,intLen,binInd,binLen);
    delete[] intInd;
    delete[] binInd;
	return smiCoreData;

}

//#############################################################################

int
SmiSmpsIO::readStochFile(SmiScnModel *smi,SmiCoreData *core, const char *c, const char *ext)
{
	
        CoinFileInput *input = 0;
        int returnCode = dealWithFileName(c,ext,input);
        if (returnCode<0) {
          return -1;
        } else if (returnCode>0) {
          delete smpsCardReader_;
          smpsCardReader_ = new SmiSmpsCardReader ( input, this);
		  if(combineRuleSet)
			  smpsCardReader_->setCoreCombineRule(combineRule_);
        }
	

	smpsCardReader_->readToNextSection();
	if ( smpsCardReader_->whichSection (  ) == COIN_NAME_SECTION ) {
		ifstoch = true;
		// check problem name and issue warning if needed
		if (strcmp(problemName_,smpsCardReader_->columnName()))
		{
			printf("Warning: Time file name %s does not match problem file name %s\n",
				smpsCardReader_->columnName(),problemName_);
		}
	} else if ( smpsCardReader_->whichSection (  ) == COIN_UNKNOWN_SECTION ) {
		handler_->message(COIN_MPS_BADFILE1,messages_)<<smpsCardReader_->card()
			<<1
			<<fileName_
			<<CoinMessageEol;
#if 0 //FIXME: SmiSmpsCardReader has no member filePointer(), also none in CoinMpsCardReader
    if (CoinFileInput::haveGzipSupport() && !smpsCardReader_->filePointer())
		  handler_->message(COIN_MPS_BADFILE2,messages_)<<CoinMessageEol;
#endif
		return -2;
	} else if ( smpsCardReader_->whichSection (  ) != COIN_EOF_SECTION ) {
		// save name of section
		free(problemName_);
		problemName_=strdup(smpsCardReader_->card());
	} else {
		handler_->message(COIN_MPS_EOF,messages_)<<fileName_
			<<CoinMessageEol;
		return -3;
	}

	if (ifstoch)
	{
		switch( smpsCardReader_->nextSmpsField() )
		{
		case SMI_INDEPENDENT_SECTION: // INDEPENDENT card
		{
			printf("Processing INDEPENDENT section\n");

			// Create discrete distribution
			SmiDiscreteDistribution *smiDD = new SmiDiscreteDistribution(core);
			SmiDiscreteRV *smiRV = NULL;

			CoinPackedVector drlo,drup,dclo,dcup,dobj;
			CoinPackedMatrix *matrix = new CoinPackedMatrix(false,0.25,0.25);
			assert(!matrix->isColOrdered());

			int nrow = core->getNumRows();
			int ncol = core->getNumCols();
			matrix->setDimensions(nrow,ncol);

			int oldi=-1;
			int oldj=-1;
					
			while( smpsCardReader_->nextSmpsField (  ) == SMI_INDEPENDENT_SECTION ) 
			{
				
				switch (smpsCardReader_->whichSmpsType())
				{
				case SMI_BL_CARD:
				{
					//this is a block 
					// not handled yet
					return -2;
				}
				case SMI_COLUMN_CARD:
				{
					int j=this->columnIndex(smpsCardReader_->columnName());
					int i=this->rowIndex(smpsCardReader_->rowName());

					double value = smpsCardReader_->value();

					
					if (j<0) // check RHS
					{
						if (!strcmp(this->getRhsName(),"")) {
							 free(rhsName_);
		 					 rhsName_=strdup(smpsCardReader_->columnName());
						} else
							assert(!strcmp(smpsCardReader_->columnName(),this->getRhsName()));
						assert(!(i<0));

						char c=this->getRowSense()[i];

						switch(c)
						{
						case 'E':
							drlo.insert(i,value);
							drup.insert(i,value);
							break;
						case 'L':
							drlo.insert(i,value);
							break;
						case 'G':
							drup.insert(i,value);
							break;
						default:
							assert(!"bad row sense: shouldn't get here");
							break;
						}
					}
					else if(i<0 || i==getNumRows()) // check OBJ
					{
						assert(!strcmp(smpsCardReader_->rowName(),this->getObjectiveName()));
						dobj.insert(j,value);
					}
					else	// add element
					{
						matrix->modifyCoefficient(i,j,value,true);
					}
					if ( !( (oldi==i)&&(oldj==j) ) ) // this is a new RV
					{
						oldi=i;
						oldj=j;
						// store the old one
						if (smiRV != NULL) 
						{
							smiDD->addDiscreteRV(smiRV);
						}
						// create a new one -- row stage dominates unless is column bound
						if (i>0)
							smiRV = new SmiDiscreteRV(core->getRowStage(i));
						else
							smiRV = new SmiDiscreteRV(core->getColStage(j));
						smiRV->addEvent(*matrix,dclo,dcup,dobj,drlo,drup,smpsCardReader_->getProb());
						matrix->clear();
						dclo.clear();
						dcup.clear();
						dobj.clear();
						drlo.clear();
						drup.clear();

					}
					else
					{
						smiRV->addEvent(*matrix,dclo,dcup,dobj,drlo,drup,smpsCardReader_->getProb());
						matrix->clear();
						dclo.clear();
						dcup.clear();
						dobj.clear();
						drlo.clear();
						drup.clear();
					}




				}
				break;
				case SMI_UNKNOWN_MPS_TYPE:
				default:
					return -2;
				}
			}
			
			if (smpsCardReader_->whichSmpsSection() == SMI_ENDATA_SECTION)
			{
				//process discrete distribution
				smi->processDiscreteDistributionIntoScenarios(smiDD);
				return 0;
			}
			else
				return -2;
			

		}

		case SMI_SCENARIOS_SECTION: // SCENARIOS card
		{
			double prob=0.0;
			int scen=0,anc=0;
			int branch=0;
			CoinPackedVector drlo,drup,dclo,dcup,dobj;
			CoinPackedMatrix *matrix = new CoinPackedMatrix(false,0.25,0.25);
			assert(!matrix->isColOrdered());

			int nrow = core->getNumRows();
			int ncol = core->getNumCols();
			matrix->setDimensions(nrow,ncol);

			while( smpsCardReader_->nextSmpsField (  ) == SMI_SCENARIOS_SECTION ) 
			{
				if (smpsCardReader_->whichSmpsType() == SMI_SC_CARD)
				{
					if (scen)
					{
						
						smi->generateScenario(core,matrix,&dclo,&dcup,&dobj,&drlo,&drup,branch,anc,prob,smpsCardReader_->getCoreCombineRule() );
						matrix->clear();
						dclo.clear();
						dcup.clear();
						dobj.clear();
						drlo.clear();
						drup.clear();
						
						matrix->setDimensions(nrow,ncol);

					}

					std::string SunStudioNeedsThis = smpsCardReader_->scenarioNew() ;
					scenarioMap_.insert(make_pair(SunStudioNeedsThis,scen++));
					prob = smpsCardReader_->value();
					if (!strncmp(smpsCardReader_->scenarioAnc(),"ROOT",4))
						anc = 0;
					else
						anc = scenarioMap_[smpsCardReader_->scenarioAnc()];
					branch = periodMap_[smpsCardReader_->periodName()];
				}
				else if(smpsCardReader_->whichSmpsType() == SMI_COLUMN_CARD)
				{
					int j=this->columnIndex(smpsCardReader_->columnName());
					int i=this->rowIndex(smpsCardReader_->rowName());
					double value = smpsCardReader_->value();

					
					if (j<0) // check RHS
					{
						if (!strcmp(this->getRhsName(),"")) {
							 free(rhsName_);
		 					 rhsName_=strdup(smpsCardReader_->columnName());
						} else
							assert(!strcmp(smpsCardReader_->columnName(),this->getRhsName()));
						assert(!(i<0));

						char c=this->getRowSense()[i];

						switch(c)
						{
						case 'E':
							drlo.insert(i,value);
							drup.insert(i,value);
							break;
						case 'L':
							//drlo.insert(i,value);
							drup.insert(i,value);
							break;
						case 'G':
							//drup.insert(i,value);
							drlo.insert(i,value);
							break;
						default:
							assert(!"bad row sense: shouldn't get here");
							break;
						}
					}
					else if(i<0 || i==getNumRows()) // check OBJ
					{
						assert(!strcmp(smpsCardReader_->rowName(),this->getObjectiveName()));
						dobj.insert(j,value);
					}
					else	// add element
					{
						matrix->modifyCoefficient(i,j,value,true);
					}
				}
				
				else // unknown type
					return -1;

			} // while

			if (smpsCardReader_->whichSmpsSection() == SMI_ENDATA_SECTION)
				{
					if (scen)
					{
						
						smi->generateScenario(core,matrix,&dclo,&dcup,&dobj,&drlo,&drup,branch,anc,prob,smpsCardReader_->getCoreCombineRule() );
						
					}
					delete matrix;
					return 0;
				}

		} 

		default:// didn't recognize section
		return -1;
		}

		

	}



	return 0;
}

//#############################################################################
//  nextSmpsField


SmiSectionType
SmiSmpsCardReader::nextSmpsField (  )
{
  
  // find next non blank character
  char *next = position_;

  static const char *blanks = (const char *) " \t";
  char valstr[COIN_MAX_FIELD_LENGTH],*after;
  
  while ( next != eol_ ) {
    if ( *next == ' ' || *next == '\t' ) {
      next++;
    } else {
      break;
    }
  }

  bool gotCard;
  
  if ( next == eol_ ) {
	  gotCard = false;
  } else {
	  gotCard = true;
  }
  while ( !gotCard ) 
  {
	  // need new image
	  
	  if ( cleanCard() ) {
		  return SMI_EOF_SECTION;
	  }
	  if ( card_[0] == ' ' ) {
		  // not a section or comment
		  position_ = card_;
		  eol_ = card_ + strlen ( card_ );


		  if ( smiSection_ == SMI_TIME_SECTION ) 
		  {
			  smiSmpsType_ = SMI_TIME_ORDERED_CORE_TYPE;
			  if (!(next = strtok(position_,blanks)))
			  {
				  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
				  return smiSection_;
			  }
			  strcpy(columnName_,next);
			  
			  if (!(next = strtok(NULL,blanks)))
			  {
				  strcpy(periodName_,columnName_);
				  smiSmpsType_ = SMI_TIME_UNORDERED_CORE_TYPE;
				  return smiSection_;
			  }
			  strcpy(rowName_,next);
			  if (!(next = strtok(NULL,blanks)))
			  {
				  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
				  return smiSection_;
			  }
			  strcpy(periodName_,next);
			  position_ = eol_;
			  
		  }

		  if ( smiSection_ == SMI_INDEPENDENT_SECTION )
		  {
			  int i;
			  if (strlen(next)==2)
			  {
				  for (i = SMI_BL_CARD; i < SMI_UNKNOWN_MPS_TYPE; ++i) {
					  if ( !strncmp ( next, smpsType[i], 2 ) ) 
					  {
						  break;
					  }
				  }
			  }
			  else
				  i = SMI_COLUMN_CARD;

			  switch(smiSmpsType_ = (SmiSmpsType) i)
			  {
			  case SMI_BL_CARD: // not handled yet
				  break;
				 
			  case SMI_COLUMN_CARD: // card info has "col,row,value,prob"
				    
				  if (!(next = strtok(position_,blanks)))
				  {
					  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
					  return smiSection_;
				  }
				  //already got the name
				  strcpy(columnName_,next);
				  
				  if (!(next = strtok(NULL,blanks)))
				  {
					  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
					  break;
				  }
				  strcpy(rowName_,next);

				  if (!(next = strtok(NULL,blanks)))
				  {
					  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
					  break;
				  }
				  strcpy(valstr,next);

				  value_ = osi_strtod(valstr,&after,0);
				  // see if error
				  assert(after>valstr);
					  
				  if (!(next = strtok(NULL,blanks)))
				  {
					  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
					  break;
				  }
				  strcpy(valstr,next);

				  prob_ = osi_strtod(valstr,&after,0);
				  // see if error
				  assert(after>valstr);

				  position_ = eol_;  //end of card -- no more information allowed
				  break;

			  case SMI_UNKNOWN_MPS_TYPE:
				  break;
			  default:
				  smiSmpsType_=SMI_UNKNOWN_MPS_TYPE;
			  }//end switch

			  return smiSection_;

		  }

		  if ( smiSection_ == SMI_SCENARIOS_SECTION )
		  {
			

			  
			  if (!(next = strtok(position_,blanks)))
			  {
				  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
				  return smiSection_;
			  }
			
			  int i;
			  if (strlen(next)==2)
			  {
				  for (i = SMI_SC_CARD; i < SMI_UNKNOWN_MPS_TYPE; ++i) {
					  if ( !strncmp ( next, smpsType[i], 2 ) ) 
					  {
						  break;
					  }
				  }
				  if (i==SMI_UNKNOWN_MPS_TYPE)
					  i = SMI_COLUMN_CARD;
			  }
			  else
				  i = SMI_COLUMN_CARD;
			
			  switch(smiSmpsType_ = (SmiSmpsType) i)
			  {
			  case SMI_SC_CARD:
				  if (!(next = strtok(NULL,blanks)))
				  {
					  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
					  break;
				  }
				  strcpy(columnName_,next);
				  
				  if (!(next = strtok(NULL,blanks)))
				  {
					  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
					  break;
				  }
				  strcpy(rowName_,next);

				  if (!(next = strtok(NULL,blanks)))
				  {
					  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
					  break;
				  }
				  strcpy(valstr,next);

				  value_ = osi_strtod(valstr,&after,0);
				  // see if error
				  assert(after>valstr);

				  if (!(next = strtok(NULL,blanks)))
				  {
					  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
					  break;
				  }
				  strcpy(periodName_,next);

				  // no futher strings allowed.
				  position_=eol_;
				  break;
			  case SMI_UNKNOWN_MPS_TYPE:
				  //this in case the column name has two characters!
			  case SMI_COLUMN_CARD:

				  //already got the name
				  strcpy(columnName_,next);
				  
				  if (!(next = strtok(NULL,blanks)))
				  {
					  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
					  break;
				  }
				  strcpy(rowName_,next);

				  if (!(next = strtok(NULL,blanks)))
				  {
					  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
					  break;
				  }
				  strcpy(valstr,next);

				  value_ = osi_strtod(valstr,&after,0);
				  // see if error
				  assert(after>valstr);
				  
				  if (!(next = strtok(NULL,blanks)))
					  position_ = eol_;
				  else
					  // strcpy(next,periodName_); // store it here until next time
					  strcpy(periodName_,next); // store it here until next time
				  break;
			  default:
				  assert(smiSmpsType_ == SMI_UNKNOWN_MPS_TYPE);
				  break;
			  }
		  }
		  
		  return smiSection_;

	  } else if ( card_[0] != '*' ) {
		  // not a comment, might be a section
		  int i;
		  
		  handler_->message(COIN_MPS_LINE,messages_)<<cardNumber_
			  <<card_<<CoinMessageEol;

		  // find the section, if there is one
		  for ( i = SMI_NAME_SECTION; i < SMI_UNKNOWN_SECTION; i++ ) {
			//printf("Comparing first 3 chars of %s with %s returns %d \n",
			//	card_, section[i], strncmp(card_,section[i],3));
			  if ( !strncmp ( card_, section[i], 3) ) {
				  break;
			  }
		  }

		  // didn't find anything so quit
		  if (i==SMI_UNKNOWN_SECTION) return (SmiSectionType) i;

		  position_ = card_;
		  eol_ = card_;
		  smiSection_ = ( SmiSectionType ) i;

		  // if its a scenario card, need to process some more info
		  if ( (smiSection_ == SMI_SCENARIOS_SECTION) || 
				(smiSection_ == SMI_INDEPENDENT_SECTION) )
		  {
			  i = SMI_SMPS_COMBINE_UNKNOWN;
			  next = strtok(position_,blanks);
			  
			  if (!(next = strtok(NULL,blanks)))
			  {
				  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
				  //break;
			  }
			  // find the section, if there is one
			  // next card should be DISCRETE
			  if (!(next = strtok(NULL,blanks)))
			  {
				  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
				  //break;
			  }
			  // find the section, if there is one
			  if (next == NULL) {
			    i = SMI_SMPS_COMBINE_REPLACE; 
			  } else {
			      for ( i = SMI_SMPS_COMBINE_ADD; i < SMI_SMPS_COMBINE_UNKNOWN; i++ ) {
				      if ( !strncmp ( next, smpsType[i], strlen ( section[i] ) ) ) {
					      break;
				      }
			      }
			  }
		   }
			  // set combine rule if it is not already set.
		   if (!combineRuleSet)
		   {
			  switch(i)
			  {
			  case SMI_SMPS_COMBINE_ADD:
				  this->setCoreCombineRule(SmiCoreCombineAdd::Instance());
				  break;
			  case SMI_SMPS_COMBINE_REPLACE:
				  this->setCoreCombineRule(SmiCoreCombineReplace::Instance());
				  break;
			  default:
				  this->setCoreCombineRule(SmiCoreCombineReplace::Instance());
				  // MESSAGE
				  printf(" Smps: setting default core combine rule to Replace\n");
			  }
			  
		  }
		  return smiSection_;
	  } else {
		  // comment
	  }
  }
  
  {
	  // if we get here, there are additional row fields on the line

				  strcpy(rowName_,periodName_); // stored the last field in periodName_
				  
				  if (!(next = strtok(NULL,blanks)))
				  {
					  smiSmpsType_ = SMI_UNKNOWN_MPS_TYPE;
					  return smiSection_;
				  }
				  strcpy(valstr,next);
				  
				  value_ = osi_strtod(valstr,&after,0);
				  // see if error
				  assert(after>valstr);

				  if (!(next = strtok(NULL,blanks)))
					  position_ = eol_;
				  else
					  strcpy(next,periodName_);
				  
				  
				  return smiSection_;
  }
}

std::string SmiSmpsIO::getModProblemName() {
    std::string name = "";
    if (strcmp(problemName_,"")==0) {
        name.append("BLANK   ");
    } else {
        if (strlen(problemName_) >= 8) {
            name.append(problemName_, 8);
        } else {
            name.append(problemName_);
            name.append(8-strlen(problemName_), ' ');
        }
    }
    return name;
}

void SmiSmpsIO::writeCoreFile(const char* filename, const char* extension, const bool strictFormat) {
    // count nels
    int nels = 0;
    for (int t = 0; t < core->getNumStages(); t++) {
        nels += core->getNode(t)->getNumMatrixElements();
    }
    
    // initialize arrays
    double* clo = new double[core->getNumCols()];
    double* cup = new double[core->getNumCols()];
    double* obj = new double[core->getNumCols()];
    double* rlo = new double[core->getNumRows()];
    double* rup = new double[core->getNumRows()];
    
    // initialize row-ordered matrix arrays
    double * dels = new double[nels];
    int * indx = new int[nels];
    int * rstrt = new int[core->getNumRows()+1];
    rstrt[0] = 0;

    int rowCount = 0, ncol = 0, nrow = 0;
    nels = 0;

    for (int t = 0; t < core->getNumStages(); t++) {

        SmiNodeData * node = core->getNode(t);
    	int stg = node->getStage();

	    core->copyRowLower(rlo+nrow,stg);
	    core->copyRowUpper(rup+nrow,stg);
	    core->copyColLower(clo+ncol,stg);
	    core->copyColUpper(cup+ncol,stg);
	    core->copyObjective(obj+ncol,stg);

	    // add rows to det. eq. matrix for current stage
	    for (int i=core->getRowStart(stg); i<core->getRowStart(stg+1) ; i++)
	    {
		    // build row explicitly into sparse arrays
		    int rowStart=rstrt[rowCount];
		    int rowNumEls=0;

		    const double *cels=node->getRowElements(i);
		    const int *cind=node->getRowIndices(i);
		    const int len=node->getRowLength(i);
		    memcpy(dels+rowStart,cels,sizeof(double)*len);
		    memcpy(indx+rowStart,cind,sizeof(int)*len);
		    rowNumEls=len;

		    //preparation for the next row
		    rowCount++;
		    nels+=rowNumEls;
		    rstrt[rowCount] = nels;
		    //done with preparations for the next row
        }
        
        ncol += core->getNumCols(stg);
	    nrow += core->getNumRows(stg);
    }
    // 
    if (solverInf_ != this->infinity_){
    // We have to transform the values from solverInfinity to infinity
        for (int i = 0; i < core->getNumCols(); i++){
            if (clo[i] == -solverInf_)
                clo[i] = -infinity_;
            if (clo[i] == solverInf_)
                clo[i] = infinity_;
            if (cup[i] == -solverInf_)
                cup[i] = -infinity_;
            if (cup[i] == solverInf_)
                cup[i] = infinity_;
            if (obj[i] == -solverInf_)
                obj[i] = -infinity_;
            if (obj[i] == solverInf_)
                obj[i] = infinity_;
        }
        for (int i = 0; i < core->getNumRows(); i++){
            if (rlo[i] == -solverInf_)
                rlo[i] = -infinity_;
            if (rlo[i] == solverInf_)
                rlo[i] = infinity_;
            if (rup[i] == -solverInf_)
                rup[i] = -infinity_;
            if (rup[i] == solverInf_)
                rup[i] = infinity_;
        }
        // We assume that there are no infinity values in the coefficient matrix..
    }
    
    CoinPackedMatrix matrix(false, ncol, nrow, nels, dels, indx, rstrt, NULL);
    
    delete [] dels;
    delete [] indx;
    delete [] rstrt;
       
    this->setMpsData(matrix, this->infinity_, clo, cup, obj, NULL, rlo, rup, this->core->getColumnNames(strictFormat), NULL);

    delete [] clo;
    delete [] cup;
    delete [] obj;
    delete [] rlo;
    delete [] rup;

    this->setProblemName(filename);

    char* integrality = new char[ncol];
    std::fill_n(integrality, ncol, 0);
    for (int i = 0; i < core->getIntegerLength(); i++) {
        integrality[core->getIntegerIndices()[i]] = 1;
    }
    this->copyInIntegerInformation(integrality);
    delete [] integrality;
    
    std::string fnWithExt = filename;
    fnWithExt.append(".").append(extension);
    
    if (strictFormat)
        this->writeMps(fnWithExt.c_str());
    else
        this->writeMps(fnWithExt.c_str(), 0, 1);
}

void SmiSmpsIO::writeTimeFile(const char* filename, const char* extension, const bool strictFormat) {
    std::string fnWithExt = filename;
    fnWithExt.append(".").append(extension);
    CoinFileOutput* output = CoinFileOutput::create(fnWithExt, CoinFileOutput::COMPRESS_NONE);

    std::stringstream line;
    line.fill(' ');
    line << "TIME          " << this->getModProblemName() << "\nPERIODS\n";
    
    for (int stg = 0; stg < core->getNumStages(); stg++) {
        line << "    ";
        if (core->getColStart(stg) == core->getNumCols())
            line << std::setw(10) << std::left << "RHS"; // if this stage has no columns, it may have some RHS values
        else if (strictFormat)
            line << std::setw(10) << std::left << this->columnName(core->getColStart(stg));
        else
            line << this->columnName(core->getColStart(stg)) << " ";
        line << std::setw(25) << std::left << this->rowName(core->getRowStart(stg));
        line << "STAGE_" << stg << "\n";
    }
    line << "ENDATA\n";
    output->puts(line.str());        
    delete output;
}

void SmiSmpsIO::writeStochFile(const char* filename, const char* extension, const bool strictFormat) {
    std::string fnWithExt = filename;
    fnWithExt.append(".").append(extension);
    CoinFileOutput* output = CoinFileOutput::create(fnWithExt, CoinFileOutput::COMPRESS_NONE);
    
    std::ostringstream line;
    line.fill(' ');

    if (strictFormat) {
        line << setprecision(11);
    }
    else {
        line << setprecision(20);
    }
    
    line << "STOCH         " << this->getModProblemName() << "\n";
    line << "SCENARIOS     DISCRETE\n";
    //line << "SCENARIOS     DISCRETE                 ";
    //if (tree->getLeaf(0)->getDataPtr()->getNode()->getCoreCombineRule() == SmiCoreCombineReplace::Instance())
    //    line << "REPLACE\n";
    //else
    //    line << "ADD\n";
        
    for (int scen = 0; scen < tree->getNumScenarios(); scen++) {
        writeScenarioToStochFile(line, tree->getLeaf(scen), scen, strictFormat);
    }
    
    line << "ENDATA\n";
    output->puts(line.str());        
    delete output;
}

void SmiSmpsIO::writeScenarioToStochFile(std::ostringstream& stream, SmiTreeNode<SmiScnNode *> * node, int scenario, bool strictFormat) {
    // first, go up the tree - if possible and the parent node is still within the same scenario and not the root node
    // if thats not the case, we reached the first node for that scenario
    if (node->hasParent() && node->getParent()->scenario() == scenario && node->getParent() != tree->getRoot()) {
        // go up the tree
        writeScenarioToStochFile(stream, node->getParent(), scenario, strictFormat);
    } else {
        // write scenario card
        stream << " SC Scen_" << std::setw(5) << std::left << node->scenario();

        if (node->getParent() != tree->getRoot()) {
            // branches from parent scenario
            stream << "Scen_" << std::setw(5) << std::left << node->getParent()->scenario();
        } else {
            // branches from root, i.e. Scen_0.
            if (scenario == 0) // First Scenario (Scen_0) branches from ROOT, other scenarios branches from Scen_0
                stream << std::setw(10) << std::left << "'ROOT'";
            else
                stream << std::setw(10) << std::left << "Scen_0";
        }

        if (strictFormat)
            stream << std::setw(15) << std::left << tree->getLeaf(scenario)->getDataPtr()->getProb(); // write the probabilty
        else
            stream << tree->getLeaf(scenario)->getDataPtr()->getProb() << " "; // write the probabilty

        stream << "STAGE_" << node->getDataPtr()->getStage() << "\n"; // write the stage number
    }
    
    // now we can print the stochastic data
    SmiNodeData* data = node->getDataPtr()->getNode();
    if (data->getNumMatrixElements() > 0) {
        // write matrix elements
        for (int i = data->getCore()->getRowStart(data->getStage()); i < data->getCore()->getRowStart(data->getStage()+1); i++) {
            for (int j = 0; j < data->getRowLength(i); j++) {
                stream << "    " << std::setw(10) << std::left << this->columnName(data->getRowIndices(i)[j]);
                stream << std::setw(10) << std::left << this->rowName(i);
                stream << data->getRowElements(i)[j] << "\n";
            }
        }
    }

    // write row lower elements
    for (int i = 0; i < data->getRowLowerLength(); i++) {
        stream << "    RHS       " << std::setw(10) << std::left << this->rowName(data->getRowLowerIndices()[i]) << data->getRowLowerElements()[i] << "\n";
    }

    // write row upper elements
    for (int i = 0; i < data->getRowUpperLength(); i++) {
        if (this->getRowSense()[data->getRowUpperIndices()[i]] != 'E')
            stream << "    RHS       " << std::setw(10) << std::left << this->rowName(data->getRowUpperIndices()[i]) << data->getRowUpperElements()[i] << "\n";
    }

    // write objective elements
    for (int i = 0; i < data->getObjectiveLength(); i++) {
        stream << "    " << std::setw(10) << std::left << this->columnName(data->getObjectiveIndices()[i]) << "OBJROW         " << data->getObjectiveElements()[i] << "\n";
    }

    // write column bounds    
    if (data->getColLowerLength() > 0 || data->getColUpperLength() > 0) {
        int icl = 0, icu = 0;
        while (icl < data->getColLowerLength() || icu < data->getColUpperLength()) {
            
            if (icl < data->getColLowerLength() &&  icu < data->getColUpperLength() && data->getColLowerIndices()[icl] == data->getColUpperIndices()[icu]) {
                // lower and upper bound given for current column
                if (data->getColLowerElements()[icl] == data->getColUpperElements()[icu]) {
                    // variable is fixed
                    stream << " FX BOUND     " << std::setw(10) << std::left << this->columnName(data->getColLowerIndices()[icl]) << data->getColLowerElements()[icl];
                } 
                else if (data->getColLowerElements()[icl] == 0.0 && data->getColUpperElements()[icu] == 1.0) {
                    // variable is binary
                    stream << " BV BOUND     " << std::setw(10) << std::left << this->columnName(data->getColLowerIndices()[icl]);
                }
                else if (data->getColLowerElements()[icl] == -this->solverInf_ && data->getColUpperElements()[icu] == this->solverInf_) {
                    // variable is free
                    stream << " FR BOUND     " << std::setw(10) << std::left << this->columnName(data->getColLowerIndices()[icl]);
                }
                icl++; icu++;
                continue;
            }
            
            if (icu >= data->getColUpperLength() || data->getColLowerIndices()[icl] < data->getColUpperIndices()[icu]) {
                // write column lower bound
                if (this->isInteger(icl)) {
                    // column lower bound for integer variable
                    stream << " LI BOUND     " << std::setw(10) << std::left << this->columnName(data->getColLowerIndices()[icl]) << data->getColLowerElements()[icl];
                }
                else {
                    if (data->getColLowerElements()[icl] == -this->solverInf_) {
                        // lower bound is -infinity
                        stream << " MI BOUND     " << std::setw(10) << std::left << this->columnName(data->getColLowerIndices()[icl]);
                    }
                    else {
                        // normal column lower bound
                        stream << " LO BOUND     " << std::setw(10) << std::left << this->columnName(data->getColLowerIndices()[icl]) << data->getColLowerElements()[icl];
                    }
                }
                icl++;
            }
            else {
                // write column upper bound
                if (this->isInteger(icu)) {
                    // column upper bound for integer variable
                    stream << " UI BOUND     " << std::setw(10) << std::left << this->columnName(data->getColUpperIndices()[icu]) << data->getColUpperElements()[icu];
                }
                else {
                    if (data->getColUpperElements()[icu] == this->solverInf_) {
                        // upper bound is infinity
                        stream << " PL BOUND     " << std::setw(10) << std::left << this->columnName(data->getColUpperIndices()[icu]);
                    }
                    else {
                        // normal column upper bound
                        stream << " UP BOUND     " << std::setw(10) << std::left << this->columnName(data->getColUpperIndices()[icu]) << data->getColUpperElements()[icu];
                    }
                }
                icu++;
            }
        }
    }

}

void SmiSmpsIO::writeSmps(const char* filename, bool winFileExtensions, bool strictFormat) {
    if (winFileExtensions) {
        //write core file
        this->writeCoreFile(filename, "cor", strictFormat);
        
        // write time file
        this->writeTimeFile(filename, "tim", strictFormat);
        
        // write stoch file
        this->writeStochFile(filename, "sto", strictFormat);
    } else {
        //write core file
        this->writeCoreFile(filename, "core", strictFormat);
        
        // write time file
        this->writeTimeFile(filename, "time", strictFormat);
        
        // write stoch file
        this->writeStochFile(filename, "stoch", strictFormat);
    }
    
    this->freeAll();
}

