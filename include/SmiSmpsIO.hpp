// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef SmiSmpsIO_HPP
#define SmiSmpsIO_HPP

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CoinMpsIO.hpp"
#include "SmiScnModel.hpp"
#include "SmiScnData.hpp"
#include "SmiScenarioTree.hpp"

#include <vector>
#include <map>
#include <string>

using namespace std;


// sections
enum SmiSectionType { SMI_NO_SECTION, SMI_NAME_SECTION,
					  SMI_ENDATA_SECTION, SMI_EOF_SECTION, 
					  SMI_TIME_SECTION, SMI_SCENARIOS_SECTION, 
					  SMI_UNKNOWN_SECTION
};
#if 1
const static char *section[] = {
  "", "NAME", "ENDATA", " ", "PERIODS", "SCENARIOS", " "
};

const static char *smpsType[] = {
	"SC","  ",""
};
#endif

enum SmiSmpsType { SMI_SC_CARD, SMI_COLUMN_CARD, SMI_TIME_UNORDERED_CORE_TYPE, 
	SMI_TIME_ORDERED_CORE_TYPE,SMI_UNKNOWN_MPS_TYPE
};

class SmiSmpsCardReader:
public CoinMpsCardReader
{
public:
	  SmiSectionType nextSmpsField (  );
	  SmiSectionType whichSmpsSection(){return smiSection_;}
	    
	  SmiSmpsType whichSmpsType() {return smiSmpsType_;}
	
	  inline const char *periodName (  ) const {return periodName_;};
	  inline const char *scenarioNew (  ) const {return columnName_;};
	  inline const char *scenarioAnc (  ) const {return rowName_;};
	  /// Constructor expects file to be open 
	  /// This one takes gzFile if fp null
	  SmiSmpsCardReader( FILE * fp, gzFile gzfp, CoinMpsIO * reader ):CoinMpsCardReader (fp, gzfp,reader ){}

	  ~SmiSmpsCardReader(){}
private:
	 /// Current third name (for SmpsIO)
	char periodName_[MAX_FIELD_LENGTH];
	float fvalue_;
	SmiSectionType smiSection_;
	SmiSmpsType smiSmpsType_;
};

class SmiSmpsIO: 
public CoinMpsIO
{
public:
	int readTimeFile(SmiScnModel *smi,const char *c,const char *ext="time");
	int readStochFile(SmiScnModel *smi,const char *c,const char *ext="stoch");
	inline int getNumStages(){ return nstag_;}
	inline int *getColumnStages(){ return cstag_;}
	inline int *getRowStages(){ return rstag_;}	
public:
	SmiSmpsIO():CoinMpsIO(),nstag_(0),cstag_(NULL),rstag_(NULL),iftime(false){}
	~SmiSmpsIO(){delete [] cstag_;delete[] rstag_;}
private:
	int nstag_;
	int *cstag_;
	int *rstag_;
	typedef 	std::map<string,int> StringIntMap;
	StringIntMap periodMap_;
	StringIntMap scenarioMap_;
	bool iftime,ifstoch;
	SmiSmpsCardReader *smpsCardReader_;
	
};









#endif //#define SmiSmpsIO_HPP
