#if ! defined (PHY_PAT_CLASS_PROB_HPP)
#define PHY_PAT_CLASS_PROB_HPP
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include "ncl/nxsallocatematrix.h"

////////////////////////////////////////////////////////////////////////////////
// header
////////////////////////////////////////////////////////////////////////////////
class NxsCharactersBlock;
class NxsSimpleNode;
class NxsSimpleTree;
class PatternSummary;
class ParsInfo;

typedef unsigned char BitField;
typedef std::vector<BitField> BitFieldRow;
typedef std::vector<BitFieldRow> BitFieldMatrix;
typedef std::map<BitField, unsigned> BitsToCount;
typedef std::pair<const NxsSimpleNode *, unsigned> NodeID;
typedef double ** TiMat;
typedef void (* TiMatFunc)(double, TiMat);


class ParsInfo {
	public:
		ParsInfo()
			:downPass(0),
			allSeen(0),
			score(0),
			numPatterns(0) {
			}
		unsigned size() const {
			return numPatterns;
		}
		
		void calculateForTip(const BitFieldRow & data);
		void calculateForInternal(const ParsInfo & leftData, const ParsInfo & rightData); 

		void write(std::ostream & o) const;
		
		
		const BitField * downPass; // alias
		const BitField * allSeen; // alias
		const unsigned * score; // alias
		
	private:
		unsigned numPatterns;
		BitFieldRow downPassOwned;
		BitFieldRow allSeenOwned;
		std::vector<unsigned> scoreOwned;
		
};
typedef std::map<NodeID, ParsInfo> NodeIDToParsInfo;


class ProbForParsScore{
    
};    
class ProbInfo {
	public:
		void calculate(const ProbInfo & leftPI, double leftEdgeLen, 
					   const ProbInfo & rightPI, double rightEdgeLen,
					   TiMatFunc fn);
		unsigned getMaxParScore() const {
		    assert(!this->byParsScore.empty());
		    return this->byParsScore.size() - 1;
		}
	protected:
	    std::vector<ProbForParsScore> byParsScore;
};

typedef std::map<NodeID, ProbInfo *> NodeIDToProbInfo;


class PatternSummary {
	public:
		void clear() {
			this->byParsScore.clear();
		}
		unsigned incrementCount(unsigned s, BitField m);
		
		void write(std::ostream & out) const;
		
	private:
		
		std::vector<BitsToCount> byParsScore; // for each parsimony score, this stores a map of the "mask" of bits and the count
		
};

void calculatePatternClassProbabilities(const NxsSimpleTree & tree, std::ostream & out );
void classifyObservedDataIntoClasses(const NxsSimpleTree & tree, const BitFieldMatrix &, std::ostream & out, PatternSummary *);
std::string convertToBitFieldMatrix(const NxsCharactersBlock & cb, BitFieldMatrix & mat, const std::set<unsigned> * toInclude=0L);
void writeBitFieldMatrix(std::ostream & out, const BitFieldMatrix & bitFieldMatrix);

void JCTiMat(double edgeLen, TiMat pMat);
	
	
#endif
