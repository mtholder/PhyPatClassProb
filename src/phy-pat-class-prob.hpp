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
class CommonInfo;

typedef unsigned char BitField;
typedef std::vector<BitField> BitFieldRow;
typedef std::vector<BitFieldRow> BitFieldMatrix;
typedef std::map<BitField, unsigned> BitsToCount;
typedef std::pair<const NxsSimpleNode *, unsigned> NodeID;
typedef double ** TiMat;
typedef double *** TiMatVec;
typedef void (* TiMatFunc)(double, TiMatVec);


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
		
		void calculateForTip(const BitFieldRow & data, const CommonInfo &);
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

typedef std::map<BitField, std::vector<double> > MaskToProbsByState;
typedef std::map<BitField, MaskToProbsByState > MaskToMaskToProbsByState;

const std::vector<double> * getProbsForStatesMask(const MaskToProbsByState *, const BitField sc);

inline const std::vector<double> * getProbsForStatesMask(const MaskToProbsByState *m, const BitField sc) {
    if (m == 0L)
        return 0L;
    MaskToProbsByState::const_iterator scIt = m->find(sc);
    if (scIt == m->end())
        return 0L;
    return &(scIt->second);
}

class ProbForParsScore{
    public:
    
        const MaskToProbsByState * getMapPtrForDownPass(const BitField sc) const {
            MaskToMaskToProbsByState::const_iterator scIt = this->byDownPass.find(sc);
            if (scIt == this->byDownPass.end())
                return 0L;
            return &(scIt->second);
        }
        const std::vector<double> * getProbsForDownPassAndMask(const BitField downPass, const BitField mask) const {
            return getProbsForStatesMask(this-> getMapPtrForDownPass(downPass), mask);
        }

    private:
        MaskToMaskToProbsByState byDownPass;
        friend class ProbInfo;
        
        
};

class ProbInfo {
	public:
	    void createForTip(const CommonInfo &);
		void calculate(const ProbInfo & leftPI, double leftEdgeLen, 
					   const ProbInfo & rightPI, double rightEdgeLen,
					   TiMatFunc fn, const CommonInfo &);
		unsigned getMaxParsScore() const {
		    assert(!this->byParsScore.empty());
		    return this->byParsScore.size() - 1;
		}
		const ProbForParsScore & getByParsScore(unsigned score) const {
		    return this->byParsScore.at(score);
		}
		unsigned getNLeavesBelow() const {
		    return nLeavesBelow;
		}
	protected:
		unsigned nLeavesBelow;
	
        void addToAncProbVec(std::vector<double> & pVec, 
            const double *** leftPMatVec, const std::vector<double> * leftProbs,
            const double *** rightPMatVec, const std::vector<double> * rightProbs,
            const CommonInfo & blob);
	    std::vector<ProbForParsScore> byParsScore;
};

typedef std::map<NodeID, ProbInfo *> NodeIDToProbInfo;


class PatternSummary {
	public:
		void clear() {
			this->byParsScore.clear();
		}
		unsigned incrementCount(unsigned s, BitField m);
		
		void write(std::ostream & out, const CommonInfo &) const;
		
	private:
		
		std::vector<BitsToCount> byParsScore; // for each parsimony score, this stores a map of the "mask" of bits and the count
		
};

void calculatePatternClassProbabilities(const NxsSimpleTree & tree, std::ostream & out, const CommonInfo &);
void classifyObservedDataIntoClasses(const NxsSimpleTree & tree, const BitFieldMatrix &, std::ostream & out, PatternSummary *, const CommonInfo &);
std::string convertToBitFieldMatrix(const NxsCharactersBlock & cb, BitFieldMatrix & mat, const std::set<unsigned> * toInclude=0L);
void writeBitFieldMatrix(std::ostream & out, const BitFieldMatrix & bitFieldMatrix);

void JCTiMat(double edgeLen, TiMat pMat);
void JCMulitCatTiMat(double edgeLen, TiMatVec pMatVec);
	
#endif
