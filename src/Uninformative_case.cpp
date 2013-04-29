#include "phy-pat-class-prob.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <libhmsbeagle/beagle.h>

#include "ncl/nxsmultiformat.h"
#include "ncl/nxsallocatematrix.h"


#include "pytbeaglehon/ccore/phylo_util.h"
#include "pytbeaglehon/ccore/asrv.h"
#include "pytbeaglehon/ccore/calc_instance.h"
#include "pytbeaglehon/ccore/discrete_state_model.h"
#include "pytbeaglehon/ccore/internal_like_calc_env.h"
NxsString errormsg;

#define DEBUGGING_OUTPUT

using namespace std;

void freeProbInfo(const std::vector<const NxsSimpleNode *> & preorderVec, NodeIDToProbInfo & nodeIDToProbInfo);

int getNextCommStSet(const int obsStSet, int i);

// Globals (all should start with g[A-Z] pattern to make it easier to find and replace them later).
const unsigned MAX_NUM_STATES = 8*sizeof(BitField);

typedef std::vector<int> stateSetContainer;

struct CommonInfo {
        bool isSymmetric;
        unsigned nStates;
        unsigned nRates;
        unsigned pVecLen; // nStates * nRates;
        BitField lastBitField;
        TiMatFunc tiMatFunc;

        std::string alphabet;
        std::vector<BitField> singleStateCodes;
        // stateIndexToStateCode is probably always going to be identical to
        //      singleStateCodes, but use stateIndexToStateCode when you need the
        //      state codes for the fundamental states in order.
        std::vector<BitField> stateIndexToStateCode;
        std::vector<BitField> multiStateCodes;
        ScopedDblThreeDMatrix firstMatVec;
        ScopedDblThreeDMatrix secondMatVec;
        std::vector<unsigned> zeroVec;
        std::vector<double> rates;
        std::vector<double> rateProb;
        std::vector<double> categStateProb;

        VMaskToVecMaskPair pairsForIntersectionForEachDownPass;
        VMaskToVecMaskPair pairsForUnionForEachDownPass;
        BitFieldMatrix statesSupersets;

        unsigned getNumStates(BitField mask) const {
            return stateCodeToNumStates[mask];
        }
        const std::string & toSymbol(BitField sc) const {
            return this->stateCodesToSymbols.find(sc)->second;
        }
        void initialize(const std::string & symbols);

        void readModel(const std::vector<std::string> & optVec);

        stateSetContainer::const_iterator stateSetBegin() const  {
            return possObsStateSet.begin();
        }

        stateSetContainer::const_iterator stateSetEnd() const  {
            return possObsStateSet.end();
        }



        DSCTModelObj * dsct_model_obj;

        std::vector<double> scaledEdgeLengths;
        std::vector<int> modelIndexVector;
        long likeCalcHandle;


        void writeModel(std::ostream & out) const;
    private:
        stateSetContainer possObsStateSet;
        void calcEigenSolution();

        std::vector<unsigned> stateCodeToNumStates;
        std::map<BitField, std::string> stateCodesToSymbols;

        ScopedDblTwoDMatrix relRateMat;
        std::vector<double> stateFreqVector;
};

CommonInfo * gBlob = 0L;

int countBits(int);

int countBits(int x)
{
    int num = 0;
    while(x > 0)
    {
        if(x& 1)
            ++num;
        x=(x>>1);
    }
    return num;
}
int convertIndexToBit(int);

int convertIndexToBit(int ind) {
    return 1 << ind; //bit shift
}

int convertBitToIndex(int);

int convertBitToIndex(int i)
{
    int ind = 0;
    while(i > 0)
    {
        if(i == 1) {  //& = bitwise intersection operator NOT reference
            return ind;
        }

        if(i& 1) {
            std::cerr << "Illegal Value to Convert Bits! \n";
            exit(1);
        }
       ind++;
       i = (i>>1);
    }
    std::cerr << "Zero to Convert Bits! \n";
    exit(1);
};

vector<int> subsetsContainingGivenState(int, int);

vector<int> subsetsContainingGivenState(int fullSet, int givenState)
{
    set<int> subsets;
    int i = 1;
    while(i<=fullSet)
    {
        std::cerr << " subsetsContainingGivenState " << fullSet << " " << givenState << " " << i << "\n";
        int j = i& fullSet;
        if(j& givenState)
            subsets.insert(j);
        i++;
    }
    return vector<int> (subsets.begin(), subsets.end());
}

vector<int> subsetsOfGivenSize(int obsStSet, int numBits);

vector<int> subsetsOfGivenSize(int obsStSet, int numBits)
{
    set<int> subsets;
    int i = 1;
    while(i<=obsStSet)
    {
        std::cerr << " subsetsOfGivenSize " << obsStSet << " " << numBits << " " << i << "\n";
        int j = i& obsStSet;
        if(countBits(j) == numBits)
            subsets.insert(j); //takes in subsets, disregards repeated values
        i++;
    }
    return vector<int> (subsets.begin(), subsets.end()); //creates iterators
}


//this function will give a -1 initially, and will return the index for the next state in the obs. ss, else -2
int getNextCommStSet(const int obsStSet, int i) {
#   if defined DEBUGGING_OUTPUT2
                     std::cerr << "from line " << __LINE__ << ":  getNextCommStSet " ; std::cerr << "obsStSet =  " << obsStSet << " i = "<< i << "\n";
#   endif
   int ind;
    int binRep;
    if(i==-1) {
        ind = 0; //index of state in normal counting sequence
        binRep = 1; // 2^i
    }
    else {
        ind = i+1;
        binRep = 1<<ind;
        if(binRep>obsStSet)
            return -2; //protects against infinite loop
    }
    while((binRep & obsStSet) == 0)  {
        binRep <<= 1;
        ind++;
    }
    return ind;

}


std::set<BitField> toElements(BitField sc) {
    std::set<BitField> ret;
    BitField curr = 1;
    while (curr <= sc) {
        if (curr & sc)
            ret.insert(curr);
        curr *= 2;
    }
    return ret;
}

void freeProbInfo(const std::vector<const NxsSimpleNode *> & preorderVec,
                  NodeIDToProbInfo & nodeIDToProbInfo) {

    for (std::vector<const NxsSimpleNode *>::const_reverse_iterator ndIt = preorderVec.rbegin();
            ndIt != preorderVec.rend();
            ++ndIt) {
        const NxsSimpleNode * nd = *ndIt;
        std::vector<NxsSimpleNode *> children = nd->GetChildren();
        const unsigned numChildren = children.size();
        NodeID currNdId(nd, 0);
        if (numChildren > 0) {
            delete nodeIDToProbInfo[currNdId];
        }
    }
}


void CommonInfo::initialize(const std::string & symbols) {
    this->isSymmetric = false;
    this->alphabet = symbols;
    this->nStates = this->alphabet.length();
    this->nRates = 1; //@TEMP no rate het
    this->pVecLen = this->nRates*this->nStates;
    int endIndex = 1<<(nStates); //end = stopping value
    for(int i=1; i<endIndex; i++) {
        possObsStateSet.push_back(i);
    }

    this->rates.assign(this->nRates, 1.0); //@TEMP no rate het
    this->rateProb.assign(this->nRates, 1.0/this->nRates); // @TEMP equally-sized rate categories
    this->categStateProb.assign(this->pVecLen, 1.0/((double)this->pVecLen)) ;

    this->singleStateCodes.clear();
    this->multiStateCodes.clear();
    this->stateCodesToSymbols.clear();


    unsigned lbfU = (1 << this->nStates) - 1;
    this->stateCodeToNumStates.assign(lbfU + 1, 0);
    this->lastBitField = BitField(lbfU);
    BitField sc = 1;
    for (;;++sc) {
        const std::set<BitField> sbf = toElements(sc);
        if (sbf.size() == 1) {
            unsigned stInd = this->singleStateCodes.size();
            this->singleStateCodes.push_back(sc);
            this->stateIndexToStateCode.push_back(sc);
            assert(this->stateIndexToStateCode[stInd] == sc);
            this->stateCodesToSymbols[sc] = this->alphabet[stInd];
            this->stateCodeToNumStates.at(sc) = 1;
        }
        else {
            this->multiStateCodes.push_back(sc);
            std::string sym;
            for (std::set<BitField>::const_iterator sbfIt = sbf.begin(); sbfIt != sbf.end(); ++sbfIt)
                sym.append(this->stateCodesToSymbols[*sbfIt]);
            this->stateCodeToNumStates.at(sc) = sbf.size();
            this->stateCodesToSymbols[sc] = sym;
        }



        if (sc == this->lastBitField)
            break;
    }

    this->pairsForUnionForEachDownPass.clear();
    this->pairsForUnionForEachDownPass.resize(this->lastBitField + 1);
    this->pairsForIntersectionForEachDownPass.clear();
    this->pairsForIntersectionForEachDownPass.resize(this->lastBitField + 1);

    for (sc = 1;;++sc) {
        if (this->getNumStates(sc) > 1) {
            VecMaskPair & forUnions = this->pairsForUnionForEachDownPass[sc];
            for (BitField leftSC = 1; leftSC < sc ; ++leftSC) {
                if ((leftSC | sc) != sc)
                    continue;
                BitField rightSC = sc - leftSC;
                assert((rightSC | leftSC) == sc);
                forUnions.push_back(MaskPair(leftSC, rightSC));
            }
        }

        VecMaskPair & forIntersections = this->pairsForIntersectionForEachDownPass[sc];
        for (BitField leftSC = 1; leftSC <= this->lastBitField ; ++leftSC) {
            for (BitField rightSC = 1; rightSC <= this->lastBitField ; ++rightSC) {
                if ((leftSC & rightSC) != sc)
                    continue;
                forIntersections.push_back(MaskPair(leftSC, rightSC));
            }
        }

        if (sc == this->lastBitField)
            break;

    }

    this->statesSupersets.clear();
    this->statesSupersets.resize(this->lastBitField + 1);
    for (sc = 1;;++sc) {
        BitFieldRow & ssRow = this->statesSupersets[sc];
        for (BitField ss = sc;; ++ss) {
            if ((ss & sc) == sc)
                ssRow.push_back(ss);
            if (ss == this->lastBitField)
                break;
        }

        if (sc == this->lastBitField)
            break;

    }



#   if defined DEBUGGING_OUTPUT
        std::cerr << "from line " << __LINE__ << ":\n";
#   endif
    for (sc = 1;;++sc) {
#       if defined DEBUGGING_OUTPUT
            std::cerr << "Combos with unions that lead to " << this->toSymbol(sc) << ":";
#       endif
        const VecMaskPair & forUnions = this->pairsForUnionForEachDownPass[sc];
#       if defined DEBUGGING_OUTPUT
            for (VecMaskPair::const_iterator fuIt = forUnions.begin(); fuIt != forUnions.end(); ++fuIt) {
                std::cerr << " (\"" << this->toSymbol(fuIt->first);
                std::cerr << "\", \"" << this->toSymbol(fuIt->second);
                std::cerr << "\")   ";
            }
            std::cerr << "\n";

            std::cerr << "Combos with intersections that lead to " << this->toSymbol(sc) << ":";
#       endif
        const VecMaskPair & forIntersections = this->pairsForIntersectionForEachDownPass[sc];
#       if defined DEBUGGING_OUTPUT
            for (VecMaskPair::const_iterator fuIt = forIntersections.begin(); fuIt != forIntersections.end(); ++fuIt) {
                std::cerr << " (\"" << this->toSymbol(fuIt->first);
                std::cerr << "\", \"" << this->toSymbol(fuIt->second);
                std::cerr << "\")   ";
            }
            std::cerr << "\n";
#       endif

        if (sc == this->lastBitField)
            break;
    }


}


void ProbInfo::createForTip(const CommonInfo & blob) {
    this->byParsScore.resize(1);
    ProbForParsScore & forZeroSteps = this->byParsScore[0];
    unsigned stateIndex = 0;
    for (std::vector<BitField>::const_iterator scIt = blob.singleStateCodes.begin();
            scIt != blob.singleStateCodes.end();
            ++scIt, ++stateIndex) {
        const BitField sc = *scIt;
        std::vector<double> & pVec = forZeroSteps.byDownPass[sc][sc];
        pVec.assign(blob.nRates*blob.nStates, 0.0);
        for (unsigned r = 0; r < blob.nRates; ++r)
            pVec[blob.nStates*r + stateIndex] = 1.0;
    }
    this->nLeavesBelow = 1;
}




void JCTiMat(double edgeLength, TiMat pMat ) {
#   if defined DEBUGGING_OUTPUT3
        std::cerr << "JCTiMat edgeLength = " << edgeLength << '\n';
#   endif
    const double exp_term = exp(-(4.0/3.0)*edgeLength);
    const double prob_change = 0.25 - 0.25*exp_term;
    const double prob_nochange = 0.25 + 0.75*exp_term;
    assert(prob_change >= 0.0);
    pMat[0][0] = prob_nochange;
    pMat[0][1] = prob_change;
    pMat[0][2] = prob_change;
    pMat[0][3] = prob_change;
    pMat[1][0] = prob_change;
    pMat[1][1] = prob_nochange;
    pMat[1][2] = prob_change;
    pMat[1][3] = prob_change;
    pMat[2][0] = prob_change;
    pMat[2][1] = prob_change;
    pMat[2][2] = prob_nochange;
    pMat[2][3] = prob_change;
    pMat[3][0] = prob_change;
    pMat[3][1] = prob_change;
    pMat[3][2] = prob_change;
    pMat[3][3] = prob_nochange;
#   if defined DEBUGGING_OUTPUT3
        std::cerr << "from line " << __LINE__ << ":\n" ;  std::cerr << "edgeLength = " << edgeLength << " probs = " << prob_nochange << ", " << prob_change << "\n";
#   endif
}

void JCMulitCatTiMat(double edgeLength, TiMatVec pMatVec) {
    assert(gBlob);
#   if defined DEBUGGING_OUTPUT3
        std::cerr << "JCMulitCatTiMat edgeLength = " << edgeLength << '\n';
#   endif
    for (unsigned i = 0; i < gBlob->rates.size(); ++i) {
#       if defined DEBUGGING_OUTPUT3
            std::cerr << "JCMulitCatTiMat gBlob->rates[" << i << "] = " << gBlob->rates[i] << '\n';
#       endif
        JCTiMat(edgeLength*(gBlob->rates[i]), pMatVec[i]); //@TEMP JC
    }
}

void GenericMulitCatTiMat(double edgeLength, TiMatVec pMatVec) {
    assert(gBlob);
    int eigenIndex = 0;
    int numToCalc = gBlob->nRates;
    gBlob->scaledEdgeLengths.resize(gBlob->nRates);
    gBlob->modelIndexVector.resize(gBlob->nRates);
#   if defined DEBUGGING_OUTPUT
        std::cerr << "GenericMulitCatTiMat edgeLength = " << edgeLength << '\n';
#   endif
    for (int i = 0; i < gBlob->nRates; ++i) {
        gBlob->scaledEdgeLengths[i] = edgeLength*gBlob->rates[i];
#       if defined DEBUGGING_OUTPUT
            std::cerr << "GenericMulitCatTiMat rates[" << i << "] = " << gBlob->rates[i] << '\n';
#       endif
        gBlob->modelIndexVector[i] = i;
    }
    int rc = calcPrMatsForHandle(gBlob->likeCalcHandle, eigenIndex, numToCalc, &(gBlob->scaledEdgeLengths[0]), &(gBlob->modelIndexVector[0]));
    if (rc != 0) {
        throw NxsException("Call to calcPrMatsForHandle failed");
    }
    for (int i = 0; i < gBlob->nRates; ++i) {
        fetchPrMat(gBlob->likeCalcHandle, i, pMatVec[i][0]);
#       if defined DEBUGGING_OUTPUT
            std::cerr << "from line " << __LINE__ << ":\n" ;  std::cerr << "scaled edgeLength = " << gBlob->scaledEdgeLengths[i] << " probs = " << pMatVec[i][0][0] << ", " << pMatVec[i][0][1] << "\n";
#       endif
    }

}



class ProbForObsStateSet{ //for each state want to set -1 to 1 and all else to 0 (will be either 1 or anything up to NumStates)
    public:
        ProbForObsStateSet(unsigned int numStates) {
            std::vector<double> initialVal(numStates, 0.0);
            noRepeatedState = initialVal;
            probVec.assign(numStates, initialVal);
        }

        std::vector<double> & getProbForCommState(int commState) {
            if(commState == -1)
                return noRepeatedState;
            return probVec.at(commState);
        }
        void writeDebug(std::ostream & o, const CommonInfo & blob) const {
            o << "ProbForObsStateSet{\n  ";
            o << "-1 ";
            for (unsigned i = 0; i < noRepeatedState.size() ; ++i) {
                o << noRepeatedState[i] << " ";
            }
            for (unsigned i = 0; i < probVec.size() ; ++i) {
                o << "\n  " << i << ' ';
                const probvec_t & pvi = probVec[i];
                for (unsigned i = 0; i < pvi.size() ; ++i) {
                    o << pvi[i] << " ";
                }
            }
            o << "}\n";
        }

    private:
        typedef std::vector<double> probvec_t;
        std::vector<probvec_t> probVec;
        probvec_t noRepeatedState;
};


class NodeDataStructure{ //data members
    public:
        NodeDataStructure(unsigned int numStates) {
            numLeaves = 0;
            int len = 1 << numStates;
            ProbForObsStateSet dummy (numStates);
            probVec.assign(len, dummy);
            std::cerr << "NodeDataStructure ctor. Address = "<< long(this) << " len = " << len <<"\n";
        }

        ProbForObsStateSet & getForObsStateSet(int obs) {
            return probVec.at(obs);
        }

        int getNumLeaves() {
            return numLeaves;
            }

        void setNumLeaves(int n) {
            numLeaves = n;
        }

        void writeDebug(std::ostream & o, const CommonInfo & blob) const {
            o << "NodeDataStructure{ numLeaves = " << this->numLeaves << "\norobVec: ";
            for (unsigned i = 0; i < probVec.size() ; ++i) {
                o << i << ' ';
                const ProbForObsStateSet & pfoss = probVec[i];
                pfoss.writeDebug(o, blob);
                o << '\n';
            }
            o << "}\n";
        }
    private:
        std::vector<ProbForObsStateSet> probVec;
        int numLeaves;
};

double calculateTransProb(int ancIndex,
                          int i,
                          double edgeLen,
                          const CommonInfo & blob) {
    TiMatFunc fn = blob.tiMatFunc;
    fn(edgeLen, blob.firstMatVec.GetAlias());
    const double *** leftPMatVec = const_cast<const double ***>(blob.firstMatVec.GetAlias());
    //@TEMP no rate heterogeneity
    const double ** tiMatrix = leftPMatVec[0];
    return tiMatrix[ancIndex][i];
}
double calcProbOfSubtreeForObsStSetAndComm(NodeDataStructure * subtreeData,
                             int ancIndex,
                             int obsBits,
                             int commonStates,
                             double edgeLen,
                             const CommonInfo & blob){
    double p = 0.0;
    ProbForObsStateSet & childProbSet = subtreeData->getForObsStateSet(obsBits);
#   if defined DEBUGGING_OUTPUT
        childProbSet.writeDebug(std::cerr, blob);
#   endif
    std::vector<double> & childProb = childProbSet.getProbForCommState(commonStates);
    for(int i = 0; i<blob.nStates; i++) {
        double transProb = calculateTransProb(ancIndex, i, edgeLen, blob); //includes lib so Mark will write this func
        double partialLike = childProb[i];
        double x = transProb * partialLike;
#       if defined DEBUGGING_OUTPUT
            std::cerr << __LINE__ << " transProb = " << transProb << "    partialLike = " << partialLike << "\n";
#       endif
        p += x; // this loop is how we sum
    }
    return p;
}


double calcProbOfSubtreeForObsStSetNoRepeated(NodeDataStructure * subtreeData,
                             int ancIndex,
                             int obsBits,
                             double edgeLen,
                             const CommonInfo & blob){
    return calcProbOfSubtreeForObsStSetAndComm(subtreeData, ancIndex, obsBits, -1, edgeLen, blob);
}


void summarizeUninformativePatternClassProbabilities(NodeDataStructure * rootData,
                                                     std::ostream & out,
                                                     const CommonInfo & blob) {
    stateSetContainer::const_iterator ssCit = blob.stateSetBegin();
    for(; ssCit!=blob.stateSetEnd(); ++ssCit){
        const int & obsStSet = *ssCit; //'dereferencing' it
        int common = -1;
        int numObsSt = countBits(obsStSet);
#       if defined DEBUGGING_OUTPUT2
            std::cerr << "from line " << __LINE__ << ":  " ; std::cerr << "obsStSet =  " << obsStSet << " numObsSt = "<< numObsSt << "\n";
#       endif
        while(common>-2) /* or for(;;)*/ { // loop over common
            ProbForObsStateSet & currNdProbSet = rootData->getForObsStateSet(obsStSet);
            std::vector<double> & currNdProbVec = currNdProbSet.getProbForCommState(-1);
#           if defined DEBUGGING_OUTPUT2
                std::cerr << "from line " << __LINE__ << ":  " ; std::cerr << " common = "<< common << "\n";
#           endif
            if(common == -1) { //no comm state
                if(rootData->getNumLeaves()==numObsSt) {
                    double p = 0.0;
                    for(int anc = 0; anc < blob.nStates; anc++) {
                        p += currNdProbVec[anc];
                    }
                    out << "Prob(obs state set = " << obsStSet << " , no repeated states) = " << p << std::endl;
                }
            }
            common = getNextCommStSet(obsStSet, common);

        }
    }
}
#if 0
                                    vector<int> leftObsStSets = subsetsOfGivenSize(obsStSet, leftNodeData->getNumLeaves());
                                    for(int j = 0; j < leftObsStSets.size(); j++) {
                                        int leftObsStSet = leftObsStSets[j];
                                        int rightObsStSet = obsStSet - leftObsStSet;
                                        std::cerr << "leftObsStSet " << leftObsStSet << '\n';
                                        std::cerr << "rightObsStSet " << rightObsStSet << '\n';
                                        double leftProb, rightProb;
                                        double leftEdgeLen = leftChild->GetEdgeToParent().GetDblEdgeLen();
                                        if(leftNodeData->getNumLeaves() == 1) {
                                            leftProb = calculateTransProb(anc, convertBitToIndex(leftObsStSet), leftEdgeLen, blob);
                                        }
                                        else {
                                            leftProb = calcProbOfSubtreeForObsStSetNoRepeated(leftNodeData, anc, leftObsStSet, leftEdgeLen, blob);
                                        }
                                        double rightEdgeLen = rightChild->GetEdgeToParent().GetDblEdgeLen();

                                        if(rightNodeData->getNumLeaves() == 1) {
                                            rightProb = calculateTransProb(anc, convertBitToIndex(rightObsStSet), rightEdgeLen, blob);
                                        }
                                        else {
                                            rightProb = calcProbOfSubtreeForObsStSetNoRepeated(rightNodeData, anc, rightObsStSet, rightEdgeLen, blob);
                                        }
#                                       if defined DEBUGGING_OUTPUT2
                                            std::cerr << "from line " << __LINE__ << ":  " ; std::cerr << " leftProb = "<< leftProb << "rightProb = " << rightProb "\n";
#                                       endif;
                                        double jointNdProb = leftProb * rightProb;
                                        currNdProbVec[anc] += jointNdProb;
                                    }
                                }
                            }
                        }
                        else {
                            int commonBits = convertIndexToBit(common);
                            for(int anc = 0; anc < blob.nStates; anc++) {
                                currNdProbVec[anc] = 0.0;
                                int leftCommSt, rightCommSt;
                                leftCommSt = common;
                                rightCommSt = common;
                                #if defined DEBUGGING_OUTPUT2
                                    std::cerr << "from line " << __LINE__ << ":  " ; std::cerr << " leftCommSt = "<< common << "rightCommSt = " << common << "\n";
                                #endif
                                vector<int> obsStSetsWithComm = subsetsContainingGivenState(obsStSet, commonBits); //prob for both
                                for(int j = 0; j < obsStSetsWithComm.size(); j++) {
                                        int leftObsStSet = obsStSetsWithComm[j];
                                        int rightObsStSet = obsStSet - leftObsStSet + commonBits;
                                        std::cerr << "leftObsStSet " << leftObsStSet << '\n';
                                        std::cerr << "rightObsStSet " << rightObsStSet << '\n';
                                        double leftProb, rightProb;
                                        double leftEdgeLen = leftChild->GetEdgeToParent().GetDblEdgeLen();
                                        leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNodeData, anc, leftObsStSet, common, leftEdgeLen, blob);
                                        double rightEdgeLen = rightChild->GetEdgeToParent().GetDblEdgeLen();
                                        rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNodeData, anc, rightObsStSet, common, rightEdgeLen, blob);
                                        double jointNdProb = leftProb * rightProb;
                                        currNdProbVec[anc] += jointNdProb;
                                    }

                                leftCommSt = -1;
                                rightCommSt = common;
                                //add probability when only right common, left not repeated
                                for(int j = 0; j < obsStSetsWithComm.size(); j++) {
                                        int rightObsStSet = obsStSetsWithComm[j];
                                        int leftObsStSet = obsStSet - rightObsStSet + commonBits;
                                        std::cerr << "rightObsStSet " << rightObsStSet << '\n';
                                        std::cerr << "leftObsStSet " << leftObsStSet << '\n';
                                        double leftProb, rightProb;
                                        double leftEdgeLen = leftChild->GetEdgeToParent().GetDblEdgeLen();
                                        leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNodeData, anc, leftObsStSet, -1, leftEdgeLen, blob);
                                        double rightEdgeLen = rightChild->GetEdgeToParent().GetDblEdgeLen();
                                        rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNodeData, anc, rightObsStSet, common, rightEdgeLen, blob);
                                        double jointNdProb = leftProb * rightProb;
                                        currNdProbVec[anc] += jointNdProb;
                                        //Now consider when the left is not displayed by commonBits as observed States
                                        leftObsStSet = obsStSet - rightObsStSet;
                                        if(leftObsStSet != 0) {
                                            leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNodeData, anc, leftObsStSet, -1, leftEdgeLen, blob);
                                            jointNdProb = leftProb * rightProb;
                                            currNdProbVec[anc] += jointNdProb;
                                        }
                                    }


                                leftCommSt = common; //do this on my own (opp from above)
                                rightCommSt = -1;
                                //add probability when only left common
                                for(int j = 0; j < obsStSetsWithComm.size(); j++) {
                                        int rightObsStSet = obsStSetsWithComm[j];
                                        int leftObsStSet = obsStSet - rightObsStSet + commonBits;
                                        std::cerr << "rightObsStSet " << rightObsStSet << '\n';
                                        std::cerr << "leftObsStSet " << leftObsStSet << '\n';
                                        double leftProb, rightProb;
                                        double leftEdgeLen = leftChild->GetEdgeToParent().GetDblEdgeLen();
                                        leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNodeData, anc, leftObsStSet, common, leftEdgeLen, blob);
                                        double rightEdgeLen = rightChild->GetEdgeToParent().GetDblEdgeLen();
                                        rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNodeData, anc, rightObsStSet, -1, rightEdgeLen, blob);
                                        double jointNdProb = leftProb * rightProb;
                                        currNdProbVec[anc] += jointNdProb;
                                        //Now consider when the right is not displayed by commonBits as observed States
                                        rightObsStSet = obsStSet - leftObsStSet;
                                        if(rightObsStSet != 0) {
                                            rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNodeData, anc, rightObsStSet, -1, rightEdgeLen, blob);
                                            jointNdProb = leftProb * rightProb;
                                            currNdProbVec[anc] += jointNdProb;
                                        }
                                    }


                                leftCommSt = -1;
                                rightCommSt = -1;
                                //add probability when neither common
                                for(int j = 0; j < obsStSetsWithComm.size(); j++) {
                                        int rightObsStSet = obsStSetsWithComm[j];
                                        int leftObsStSet = obsStSet - rightObsStSet + commonBits;
                                        std::cerr << "rightObsStSet " << rightObsStSet << '\n';
                                        std::cerr << "leftObsStSet " << leftObsStSet << '\n';
                                        double leftProb, rightProb;
                                        double leftEdgeLen = leftChild->GetEdgeToParent().GetDblEdgeLen();
                                        leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNodeData, anc, leftObsStSet, -1, leftEdgeLen, blob);
                                        double rightEdgeLen = rightChild->GetEdgeToParent().GetDblEdgeLen();
                                        rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNodeData, anc, rightObsStSet, -1, rightEdgeLen, blob);
                                        double jointNdProb = leftProb * rightProb;
                                        currNdProbVec[anc] += jointNdProb;
                                }

                            }

                        }
                        common = getNextCommStSet(obsStSet, common);
                    }
                }
            }
        }
    }
    catch (...) {
        throw;
    }
    //std::vector<const NxsSimplenode.....
    //
    //if(tipProbInfo = constant)
    //else if(tipProbInfo changes) ....
    //return false;
    //
    //bool needToDelRootProbInfo = true;    <--should this be true since it's ininformative?
    //Do we not need 'tipProbInfo.createForTip(blob); (since we don't need to know the tip info?)
    //
    //
    //if (blob.isSymmetric) {
    //              currProbInfo->calculateSymmetric(*leftPI, leftNd->GetEdgeToParent().GetDblEdgeLen(),
    //                                      *rightPI, rightNd->GetEdgeToParent().GetDblEdgeLen(),
    //                                      tiMatFunc, blob);
    //          }
    //          else {
    //              currProbInfo->calculate(*leftPI, leftNd->GetEdgeToParent().GetDblEdgeLen(),
    //                                      *rightPI, rightNd->GetEdgeToParent().GetDblEdgeLen(),
    //                                      tiMatFunc, blob);
    //will the above stuff be the same?
    //
    //
    //need to define Z, t,c,a (t observed states) (c common state - for c=-1 we will have no repeated states)
    // for(c=-1)
    //      {nStates = -1}
    //
    //
}
#endif

 NodeDataStructure * calculateUninformativePatternClassProbabilities(const NxsSimpleTree & tree,
                                                     std::ostream & out,
                                                     const CommonInfo & blob) {
    cout << "blah\n";
    std::vector<const NxsSimpleNode *> preorderVec = tree.GetPreorderTraversal();
    std::map<const NxsSimpleNode *, NodeDataStructure *> node2dataMap;
#   if defined DEBUGGING_OUTPUT
        std::cerr << "blob.nStates = " << blob.nStates << '\n';
#   endif
    NodeDataStructure * currNdData = 0L;
    try {
        int ndInd = preorderVec.size() - 1;
        for (; ndInd >= 0; --ndInd) {
            const NxsSimpleNode * nd = preorderVec[ndInd];
            std::vector<NxsSimpleNode *> children = nd->GetChildren();
            const unsigned numChildren = children.size();
#           if defined DEBUGGING_OUTPUT
                std::cerr << "from line " << __LINE__ << ":  " ; std::cerr << "In calculateUninformativePatternClassProbabilities at node " << nd->GetTaxonIndex() << ", #children = " << numChildren << "\n";
#           endif
            NodeID currNdId(nd, 0);
            currNdData = new NodeDataStructure(blob.nStates); //curNdData = ancestor (in lower loop)
            node2dataMap[nd] = currNdData;
            if (numChildren == 0) {
                for(int i=0; i<blob.nStates; i++)
                {
                    int ss=1 << i;
                    ProbForObsStateSet & p = currNdData->getForObsStateSet(ss);
                    std::vector<double> & v = p.getProbForCommState(-1);
                    v[i] = 1.0;
                    currNdData->setNumLeaves(1);
                }
            }
            else {
                if (numChildren != 2)
                    throw NxsException("Trees must be of degree 2\n");
                NxsSimpleNode * leftChild = children[0];
                NodeDataStructure * leftNodeData = node2dataMap[leftChild];

                NxsSimpleNode * rightChild = children[1];
                NodeDataStructure * rightNodeData = node2dataMap[rightChild];
                currNdData->setNumLeaves(leftNodeData->getNumLeaves()+rightNodeData->getNumLeaves());
#               if defined DEBUGGING_OUTPUT
                     std::cerr << "from line " << __LINE__ << ":  " ; std::cerr << "currNdData->setNumLeaves set to " << leftNodeData->getNumLeaves()+rightNodeData->getNumLeaves() << "\n";
#               endif

                stateSetContainer::const_iterator ssCit = blob.stateSetBegin();
                for(; ssCit!=blob.stateSetEnd(); ++ssCit){
                    const int & obsStSet = *ssCit; //'dereferencing' it
                    int common = -1;
                    int numObsSt = countBits(obsStSet);
#                   if defined DEBUGGING_OUTPUT
                        std::cerr << "from line " << __LINE__ << ":  " ; std::cerr << "obsStSet =  " << obsStSet << " numObsSt = "<< numObsSt << "\n";
#                   endif
                    while(common>-2) /* or for(;;)*/ { // loop over common
                        ProbForObsStateSet & currNdProbSet = currNdData->getForObsStateSet(obsStSet);
                        std::vector<double> & currNdProbVec = currNdProbSet.getProbForCommState(-1);
#                       if defined DEBUGGING_OUTPUT
                                std::cerr << "from line " << __LINE__ << ":  " ; std::cerr << " common = "<< common << "\n";
#                       endif

                        if(common == -1) { //no comm state
 #                          if defined DEBUGGING_OUTPUT
                                std::cerr << __LINE__ <<  " currNdData->getNumLeaves() = " << currNdData->getNumLeaves() << "\n";
                                std::cerr << __LINE__ <<  " numObsSt = " << numObsSt << "\n";
 #                          endif
                            if(currNdData->getNumLeaves() == numObsSt) {
                                for(int anc = 0; anc < blob.nStates; anc++) {
                                    //vector<int> subsetsOfGivenSize(int, int);
                                    std::cerr << "ObsStSet " << obsStSet << '\n';
                                    currNdProbVec[anc] = 0.0;
                                    vector<int> leftObsStSets = subsetsOfGivenSize(obsStSet, leftNodeData->getNumLeaves());
                                    for(int j = 0; j < leftObsStSets.size(); j++) {
                                        int leftObsStSet = leftObsStSets[j];
                                        int rightObsStSet = obsStSet - leftObsStSet;
#                                       if defined DEBUGGING_OUTPUT2
                                            std::cerr << "leftObsStSet " << leftObsStSet << '\n';
                                            std::cerr << "rightObsStSet " << rightObsStSet << '\n';
#                                       endif
                                        double leftProb, rightProb;
                                        double leftEdgeLen = leftChild->GetEdgeToParent().GetDblEdgeLen();
                                        if(leftNodeData->getNumLeaves() == 1) {
                                            leftProb = calculateTransProb(anc, convertBitToIndex(leftObsStSet), leftEdgeLen, blob);
                                        }
                                        else {
                                            leftProb = calcProbOfSubtreeForObsStSetNoRepeated(leftNodeData, anc, leftObsStSet, leftEdgeLen, blob);
                                        }
                                        double rightEdgeLen = rightChild->GetEdgeToParent().GetDblEdgeLen();

                                        if(rightNodeData->getNumLeaves() == 1) {
                                            rightProb = calculateTransProb(anc, convertBitToIndex(rightObsStSet), rightEdgeLen, blob);
                                        }
                                        else {
                                            rightProb = calcProbOfSubtreeForObsStSetNoRepeated(rightNodeData, anc, rightObsStSet, rightEdgeLen, blob);
                                        }
                                        double jointNdProb = leftProb * rightProb;
                                        currNdProbVec[anc] += jointNdProb;
                                    }
                                }
                            }
                        }
                        else {
                            int commonBits = convertIndexToBit(common);
#                           if defined DEBUGGING_OUTPUT
                                std::cerr << "commonBits = " << commonBits << '\n';
#                           endif
                            for(int anc = 0; anc < blob.nStates; anc++) {
                                currNdProbVec[anc] = 0.0;
                                int leftCommSt, rightCommSt;
                                leftCommSt = common;
                                rightCommSt = common;
                                vector<int> obsStSetsWithComm = subsetsContainingGivenState(obsStSet, commonBits); //prob for both
#                               if defined DEBUGGING_OUTPUT
                                    std::cerr << "leftCommSt, rightCommSt " << leftCommSt << ',' << rightCommSt << '\n';
                                    std::cerr << "obsStSetsWithComm.size() = " << obsStSetsWithComm.size() << '\n';
#                               endif
                                for(int j = 0; j < obsStSetsWithComm.size(); j++) {
                                        int leftObsStSet = obsStSetsWithComm[j];
                                        int rightObsStSet = obsStSet - leftObsStSet + commonBits;
#                                       if defined DEBUGGING_OUTPUT
                                            std::cerr << "leftObsStSet " << leftObsStSet << '\n';
                                            std::cerr << "rightObsStSet " << rightObsStSet << '\n';
#                                       endif
                                        double leftProb, rightProb;
                                        double leftEdgeLen = leftChild->GetEdgeToParent().GetDblEdgeLen();
                                        leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNodeData, anc, leftObsStSet, common, leftEdgeLen, blob);
                                        double rightEdgeLen = rightChild->GetEdgeToParent().GetDblEdgeLen();
                                        rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNodeData, anc, rightObsStSet, common, rightEdgeLen, blob);
                                        double jointNdProb = leftProb * rightProb;
#                                       if defined DEBUGGING_OUTPUT
                                            std::cerr << "leftProb " << leftProb << '\n';
                                            std::cerr << "rightProb " << rightProb << '\n';
#                                       endif
                                        currNdProbVec[anc] += jointNdProb;
                                    }

                                leftCommSt = -1;
                                rightCommSt = common;
                                //add probability when only right common, left not repeated
                                for(int j = 0; j < obsStSetsWithComm.size(); j++) {
                                        int rightObsStSet = obsStSetsWithComm[j];
                                        int leftObsStSet = obsStSet - rightObsStSet + commonBits;
#                                       if defined DEBUGGING_OUTPUT
                                            std::cerr << "leftCommSt, rightCommSt" << leftCommSt << ',' << rightCommSt << '\n';
                                            std::cerr << "rightObsStSet " << rightObsStSet << '\n';
                                            std::cerr << "leftObsStSet " << leftObsStSet << '\n';
#                                       endif
                                        double leftProb, rightProb;
                                        double leftEdgeLen = leftChild->GetEdgeToParent().GetDblEdgeLen();
                                        leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNodeData, anc, leftObsStSet, -1, leftEdgeLen, blob);
                                        double rightEdgeLen = rightChild->GetEdgeToParent().GetDblEdgeLen();
                                        rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNodeData, anc, rightObsStSet, common, rightEdgeLen, blob);
                                        double jointNdProb = leftProb * rightProb;
                                        currNdProbVec[anc] += jointNdProb;
                                        //Now consider when the left is not displayed by commonBits as observed States
                                        leftObsStSet = obsStSet - rightObsStSet;
                                        if(leftObsStSet != 0) {
                                            leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNodeData, anc, leftObsStSet, -1, leftEdgeLen, blob);
                                            jointNdProb = leftProb * rightProb;
#                                       if defined DEBUGGING_OUTPUT
                                            std::cerr << "leftProb " << leftProb << '\n';
                                            std::cerr << "rightProb " << rightProb << '\n';
#                                       endif
                                            currNdProbVec[anc] += jointNdProb;
                                        }
                                    }


                                leftCommSt = common; //do this on my own (opp from above)
                                rightCommSt = -1;
                                //add probability when only left common
                                for(int j = 0; j < obsStSetsWithComm.size(); j++) {
                                        int rightObsStSet = obsStSetsWithComm[j];
                                        int leftObsStSet = obsStSet - rightObsStSet + commonBits;
#                                       if defined DEBUGGING_OUTPUT
                                            std::cerr << "leftCommSt, rightCommSt" << leftCommSt << ',' << rightCommSt << '\n';
                                            std::cerr << "rightObsStSet " << rightObsStSet << '\n';
                                            std::cerr << "leftObsSmatSet " << leftObsStSet << '\n';
#                                       endif
                                        double leftProb, rightProb;
                                        double leftEdgeLen = leftChild->GetEdgeToParent().GetDblEdgeLen();
                                        leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNodeData, anc, leftObsStSet, common, leftEdgeLen, blob);
                                        double rightEdgeLen = rightChild->GetEdgeToParent().GetDblEdgeLen();
                                        rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNodeData, anc, rightObsStSet, -1, rightEdgeLen, blob);
                                        double jointNdProb = leftProb * rightProb;
                                        currNdProbVec[anc] += jointNdProb;
                                        //Now consider when the right is not displayed by commonBits as observed States
                                        rightObsStSet = obsStSet - leftObsStSet;
                                        if(rightObsStSet != 0) {
                                            rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNodeData, anc, rightObsStSet, -1, rightEdgeLen, blob);
                                            jointNdProb = leftProb * rightProb;
#                                       if defined DEBUGGING_OUTPUT
                                            std::cerr << "leftProb " << leftProb << '\n';
                                            std::cerr << "rightProb " << rightProb << '\n';
#                                       endif
                                            currNdProbVec[anc] += jointNdProb;
                                        }
                                    }


                                leftCommSt = -1;
                                rightCommSt = -1;
                                //add probability when neither common
                                for(int j = 0; j < obsStSetsWithComm.size(); j++) {
                                        int rightObsStSet = obsStSetsWithComm[j];
                                        int leftObsStSet = obsStSet - rightObsStSet + commonBits;
#                                       if defined DEBUGGING_OUTPUT
                                            std::cerr << "leftCommSt, rightCommSt" << leftCommSt << ',' << rightCommSt << '\n';
                                            std::cerr << "rightObsStSet " << rightObsStSet << '\n';
                                            std::cerr << "leftObsStSet " << leftObsStSet << '\n';
#                                       endif
                                        double leftProb, rightProb;
                                        double leftEdgeLen = leftChild->GetEdgeToParent().GetDblEdgeLen();
                                        leftProb = calcProbOfSubtreeForObsStSetAndComm(leftNodeData, anc, leftObsStSet, -1, leftEdgeLen, blob);
                                        double rightEdgeLen = rightChild->GetEdgeToParent().GetDblEdgeLen();
                                        rightProb = calcProbOfSubtreeForObsStSetAndComm(rightNodeData, anc, rightObsStSet, -1, rightEdgeLen, blob);
                                        double jointNdProb = leftProb * rightProb;
#                                       if defined DEBUGGING_OUTPUT
                                            std::cerr << "leftProb " << leftProb << '\n';
                                            std::cerr << "rightProb " << rightProb << '\n';
#                                       endif
                                        currNdProbVec[anc] += jointNdProb;
                                }
#                               if defined DEBUGGING_OUTPUT
                                    std::cerr << "currNdProbVec[" << anc << "] " << currNdProbVec[anc] << '\n';
#                               endif
                            }
                        }
#                       if defined DEBUGGING_OUTPUT
                            //std::cerr << "EARLY EXIT!!!\n";
                            //exit(1);
#                       endif
                        common = getNextCommStSet(obsStSet, common);
                    }
                }
            }
        }
    }
    catch (...) {
        throw;
    }
    return currNdData;
    //std::vector<const NxsSimplenode.....
    //
    //if(tipProbInfo = constant)
    //else if(tipProbInfo changes) ....
    //return false;
    //
    //bool needToDelRootProbInfo = true;    <--should this be true since it's ininformative?
    //Do we not need 'tipProbInfo.createForTip(blob); (since we don't need to know the tip info?)
    //
    //
    //if (blob.isSymmetric) {
    //              currProbInfo->calculateSymmetric(*leftPI, leftNd->GetEdgeToParent().GetDblEdgeLen(),
    //                                      *rightPI, rightNd->GetEdgeToParent().GetDblEdgeLen(),
    //                                      tiMatFunc, blob);
    //          }
    //          else {
    //              currProbInfo->calculate(*leftPI, leftNd->GetEdgeToParent().GetDblEdgeLen(),
    //                                      *rightPI, rightNd->GetEdgeToParent().GetDblEdgeLen(),
    //                                      tiMatFunc, blob);
    //will the above stuff be the same?
    //
    //
    //need to define Z, t,c,a (t observed states) (c common state - for c=-1 we will have no repeated states)
    // for(c=-1)
    //      {nStates = -1}
    //
    //
}

/// \ NxsException for gaps or ambiguous cells
/// \returns a string of the symbols for each state (length == nStates).
std::string convertToBitFieldMatrix(const NxsCharactersBlock & cb, BitFieldMatrix & mat, const NxsUnsignedSet * toInclude) {
    std::vector<const NxsDiscreteDatatypeMapper *> mappers = cb.GetAllDatatypeMappers();
    if (mappers.empty() || mappers[0] == NULL)
        throw NxsException("no mappers");

    if (mappers.size() != 1)
        throw NxsException("Expecting an unmixed characters block, but found a matrix with datatype = mixed or a datatype with augmented symbols\n");
    NxsUnsignedSet scratchSet;
    if (toInclude == 0L) {
        for (unsigned i = 0; i < cb.GetNChar(); ++i)
            scratchSet.insert(i);
        toInclude = & scratchSet;
     }

    std::set <const NxsDiscreteDatatypeMapper * > usedMappers;
    for (NxsUnsignedSet::const_iterator indIt = toInclude->begin(); indIt != toInclude->end(); ++indIt) {
        unsigned charIndex = *indIt;
        usedMappers.insert(cb.GetDatatypeMapperForChar(charIndex));
    }

    if (usedMappers.size() > 1)
        throw NxsException("too many mappers");
    if (usedMappers.empty())
        throw NxsException("no mappers - or empty charset");


    const NxsDiscreteDatatypeMapper & mapper = **usedMappers.begin();

    NxsCharactersBlock::DataTypesEnum inDatatype = mapper.GetDatatype();
    const unsigned nStates =  mapper.GetNumStates();
    if (nStates > MAX_NUM_STATES)
        throw NxsException("MAX_NUM_STATES exceeded.  Recompile with a larger datatype for BitField.\n");

    for (NxsDiscreteStateCell i = 0; i < (NxsDiscreteStateCell)nStates; ++i) {
        if (mapper.GetStateSetForCode(i).size() != 1)
            throw NxsException("Expecting the first states to correspond to state sets with size 1");
    }
    const std::string fundamentalSymbols = mapper.GetSymbols();
    if (fundamentalSymbols.length() != nStates) {
        errormsg << "Expecting the fundamental symbols (" << fundamentalSymbols << ") to have length " << nStates;
        throw NxsException(errormsg);
    }
    NCL_ASSERT((int)NXS_MISSING_CODE < 0);
    NCL_ASSERT((int)NXS_GAP_STATE_CODE < 0);


    const unsigned nTaxa = cb.GetNTax();
    const unsigned includedNChar = toInclude->size();
    mat.resize(nTaxa);
    for (unsigned i = 0; i < nTaxa; ++i) {
        BitFieldRow & bfRow = mat[i];
        bfRow.resize(includedNChar);
        const NxsDiscreteStateRow & row = cb.GetDiscreteMatrixRow(i);
        if (row.empty()) {
            errormsg << "Empty row encountered for taxon " << (1+i) << ". Missing data is not supported.\n";
            throw NxsException(errormsg);
        }
        unsigned j = 0;
        for (NxsUnsignedSet::const_iterator tIncIt = toInclude->begin(); tIncIt != toInclude->end(); ++tIncIt, ++j) {
            const NxsDiscreteStateCell & cell = row.at(*tIncIt);
            if (cell < 0 || cell >= (NxsDiscreteStateCell) nStates) {
                errormsg << "Ambiguous/missing/gap data found for taxon " << (1+i) << " at site " << *tIncIt << ". Missing data is not supported.\n";
                throw NxsException(errormsg);
            }
            const int bfi = 1 << (int) cell;
            bfRow[j] = BitField(bfi);
        }
    }
    return fundamentalSymbols;
}

void writeBitFieldMatrix(std::ostream & out, const BitFieldMatrix & bitFieldMatrix) {
    for (unsigned i = 0; i < bitFieldMatrix.size(); ++i) {
        out << "taxon " << (i + 1) << ":\n";
        const BitFieldRow & row = bitFieldMatrix[i];
        for (unsigned j = 0; j < row.size(); ++j)
            out << ' ' << (int) row[j];
        out << '\n';
    }
}



void CommonInfo::calcEigenSolution() {
    const double * stateFreqArray = &(this->stateFreqVector[0]);
    const double ** relRateMatAlias = const_cast<const double **>(this->relRateMat.GetAlias());
    ScopedDblTwoDMatrix qMat(this->nStates, this->nStates);
    double ** qMatAlias = qMat.GetAlias();
    double diagWeightedSum = 0.0;
    for (unsigned a = 0; a < this->nStates; ++a) {
        double offDiagSum = 0.0;
        for (unsigned d = 0; d < this->nStates; ++d) {
            const double unNormRate = stateFreqArray[d]*relRateMatAlias[a][d];
            qMatAlias[a][d] = unNormRate;
            if (a != d)
                offDiagSum += unNormRate;
        }
        qMatAlias[a][a] = -offDiagSum;
        diagWeightedSum += stateFreqArray[a]*offDiagSum;
    }

    for (unsigned a = 0; a < this->nStates; ++a) {
        for (unsigned d = 0; d < this->nStates; ++d) {
            qMatAlias[a][d] /= diagWeightedSum;
        }
    }

    setQMatForHandle(this->likeCalcHandle, 0, const_cast<const double **>(qMatAlias));

    int rc = calc_eigen_mat(dsct_model_obj, 0);
    if (rc == 0) {
        errormsg << "recalc_eigen_mat failed";
        throw NxsException(errormsg);
    }
}

std::vector<double> toVecDouble(const std::string &s, const char * optName) {
    try {
        size_t begInd = 0;
        std::vector<double> v;
        for (;;) {
            size_t endInd = s.find(',', begInd);
            if (endInd == std::string::npos) {
                const NxsString n(s.substr(begInd).c_str());
                v.push_back(n.ConvertToDouble());
                return v;
            }
            else {
                const NxsString n(s.substr(begInd, endInd).c_str());
                v.push_back(n.ConvertToDouble());
                begInd = endInd + 1;
                if (begInd == s.length()) {
                    errormsg << "Error parsing " << optName << ": not expecting a series of numbers to end with a comma.";
                    throw NxsException(errormsg);
                }
            }
        }
    }
    catch (const NxsString::NxsX_NotANumber &) {
        errormsg << "Error parsing " << optName << ": expecting a series of comma-seperated numbers.";
        throw NxsException(errormsg);
    }
}

void CommonInfo::readModel(const std::vector<std::string> &optVec) {
    this->relRateMat.Initialize(this->nStates, this->nStates);
    for (unsigned a = 0; a < this->nStates; ++a) {
        for (unsigned d = 0; d < this->nStates; ++d) {
            this->relRateMat.GetAlias()[a][d] = 1.0;
        }
    }
    this->stateFreqVector.assign(this->nStates, 1.0/this->nStates);
    for (std::vector<std::string>::const_iterator ovIt = optVec.begin(); ovIt != optVec.end(); ++ovIt) {
        const std::string opt = *ovIt;
        assert(opt[0] == '-');
        const std::string flagWithDash = opt.substr(0,2);
        const char flag = opt[1];
        const std::string val(opt.c_str() + 2);
        const bool hasValues = (flag != 's');
        std::vector<double> v;
        if (hasValues) {
            v = toVecDouble(val, flagWithDash.c_str());
        }
        double x, y;
        int el = 0;
        const int numRelRates = (((this->nStates - 1)*(this->nStates))/2) - 1;
#       if defined DEBUGGING_OUTPUT
            std::cerr << "flagWithDash = " << flagWithDash << "   val = " << val << '\n';
            std::cerr << "parsedToDoubles v =";
            for (std::vector<double>::const_iterator vIt = v.begin(); vIt != v.end(); ++vIt) {
                std::cerr << ' ' << *vIt;
            }
            std::cerr << '\n';
#       endif
        switch (flag) {
            case 'f' :
                if (v.size() != this->nStates - 1) {
                    errormsg << "Expecting " << this->nStates - 1 << " state frequencies (the last will be obtained by subtraction)\n";
                    throw NxsException(errormsg);
                }
                x = 1.0;
                for (unsigned i = 0; i < this->nStates - 1; ++i) {
                    if (v[i] <= 0.0 || v[i] >= 1.0) {
                        throw NxsException("State frequencies must be between 0 and 1");
                    }
                    x -= v[i];
                }
                if (x < 0.0) {
                        throw NxsException("The sum of the state frequencies must be less than 1");
                }

                v.push_back(x);
                this->stateFreqVector = v;
                break;
            case 'r' :
                if (v.size() != numRelRates) {
                    errormsg << "Expecting " << numRelRates << " relative rates frequenies (the last will be set to 1.0)\n";
                    throw NxsException(errormsg);
                }
                v.push_back(1.0);
                el = 0;
                for (unsigned a = 0; a < this->nStates; ++a) {
                    for (unsigned d = a + 1; d < this->nStates; ++d) {
                        if (v[el] <= 0.0) {
                            throw NxsException("relative rates must be positive");
                        }
                        this->relRateMat.GetAlias()[a][d] = v[el];
                        this->relRateMat.GetAlias()[d][a] = v[el];
                        el++;
                    }
                }
                assert(el == 1 + numRelRates);
                break;
            case 'm' :
                this->rates = v;
                this->nRates = this->rates.size();

                for (unsigned i = 0; i < this->nRates; ++i) {
                    if (v[i] < 0.0) {
                        throw NxsException("rate multipliers must be positive");
                    }
                }
                break;
            case 'p' :
                this->rateProb = v;
                x = 1.0;
                for (unsigned i = 0; i < v.size(); ++i) {
                    if (v[i] <= 0.0 || v[i] >= 1.0) {
                        throw NxsException("rate category probabilities must be between 0 and 1");
                    }
                    x -= v[i];
                }
                if (x < 0.0) {
                    throw NxsException("The sum of rate category probabilities must be less than 1");
                }
                this->rateProb.push_back(x);
                break;
            case 's' :
                this->isSymmetric = true;
                break;
        }
    }

    this->nRates = this->rates.size();
    if (this->nRates != this->rateProb.size()) {
        throw NxsException("The number of rates and the number of rate category probabilities must agree");
    }

    this->pVecLen = this->nRates*this->nStates;
    this->categStateProb.resize(this->pVecLen) ;

    for (unsigned i = 0; i < this->nRates; ++i) {
        const double rp = this->rateProb[i];
        for (unsigned j = 0; j < this->nStates; ++j) {
            this->categStateProb[i*this->nStates + j] = rp*this->stateFreqVector[j];
        }
    }

    this->writeModel(std::cerr);


    this->firstMatVec.Free();
    this->firstMatVec.Initialize(this->nRates, this->nStates, this->nStates);
    this->secondMatVec.Free();
    this->secondMatVec.Initialize(this->nRates, this->nStates, this->nStates);



    int numLeaves = 2; // we don't need any, but I'm afraid of using 0 as an arg to beagle
    long numPatterns = 1 ; // we don't need any, but I'm afraid of using 0 as an arg to beagle
    const double * patternWeights=0L;
    int numStateCodeArrays = 2; // we don't need any, but I'm afraid of using 0 as an arg to beagle
    int numPartialStructs = 1;  // we don't need any, but I'm afraid of using 0 as an arg to beagle
    int numInstRateModels = 1 ; // we need a single q-matrix
    const ASRVObj ** asrvObjectArray=0L;  // we won't use pytbeaglehon's asrv
    int numProbMats= 2*(this->nRates); // we need a left and right for each rate category
    int numEigenStorage=1;
    int numRescalingsMultipliers=0;
    int resourceIndex=-1;
    long resourcePref=64;
    long resourceReq=64;
    int eigenIndex, numToCalc ;
    this->dsct_model_obj = 0L;
    this->likeCalcHandle = createLikelihoodCalcInstance(numLeaves,
                                                        numPatterns,
                                                        patternWeights,
                                                        this->nStates,
                                                        numStateCodeArrays,
                                                        numPartialStructs,
                                                        numInstRateModels,
                                                        asrvObjectArray,
                                                        numProbMats,
                                                        numEigenStorage,
                                                        numRescalingsMultipliers,
                                                        resourceIndex,
                                                        resourcePref,
                                                        resourceReq);

    if (this->likeCalcHandle < 0) {
        throw NxsException("Could not initialize a LikelihoodCalcInstance!");
    }
    const struct LikeCalculatorInstance * LCI = getLikeCalculatorInstance(this->likeCalcHandle);
    if (LCI == 0L) {
        throw NxsException("Call to getLikeCalculatorInstance failed");
    }
    this->dsct_model_obj = LCI->probModelArray[0];


    calcEigenSolution();

}

void CommonInfo::writeModel(std::ostream & out) const {
    int a, d;
    out << "Model state_freqs = (";
    for (unsigned i = 0; i < this->nStates; ++i) {
        out << ' ' << this->stateFreqVector.at(i);
    }
    out << ")\n";
    out << "Model rel_rates = (";
    for (a = 0; a < this->nStates; ++a) {
        out << '(';
        for (d = a + 1; d < this->nStates; ++d)
            out << ' ' << this->relRateMat.GetAlias()[a][d];
        out << ") ";
    }
    out << ")\n";
    for (a = 0; a < this->nRates; ++a) {
        out << "Model rate_multiplier_category" << a << " = (rate=" << this->rates.at(a) << ", prob=" << this->rateProb.at(a) << ")\n";
    }

    out << "Model_arg -f";
    for (unsigned i = 0; i < this->nStates; ++i) {
        if (i > 0)
            out << ',';
        out<< this->stateFreqVector.at(i);
    }
    out << "\nModel_arg -r";
    for (a = 0; a < this->nStates - 1; ++a) {
        for (d = a + 1; d < this->nStates; ++d) {
            if (d > 1)
                out << ',';
            out << this->relRateMat.GetAlias()[a][d];
        }
    }
    out << "\nModel_arg -m";
    for (a = 0; a < this->nRates; ++a) {
        if (a > 0)
            out << ',';
        out << this->rates.at(a);
    }
    out << "\nModel_arg -p";
    for (a = 0; a < this->nRates; ++a) {
        if (a > 0)
            out << ',';
        out << this->rateProb.at(a);
    }
    if (this->isSymmetric) {
        out << "\nModel_arg -s";
    }
    else {
        out << "\nModel not fully symmetric (not an Mk model)";
    }
    out << "\n";

    out << "[PAUP] begin paup;    lset nst = 6 BaseFreq = (";
    for (unsigned i = 0; i < this->nStates - 1; ++i) {
        out << ' ' << this->stateFreqVector.at(i);
    }
    out << ") rMat = (";
    for (a = 0; a < this->nStates - 2; ++a) {
        for (d = a + 1; d < this->nStates; ++d) {
            out << ' ' << this->relRateMat.GetAlias()[a][d];
        }
    }
    out << ");  end;\n";
}


void printHelp(std::ostream & o) {
    o << "Takes the path to a NEXUS file with a single characters block as an argument\n";
    o << "\nOptions:\n";
    o << "   -f#,#,#        state frequencies (last is calculated by subtraction).\n";
    o << "   -r#,#,#,#,#    relative rates (the last is defined to be 1.0).\n";
    o << "   -m#,#,#,...    rate multipliers (one for each rate category).\n";
    o << "   -p#,#,#,...    rate category probabilities (one fewer than the number of rates - the last is calculated by subtraction).\n";
    o << "   -s             use algorithms specifically designed for the Mk-type models\n";
}


int main(int argc, char * argv[]) {
    /* Usually it is easiest to surround interactions with NCL in a try block
        that catches NxsException instances.
        Reading files can certainly generate these errors, but even queries
        after the parse can result in these exceptions.
    */
    std::cout << "argc = " << argc << std::endl;
    std::cout << "argv[0] = " << argv[0] << std::endl;
    std::cout << "argv[1] = " << argv[1] << std::endl;
    std::string filename;
    std::vector<std::string> optVec;
    for (int argi = 1; argi < argc; ++argi) {
        if (strlen(argv[argi]) > 0) {
            std::string argstr(argv[argi]);
            if (argstr[0] == '-' && strlen(argv[argi]) > 1) {
                if (argstr[1] == 'h' ) {
                    printHelp(std::cerr);
                    return 0;
                }
                if (argstr[1] == 's' ) {
                    optVec.push_back(argstr);
                }
                else if (strlen(argv[argi]) > 2)
                    optVec.push_back(argstr);
                else {
                    std::cerr << "Expecting a numerical argument after " << argstr << '\n';
                    return 1;
                }

            }
            else if (filename.empty()) {
                filename = argstr;
            }
            else {
                std::cerr << "Expecting only one filename argument\n";
                return 1;
            }
        }
    }




    if (filename.empty()) {
        std::cerr << "Expecting a filename as an argument.\n";
        return 1;
    }


    try {
        int blocksToRead =  (PublicNexusReader::NEXUS_TAXA_BLOCK_BIT
                            | PublicNexusReader::NEXUS_CHARACTERS_BLOCK_BIT
                            | PublicNexusReader::NEXUS_ASSUMPTIONS_BLOCK_BIT
                            | PublicNexusReader::NEXUS_SETS_BLOCK_BIT
                            | PublicNexusReader::NEXUS_TREES_BLOCK_BIT
                            );
        MultiFormatReader nexusReader(blocksToRead, NxsReader::WARNINGS_TO_STDERR);

        std::cerr << "Reading " << filename << "\n";
        try {
            nexusReader.ReadFilepath(filename.c_str(), MultiFormatReader::NEXUS_FORMAT);
        }
        catch(const NxsException &x) {
            std::cerr << "Error:\n " << x.msg << std::endl;
            if (x.line > 0 || x.pos > 0)
                std::cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << std::endl;
            return 2;
        }
        catch(...) {
            nexusReader.DeleteBlocksFromFactories();
            std::cerr << "Exiting with an unknown error" << std::endl;
            return 1;
        }
        const unsigned numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
        if (numTaxaBlocks != 1) {
            std::cerr << "Expecting a file with exactly 1 TAXA block, but found " << numTaxaBlocks << " in the file " << filename << ".\n";
            return 2;
        }
        NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(0);
        const unsigned nCharsBlocks = nexusReader.GetNumCharactersBlocks(taxaBlock);
        if (nCharsBlocks != 1) {
            std::cerr << "Expecting a file with exactly 1 CHARACTERS/DATA block, but found " << nCharsBlocks << " in the file " << filename << ".\n";
            return 3;
        }
        NCL_COULD_BE_CONST  NxsCharactersBlock * charsBlock = nexusReader.GetCharactersBlock(taxaBlock, 0);
        assert(charsBlock);

        const NxsTransformationManager &tm = charsBlock->GetNxsTransformationManagerRef();
        std::vector<int> patternWeights = tm.GetDefaultIntWeights();
        const int * pwPtr = (patternWeights.empty() ? 0L : &patternWeights[0]);
        BitFieldMatrix bitFieldMatrix;
        std::string symbols = convertToBitFieldMatrix(*charsBlock, bitFieldMatrix);
        CommonInfo blob;
        blob.initialize(symbols);

        blob.readModel(optVec);


        blob.zeroVec.assign(bitFieldMatrix[0].size(), 0);
        gBlob = &blob;
        std::cerr << "from line " << __LINE__ << ":\n" ; writeBitFieldMatrix(std::cerr, bitFieldMatrix);


        const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(taxaBlock);
        if (nTreesBlocks == 0) {
            std::cerr << "Expecting a file with at least 1 TREES block.\n";
            return 3;
        }
        unsigned totalTreeIndex = 0;
        for (unsigned treesBlockInd = 0; treesBlockInd < nTreesBlocks; ++treesBlockInd) {
            NCL_COULD_BE_CONST  NxsTreesBlock * treesBlock = nexusReader.GetTreesBlock(taxaBlock, treesBlockInd);
            assert(treesBlock);
            treesBlock->ProcessAllTrees();
            const unsigned nTreesThisBlock = treesBlock->GetNumTrees();
            for (unsigned treeInd = 0; treeInd < nTreesThisBlock; ++treeInd) {
                const NxsFullTreeDescription & ftd = treesBlock->GetFullTreeDescription(treeInd);
                if (ftd.AllEdgesHaveLengths()) {
                    NxsSimpleTree nclTree(ftd, 0, 0.0);
                    blob.tiMatFunc = JCMulitCatTiMat; //@TEMP JC GenericMulitCatTiMat;
                    NodeDataStructure * rootData =  calculateUninformativePatternClassProbabilities(nclTree,
                                                                                                    std::cout,
                                                                                                    blob);
                    summarizeUninformativePatternClassProbabilities(rootData, std::cout, blob);

/*
                    if (false) {
                        PatternSummary observed;
                        classifyObservedDataIntoClasses(nclTree, bitFieldMatrix, pwPtr, std::cout, &observed, blob);
                    }
*/
                }
                else {
                    std::cerr << "Tree " << (1 + treeInd) << " of TREES block " << (1 + treesBlockInd) << " does not lengths for all of the edges. Skipping this tree.\n";
                }
            }
        }
        nexusReader.DeleteBlocksFromFactories();
    }
    catch(const NxsException &x) {
        std::cerr << "Error:\n " << x.msg << std::endl;
        if (x.line > 0 || x.pos > 0)
            std::cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << std::endl;
        return 2;
    }
    catch(...) {
        std::cerr << "Exiting with an unknown error" << std::endl;
        return 1;
    }
    return 0;
}

