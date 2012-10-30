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

using namespace std;

void freeProbInfo(const std::vector<const NxsSimpleNode *> & preorderVec, NodeIDToProbInfo & nodeIDToProbInfo);


// Globals (all should start with g[A-Z] pattern to make it easier to find and replace them later).
const unsigned MAX_NUM_STATES = 8*sizeof(BitField);


struct CommonInfo {
		bool isSymmetric;
		unsigned nStates;
		unsigned nRates;
		unsigned pVecLen; // nStates * nRates;
		BitField lastBitField;

		std::string alphabet;
		std::vector<BitField> singleStateCodes;
		// stateIndexToStateCode is probably always going to be identical to
		//		singleStateCodes, but use stateIndexToStateCode when you need the
		//		state codes for the fundamental states in order.
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

		DSCTModelObj * dsct_model_obj;

		std::vector<double> scaledEdgeLengths;
		std::vector<int> modelIndexVector;
		long likeCalcHandle;


		void writeModel(std::ostream & out) const;
	private:
		void calcEigenSolution();

		std::vector<unsigned> stateCodeToNumStates;
		std::map<BitField, std::string> stateCodesToSymbols;

		ScopedDblTwoDMatrix relRateMat;
		std::vector<double> stateFreqVector;
};

CommonInfo * gBlob = 0L;



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



#	if defined DEBUGGING_OUTPUT
		std::cerr << "from line " << __LINE__ << ":\n";
#	endif
	for (sc = 1;;++sc) {
#		if defined DEBUGGING_OUTPUT
			std::cerr << "Combos with unions that lead to " << this->toSymbol(sc) << ":";
#		endif
		const VecMaskPair & forUnions = this->pairsForUnionForEachDownPass[sc];
#		if defined DEBUGGING_OUTPUT
			for (VecMaskPair::const_iterator fuIt = forUnions.begin(); fuIt != forUnions.end(); ++fuIt) {
				std::cerr << " (\"" << this->toSymbol(fuIt->first);
				std::cerr << "\", \"" << this->toSymbol(fuIt->second);
				std::cerr << "\")	";
			}
			std::cerr << "\n";

			std::cerr << "Combos with intersections that lead to " << this->toSymbol(sc) << ":";
#		endif
		const VecMaskPair & forIntersections = this->pairsForIntersectionForEachDownPass[sc];
#		if defined DEBUGGING_OUTPUT
			for (VecMaskPair::const_iterator fuIt = forIntersections.begin(); fuIt != forIntersections.end(); ++fuIt) {
				std::cerr << " (\"" << this->toSymbol(fuIt->first);
				std::cerr << "\", \"" << this->toSymbol(fuIt->second);
				std::cerr << "\")	";
			}
			std::cerr << "\n";
#		endif

		if (sc == this->lastBitField)
			break;
	}


}


void ParsInfo::calculateForTip(const BitFieldRow & data, const CommonInfo & blob) {
	assert(blob.zeroVec.size() >= data.size());
	this->numPatterns = data.size();
	this->downPass = &data[0];
	this->allSeen = &data[0];
	this->score = &blob.zeroVec[0];
}
void ParsInfo::calculateForInternal(const ParsInfo & leftData, const ParsInfo & rightData) {
	this->numPatterns = leftData.size();
	assert(numPatterns == rightData.size());
	this->downPassOwned.resize(numPatterns);
	this->downPass = &this->downPassOwned[0];
	this->allSeenOwned.resize(numPatterns);
	this->allSeen = &this->allSeenOwned[0];
	this->scoreOwned.resize(numPatterns);
	this->score = &this->scoreOwned[0];
	const BitField * leftDown = leftData.downPass;
	const BitField * leftAll = leftData.allSeen;
	const unsigned * leftScore = leftData.score;
	const BitField * rightDown = rightData.downPass;
	const BitField * rightAll = rightData.allSeen;
	const unsigned * rightScore = rightData.score;
	for (unsigned i = 0; i < numPatterns; ++i, ++leftDown, ++leftAll, ++leftScore, ++rightDown, ++rightAll, ++rightScore) {
		const BitField lrIntersection = (*leftDown) & (*rightDown);
		const unsigned accumulScore = *leftScore + *rightScore;
		this->allSeenOwned[i] = (*leftAll) | (*rightAll);
		if (lrIntersection == BitField(0)) {
			const BitField lrUnion = (*leftDown) | (*rightDown);
			this->downPassOwned[i] = lrUnion;
			this->scoreOwned[i] = 1 + accumulScore;
		}
		else {
			this->downPassOwned[i] = lrIntersection;
			this->scoreOwned[i] = accumulScore;
		}
	}
}

void ParsInfo::write(std::ostream & o) const {
	unsigned totalScore = 0;
	for (unsigned i = 0; i < numPatterns; ++i) {
		o << "pattern = " << i << " score = " << (int)score[i] << " downpass = " << (int)downPass[i] << " states = " << (int)allSeen[i] << '\n';
		totalScore += score[i];
	}
	o << "totalScore = " << totalScore << '\n';
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

void ProbInfo::calculateSymmetric(const ProbInfo & leftPI, double leftEdgeLen,
					   const ProbInfo & rightPI, double rightEdgeLen,
					   TiMatFunc fn,
					   const CommonInfo & blob) {
#	if defined DEBUGGING_OUTPUT
		std::cerr << "from line " << __LINE__ << ":\n" ; std::cerr << "calculateSymmetric \n" ;
#	endif

	// fn will fill the prob matrix (blob.firstMatVec) for the appropriate branch length
	fn(leftEdgeLen, blob.firstMatVec.GetAlias());
	const double *** leftPMatVec = const_cast<const double ***>(blob.firstMatVec.GetAlias());
	fn(rightEdgeLen, blob.secondMatVec.GetAlias());
	const double *** rightPMatVec = const_cast<const double ***>(blob.secondMatVec.GetAlias());

	// now leftPMatVec and rightPMatVec correctly populated with transition probabilities


	const unsigned int leftMaxP = leftPI.getMaxParsScore();
	const unsigned int rightMaxP = rightPI.getMaxParsScore();
	const unsigned int maxParsScore = 1 + leftMaxP + rightMaxP;

	this->byParsScore.clear();
	this->byParsScore.resize(maxParsScore + 1); // add one to account zero

	this->nLeavesBelow = leftPI.getNLeavesBelow() + rightPI.getNLeavesBelow();

	unsigned obsMaxParsScore = 0;
	// do the calculations for staying in the constant patterns, these are more simple than the general calcs...
	if (true) { //@TEMP if true so that variables are scoped.
		const ProbForParsScore & leftFPS = leftPI.getByParsScore(0);
		const ProbForParsScore & rightFPS = rightPI.getByParsScore(0);
		ProbForParsScore & forZeroSteps = this->byParsScore[0];

		unsigned stateIndex = 0;
		const BitField sc = 1;

		const std::vector<double> * leftProbs = leftFPS.getProbsForDownPassAndObsMask(sc, sc);
		assert(leftProbs != 0L);
		const std::vector<double> * rightProbs = rightFPS.getProbsForDownPassAndObsMask(sc, sc);
		assert(rightProbs != 0L);

		std::vector<double> & toFillVec = forZeroSteps.byDownPass[sc][sc];
		toFillVec.assign(blob.nRates*blob.nStates, 0.0);

		addToAncProbVec(toFillVec, leftPMatVec, leftProbs, rightPMatVec, rightProbs, blob);
	}

	unsigned currScore = 1;
	std::cerr << "from line: " << __LINE__<< ": currScore = " << currScore << ".\n";
	bool scObserved = true;
	ProbForParsScore & forCurrScore = this->byParsScore[currScore];
	// if the ancestor has state set of {0}, and currScore = 1 then
	// 	one child must display a single state, and the other must have 2 states observed
	BitField downPass = 1;
	if (true) { // braces for reducing the scoping of variables, only
		unsigned int leftAccum = 0;
		const ProbForParsScore * leftFPS = &(leftPI.getByParsScore(leftAccum));
		unsigned int rightAccum = 1; // one step on the right side
		const ProbForParsScore * rightFPS = &(rightPI.getByParsScore(rightAccum));

		for (int leftright = 0; leftright < 2; ++leftright) {
			BitField leftDown = 1; // set {0}  is 1 in our bitfield notation
			const MaskToProbsByState * leftM2PBS = leftFPS->getMapPtrForDownPass(leftDown);
			assert(leftM2PBS != 0L);
			const BitField leftAllStates = 1; // set {0}  is 1 in our bitfield notation
			const std::vector<double> * leftProbs = getProbsForStatesMask(leftM2PBS, leftAllStates);
			assert(leftProbs != 0L);
			const BitField rightAllStates = (1|2); // set {0}  is 1 in our bitfield notation
			BitField rightDownArray[] = {1, 1|2};
			for (int rdi = 0 ; rdi < 2; ++rdi) {
				BitField rightDown = rightDownArray[rdi];
				const MaskToProbsByState * rightM2PBS = rightFPS->getMapPtrForDownPass(rightDown);
				assert(rightM2PBS != 0L);
				const std::vector<double> * rightProbs = getProbsForStatesMask(rightM2PBS, rightAllStates);
				assert(rightProbs != 0L);

				const BitField ancAllField = (1|2);
				MaskToProbsByState & forCurrScoreDownPass = forCurrScore.byDownPass[downPass];
				std::vector<double> * ancVec = getMutableProbsForStatesMask(&forCurrScoreDownPass, ancAllField);
				if (ancVec == 0L) { // if we have not visited this set of probabilities, start with a vector of 0's
					ancVec = &(forCurrScoreDownPass[ancAllField]); // get the memory
					ancVec->assign(blob.nRates*blob.nStates, 0.0); // set it to 0.0
				}
				addToAncProbVec(*ancVec, leftPMatVec, leftProbs, rightPMatVec, rightProbs, blob);
			}
			std::swap(leftFPS, rightFPS);
		}
	}
	// if the ancestor has state set of {0,1}, and currScore = 1 then
	// 	both children must display a single state.
	downPass = 3;
	unsigned int leftAccum = 0;
	BitField leftDown = 1; // set {0}  is 1 in our bitfield notation
	unsigned int rightAccum = 0;
	BitField rightDown = 1; // set {0}  is 1 in our bitfield notation
	const ProbForParsScore & leftFPS = leftPI.getByParsScore(leftAccum);
	const MaskToProbsByState * leftM2PBS = leftFPS.getMapPtrForDownPass(leftDown);
	assert(leftM2PBS != 0L);
	const ProbForParsScore & rightFPS = rightPI.getByParsScore(rightAccum);
	const MaskToProbsByState * rightM2PBS = rightFPS.getMapPtrForDownPass(rightDown);
	assert(rightM2PBS != 0L);
	const BitField leftAllStates = 1; // set {0}  is 1 in our bitfield notation
	const BitField rightAllStates = 1; // set {0}  is 1 in our bitfield notation
	const std::vector<double> * leftProbs = getProbsForStatesMask(leftM2PBS, leftAllStates);
	assert(leftProbs != 0L);
	const std::vector<double> * rightProbs = getProbsForStatesMask(rightM2PBS, rightAllStates);
	assert(rightProbs != 0L);

	// even though, we are accessing state set {0} for each child, we are using
	// symmetry to treat one of the children as having state {1}.
	// thus, our ancestor has observed state set {0, 1} (or 3 in BitField notation).
	const BitField ancAllField = 3;
	MaskToProbsByState & forCurrScoreDownPass = forCurrScore.byDownPass[0];
	std::vector<double> * ancVec = getMutableProbsForStatesMask(&forCurrScoreDownPass, ancAllField);
	if (ancVec == 0L) { // if we have not visited this set of probabilities, start with a vector of 0's
		ancVec = &(forCurrScoreDownPass[ancAllField]); // get the memory
		ancVec->assign(blob.nRates*blob.nStates, 0.0); // set it to 0.0
	}
	std::vector<unsigned int> stateCodeTranslationVec;
	for (unsigned i = 0; i < blob.nStates; ++i) {
		stateCodeTranslationVec.push_back(i);
	}
	stateCodeTranslationVec[0] = 1;
    stateCodeTranslationVec[1] = 0;
	addToAncProbVecSymmetric(*ancVec, leftPMatVec, leftProbs, rightPMatVec, rightProbs, stateCodeTranslationVec, blob);
	addToAncProbVecSymmetric(*ancVec, rightPMatVec, rightProbs, leftPMatVec, leftProbs, stateCodeTranslationVec, blob);

/*
	for (BitField downPass = 1; ; ++downPass) {
		const unsigned numStatesInMask = blob.getNumStates(downPass);
		if (numStatesInMask - 1 <= currScore) { // we cannot demand 3 states seen, but only 1 parsimony change... (all the probs will be zero, so we can skip them)

			MaskToProbsByState & forCurrScoreDownPass = forCurrScore.byDownPass[downPass];

			if (blob.getNumStates(downPass) > 1) {
				std::cerr << "from line: " << __LINE__<< ": downPass = " << (int)downPass << " " << blob.toSymbol(downPass) << " UNIONS:\n";
				const VecMaskPair & forUnions = blob.pairsForUnionForEachDownPass[downPass];
				bool didCalculations = this->allCalcsForAllPairs(forCurrScoreDownPass,
										  forUnions,
										  leftPI,
										  leftPMatVec,
										  rightPI,
										  rightPMatVec,
										  currScore - 1,
										  false,
										  blob);
				if (didCalculations)
					scObserved = true;
			}
			if (leftMaxP + rightMaxP >= currScore) {
				std::cerr << "from line: " << __LINE__<< ": downPass = " << blob.toSymbol(downPass) << " INTERSECTIONS:\n";
				const VecMaskPair & forIntersections = blob.pairsForIntersectionForEachDownPass[downPass];
				bool didCalculations = this->allCalcsForAllPairs(forCurrScoreDownPass,
										  forIntersections,
										  leftPI,
										  leftPMatVec,
										  rightPI,
										  rightPMatVec,
										  currScore,
										  true,
										  blob) || scObserved;
				if (didCalculations)
					scObserved = true;
			}
		}
		if (downPass == blob.lastBitField)
			break;
		assert(downPass < blob.lastBitField);
	}
	*/
	if (scObserved)
		obsMaxParsScore = currScore;

	for (currScore = 2;	 currScore <= maxParsScore; ++currScore) {
#if defined DEBUGGING_OUTPUT
		std::cerr << "from line: " << __LINE__<< ": currScore = " << currScore << ".\n";
#endif
		bool scObserved = false;
		ProbForParsScore & forCurrScore = this->byParsScore[currScore];
		for (BitField downPass = 1; ; ++downPass) {
			const unsigned numStatesInMask = blob.getNumStates(downPass);
			if (numStatesInMask - 1 <= currScore) { // we cannot demand 3 states seen, but only 1 parsimony change... (all the probs will be zero, so we can skip them)

				MaskToProbsByState & forCurrScoreDownPass = forCurrScore.byDownPass[downPass];

				if (blob.getNumStates(downPass) > 1) {
#					if defined DEBUGGING_OUTPUT
						std::cerr << "from line: " << __LINE__<< ": downPass = " << (int)downPass << " " << blob.toSymbol(downPass) << " UNIONS:\n";
#					endif
					const VecMaskPair & forUnions = blob.pairsForUnionForEachDownPass[downPass];
					scObserved = this->allCalcsForAllPairs(forCurrScoreDownPass,
											  forUnions,
											  leftPI,
											  leftPMatVec,
											  rightPI,
											  rightPMatVec,
											  currScore - 1,
											  false,
											  blob) || scObserved;
				}
				if (leftMaxP + rightMaxP >= currScore) {
#					if defined DEBUGGING_OUTPUT
						std::cerr << "from line: " << __LINE__<< ": downPass = " << blob.toSymbol(downPass) << " INTERSECTIONS:\n";
#					endif
					const VecMaskPair & forIntersections = blob.pairsForIntersectionForEachDownPass[downPass];
					scObserved = this->allCalcsForAllPairs(forCurrScoreDownPass,
											  forIntersections,
											  leftPI,
											  leftPMatVec,
											  rightPI,
											  rightPMatVec,
											  currScore,
											  true,
											  blob) || scObserved;
				}
			}
			if (downPass == blob.lastBitField)
				break;
			assert(downPass < blob.lastBitField);
		}
		if (scObserved)
			obsMaxParsScore = currScore;
	}
	if (obsMaxParsScore < maxParsScore) {
#		if defined DEBUGGING_OUTPUT
			std::cerr << "from line: " << __LINE__<< '\n'; std::cerr << "maxParsScore  = " << maxParsScore << " obsMaxParsScore = " << obsMaxParsScore << "\n";
#		endif
		this->byParsScore.resize(obsMaxParsScore + 1);
	}
}

void ProbInfo::calculate(const ProbInfo & leftPI, double leftEdgeLen,
					   const ProbInfo & rightPI, double rightEdgeLen,
					   TiMatFunc fn,
					   const CommonInfo & blob) {
#	if defined DEBUGGING_OUTPUT
		std::cerr << "from line " << __LINE__ << ":\n" ; std::cerr << "calculate \n" ;
#	endif
	fn(leftEdgeLen, blob.firstMatVec.GetAlias());
	const double *** leftPMatVec = const_cast<const double ***>(blob.firstMatVec.GetAlias());
	fn(rightEdgeLen, blob.secondMatVec.GetAlias());
	const double *** rightPMatVec = const_cast<const double ***>(blob.secondMatVec.GetAlias());

	const unsigned int leftMaxP = leftPI.getMaxParsScore();
	const unsigned rightMaxP = rightPI.getMaxParsScore();
	const unsigned maxParsScore = 1 + leftMaxP + rightMaxP;
	this->byParsScore.clear();
	this->byParsScore.resize(maxParsScore + 1); // add one to account zero
	this->nLeavesBelow = leftPI.getNLeavesBelow() + rightPI.getNLeavesBelow();
	unsigned obsMaxParsScore = 0;
	// do the calculations for staying in the constant patterns, these are more simple than the general calcs...
	if (true) { //@TEMP if true so that variables are scoped.
		const ProbForParsScore & leftFPS = leftPI.getByParsScore(0);
		const ProbForParsScore & rightFPS = rightPI.getByParsScore(0);
		ProbForParsScore & forZeroSteps = this->byParsScore[0];


		unsigned stateIndex = 0;
		for (std::vector<BitField>::const_iterator scIt = blob.singleStateCodes.begin();
				scIt != blob.singleStateCodes.end();
				++scIt, ++stateIndex) {
			const BitField sc = *scIt;
			const std::vector<double> * leftProbs = leftFPS.getProbsForDownPassAndObsMask(sc, sc);
			assert(leftProbs != 0L);
			const std::vector<double> * rightProbs = rightFPS.getProbsForDownPassAndObsMask(sc, sc);
			assert(rightProbs != 0L);

			std::vector<double> & toFillVec = forZeroSteps.byDownPass[sc][sc];
			toFillVec.assign(blob.nRates*blob.nStates, 0.0);

			addToAncProbVec(toFillVec, leftPMatVec, leftProbs, rightPMatVec, rightProbs, blob);
		}
	}

	// order N
	for (unsigned currScore = 1; currScore <= maxParsScore; ++currScore) {
#		if defined DEBUGGING_OUTPUT
			std::cerr << "from line: " << __LINE__<< ": currScore = " << currScore << ".\n";
#		endif
		bool scObserved = false;
		ProbForParsScore & forCurrScore = this->byParsScore[currScore];
		for (BitField downPass = 1; ; ++downPass) {
			const unsigned numStatesInMask = blob.getNumStates(downPass);
			if (numStatesInMask - 1 <= currScore) { // we cannot demand 3 states seen, but only 1 parsimony change... (all the probs will be zero, so we can skip them)

				MaskToProbsByState & forCurrScoreDownPass = forCurrScore.byDownPass[downPass];

				if (blob.getNumStates(downPass) > 1) {
#					if defined DEBUGGING_OUTPUT
						std::cerr << "from line: " << __LINE__<< ": downPass = " << (int)downPass << " " << blob.toSymbol(downPass) << " UNIONS:\n";
#					endif
					const VecMaskPair & forUnions = blob.pairsForUnionForEachDownPass[downPass];
					scObserved = this->allCalcsForAllPairs(forCurrScoreDownPass,
											  forUnions,
											  leftPI,
											  leftPMatVec,
											  rightPI,
											  rightPMatVec,
											  currScore - 1,
											  false,
											  blob) || scObserved;
				}
				if (leftMaxP + rightMaxP >= currScore) {
#					if defined DEBUGGING_OUTPUT
						std::cerr << "from line: " << __LINE__<< ": downPass = " << blob.toSymbol(downPass) << " INTERSECTIONS:\n";
#					endif
					const VecMaskPair & forIntersections = blob.pairsForIntersectionForEachDownPass[downPass];
					scObserved = this->allCalcsForAllPairs(forCurrScoreDownPass,
											  forIntersections,
											  leftPI,
											  leftPMatVec,
											  rightPI,
											  rightPMatVec,
											  currScore,
											  true,
											  blob) || scObserved;
				}
			}
			if (downPass == blob.lastBitField)
				break;
			assert(downPass < blob.lastBitField);
		}
		if (scObserved)
			obsMaxParsScore = currScore;
	}
	if (obsMaxParsScore < maxParsScore) {
#		if defined DEBUGGING_OUTPUT
			std::cerr << "from line: " << __LINE__<< '\n'; std::cerr << "maxParsScore  = " << maxParsScore << " obsMaxParsScore = " << obsMaxParsScore << "\n";
#		endif
		this->byParsScore.resize(obsMaxParsScore + 1);
	}
}

// \returns true if there were probabalities that were summed
bool ProbInfo::allCalcsForAllPairs(
			MaskToProbsByState & forCurrScoreDownPass,
			const VecMaskPair & pairVec,
			const ProbInfo & leftPI,
			const double *** leftPMatVec,
			const ProbInfo & rightPI,
			const double *** rightPMatVec,
			const unsigned accumScore,
			const bool doingIntersection,
			const CommonInfo & blob)
{	// order (2^k)^2
	const unsigned leftMaxP = leftPI.getMaxParsScore();
	const unsigned rightMaxP = rightPI.getMaxParsScore();
	bool probsAdded = false;
	for (VecMaskPair::const_iterator fuIt = pairVec.begin(); fuIt != pairVec.end(); ++fuIt) {
		const BitField leftDown = fuIt->first;
		const BitField rightDown = fuIt->second;
#		if defined DEBUGGING_OUTPUT
			std::cerr << "from line: " << __LINE__<< ": "; std::cerr << "leftDown = " << blob.toSymbol(leftDown) << " rightDown = " << blob.toSymbol(rightDown) << '\n';
#		endif
		assert(leftDown > 0);
		assert(rightDown > 0);
		if (doingIntersection) {
			assert((leftDown & rightDown) != 0);
		}
		else {
			assert((leftDown & rightDown) == 0);
		}
		const unsigned leftMinAccum = blob.getNumStates(leftDown) - 1;
		const unsigned rightMinAccum = blob.getNumStates(rightDown) - 1;
		if (leftMinAccum + rightMinAccum > accumScore) {
#			if defined DEBUGGING_OUTPUT
				std::cerr << "from line: " << __LINE__<< ": minScore exceed required score\n";
#			endif
			continue;
		}
		const unsigned acMaxLeftAccum = std::min(accumScore - rightMinAccum, leftMaxP);
		const unsigned acMaxRightAccum = std::min(accumScore - leftMinAccum, rightMaxP);
		const unsigned acMinLeftAccum = std::max(leftMinAccum, accumScore - acMaxRightAccum);
		if (false) {
#			if defined DEBUGGING_OUTPUT
				std::cerr << "accumScore = " << accumScore << '\n';
				std::cerr << "accumScore = " << accumScore << '\n';
				std::cerr << "leftMinAccum = " << leftMinAccum << '\n';
				std::cerr << "rightMinAccum = " << rightMinAccum << '\n';
				std::cerr << "acMaxLeftAccum = " << acMaxLeftAccum << '\n';
				std::cerr << "acMaxRightAccum = " << acMaxRightAccum << '\n';
				std::cerr << "acMinLeftAccum = " << acMinLeftAccum << '\n';
#			endif
		}
		// order N
		for (unsigned leftAccum = acMinLeftAccum; leftAccum <= acMaxLeftAccum; ++leftAccum) {
			const unsigned rightAccum = accumScore - leftAccum;
			assert(rightAccum <= rightMaxP);
			const ProbForParsScore & leftFPS = leftPI.getByParsScore(leftAccum);
			const MaskToProbsByState * leftM2PBS = leftFPS.getMapPtrForDownPass(leftDown);
			if (leftM2PBS == 0L) {
#				if defined DEBUGGING_OUTPUT
					std::cerr << "from line: " << __LINE__<< ": left child empty row. Skipping...\n";
#				endif
				continue;
			}
			const ProbForParsScore & rightFPS = rightPI.getByParsScore(rightAccum);
			const MaskToProbsByState * rightM2PBS = rightFPS.getMapPtrForDownPass(rightDown);
			if (rightM2PBS == 0L) {
#				if defined DEBUGGING_OUTPUT
					std::cerr << "from line: " << __LINE__<< ": right child empty row. Skipping...\n";
#				endif
				continue;
			}
			const BitFieldRow & leftSSRow = blob.statesSupersets[leftDown];
			// order (2^k)
			for (BitFieldRow::const_iterator lasIt = leftSSRow.begin(); lasIt != leftSSRow.end(); ++lasIt) {
				const BitField leftAllStates = *lasIt;
#				if defined DEBUGGING_OUTPUT
					std::cerr << "from line: " << __LINE__<< ":	 leftAllStates="  << blob.toSymbol(leftAllStates) << '\n';
#				endif
				const std::vector<double> * leftProbs = getProbsForStatesMask(leftM2PBS, leftAllStates);
				if (leftProbs == 0L) {
#					if defined DEBUGGING_OUTPUT
						std::cerr << "from line: " << __LINE__<< ": left child empty bin. Skipping...\n";
#					endif
					continue;
				}
				const BitFieldRow & rightSSRow = blob.statesSupersets[rightDown];
				// order (2^k)
				for (BitFieldRow::const_iterator rasIt = rightSSRow.begin(); rasIt != rightSSRow.end(); ++rasIt) {
					const BitField rightAllStates = *rasIt;
//					std::cerr << "from line: " << __LINE__<< ":  rightAllStates="  << blob.toSymbol(rightAllStates) << '\n';
					const std::vector<double> * rightProbs = getProbsForStatesMask(rightM2PBS, rightAllStates);
					if (rightProbs == 0L) {
//						std::cerr << "from line: " << __LINE__<< ": right child empty bin. Skipping...\n";
						continue;
					}
					const BitField ancAllField = (rightAllStates|leftAllStates);

					std::vector<double> * ancVec = getMutableProbsForStatesMask(&forCurrScoreDownPass, ancAllField);
					if (ancVec == 0L) {
						ancVec = &(forCurrScoreDownPass[ancAllField]);
						ancVec->assign(blob.nRates*blob.nStates, 0.0);
					}
					probsAdded = true;
//					std::cerr << __LINE__ << " adding:";
//					std::cerr << " leftDown="  << blob.toSymbol(leftDown)  << " leftAccum="  << leftAccum  << " leftAllStates="  << blob.toSymbol(leftAllStates);
//					std::cerr << " rightDown=" << blob.toSymbol(rightDown) << " rightAccum=" << rightAccum << " rightAllStates=" << blob.toSymbol(rightAllStates) << " \n";
					addToAncProbVec(*ancVec, leftPMatVec, leftProbs, rightPMatVec, rightProbs, blob);
				}
			}
		}
	}
	return probsAdded;
}


void ProbInfo::addToAncProbVecSymmetric(std::vector<double> & pVec,
		const double *** leftPMatVec, const std::vector<double> * leftProbs,
		const double *** rightPMatVec, const std::vector<double> * rightProbs,
		const std::vector<unsigned int> & rightChildStateCodeTranslation,
		const CommonInfo & blob) {
	if (leftProbs == 0L || rightProbs == 0L)
		return;
	unsigned rOffset = 0;
	// ignore loop over rates blob.nRates = 1
	for (unsigned r = 0; r < blob.nRates; ++r) {
		// this code looks up the correct transition prob matrix
		const double ** leftPMat = leftPMatVec[r];
		const double ** rightPMat = rightPMatVec[r];


		for (unsigned ancState = 0; ancState < blob.nStates; ++ancState) {
			double leftProb = 0.0;
			double rightProb = 0.0;
			for (unsigned desState = 0; desState < blob.nStates; ++desState) {

				const double leftTiProb = leftPMat[ancState][desState];
				const double leftAccumProb = (*leftProbs)[rOffset + desState];
				leftProb += leftTiProb*leftAccumProb;


				const unsigned int translatedState = rightChildStateCodeTranslation[desState];
				const double rightTiProb = rightPMat[ancState][translatedState];
				const double rightAccumProb = (*rightProbs)[rOffset + desState];
#				if defined DEBUGGING_OUTPUT
					std::cerr << "addToAncProbVecSymmetric tr = " << translatedState << " des " << desState << " rightTiProb= " << rightTiProb << " rightAccumProb = " << rightAccumProb <<'\n';
#				endif
				rightProb += rightTiProb*rightAccumProb;
			}
			pVec[rOffset + ancState] += leftProb*rightProb;
		}
		rOffset += blob.nStates;
	}
}

void ProbInfo::addToAncProbVec(std::vector<double> & pVec,
		const double *** leftPMatVec, const std::vector<double> * leftProbs,
		const double *** rightPMatVec, const std::vector<double> * rightProbs,
		const CommonInfo & blob) {
	if (leftProbs == 0L || rightProbs == 0L)
		return;
	unsigned rOffset = 0;
	// ignore loop over rates blob.nRates = 1
	for (unsigned r = 0; r < blob.nRates; ++r) {
		// this code looks up the correct transition prob matrix
		const double ** leftPMat = leftPMatVec[r];
		const double ** rightPMat = rightPMatVec[r];


		for (unsigned ancState = 0; ancState < blob.nStates; ++ancState) {
			double leftProb = 0.0;
			double rightProb = 0.0;
			for (unsigned desState = 0; desState < blob.nStates; ++desState) {
				const double leftTiProb = leftPMat[ancState][desState];
				const double leftAccumProb = (*leftProbs)[rOffset + desState];
				leftProb += leftTiProb*leftAccumProb;



				const double rightTiProb = rightPMat[ancState][desState];
				const double rightAccumProb = (*rightProbs)[rOffset + desState];
				rightProb += rightTiProb*rightAccumProb;
			}
			pVec[rOffset + ancState] += leftProb*rightProb;
		}
		rOffset += blob.nStates;
	}
}

void JCTiMat(double edgeLength, TiMat pMat ) {
#	if defined DEBUGGING_OUTPUT
		std::cerr << "JCTiMat edgeLength = " << edgeLength << '\n';
#	endif
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
#	if defined DEBUGGING_OUTPUT
		std::cerr << "from line " << __LINE__ << ":\n" ;  std::cerr << "edgeLength = " << edgeLength << " probs = " << prob_nochange << ", " << prob_change << "\n";
#	endif
}

void JCMulitCatTiMat(double edgeLength, TiMatVec pMatVec) {
	assert(gBlob);
#	if defined DEBUGGING_OUTPUT
		std::cerr << "JCMulitCatTiMat edgeLength = " << edgeLength << '\n';
#	endif
	for (unsigned i = 0; i < gBlob->rates.size(); ++i) {
#		if defined DEBUGGING_OUTPUT
			std::cerr << "JCMulitCatTiMat gBlob->rates[" << i << "] = " << gBlob->rates[i] << '\n';
#		endif
		JCTiMat(edgeLength*(gBlob->rates[i]), pMatVec[i]); //@TEMP JC
	}
}

void GenericMulitCatTiMat(double edgeLength, TiMatVec pMatVec) {
	assert(gBlob);
	int eigenIndex = 0;
	int numToCalc = gBlob->nRates;
	gBlob->scaledEdgeLengths.resize(gBlob->nRates);
	gBlob->modelIndexVector.resize(gBlob->nRates);
#	if defined DEBUGGING_OUTPUT
		std::cerr << "JCMulitCatTiMat edgeLength = " << edgeLength << '\n';
#	endif
	for (int i = 0; i < gBlob->nRates; ++i) {
		gBlob->scaledEdgeLengths[i] = edgeLength*gBlob->rates[i];
#		if defined DEBUGGING_OUTPUT
			std::cerr << "GenericMulitCatTiMat rates[" << i << "] = " << gBlob->rates[i] << '\n';
#		endif
		gBlob->modelIndexVector[i] = i;
	}
	int rc = calcPrMatsForHandle(gBlob->likeCalcHandle, eigenIndex, numToCalc, &(gBlob->scaledEdgeLengths[0]), &(gBlob->modelIndexVector[0]));
	if (rc != 0) {
		throw NxsException("Call to calcPrMatsForHandle failed");
	}
	for (int i = 0; i < gBlob->nRates; ++i) {
		fetchPrMat(gBlob->likeCalcHandle, i, pMatVec[i][0]);
#		if defined DEBUGGING_OUTPUT
			std::cerr << "from line " << __LINE__ << ":\n" ;  std::cerr << "scaled edgeLength = " << gBlob->scaledEdgeLengths[i] << " probs = " << pMatVec[i][0][0] << ", " << pMatVec[i][0][1] << "\n";
#		endif
	}

}

unsigned PatternSummary::incrementCount(unsigned s, BitField m, unsigned toAdd) {
	if (s >= this->byParsScore.size())
		this->byParsScore.resize(s + 1);
	BitsToCount & mapper = this->byParsScore[s];
	BitsToCount::iterator mIt = mapper.find(m);
	if (mIt == mapper.end()) {
		mapper[m] = toAdd;
		return toAdd;
	}
	unsigned prev = mIt->second;
	mIt->second = toAdd + prev;
	return toAdd + prev;
}

void PatternSummary::write(std::ostream & out, const CommonInfo & blob) const {
	const BitsToCount & constPatsMap = this->byParsScore[0];
	for (std::vector<BitField>::const_iterator cIt = blob.singleStateCodes.begin(); cIt != blob.singleStateCodes.end(); ++cIt) {
		BitsToCount::const_iterator queryIt = constPatsMap.find(*cIt);
		unsigned obsCount = (queryIt == constPatsMap.end() ? 0 : queryIt->second);
		out << "Observed steps = 0 states = " << blob.toSymbol(*cIt) << " count = " << obsCount	 << '\n';
	}
	for (unsigned score = 1; score < this->byParsScore.size(); ++score) {
		const BitsToCount & currPatsMap = this->byParsScore[score];
		for (std::vector<BitField>::const_iterator cIt = blob.multiStateCodes.begin(); cIt != blob.multiStateCodes.end(); ++cIt) {
			BitsToCount::const_iterator queryIt = currPatsMap.find(*cIt);
			unsigned obsCount = (queryIt == currPatsMap.end() ? 0 : queryIt->second);
			const unsigned ns = blob.getNumStates(*cIt);
			if (ns > 1 && (ns - 1) <= score) {
				out << "Observed steps = " << score << " states = " << blob.toSymbol(*cIt) << " count = " << obsCount  << '\n';
			}
			else {
				assert(obsCount == 0);
			}
		}
	}
}

void calculateUninformativePatternClassProbabilities(const NxsSimpleTree & tree, std::ostream & out, TiMatFunc tiMatFunc, const CommonInfo & blob) {
    cout << "blah\n";
}
void calculatePatternClassProbabilities(const NxsSimpleTree & tree, std::ostream & out, TiMatFunc tiMatFunc, const CommonInfo & blob) {
	std::vector<const NxsSimpleNode *> preorderVec = tree.GetPreorderTraversal();
	NodeIDToProbInfo nodeIDToProbInfo;
	ProbInfo * rootProbInfo = 0L;
	bool needToDelRootProbInfo = false;
	ProbInfo tipProbInfo;
	tipProbInfo.createForTip(blob);
	try {
		int ndInd = preorderVec.size() - 1;
		for (; ndInd >= 0; --ndInd) {
			const NxsSimpleNode * nd = preorderVec[ndInd];
			std::vector<NxsSimpleNode *> children = nd->GetChildren();
			const unsigned numChildren = children.size();
#			if defined DEBUGGING_OUTPUT
				std::cerr << "from line " << __LINE__ << ":\n" ; std::cerr << "In calculatePatternClassProbabilities at node " << nd->GetTaxonIndex() << ", #children = " << numChildren << "\n";
#			endif
			NodeID currNdId(nd, 0);
			if (numChildren == 0) {
				nodeIDToProbInfo[currNdId] = &tipProbInfo;
			}
			else {
				if (numChildren == 1)
					throw NxsException("Trees of degree 2 are not supported, yet\n");
				ProbInfo * currProbInfo = new ProbInfo(); // allocation
				nodeIDToProbInfo[currNdId] = currProbInfo;

				const NxsSimpleNode * leftNd = children[0];
				const NxsSimpleNode * rightNd = children[1];
				NodeIDToProbInfo::const_iterator leftPIIt= nodeIDToProbInfo.find(NodeID(leftNd, 0));
				assert(leftPIIt != nodeIDToProbInfo.end());
				ProbInfo * leftPI = leftPIIt->second;
				assert(leftPI);
				NodeIDToProbInfo::const_iterator rightPIIt= nodeIDToProbInfo.find(NodeID(rightNd, 0));
				assert(rightPIIt != nodeIDToProbInfo.end());
				ProbInfo * rightPI = rightPIIt->second;
				assert(rightPI != 0L);
				ProbInfo lt, rt;
				if (blob.isSymmetric) {
					currProbInfo->calculateSymmetric(*leftPI, leftNd->GetEdgeToParent().GetDblEdgeLen(),
											*rightPI, rightNd->GetEdgeToParent().GetDblEdgeLen(),
											tiMatFunc, blob);
				}
				else {
					currProbInfo->calculate(*leftPI, leftNd->GetEdgeToParent().GetDblEdgeLen(),
											*rightPI, rightNd->GetEdgeToParent().GetDblEdgeLen(),
											tiMatFunc, blob);
				}
				if (leftPI->getNLeavesBelow() > 1) {
					delete leftPI;
					nodeIDToProbInfo[NodeID(leftNd, 0)] = 0L;
				}
				if (rightPI->getNLeavesBelow() > 1) {
					delete rightPI;
					nodeIDToProbInfo[NodeID(rightNd, 0)] = 0L;
				}


				if (numChildren > 2) {
					if (nd != preorderVec[0] || numChildren > 3)
						throw NxsException("Parsimony scoring on non-binary trees is not supported, yet\n");
					const NxsSimpleNode * lastNd = children.at(2);
					NodeIDToProbInfo::const_iterator lastPIIt= nodeIDToProbInfo.find(NodeID(lastNd, 0));
					assert(lastPIIt != nodeIDToProbInfo.end());
					ProbInfo * lastPI = lastPIIt->second;
					assert(lastPI);
					rootProbInfo = new ProbInfo();
					needToDelRootProbInfo = true;
					if (blob.isSymmetric) {
						rootProbInfo->calculateSymmetric(*currProbInfo, 0.0,
												*lastPI, lastNd->GetEdgeToParent().GetDblEdgeLen(),
												tiMatFunc, blob);
					}
					else {
						rootProbInfo->calculate(*currProbInfo, 0.0,
												*lastPI, lastNd->GetEdgeToParent().GetDblEdgeLen(),
												tiMatFunc, blob);
					}
				}
				else if (nd == preorderVec[0])
					rootProbInfo = currProbInfo;
			}
		}
		assert(rootProbInfo != 0L);

		const ExpectedPatternSummary eps(*rootProbInfo, blob);
		eps.write(std::cout, blob);

	}
	catch (...) {
		freeProbInfo(preorderVec, nodeIDToProbInfo);
		if (needToDelRootProbInfo)
			delete rootProbInfo;
		throw;
	}
	freeProbInfo(preorderVec, nodeIDToProbInfo);
	if (needToDelRootProbInfo)
		delete rootProbInfo;
}
// marginalize over downpass set, rate categories, and anc states
ExpectedPatternSummary::ExpectedPatternSummary(const ProbInfo & rootProbInfo, const CommonInfo & blob) {
	if (blob.isSymmetric) {

		const unsigned maxNumSteps = rootProbInfo.getMaxParsScore();
		std::vector<double> emptyRow(blob.lastBitField + 1, 0.0);
		this->probsByStepsThenObsStates.resize(maxNumSteps + 1, emptyRow);
		// special case for constant patterns
		const ProbForParsScore & constFPS = rootProbInfo.getByParsScore(0);
		const std::vector<double> * pVec = constFPS.getProbsForDownPassAndObsMask(1, 1);
		assert(pVec);
		double patClassProb = 0.0;
		std::vector<double>::const_iterator wtIt = blob.categStateProb.begin();
		std::vector<double>::const_iterator pIt = pVec->begin();
		for (; wtIt != blob.categStateProb.end() ; ++wtIt, ++pIt) {
			assert(pIt != pVec->end());
			patClassProb += (*wtIt) * (*pIt);
		}
		for (std::vector<BitField>::const_iterator scIt = blob.singleStateCodes.begin();
				scIt != blob.singleStateCodes.end();
				++scIt) {
			 const BitField downPass = *scIt; // for the constant pattern, there will only be one downpass and allstates
			 this->probsByStepsThenObsStates[0][downPass] = patClassProb;
		}

		// all other number of states
		for (unsigned i = 1; i <= maxNumSteps; ++i) {
			const ProbForParsScore & fps = rootProbInfo.getByParsScore(i);
			 for (BitField downPass= 1; ;++downPass) {
				for (BitField obsStates = 1; ; ++obsStates ) {
					double patClassProb = 0.0;
					 std::vector<double>::const_iterator wtIt = blob.categStateProb.begin();
					 const std::vector<double> * pVec = fps.getProbsForDownPassAndObsMask(downPass, obsStates);
					 if (pVec) {
						 std::vector<double>::const_iterator pIt = pVec->begin();
						 for (; wtIt != blob.categStateProb.end() ; ++wtIt, ++pIt) {
							assert(pIt != pVec->end());
							patClassProb += (*wtIt) * (*pIt);
						 }
					}
					this->probsByStepsThenObsStates[i][obsStates] += patClassProb;
					if (obsStates == blob.lastBitField)
						break;
				}
				if (downPass == blob.lastBitField)
					break;
			 }
		}
	}
	else {
		const unsigned maxNumSteps = rootProbInfo.getMaxParsScore();
		std::vector<double> emptyRow(blob.lastBitField + 1, 0.0);
		this->probsByStepsThenObsStates.resize(maxNumSteps + 1, emptyRow);
		const ProbForParsScore & constFPS = rootProbInfo.getByParsScore(0);
		for (std::vector<BitField>::const_iterator scIt = blob.singleStateCodes.begin();
				scIt != blob.singleStateCodes.end();
				++scIt) {
			 const BitField downPass = *scIt; // for the constant pattern, there will only be one downpass and allstates
			 const std::vector<double> * pVec = constFPS.getProbsForDownPassAndObsMask(downPass, downPass);
			 assert(pVec);
			 double patClassProb = 0.0;
			 std::vector<double>::const_iterator wtIt = blob.categStateProb.begin();
			 std::vector<double>::const_iterator pIt = pVec->begin();
			 for (; wtIt != blob.categStateProb.end() ; ++wtIt, ++pIt) {
				assert(pIt != pVec->end());
				patClassProb += (*wtIt) * (*pIt);
			 }
			 this->probsByStepsThenObsStates[0][downPass] = patClassProb;
		}

		for (unsigned i = 1; i <= maxNumSteps; ++i) {
			const ProbForParsScore & fps = rootProbInfo.getByParsScore(i);
			 for (BitField downPass= 1; ;++downPass) {
				for (BitField obsStates = 1; ; ++obsStates ) {
					double patClassProb = 0.0;
					 std::vector<double>::const_iterator wtIt = blob.categStateProb.begin();
					 const std::vector<double> * pVec = fps.getProbsForDownPassAndObsMask(downPass, obsStates);
					 if (pVec) {
						 std::vector<double>::const_iterator pIt = pVec->begin();
						 for (; wtIt != blob.categStateProb.end() ; ++wtIt, ++pIt) {
							assert(pIt != pVec->end());
							patClassProb += (*wtIt) * (*pIt);
						 }
					}
					this->probsByStepsThenObsStates[i][obsStates] += patClassProb;
					if (obsStates == blob.lastBitField)
						break;
				}
				if (downPass == blob.lastBitField)
					break;
			 }
		}
	}
}
void ExpectedPatternSummary::write(std::ostream & out, const CommonInfo & blob) const {

	const unsigned maxNumSteps = this->probsByStepsThenObsStates.size() - 1;
	const std::vector<double> & constFPS = this->probsByStepsThenObsStates[0];
	double totalProb = 0.0;
	for (std::vector<BitField>::const_iterator scIt = blob.singleStateCodes.begin();
		 	scIt != blob.singleStateCodes.end();
		 	++scIt) {
		 out << "Expected steps = 0 states = " << blob.toSymbol(*scIt) << " prob = " << constFPS[*scIt] << '\n';
		 totalProb += constFPS[*scIt];
	}

	for (unsigned i = 1; i <= maxNumSteps; ++i) {
		const std::vector<double> & currFPS = this->probsByStepsThenObsStates[i];
		for (BitField obsStates= 1; ;++obsStates) {
			const double p = currFPS[obsStates];
			const unsigned ns = blob.getNumStates(obsStates);
			if (ns > 1 && (ns - 1) <= i) {
				out << "Expected steps = " << i << " states = " << blob.toSymbol(obsStates) << " prob = " << p << '\n';
				totalProb += p;
			}
			else {
				assert(p == 0.0);
			}
			if (obsStates == blob.lastBitField)
			 	break;
		 }
	}
	out << "totalprob = " << totalProb << '\n';
}

void classifyObservedDataIntoClasses(
		const NxsSimpleTree & tree,
		const BitFieldMatrix &mat,
		const int * pwPtr,
		std::ostream & out,
		PatternSummary *summary,
		const CommonInfo & blob) {
	std::vector<const NxsSimpleNode *> preorderVec = tree.GetPreorderTraversal();
	NodeIDToParsInfo nodeIDToParsInfo;
	const ParsInfo * rootParsInfo = 0L;
	assert(blob.zeroVec.size() == mat[0].size());
	int ndInd = preorderVec.size() - 1;
	for (; ndInd >= 0; ndInd--) {
		const NxsSimpleNode * nd = preorderVec[ndInd];
		std::vector<NxsSimpleNode *> children = nd->GetChildren();
		const unsigned numChildren = children.size();
#		if defined DEBUGGING_OUTPUT
			std::cerr << "from line " << __LINE__ << ":\n" ; std::cerr << "In classifyObservedDataIntoClasses at node " << nd->GetTaxonIndex() << ", #children = " << numChildren << "\n";
#		endif
		NodeID currNdId(nd, 0);
		ParsInfo & currParsInfo = nodeIDToParsInfo[currNdId];
		if (numChildren == 0) {
			const unsigned taxInd = nd->GetTaxonIndex();
			currParsInfo.calculateForTip(mat.at(taxInd), blob);
		}
		else {
			if (numChildren == 1)
				throw NxsException("Trees of degree 2 are not supported, yet\n");
			const NxsSimpleNode * leftNd = children[0];
			const NxsSimpleNode * rightNd = children[1];
			NodeIDToParsInfo::const_iterator leftPIIt= nodeIDToParsInfo.find(NodeID(leftNd, 0));
			assert(leftPIIt != nodeIDToParsInfo.end());
			NodeIDToParsInfo::const_iterator rightPIIt= nodeIDToParsInfo.find(NodeID(rightNd, 0));
			assert(rightPIIt != nodeIDToParsInfo.end());
			currParsInfo.calculateForInternal(leftPIIt->second, rightPIIt->second);

			if (numChildren > 2) {
				if (nd != preorderVec[0] || numChildren > 3)
					throw NxsException("Parsimony scoring on non-binary trees is not supported, yet\n");
				NodeID lastNdId(nd, 1);
				ParsInfo & lastParsInfo = nodeIDToParsInfo[currNdId];
				const NxsSimpleNode * lastNd = children.at(2);
				NodeIDToParsInfo::const_iterator lastPIIt= nodeIDToParsInfo.find(NodeID(lastNd, 0));
				assert(lastPIIt != nodeIDToParsInfo.end());
				lastParsInfo.calculateForInternal(currParsInfo, lastPIIt->second);
				rootParsInfo = &lastParsInfo;
			}
			else if (nd == preorderVec[0])
				rootParsInfo = &currParsInfo;
		}
	}
	assert(rootParsInfo != 0L);
#	if defined DEBUGGING_OUTPUT
		std::cerr << "from line " << __LINE__ << ":\n" ; rootParsInfo->write(std::cerr);
#	endif
	if (summary) {
		summary->clear();
		for (unsigned p = 0; p < rootParsInfo->size(); ++p) {
			unsigned toAdd = (pwPtr != 0L ? (unsigned)pwPtr[p] : 1);
			summary->incrementCount(rootParsInfo->score[p], rootParsInfo->allSeen[p], toAdd);
		}
		summary->write(std::cout, blob);
	}
}

/// \throws NxsException for gaps or ambiguous cells
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
		std::vector<double>	v;
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
#		if defined DEBUGGING_OUTPUT
			std::cerr << "flagWithDash = " << flagWithDash << "	  val = " << val << '\n';
			std::cerr << "parsedToDoubles v =";
			for (std::vector<double>::const_iterator vIt = v.begin(); vIt != v.end(); ++vIt) {
				std::cerr << ' ' << *vIt;
			}
			std::cerr << '\n';
#		endif
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
		int blocksToRead =	(PublicNexusReader::NEXUS_TAXA_BLOCK_BIT
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
		NCL_COULD_BE_CONST	NxsCharactersBlock * charsBlock = nexusReader.GetCharactersBlock(taxaBlock, 0);
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
			NCL_COULD_BE_CONST	NxsTreesBlock * treesBlock = nexusReader.GetTreesBlock(taxaBlock, treesBlockInd);
			assert(treesBlock);
			treesBlock->ProcessAllTrees();
            const unsigned nTreesThisBlock = treesBlock->GetNumTrees();
			for (unsigned treeInd = 0; treeInd < nTreesThisBlock; ++treeInd) {
				const NxsFullTreeDescription & ftd = treesBlock->GetFullTreeDescription(treeInd);
				if (ftd.AllEdgesHaveLengths()) {
					NxsSimpleTree nclTree(ftd, 0, 0.0);
//					calculatePatternClassProbabilities(nclTree, std::cout, JCMulitCatTiMat, blob); //@TEMP JC
					calculateUninformativePatternClassProbabilities(nclTree, std::cout, GenericMulitCatTiMat, blob); //@TEMP JC

					if (false) {
						PatternSummary observed;
						classifyObservedDataIntoClasses(nclTree, bitFieldMatrix, pwPtr, std::cout, &observed, blob);
					}
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

