#include "phy-pat-class-prob.hpp"
#include "ncl/nxsmultiformat.h"
#include "ncl/nxsallocatematrix.h"



void freeProbInfo(const std::vector<const NxsSimpleNode *> & preorderVec, NodeIDToProbInfo & nodeIDToProbInfo);


// Globals (all should start with g[A-Z] pattern to make it easier to find and replace them later).
const unsigned MAX_NUM_STATES = 8*sizeof(BitField);

struct CommonInfo {
	unsigned nStates;
	unsigned nRates;
	unsigned pVecLen; // nStates * nRates;
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
	BitField lastBitField;
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

	private:
		std::vector<unsigned> stateCodeToNumStates;
		std::map<BitField, std::string> stateCodesToSymbols;
		
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

void freeProbInfo(const std::vector<const NxsSimpleNode *> & preorderVec, NodeIDToProbInfo & nodeIDToProbInfo) {

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
	this->alphabet = symbols;
	this->nStates = this->alphabet.length();
	this->nRates = 1; //@TEMP no rate het
	this->pVecLen = this->nRates*this->nStates;
	
	this->rates.assign(this->nRates, 1.0); //@TEMP no rate het
	this->categStateProb.assign(this->pVecLen, 1.0/((double)this->pVecLen)) ; //@TEMP no rate het //@TEMP JC
	
	this->singleStateCodes.clear();
	this->multiStateCodes.clear();
	this->stateCodesToSymbols.clear();
	
	this->firstMatVec.Free();
	this->firstMatVec.Initialize(this->nRates, this->nStates, this->nStates);
	this->secondMatVec.Free();
	this->secondMatVec.Initialize(this->nRates, this->nStates, this->nStates);
	
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


	
	std::cerr << "from line " << __LINE__ << ":\n";
	for (sc = 1;;++sc) {
		std::cerr << "Combos with unions that lead to " << this->toSymbol(sc) << ":";
		const VecMaskPair & forUnions = this->pairsForUnionForEachDownPass[sc];
		for (VecMaskPair::const_iterator fuIt = forUnions.begin(); fuIt != forUnions.end(); ++fuIt) {
			std::cerr << " (\"" << this->toSymbol(fuIt->first);
			std::cerr << "\", \"" << this->toSymbol(fuIt->second);
			std::cerr << "\")   ";
		}
		std::cerr << "\n";

		std::cerr << "Combos with intersections that lead to " << this->toSymbol(sc) << ":";
		const VecMaskPair & forIntersections = this->pairsForIntersectionForEachDownPass[sc];
		for (VecMaskPair::const_iterator fuIt = forIntersections.begin(); fuIt != forIntersections.end(); ++fuIt) {
			std::cerr << " (\"" << this->toSymbol(fuIt->first);
			std::cerr << "\", \"" << this->toSymbol(fuIt->second);
			std::cerr << "\")   ";
		}
		std::cerr << "\n";
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

void ProbInfo::calculate(const ProbInfo & leftPI, double leftEdgeLen, 
					   const ProbInfo & rightPI, double rightEdgeLen,
					   TiMatFunc fn,
					   const CommonInfo & blob) {
	std::cerr << "from line " << __LINE__ << ":\n" ; std::cerr << "calculate \n" ;
	fn(leftEdgeLen, blob.firstMatVec.GetAlias());
	const double *** leftPMatVec = const_cast<const double ***>(blob.firstMatVec.GetAlias());
	fn(rightEdgeLen, blob.secondMatVec.GetAlias());
	const double *** rightPMatVec = const_cast<const double ***>(blob.secondMatVec.GetAlias());
	
	const unsigned leftMaxP = leftPI.getMaxParsScore();
	const unsigned rightMaxP = rightPI.getMaxParsScore();
	const unsigned maxParsScore = 1 + leftMaxP + rightMaxP;
	this->byParsScore.clear();
	this->byParsScore.resize(maxParsScore + 1);
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
			const std::vector<double> * leftProbs = leftFPS.getProbsForDownPassAndMask(sc, sc);
			assert(leftProbs);
			const std::vector<double> * rightProbs = rightFPS.getProbsForDownPassAndMask(sc, sc);
			assert(rightProbs);
	
			std::vector<double> & pVec = forZeroSteps.byDownPass[sc][sc];
			pVec.assign(blob.nRates*blob.nStates, 0.0);
			
			addToAncProbVec(pVec, leftPMatVec, leftProbs, rightPMatVec, rightProbs, blob);
		}
	}
	
	// order N
	for (unsigned currScore = 1; currScore <= maxParsScore; ++currScore) {
		std::cerr << "from line: " << __LINE__<< ": currScore = " << currScore << ".\n";
		bool scObserved = false;
		ProbForParsScore & forCurrScore = this->byParsScore[currScore];
		for (BitField downPass = 1; ; ++downPass) {
			const unsigned numStatesInMask = blob.getNumStates(downPass);
			if (numStatesInMask - 1 <= currScore) { // we cannot demand 3 states seen, but only 1 parsimony change... (all the probs will be zero, so we can skip them)
			
				MaskToProbsByState & forCurrScoreDownPass = forCurrScore.byDownPass[downPass];
				
				if (blob.getNumStates(downPass) > 1) {
					std::cerr << "from line: " << __LINE__<< ": downPass = " << (int)downPass << " " << blob.toSymbol(downPass) << " UNIONS:\n";
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
					std::cerr << "from line: " << __LINE__<< ": downPass = " << blob.toSymbol(downPass) << " INTERSECTIONS:\n";
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
		std::cerr << "from line: " << __LINE__<< '\n'; std::cerr << "maxParsScore  = " << maxParsScore << " obsMaxParsScore = " << obsMaxParsScore << "\n";
		this->byParsScore.resize(obsMaxParsScore + 1);
	}
}

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
		std::cerr << "from line: " << __LINE__<< ": "; std::cerr << "leftDown = " << blob.toSymbol(leftDown) << " rightDown = " << blob.toSymbol(rightDown) << '\n';
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
			std::cerr << "from line: " << __LINE__<< ": minScore exceed required score\n";
			continue;
		}
		const unsigned acMaxLeftAccum = std::min(accumScore - rightMinAccum, leftMaxP);
		const unsigned acMaxRightAccum = std::min(accumScore - leftMinAccum, rightMaxP);
		const unsigned acMinLeftAccum = std::max(leftMinAccum, accumScore - acMaxRightAccum);
		if (false) {
			std::cerr << "accumScore = " << accumScore << '\n';
			std::cerr << "accumScore = " << accumScore << '\n';
			std::cerr << "leftMinAccum = " << leftMinAccum << '\n';
			std::cerr << "rightMinAccum = " << rightMinAccum << '\n';
			std::cerr << "acMaxLeftAccum = " << acMaxLeftAccum << '\n';
			std::cerr << "acMaxRightAccum = " << acMaxRightAccum << '\n';
			std::cerr << "acMinLeftAccum = " << acMinLeftAccum << '\n';
		}
		// order N
		for (unsigned leftAccum = acMinLeftAccum; leftAccum <= acMaxLeftAccum; ++leftAccum) {
			const unsigned rightAccum = accumScore - leftAccum;
			assert(rightAccum <= rightMaxP);
			const ProbForParsScore & leftFPS = leftPI.getByParsScore(leftAccum);
			const MaskToProbsByState * leftM2PBS = leftFPS.getMapPtrForDownPass(leftDown);
			if (leftM2PBS == 0L) {
				std::cerr << "from line: " << __LINE__<< ": left child empty row. Skipping...\n";
				continue;
			}
			const ProbForParsScore & rightFPS = rightPI.getByParsScore(rightAccum);
			const MaskToProbsByState * rightM2PBS = rightFPS.getMapPtrForDownPass(rightDown);
			if (rightM2PBS == 0L) {
				std::cerr << "from line: " << __LINE__<< ": right child empty row. Skipping...\n";
				continue;
			}
			const BitFieldRow & leftSSRow = blob.statesSupersets[leftDown];
			// order (2^k)
			for (BitFieldRow::const_iterator lasIt = leftSSRow.begin(); lasIt != leftSSRow.end(); ++lasIt) {
				const BitField leftAllStates = *lasIt;
				std::cerr << "from line: " << __LINE__<< ":  leftAllStates="  << blob.toSymbol(leftAllStates) << '\n';
				const std::vector<double> * leftProbs = getProbsForStatesMask(leftM2PBS, leftAllStates);
				if (leftProbs == 0L) {
					std::cerr << "from line: " << __LINE__<< ": left child empty bin. Skipping...\n";
					continue; 
				}
				const BitFieldRow & rightSSRow = blob.statesSupersets[rightDown];
				// order (2^k)
				for (BitFieldRow::const_iterator rasIt = rightSSRow.begin(); rasIt != rightSSRow.end(); ++rasIt) {
					const BitField rightAllStates = *rasIt;
					std::cerr << "from line: " << __LINE__<< ":  rightAllStates="  << blob.toSymbol(rightAllStates) << '\n';
					const std::vector<double> * rightProbs = getProbsForStatesMask(rightM2PBS, rightAllStates);
					if (rightProbs == 0L) {
						std::cerr << "from line: " << __LINE__<< ": right child empty bin. Skipping...\n";
						continue; 
					}
					const BitField ancAllField = (rightAllStates|leftAllStates);
					
					std::vector<double> * ancVec = getMutableProbsForStatesMask(&forCurrScoreDownPass, ancAllField);
					if (ancVec == 0L) {
						ancVec = &(forCurrScoreDownPass[ancAllField]);
						ancVec->assign(blob.nRates*blob.nStates, 0.0);
					}
					probsAdded = true;
					std::cerr << __LINE__ << " adding:";
					std::cerr << " leftDown="  << blob.toSymbol(leftDown)  << " leftAccum="  << leftAccum  << " leftAllStates="  << blob.toSymbol(leftAllStates);
					std::cerr << " rightDown=" << blob.toSymbol(rightDown) << " rightAccum=" << rightAccum << " rightAllStates=" << blob.toSymbol(rightAllStates) << " \n";
					addToAncProbVec(*ancVec, leftPMatVec, leftProbs, rightPMatVec, rightProbs, blob);
				}
			}
		}
	}
	return probsAdded;
}


void ProbInfo::addToAncProbVec(std::vector<double> & pVec, 
		const double *** leftPMatVec, const std::vector<double> * leftProbs,
		const double *** rightPMatVec, const std::vector<double> * rightProbs,
		const CommonInfo & blob) {
	if (leftProbs == 0L || rightProbs == 0L)
		return;
	unsigned rOffset = 0;
	for (unsigned r = 0; r < blob.nRates; ++r) {
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

void JCTiMat(double edgeLength, TiMat pMat) {
	std::cerr << "JCTiMat edgeLength = " << edgeLength << '\n';
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
	std::cerr << "from line " << __LINE__ << ":\n" ;  std::cerr << "edgeLength = " << edgeLength << " probs = " << prob_nochange << ", " << prob_change << "\n";
}

void JCMulitCatTiMat(double edgeLen, TiMatVec pMatVec) {
	assert(gBlob);
	std::cerr << "JCMulitCatTiMat edgeLen = " << edgeLen << '\n';
	for (unsigned i = 0; i < gBlob->rates.size(); ++i) {
		std::cerr << "JCMulitCatTiMat gBlob->rates[" << i << "] = " << gBlob->rates[i] << '\n';
		JCTiMat(edgeLen*(gBlob->rates[i]), pMatVec[i]); //@TEMP JC
	}
}

unsigned PatternSummary::incrementCount(unsigned s, BitField m) {
	if (s >= this->byParsScore.size())
		this->byParsScore.resize(s + 1);
	BitsToCount & mapper = this->byParsScore[s];
	BitsToCount::iterator mIt = mapper.find(m);
	if (mIt == mapper.end()) {
		mapper[m] = 1;
		return 1;
	}
	unsigned prev = mIt->second;
	mIt->second = 1 + prev;
	return 1 + prev;
}
		
void PatternSummary::write(std::ostream & out, const CommonInfo & blob) const {
	const BitsToCount & constPatsMap = this->byParsScore[0];
	for (std::vector<BitField>::const_iterator cIt = blob.singleStateCodes.begin(); cIt != blob.singleStateCodes.end(); ++cIt) {
		BitsToCount::const_iterator queryIt = constPatsMap.find(*cIt);
		unsigned obsCount = (queryIt == constPatsMap.end() ? 0 : queryIt->second);
		out << "Observed steps = 0 states = " << blob.toSymbol(*cIt) << " count = " << obsCount  << '\n';
	}
	for (unsigned score = 1; score < this->byParsScore.size(); ++score) {
		const BitsToCount & currPatsMap = this->byParsScore[score];
		for (std::vector<BitField>::const_iterator cIt = blob.multiStateCodes.begin(); cIt != blob.multiStateCodes.end(); ++cIt) {
			BitsToCount::const_iterator queryIt = currPatsMap.find(*cIt);
			unsigned obsCount = (queryIt == currPatsMap.end() ? 0 : queryIt->second);
			out << "Observed steps = " << score << " states = " << blob.toSymbol(*cIt) << " count = " << obsCount  << '\n';
		}
	}
}


void calculatePatternClassProbabilities(const NxsSimpleTree & tree, std::ostream & out, TiMatFunc tiMatFunc, const CommonInfo & blob) {
	std::vector<const NxsSimpleNode *> preorderVec = tree.GetPreorderTraversal();
	NodeIDToProbInfo nodeIDToProbInfo;
	ProbInfo * rootProbInfo = 0L;
	bool needToDelRootProbInfo = false;
	ProbInfo tipProbInfo;
	tipProbInfo.createForTip(blob);
	try {
		for (std::vector<const NxsSimpleNode *>::const_reverse_iterator ndIt = preorderVec.rbegin();
				ndIt != preorderVec.rend();
				++ndIt) {
			const NxsSimpleNode * nd = *ndIt;
			std::vector<NxsSimpleNode *> children = nd->GetChildren();
			const unsigned numChildren = children.size();
			std::cerr << "from line " << __LINE__ << ":\n" ; std::cerr << "In calculatePatternClassProbabilities at node " << nd->GetTaxonIndex() << ", #children = " << numChildren << "\n";
			NodeID currNdId(nd, 0);
			if (numChildren == 0) {
				nodeIDToProbInfo[currNdId] = &tipProbInfo;
			}
			else {
				if (numChildren == 1)
					throw NxsException("Trees of degree 2 are not supported, yet\n");
				ProbInfo * currProbInfo = new ProbInfo();
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
				assert(rightPI);
				ProbInfo lt, rt;
				currProbInfo->calculate(*leftPI, leftNd->GetEdgeToParent().GetDblEdgeLen(),
										*rightPI, rightNd->GetEdgeToParent().GetDblEdgeLen(),
										tiMatFunc, blob);
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
					rootProbInfo->calculate(*currProbInfo, 0.0,
											*lastPI, lastNd->GetEdgeToParent().GetDblEdgeLen(),
											tiMatFunc, blob);
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
	
	const unsigned maxNumSteps = rootProbInfo.getMaxParsScore();
	std::vector<double> emptyRow(blob.lastBitField + 1, 0.0);
	this->probsByStepsThenObsStates.resize(maxNumSteps + 1, emptyRow);
	const ProbForParsScore & constFPS = rootProbInfo.getByParsScore(0);
	for (std::vector<BitField>::const_iterator scIt = blob.singleStateCodes.begin(); 
		 	scIt != blob.singleStateCodes.end();
		 	++scIt) {
		 const BitField downPass = *scIt; // for the constant pattern, there will only be one downpass and allstates
		 const std::vector<double> * pVec = constFPS.getProbsForDownPassAndMask(downPass, downPass);
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
				 const std::vector<double> * pVec = fps.getProbsForDownPassAndMask(downPass, obsStates);
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
			out << "Expected steps = " << i << " states = " << blob.toSymbol(obsStates) << " prob = " << p << '\n';
			totalProb += p;
			if (obsStates == blob.lastBitField)
			 	break;
		 }
	}
	out << "totalprob = " << totalProb << '\n';
}

void classifyObservedDataIntoClasses(const NxsSimpleTree & tree, const BitFieldMatrix &mat, std::ostream & out, PatternSummary *summary, const CommonInfo & blob) {
	std::vector<const NxsSimpleNode *> preorderVec = tree.GetPreorderTraversal();
	NodeIDToParsInfo nodeIDToParsInfo;
	const ParsInfo * rootParsInfo = 0L;
	assert(blob.zeroVec.size() == mat[0].size());
	for (std::vector<const NxsSimpleNode *>::const_reverse_iterator ndIt = preorderVec.rbegin();
			ndIt != preorderVec.rend();
			++ndIt) {
		const NxsSimpleNode * nd = *ndIt;
		std::vector<NxsSimpleNode *> children = nd->GetChildren();
		const unsigned numChildren = children.size();
		std::cerr << "from line " << __LINE__ << ":\n" ; std::cerr << "In classifyObservedDataIntoClasses at node " << nd->GetTaxonIndex() << ", #children = " << numChildren << "\n";
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
	std::cerr << "from line " << __LINE__ << ":\n" ; rootParsInfo->write(std::cerr);
	if (summary) {
		summary->clear();
		for (unsigned p = 0; p < rootParsInfo->size(); ++p) {
			summary->incrementCount(rootParsInfo->score[p], rootParsInfo->allSeen[p]);
		}
		summary->write(std::cout, blob);
	}
}

/// \throws NxsException for gaps or ambiguous cells
/// \returns a string of the symbols for each state (length == nStates).
std::string convertToBitFieldMatrix(const NxsCharactersBlock & cb, BitFieldMatrix & mat, const NxsUnsignedSet * toInclude) {
	NxsString errormsg;
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
	
int main(int argc, char * argv[]) {
	/* Usually it is easiest to surround interactions with NCL in a try block
		that catches NxsException instances.
		Reading files can certainly generate these errors, but even queries
		after the parse can result in these exceptions.
	*/
	if (argc < 2) {
		std::cerr << "Expecting one arguments: a file name\n";
		return 1;
	}
	if (argv[1][0] == '-' &&  argv[1][1] == 'h' && argv[1][2] == '\0' ) {
		std::cerr << "Takes a arguments: The path to a NEXUS file with a single characters block\n";
		return 0;
	}

	std::string filename(argv[1]);
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
		BitFieldMatrix bitFieldMatrix;
		std::string symbols = convertToBitFieldMatrix(*charsBlock, bitFieldMatrix);
		CommonInfo blob;
		blob.initialize(symbols);
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
					calculatePatternClassProbabilities(nclTree, std::cout, JCMulitCatTiMat, blob); //@TEMP JC
					PatternSummary observed;
					classifyObservedDataIntoClasses(nclTree, bitFieldMatrix, std::cout, &observed, blob);
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
