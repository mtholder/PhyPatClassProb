
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
	//	one child must display a single state, and the other must have 2 states observed
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
	//	both children must display a single state.
	downPass = 3;
	unsigned int leftAccum = 0;
	BitField leftDown = 1; // set {0}  is 1 in our bitfield notation
	unsigned int rightAccum = 0;
	BitField rightDown = 1; // set {0}	is 1 in our bitfield notation
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
//					std::cerr << "from line: " << __LINE__<< ":	 rightAllStates="  << blob.toSymbol(rightAllStates) << '\n';
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
//					std::cerr << " leftDown="  << blob.toSymbol(leftDown)  << " leftAccum="	 << leftAccum  << " leftAllStates="	 << blob.toSymbol(leftAllStates);
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






