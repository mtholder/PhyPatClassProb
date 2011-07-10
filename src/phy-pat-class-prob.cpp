#include <iostream>
#include <cassert>
#include <vector>
#include "ncl/nxsmultiformat.h"

////////////////////////////////////////////////////////////////////////////////
// header
////////////////////////////////////////////////////////////////////////////////

typedef unsigned char BitField;
const unsigned MAX_NUM_STATES = 8*sizeof(BitField);
typedef std::vector<BitField> BitFieldRow;
typedef std::vector<BitFieldRow> BitFieldMatrix;


/// \throws NxsException for gaps or ambiguous cells
/// \returns a string of the symbols for each state (length == nStates).
std::string ConvertToBitFieldMatrix(const NxsCharactersBlock & cb, BitFieldMatrix & mat, const NxsUnsignedSet * toInclude=0L);
void writeBitFieldMatrix(std::ostream & out, const BitFieldMatrix & bitFieldMatrix);


////////////////////////////////////////////////////////////////////////////////
// cpp
////////////////////////////////////////////////////////////////////////////////

std::string ConvertToBitFieldMatrix(const NxsCharactersBlock & cb, BitFieldMatrix & mat, const NxsUnsignedSet * toInclude)
    {
    NxsString errormsg;
	std::vector<const NxsDiscreteDatatypeMapper *> mappers = cb.GetAllDatatypeMappers();
	if (mappers.empty() || mappers[0] == NULL)
		throw NxsException("no mappers");

    if (mappers.size() != 1)
        throw NxsException("Expecting an unmixed characters block, but found a matrix with datatype = mixed or a datatype with augmented symbols\n");
    NxsUnsignedSet scratchSet;
	if (toInclude == 0L)
		{
		for (unsigned i = 0; i < cb.GetNChar(); ++i)
			scratchSet.insert(i);
		toInclude = & scratchSet;
	 	}
	
    std::set <const NxsDiscreteDatatypeMapper * > usedMappers;
	for (NxsUnsignedSet::const_iterator indIt = toInclude->begin(); indIt != toInclude->end(); ++indIt)
		{
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
    
    for (NxsDiscreteStateCell i = 0; i < (NxsDiscreteStateCell)nStates; ++i)
        {
        if (mapper.GetStateSetForCode(i).size() != 1)
            throw NxsException("Expecting the first states to correspond to state sets with size 1");
        }
	const std::string fundamentalSymbols = mapper.GetSymbols();
	if (fundamentalSymbols.length() != nStates)
	    {
	    errormsg << "Expecting the fundamental symbols (" << fundamentalSymbols << ") to have length " << nStates;
	    throw NxsException(errormsg);
        }
	NCL_ASSERT((int)NXS_MISSING_CODE < 0);
	NCL_ASSERT((int)NXS_GAP_STATE_CODE < 0);

    
    const unsigned nTaxa = cb.GetNTax();
    const unsigned includedNChar = toInclude->size();
    mat.resize(nTaxa);
    for (unsigned i = 0; i < nTaxa; ++i)
        {
        BitFieldRow & bfRow = mat[i];
        bfRow.resize(includedNChar);
    	const NxsDiscreteStateRow & row = cb.GetDiscreteMatrixRow(i);
    	if (row.empty())
		    {
		    errormsg << "Empty row encountered for taxon " << (1+i) << ". Missing data is not supported.\n";
			throw NxsException(errormsg);
		    }
		unsigned j = 0; 
		for (NxsUnsignedSet::const_iterator tIncIt = toInclude->begin(); tIncIt != toInclude->end(); ++tIncIt, ++j)
		    {
		    const NxsDiscreteStateCell & cell = row.at(*tIncIt);
		    if (cell < 0 || cell >= (NxsDiscreteStateCell) nStates)
		        {
	    	    errormsg << "Ambiguous/missing/gap data found for taxon " << (1+i) << " at site " << *tIncIt << ". Missing data is not supported.\n";
    			throw NxsException(errormsg);
    			}
    		const int bfi = 1 << (int) cell;
    		bfRow[j] = BitField(bfi);
		    }
		}
	return fundamentalSymbols;
    }   

void writeBitFieldMatrix(std::ostream & out, const BitFieldMatrix & bitFieldMatrix)
    {
    for (unsigned i = 0; i < bitFieldMatrix.size(); ++i)
        {
        out << "taxon " << (i + 1) << ":\n";
        const BitFieldRow & row = bitFieldMatrix[i];
        for (unsigned j = 0; j < row.size(); ++j)
            out << ' ' << (int) row[j];
        out << '\n';
        }
    }
    
int main(int argc, char * argv[])
{
	/* Usually it is easiest to surround interactions with NCL in a try block
		that catches NxsException instances.
		Reading files can certainly generate these errors, but even queries
		after the parse can result in these exceptions.
	*/
	if (argc < 2)
		{
		std::cerr << "Expecting one arguments: a file name\n";
		return 1;
		}
	if (argv[1][0] == '-' &&  argv[1][1] == 'h' && argv[1][2] == '\0' )
		{
		std::cerr << "Takes a arguments: The path to a NEXUS file with a single characters block\n";
		return 0;
		}

	std::string filename(argv[1]);
	try {
		int blocksToRead =  (PublicNexusReader::NEXUS_TAXA_BLOCK_BIT
							| PublicNexusReader::NEXUS_CHARACTERS_BLOCK_BIT
							| PublicNexusReader::NEXUS_ASSUMPTIONS_BLOCK_BIT
							| PublicNexusReader::NEXUS_SETS_BLOCK_BIT
							);
		MultiFormatReader nexusReader(blocksToRead, NxsReader::WARNINGS_TO_STDERR);
        
		std::cerr << "Reading " << filename << "\n";
		try {
			nexusReader.ReadFilepath(filename.c_str(), MultiFormatReader::NEXUS_FORMAT);
			}
		catch(const NxsException &x)
			{
			std::cerr << "Error:\n " << x.msg << std::endl;
			if (x.line > 0 || x.pos > 0)
				std::cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << std::endl;
			return 2;
			}
		catch(...)
			{
			nexusReader.DeleteBlocksFromFactories();
			std::cerr << "Exiting with an unknown error" << std::endl;
			return 1;
			}

		const unsigned numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
		if (numTaxaBlocks != 1)
			{
			std::cerr << "Expecting a file with exactly 1 TAXA block, but found " << numTaxaBlocks << " in the file " << filename << ".\n";
			return 2;
			}
		NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(0);
		const unsigned nCharsBlocks = nexusReader.GetNumCharactersBlocks(taxaBlock);
		if (nCharsBlocks != 1)
			{
			std::cerr << "Expecting a file with exactly 1 CHARACTERS/DATA block, but found " << nCharsBlocks << " in the file " << filename << ".\n";
			return 3;
			}
		NCL_COULD_BE_CONST  NxsCharactersBlock * charsBlock = nexusReader.GetCharactersBlock(taxaBlock, 0);
		assert(charsBlock);
		BitFieldMatrix bitFieldMatrix;
		std::string symbols = ConvertToBitFieldMatrix(*charsBlock, bitFieldMatrix);

        writeBitFieldMatrix(std::cerr, bitFieldMatrix);

		const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(taxaBlock);
		if (nTreesBlocks == 0)
			{
			std::cerr << "Expecting a file with at least 1 TREES block.\n";
			return 3;
			}
		unsigned totalTreeIndex = 0;
		for (unsigned treesBlockInd = 0; treesBlockInd < nTreesBlocks; ++treesBlockInd) 
		    {
            NCL_COULD_BE_CONST  NxsTreesBlock * treesBlock = nexusReader.GetTreesBlock(taxaBlock, treesBlockInd);
            assert(treesBlock);
            treesBlock->ProcessAllTrees();
            const unsigned nTreesThisBlock = treesBlock->GetNumTrees();
            for (unsigned treeInd = 0; treeInd < nTreesThisBlock; ++treeInd)
                {
        		const NxsFullTreeDescription & ftd = treesBlock->GetFullTreeDescription(treeInd);
        		if (ftd.AllEdgesHaveLengths())
        		    {
        		    NxsSimpleTree nclTree(ftd, 0, 0.0);
        		    }
        		else
        		    {
        		    std::cerr << "Tree " << (1 + treeInd) << " of TREES block " << (1 + treesBlockInd) << " does not lengths for all of the edges. Skipping this tree.\n";
        		    }
                    

                }
            }		
		nexusReader.DeleteBlocksFromFactories();
	}
	catch(const NxsException &x)
		{
		std::cerr << "Error:\n " << x.msg << std::endl;
		if (x.line > 0 || x.pos > 0)
			std::cerr << "at line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << std::endl;
		return 2;
		}
	catch(...)
		{
		std::cerr << "Exiting with an unknown error" << std::endl;
		return 1;
		}
	return 0;
}
