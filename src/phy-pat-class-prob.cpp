#include<iostream>
#include "ncl/nxsmultiformat.h"



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
		const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(taxaBlock);
		if (nTreesBlocks != 1)
			{
			std::cerr << "Expecting a file with exactly 1 CHARACTERS/DATA block, but found " << nTreesBlocks << " in the file " << filename << ".\n";
			return 3;
			}
		NCL_COULD_BE_CONST  NxsTreesBlock * treesBlock = nexusReader.GetTreesBlock(taxaBlock, 0);
		
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
