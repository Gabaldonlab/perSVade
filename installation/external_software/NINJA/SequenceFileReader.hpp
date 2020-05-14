/*
 * SequenceFIleReader.hpp
 *
 *  Created on: Jan 28, 2016
 *      Author: michel
 */

#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "ExceptionHandler.hpp"


#ifndef Included_SequenceFileReader_H
#define Included_SequenceFileReader_H

class SequenceFileReader{
	public:
		enum AlphabetType {dna, amino, null};
		std::string** seqs;
		std::string** names;
		AlphabetType alphType;
		int numSeqs;
		enum fileType {fasta, stockholm};

		SequenceFileReader (std::string *filename, AlphabetType alphTypeIn);
		~SequenceFileReader();

		fileType filetype;

		std::string **getSeqs();
		std::string **getNames();
		AlphabetType getAlphType();
};

#endif
