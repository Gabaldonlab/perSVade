/*
 * SequenceFileReader.cpp
 *
 *  Created on: Jan 28, 2016
 *      Author: michel
 */

#include "SequenceFileReader.hpp"
#include <errno.h>
#include <algorithm>

SequenceFileReader::SequenceFileReader(std::string *filename, AlphabetType alphTypeIn){
	//char dna_chars[] = {'A','G','C','T','U'};
	FILE* inFile = NULL;
	long int size = 0;
	if(filename == NULL or filename->empty()){
		inFile = stdin;
		fseek(inFile,0,SEEK_END);
		size = ftell(inFile);
		fseek(inFile,0,SEEK_SET);
	}else{
		inFile = fopen(filename->c_str(),"r");
		if( inFile == NULL){
			Exception::criticalErrno(filename->c_str());
		}
		fseek(inFile,0,SEEK_END);
		size = ftell(inFile);
		fseek(inFile,0,SEEK_SET);
	}
	if (size==0) {
		fprintf(stderr,"No sequences to align. Quitting\n");
		Exception::critical();
	}
	char *x = new char[size];
	long int numRead = (long int)fread(x,1,size,inFile);
	if (size != numRead){
		Exception::criticalErrno((const char*)NULL);
	}
	fclose(inFile);

	std::vector<std::string> names;
	std::vector<std::string> seq;
	int begin = 0;
	int end = 0;
	size_t charSize = sizeof(char);
	this->filetype = fasta;
	//TODO: If I assume a constant length throughout the sequences, it`s easy to increase the performance of this read
	//reads x and format the names and sequences
	for(int i=0;i<size;i++){ //try to read as a fasta format file
		if (x[i] == '>'){
			i++;
			begin = i;
			while(x[i] != ' ' && x[i] != '\t' && x[i] != '\n')
				i++;
			end = i;
			std::string auxName;
			auxName.assign(x + begin*charSize,(end-begin)*charSize);
			names.push_back(auxName);
			i++;
			begin = i;
			while(i<size && x[i] != '>')
				i++;
			end = --i;
			auxName.assign(x + begin*charSize,(end-begin)*charSize);

			auxName.erase(remove_if(auxName.begin(),auxName.end(),isspace),auxName.end());
			seq.push_back(auxName);
		}else{
			this->filetype = stockholm;
			break;
			fprintf(stderr,"Wrong formatting. Only fasta files are currently accepted. Quitting\n");
			Exception::critical();
		}
	}
//	if(this->filetype == stockholm){ //TODO: fix and optimize
//		//try to read as stockholm format file
//		for(int i=0;i<size;i++){
//			if (x[i] == '#' || x[i] == '\n' || x[i] == '\t'){
//				if (x[i] != '\n' || x[i+1] == '\n' || x[i+1] == '#'){
//					while(x[i]!='\n'){
//						i++;
//					}
//				}
//			}else{
//				if(x[i]==' ' || x[i] == '\t' || x[i] == '\n') i++;
//				begin = i;
//				while(x[i] != ' ' && x[i] != '\t' && x[i] != '\n')
//					i++;
//				end = i;
//				std::string auxName;
//				auxName.assign(x + begin*charSize,(end-begin)*charSize);
//				if(auxName == "//")
//					break;
//				names.push_back(auxName);
//				while(x[i]!=' '&& x[i]!= '\t' && x[i] != '\n'){
//					i++;
//				}
//				i++;
//				begin = i;
//				while(i<size && x[i] != '\n')
//					i++;
//				end = --i;
//				auxName.assign(x + begin*charSize,(end-begin)*charSize);
//
//				auxName.erase(remove_if(auxName.begin(),auxName.end(),isspace),auxName.end());
//				seq.push_back(auxName);
//				i++;
//			}
//		}
//	}
	delete[](x);
	int numSeqs = names.size();
	this->numSeqs = numSeqs;
	std::string** seqNames = new std::string*[numSeqs];
	std::string** sequences = new std::string*[numSeqs];
	for (int i=0;i<numSeqs;i++){ //get strings into a strings array instead of a vector
		seqNames[i] = new std::string();
		seqNames[i]->assign(names[i]);
		sequences[i] = new std::string();
		sequences[i]->assign(seq[i]);
	}

	//int k = 0;
	if (alphTypeIn != null)
		this->alphType = alphTypeIn;
	else
		this->alphType = dna;


	//Detect if it is DNA or AMINO ACID
	if (this->alphType == null){
		for (int i=0; i<numSeqs; i++){
			for (int j=0; j<(signed)sequences[i]->size(); j++) {
				if (sequences[i]->at(j) != 'A' && sequences[i]->at(j) != 'G' && sequences[i]->at(j) != 'C' && sequences[i]->at(j) != 'T'){
					this->alphType = amino;
					break;
				}
			}
			if (this->alphType != null)
				break;
		}
		if (this->alphType == null)
			this->alphType = dna;
	}

	for (int i=0; i<numSeqs; i++){
		//bool good = true;
		std::transform(sequences[i]->begin(), sequences[i]->end(),sequences[i]->begin(), ::toupper);
		for (int j=0; j<(signed)sequences[i]->size(); j++) {
			if (sequences[i]->at(j) == '.'){
				(*sequences[i])[j] = '-';
				//good = false;
			}else{
				if (this->alphType == dna && sequences[i]->at(j) == 'U') //change U to T
					(*sequences[i])[j] = 'T';
			}
		}
	}
	this->names = seqNames;
	this->seqs = sequences;
}
SequenceFileReader::~SequenceFileReader(){
	/*
	if (this->names != NULL){
		for(int i=0;i<this->seqSize;i++){
			if (this->names[i] != NULL)
				delete(this->names[i]);
		}
		delete[](this->names);
	}
	*/

/*	for(int i=0;i<this->seqSize;i++)
		delete(this->seqs[i]);
	delete[](this->seqs);*/
}

std::string **SequenceFileReader::getSeqs(){
	return this->seqs;
}
std::string **SequenceFileReader::getNames(){
	return this->names;
}
SequenceFileReader::AlphabetType SequenceFileReader::getAlphType(){
	return this->alphType;
}
