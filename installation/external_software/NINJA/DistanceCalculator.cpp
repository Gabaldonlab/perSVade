/*
 * DistanceCalculator.cpp
 *
 *  Created on: Feb 13, 2016
 *      Author: michel
 */
#include "DistanceCalculator.hpp"
#include <iostream>

DistanceCalculator::DistanceCalculator (std::string** A /*alignment*/, AlphabetType alphType, CorrectionType corrType, int numberOfSequences, bool useSSE) :
bl45{
				  {0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
			      {0.0, 0, 1.31097856, 1.06573001, 1.26827829, 0.90471293, 1.05855446, 1.05232790, 0.76957444, 1.27579668, 0.96460409, 0.98717819, 1.05007594, 1.05464162, 1.19859874, 0.96740447, 0.70049019, 0.88006018, 1.09748548, 1.28141710, 0.80003850},
			      {0.0, 1.31097856, 0, 0.80108902, 0.95334071, 1.36011107, 0.63154377, 0.79101490, 1.15694899, 0.76115257, 1.45014917, 1.17792001, 0.39466107, 0.99880755, 1.13514340, 1.15432562, 1.05309036, 1.05010474, 1.03938321, 0.96321690, 1.20274751},
			      {0.0, 1.06573001, 0.80108902, 0, 0.48821721, 1.10567116, 0.81497020, 0.81017644, 0.74648741, 0.61876156, 1.17886558, 1.52003670, 0.80844267, 1.28890258, 1.16264109, 1.18228799, 0.67947568, 0.85365861, 1.68988558, 1.24297493, 1.55207513},
			      {0.0, 1.26827829, 0.95334071, 0.48821721, 0, 1.31581050, 0.76977847, 0.48207762, 0.88836175, 0.73636084, 1.76756333, 1.43574761, 0.76361291, 1.53386612, 1.74323672, 0.88634740, 0.80861404, 1.01590147, 1.59617804, 1.17404948, 1.46600946},
			      {0.0, 0.90471293, 1.36011107, 1.10567116, 1.31581050, 0, 1.38367893, 1.37553994, 1.26740695, 1.32361065, 1.26087264, 1.02417540, 1.37259631, 1.09416720, 0.98698208, 1.59321190, 0.91563878, 0.91304285, 1.80744143, 1.32944171, 0.83002214},
			      {0.0, 1.05855446, 0.63154377, 0.81497020, 0.76977847, 1.38367893, 0, 0.50694279, 1.17699648, 0.61459544, 1.17092829, 1.19833088, 0.63734107, 0.80649084, 1.83315144, 0.93206447, 0.85032169, 1.06830084, 1.05739353, 0.97990742, 1.54162503},
			      {0.0, 1.05232790, 0.79101490, 0.81017644, 0.48207762, 1.37553994, 0.50694279, 0, 1.17007322, 0.76978695, 1.46659942, 1.19128214, 0.63359215, 1.27269395, 1.44641491, 0.73542857, 0.84531998, 1.06201695, 1.32439599, 1.22734387, 1.53255698},
			      {0.0, 0.76957444, 1.15694899, 0.74648741, 0.88836175, 1.26740695, 1.17699648, 1.17007322, 0, 1.12590070, 1.70254155, 1.38293205, 1.16756929, 1.17264582, 1.33271035, 1.07564768, 0.77886828, 1.23287107, 0.96853965, 1.42479529, 1.41208067},
			      {0.0, 1.27579668, 0.76115257, 0.61876156, 0.73636084, 1.32361065, 0.61459544, 0.76978695, 1.12590070, 0, 1.41123246, 1.14630894, 0.96779528, 0.77147945, 1.10468029, 1.12334774, 1.02482926, 1.28754326, 1.27439749, 0.46868384, 1.47469999},
			      {0.0, 0.96460409, 1.45014917, 1.17886558, 1.76756333, 1.26087264, 1.17092829, 1.46659942, 1.70254155, 1.41123246, 0, 0.43335051, 1.46346092, 0.46296554, 0.66291968, 1.07010201, 1.23000200, 0.97348545, 0.96354620, 0.70872476, 0.35120011},
			      {0.0, 0.98717819, 1.17792001, 1.52003670, 1.43574761, 1.02417540, 1.19833088, 1.19128214, 1.38293205, 1.14630894, 0.43335051, 0, 1.49770950, 0.47380007, 0.53847312, 1.37979627, 1.58597231, 0.99626739, 0.98609554, 0.72531066, 0.57054219},
			      {0.0, 1.05007594, 0.39466107, 0.80844267, 0.76361291, 1.37259631, 0.63734107, 0.63359215, 1.16756929, 0.96779528, 1.46346092, 1.49770950, 0, 1.00797618, 1.44331961, 0.92459908, 1.06275728, 1.05974425, 1.04892430, 0.97205882, 1.21378822},
			      {0.0, 1.05464162, 0.99880755, 1.28890258, 1.53386612, 1.09416720, 0.80649084, 1.27269395, 1.17264582, 0.77147945, 0.46296554, 0.47380007, 1.00797618, 0, 0.72479754, 1.16998686, 1.34481214, 1.06435197, 1.05348497, 0.77487815, 0.60953285},
			      {0.0, 1.19859874, 1.13514340, 1.16264109, 1.74323672, 0.98698208, 1.83315144, 1.44641491, 1.33271035, 1.10468029, 0.66291968, 0.53847312, 1.44331961, 0.72479754, 0, 1.32968844, 1.21307373, 0.96008757, 0.47514255, 0.34948536, 0.69273324},
			      {0.0, 0.96740447, 1.15432562, 1.18228799, 0.88634740, 1.59321190, 0.93206447, 0.73542857, 1.07564768, 1.12334774, 1.07010201, 1.37979627, 0.92459908, 1.16998686, 1.32968844, 0, 0.97908742, 0.97631161, 1.21751652, 1.42156458, 1.40887880},
			      {0.0, 0.70049019, 1.05309036, 0.67947568, 0.80861404, 0.91563878, 0.85032169, 0.84531998, 0.77886828, 1.02482926, 1.23000200, 1.58597231, 1.06275728, 1.34481214, 1.21307373, 0.97908742, 0, 0.56109848, 1.76318885, 1.29689226, 1.02015839},
			      {0.0, 0.88006018, 1.05010474, 0.85365861, 1.01590147, 0.91304285, 1.06830084, 1.06201695, 1.23287107, 1.28754326, 0.97348545, 0.99626739, 1.05974425, 1.06435197, 0.96008757, 0.97631161, 0.56109848, 0, 1.39547634, 1.02642577, 0.80740466},
			      {0.0, 1.09748548, 1.03938321, 1.68988558, 1.59617804, 1.80744143, 1.05739353, 1.32439599, 0.96853965, 1.27439749, 0.96354620, 0.98609554, 1.04892430, 1.05348497, 0.47514255, 1.21751652, 1.76318885, 1.39547634, 0, 0.32000293, 1.26858915},
			      {0.0, 1.28141710, 0.96321690, 1.24297493, 1.17404948, 1.32944171, 0.97990742, 1.22734387, 1.42479529, 0.46868384, 0.70872476, 0.72531066, 0.97205882, 0.77487815, 0.34948536, 1.42156458, 1.29689226, 1.02642577, 0.32000293, 0, 0.93309543},
			      {0.0, 0.80003850, 1.20274751, 1.55207513, 1.46600946, 0.83002214, 1.54162503, 1.53255698, 1.41208067, 1.47469999, 0.35120011, 0.57054219, 1.21378822, 0.60953285, 0.69273324, 1.40887880, 1.02015839, 0.80740466, 1.26858915, 0.93309543, 0}
				}
//				      {0, 1.31097856, 1.06573001, 1.26827829, 0.90471293, 1.05855446, 1.05232790, 0.76957444, 1.27579668, 0.96460409, 0.98717819, 1.05007594, 1.05464162, 1.19859874, 0.96740447, 0.70049019, 0.88006018, 1.09748548, 1.28141710, 0.80003850},
//				      {1.31097856, 0, 0.80108902, 0.95334071, 1.36011107, 0.63154377, 0.79101490, 1.15694899, 0.76115257, 1.45014917, 1.17792001, 0.39466107, 0.99880755, 1.13514340, 1.15432562, 1.05309036, 1.05010474, 1.03938321, 0.96321690, 1.20274751},
//				      {1.06573001, 0.80108902, 0, 0.48821721, 1.10567116, 0.81497020, 0.81017644, 0.74648741, 0.61876156, 1.17886558, 1.52003670, 0.80844267, 1.28890258, 1.16264109, 1.18228799, 0.67947568, 0.85365861, 1.68988558, 1.24297493, 1.55207513},
//				      {1.26827829, 0.95334071, 0.48821721, 0, 1.31581050, 0.76977847, 0.48207762, 0.88836175, 0.73636084, 1.76756333, 1.43574761, 0.76361291, 1.53386612, 1.74323672, 0.88634740, 0.80861404, 1.01590147, 1.59617804, 1.17404948, 1.46600946},
//				      {0.90471293, 1.36011107, 1.10567116, 1.31581050, 0, 1.38367893, 1.37553994, 1.26740695, 1.32361065, 1.26087264, 1.02417540, 1.37259631, 1.09416720, 0.98698208, 1.59321190, 0.91563878, 0.91304285, 1.80744143, 1.32944171, 0.83002214},
//				      {1.05855446, 0.63154377, 0.81497020, 0.76977847, 1.38367893, 0, 0.50694279, 1.17699648, 0.61459544, 1.17092829, 1.19833088, 0.63734107, 0.80649084, 1.83315144, 0.93206447, 0.85032169, 1.06830084, 1.05739353, 0.97990742, 1.54162503},
//				      {1.05232790, 0.79101490, 0.81017644, 0.48207762, 1.37553994, 0.50694279, 0, 1.17007322, 0.76978695, 1.46659942, 1.19128214, 0.63359215, 1.27269395, 1.44641491, 0.73542857, 0.84531998, 1.06201695, 1.32439599, 1.22734387, 1.53255698},
//				      {0.76957444, 1.15694899, 0.74648741, 0.88836175, 1.26740695, 1.17699648, 1.17007322, 0, 1.12590070, 1.70254155, 1.38293205, 1.16756929, 1.17264582, 1.33271035, 1.07564768, 0.77886828, 1.23287107, 0.96853965, 1.42479529, 1.41208067},
//				      {1.27579668, 0.76115257, 0.61876156, 0.73636084, 1.32361065, 0.61459544, 0.76978695, 1.12590070, 0, 1.41123246, 1.14630894, 0.96779528, 0.77147945, 1.10468029, 1.12334774, 1.02482926, 1.28754326, 1.27439749, 0.46868384, 1.47469999},
//				      {0.96460409, 1.45014917, 1.17886558, 1.76756333, 1.26087264, 1.17092829, 1.46659942, 1.70254155, 1.41123246, 0, 0.43335051, 1.46346092, 0.46296554, 0.66291968, 1.07010201, 1.23000200, 0.97348545, 0.96354620, 0.70872476, 0.35120011},
//				      {0.98717819, 1.17792001, 1.52003670, 1.43574761, 1.02417540, 1.19833088, 1.19128214, 1.38293205, 1.14630894, 0.43335051, 0, 1.49770950, 0.47380007, 0.53847312, 1.37979627, 1.58597231, 0.99626739, 0.98609554, 0.72531066, 0.57054219},
//				      {1.05007594, 0.39466107, 0.80844267, 0.76361291, 1.37259631, 0.63734107, 0.63359215, 1.16756929, 0.96779528, 1.46346092, 1.49770950, 0, 1.00797618, 1.44331961, 0.92459908, 1.06275728, 1.05974425, 1.04892430, 0.97205882, 1.21378822},
//				      {1.05464162, 0.99880755, 1.28890258, 1.53386612, 1.09416720, 0.80649084, 1.27269395, 1.17264582, 0.77147945, 0.46296554, 0.47380007, 1.00797618, 0, 0.72479754, 1.16998686, 1.34481214, 1.06435197, 1.05348497, 0.77487815, 0.60953285},
//				      {1.19859874, 1.13514340, 1.16264109, 1.74323672, 0.98698208, 1.83315144, 1.44641491, 1.33271035, 1.10468029, 0.66291968, 0.53847312, 1.44331961, 0.72479754, 0, 1.32968844, 1.21307373, 0.96008757, 0.47514255, 0.34948536, 0.69273324},
//				      {0.96740447, 1.15432562, 1.18228799, 0.88634740, 1.59321190, 0.93206447, 0.73542857, 1.07564768, 1.12334774, 1.07010201, 1.37979627, 0.92459908, 1.16998686, 1.32968844, 0, 0.97908742, 0.97631161, 1.21751652, 1.42156458, 1.40887880},
//				      {0.70049019, 1.05309036, 0.67947568, 0.80861404, 0.91563878, 0.85032169, 0.84531998, 0.77886828, 1.02482926, 1.23000200, 1.58597231, 1.06275728, 1.34481214, 1.21307373, 0.97908742, 0, 0.56109848, 1.76318885, 1.29689226, 1.02015839},
//				      {0.88006018, 1.05010474, 0.85365861, 1.01590147, 0.91304285, 1.06830084, 1.06201695, 1.23287107, 1.28754326, 0.97348545, 0.99626739, 1.05974425, 1.06435197, 0.96008757, 0.97631161, 0.56109848, 0, 1.39547634, 1.02642577, 0.80740466},
//				      {1.09748548, 1.03938321, 1.68988558, 1.59617804, 1.80744143, 1.05739353, 1.32439599, 0.96853965, 1.27439749, 0.96354620, 0.98609554, 1.04892430, 1.05348497, 0.47514255, 1.21751652, 1.76318885, 1.39547634, 0, 0.32000293, 1.26858915},
//				      {1.28141710, 0.96321690, 1.24297493, 1.17404948, 1.32944171, 0.97990742, 1.22734387, 1.42479529, 0.46868384, 0.70872476, 0.72531066, 0.97205882, 0.77487815, 0.34948536, 1.42156458, 1.29689226, 1.02642577, 0.32000293, 0, 0.93309543},
//				      {0.80003850, 1.20274751, 1.55207513, 1.46600946, 0.83002214, 1.54162503, 1.53255698, 1.41208067, 1.47469999, 0.35120011, 0.57054219, 1.21378822, 0.60953285, 0.69273324, 1.40887880, 1.02015839, 0.80740466, 1.26858915, 0.93309543, 0}
//					}
{
	this->A  = A;
	this->corr_type = corrType;
	this->alph_type = alphType;
	this->dna_chars = "AGCT";
	this->aa_chars = "ARNDCQEGHILKMFPSTWYV";
	this->numberOfSequences = numberOfSequences;
	this->lengthOfSequences = A[0]->size();
	this->newCalculation = useSSE;

	if(this->newCalculation && this->alph_type == this->dna)
		convertAllDNA();
	else if(this->alph_type == this->amino){
		convertAllProtein();
		#ifdef TEST_DIFF
		generateProteinOriginalDict(this->protein_dict_original);
		#endif
	}

	if (this->corr_type == not_assigned) {
		if (this->alph_type == amino) {
			this->corr_type = FastTree;
		} else {
			this->corr_type = Kimura2;
		}
	}

	this->inv_alph = DistanceCalculator::getInverseAlphabet( this->alph_type==dna ? this->dna_chars : this->aa_chars, this->alph_type==dna ? 4 : 20 );
}


int* DistanceCalculator::getInverseAlphabet (std::string alph, int length) {
	int* inv_alph = new int[256];
	for (int i=0; i<256; i++)
		inv_alph[i]=-1;
	for (int i=0; i<length; i++)
		inv_alph[(int)alph[i]] = i;
	return inv_alph;
}
inline void DistanceCalculator::count128(register __m128i &seq1, register __m128i &seq2, register __m128i &gap1, register __m128i &gap2, register __m128i &tmp, register __m128i &tmp2, register __m128i &tmp3, register __m128i &count_transversions, register __m128i &count_transitions, register __m128i &count_gaps){
	/*
	 * Maps and their description:
	 *
	 * GAPS_COUNT_MASK (count for gaps, basically how many 0`s there are)
	 * 	0000	4
	 * 	0001	3
	 * 	0010	3
	 * 	0011	2
	 * 	0100	3
	 * 	0101	2
	 * 	0110	2
	 * 	0111	1
	 * 	1000	3
	 * 	1001	2
	 * 	1010	2
	 * 	1011	1
	 * 	1100	2
	 * 	1101	1
	 * 	1110	1
	 * 	1111	0
	 *
	 * 	DECOMPRESSED_GAPS (transform 1-bit representation to 2-bits)
	 *
	 * 	0000	00 00 00 00
	 * 	0001	00 00 00 11
	 * 	0010	00 00 11 00
	 * 	0011	00 00 11 11
	 * 	0100	00 11 00 00
	 * 	0101	00 11 00 11
	 * 	0110	00 11 11 00
	 * 	0111	00 11 11 11
	 * 	1000	11 00 00 00
	 * 	1001	11 00 00 11
	 * 	1010	11 00 11 00
	 * 	1011	11 00 11 11
	 * 	1100	11 11 00 00
	 * 	1101	11 11 00 11
	 * 	1110	11 11 11 00
	 * 	1111	11 11 11 11
	 *
	 * 	COUNTS_MASK (clear upper 4 bits)
	 *
	 * 	0000 1111
	 *
	 * 	TRANSITIONS_MASK (count 10`s)
	 *
	 * 	0000	0
	 * 	0001	0
	 * 	0010	1
	 * 	0011	0
	 * 	0100	0
	 * 	0101	0
	 * 	0110	1
	 * 	0111	0
	 * 	1000	1
	 * 	1001	1
	 * 	1010	2
	 * 	1011	1
	 * 	1100	0
	 * 	1101	0
	 * 	1110	1
	 * 	1111	0
	 *
	 * 	TRANSVERSIONS_MASK (count 01`s and 11`s)
	 *
	 * 	0000	0
	 * 	0001	1
	 * 	0010	0
	 * 	0011	1
	 * 	0100	1
	 * 	0101	2
	 * 	0110	1
	 * 	0111	2
	 * 	1000	0
	 * 	1001	1
	 * 	1010	0
	 * 	1011	1
	 * 	1100	1
	 * 	1101	2
	 * 	1110	1
	 * 	1111	2
	 */

	tmp = _mm_xor_si128(seq1, seq2);
	tmp2  = _mm_and_si128(gap1,gap2);

	tmp3 = _mm_shuffle_epi8(GAPS_COUNT_MASK, tmp2);

	count_gaps = _mm_add_epi8(count_gaps, tmp3);

	tmp3 = _mm_shuffle_epi8(DECOMPRESSED_GAPS, tmp2);



	tmp  = _mm_and_si128(tmp,tmp3);


	tmp2  = _mm_and_si128(tmp,COUNTS_MASK);

	tmp3 = _mm_srli_epi64(tmp,4);

	tmp3 = _mm_and_si128(tmp3,COUNTS_MASK);


	tmp = _mm_shuffle_epi8(TRANSITIONS_MASK, tmp2);

	count_transitions = _mm_add_epi8(count_transitions, tmp);

	tmp = _mm_shuffle_epi8(TRANSITIONS_MASK, tmp3);

	count_transitions = _mm_add_epi8(count_transitions, tmp);



	tmp = _mm_shuffle_epi8(TRANSVERSIONS_MASK, tmp2);

	count_transversions = _mm_add_epi8(count_transversions, tmp);

	tmp = _mm_shuffle_epi8(TRANSVERSIONS_MASK, tmp3);

	count_transversions = _mm_add_epi8(count_transversions, tmp);


}

// RMH
// I got sick and tired of the inconsistent results
// when casting m128i to int arrays so I bit the bullet
// and just did it the right(?) way. 
//
// There are two functions here.  The first one,
// print_m128i_bits() prints the register in big
// endian order (MSB to LSB).  This not the natural
// order of the register but makes more sense for 
// standered bit shifting functions ( e.g. right shift
// does actually shift things right ).
//
// The print_m128i_bits_le() function prints the 
// register in little endian order (LSB to MSB).
//
void print_m128i_bits( __m128i value, char * prefix ) {
  alignas(16) uint32_t arr[4];

  fprintf(stderr, prefix);
  fprintf(stderr, " msb[127] ");
  _mm_store_si128((__m128i*)arr, value);
  for ( int m=3; m >= 0; m-- ) 
      for ( int k = 0; k < 32; k++ ) 
        fprintf(stderr,"%d",(arr[m] & (1<<(31-k)))>>(31-k));
  fprintf(stderr," lsb[ 0 ]\n");
}

void print_m128i_bits_le( __m128i value, char * prefix ) {
  alignas(16) uint32_t arr[4];

  fprintf(stderr, prefix);
  fprintf(stderr, " lsb[ 0 ] ");
  _mm_store_si128((__m128i*)arr, value);
  // One little, two little, three little endians
  for ( int m=0; m < 4; m++ ) 
    for ( int k = 0; k < 32; k++ ) 
      fprintf(stderr,"%d",(arr[m] & (1<<(31-k)))>>(31-k));
  fprintf(stderr," msb[127]\n");
}

// Bit Counting
// Borrowed quickly from https://stackoverflow.com/questions/17354971/fast-counting-the-number-of-set-bits-in-m128i-register
// to make some progress.  This could be streamlined for SSSE3/AMD if there is a desire to support specific processors or
// SSE versions.  Also, I am not sure if somewhere in this codebase a basic bit counter is already coded.
static const __m128i popcount_mask = _mm_set1_epi8(0x0F);
static const __m128i popcount_table = _mm_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
static inline __m128i popcnt8(__m128i n) {
  const __m128i pcnt0 = _mm_shuffle_epi8(popcount_table, _mm_and_si128(n, popcount_mask));
  const __m128i pcnt1 = _mm_shuffle_epi8(popcount_table, _mm_and_si128(_mm_srli_epi16(n, 4), popcount_mask));
  return _mm_add_epi8(pcnt0, pcnt1);
}
static inline __m128i popcnt64(__m128i n) {
  const __m128i cnt8 = popcnt8(n);
  return _mm_sad_epu8(cnt8, _mm_setzero_si128());
}
static inline int popcnt128(__m128i n) {
  const __m128i cnt64 = popcnt64(n);
  const __m128i cnt64_hi = _mm_unpackhi_epi64(cnt64, cnt64);
  const __m128i cnt128 = _mm_add_epi32(cnt64, cnt64_hi);
  return _mm_cvtsi128_si32(cnt128);
}
// END RMH



double DistanceCalculator::newCalcDNA(int a, int b){
	/*
	 * A = 00
	 * C = 01
	 * G = 10
	 * T = 11
	 *
	 *
	 * Explanation of what`s coded below:
	 *
	 *
	 *
	 * TMP = XOR (Seq1, Seq2)
	 *
	 * After the XOR we can see the following pattern
	 *
	 * Transition = 10
	 * Transversion = 11,01
	 *
	 * However, there might be some misleading counts because of gaps (currently encoded the same way as A, which is 00), which forces us to use
	 * a mask which accounts for place where there are gaps and places where there are not. To account for that we do the following:
	 *
	 * TMP2 = AND(Gap1, Gap2)
	 *
	 * We can now use a count process based on SHUFFLE to figure out how many gaps there are:
	 *
	 * TMP3 = SHUFFLE(GAPS_COUNT_MASK, TMP2)
	 *
	 * COUNTS_GAP = ADD(COUNTS_GAP, TMP3)
	 *
	 * Now we have added the 8-bit count bins to COUNTS_GAP. We can then proceed to calculate the transitions and transversions:
	 *
	 * TMP3 = SHUFFLE(DECOMPRESSED_GAPS, TPM2) //we use the 0000 - 1 bit representation per gap, and map it to a 2 bit representation of the gaps, with 11 as gap, and 00 as not-gap
	 *
	 * Now that we have the right mask we AND the mask and the result of the XOR:
	 *
	 * TMP = AND(TMP, TMP3)
	 *
	 * With the right vector in hands we will calculate transitions using the SHUFFLE count approach:
	 *
	 * TMP2 = AND(TMP, COUNTS_MASK) //the mask is 0000 1111 , and it basically removes the higher 4 bits
	 *
	 * TMP3 = SHIFT_RIGHT(TMP, 4)
	 *
	 * TMP3 = AND(TMP3, COUNTS_MASK) // same here, but then we have the higher 4 bits, originally, in the lower 4 bits now
	 *
	 * Finally get the transitions count:
	 *
	 * TMP = SHUFFLE(TRANSITIONS_MASK, TMP2) //get the counts
	 *
	 * TRANSITIONS_COUNT = ADD(TRANSITIONS_COUNT, TMP) //add to the transitions vector
	 *
	 * TMP = SHUFFLE(TRANSITIONS_MASK, TMP3) //get the counts
	 *
	 * TRANSITIONS_COUNT = ADD(TRANSITIONS_COUNT, TMP) //add to the transitions vector
	 *
	 *
	 * And then the same happens to transversions:
	 *
	 * TMP = SHUFFLE(TRANSVERSIONS_MASK, TMP2) //get the counts
	 *
	 * TRANSITIONS_COUNT = ADD(TRANSITIONS_COUNT, TMP) //add to the transitions vector
	 *
	 * TMP = SHUFFLE(TRANSVERSIONS_MASK, TMP3) //get the counts
	 *
	 * TRANSITIONS_COUNT = ADD(TRANSITIONS_COUNT, TMP) //add to the transitions vector
	 *
	 *
	 * And that`s it. We can calculate 64 characters in 17 operations. There is also a gather/extract part after all vectors which is 5 instructions
	 * for each count, which equals 15.
	 *
	 */

	register __m128i seq1;
	register __m128i seq2;
	register __m128i gap1;
	register __m128i gap2;
	register __m128i tmp;
	register __m128i tmp2;
	register __m128i tmp3;
	register __m128i counts_transversions;
	register __m128i counts_gaps;
	register __m128i counts_transitions;

	int numOfInts = ceil((float)this->lengthOfSequences/16.0);
	if(numOfInts % 4 != 0)
		numOfInts += 4 - (numOfInts % 4);

	int length = this->lengthOfSequences;

	const unsigned int* Achar = this->convertedSequences[a];
	const unsigned int* Bchar = this->convertedSequences[b];

	const unsigned int* Agap = this->gapInTheSequences[a];
	const unsigned int* Bgap = this->gapInTheSequences[b];

	int num_transversions = 0;

	int num_transitions = 0;

	int gaps = 0;

	int i = 0;

        // RMH: For gap initiation counting
        int inGap = 0;
        int gap1Cnt = 0;
        int gap2Cnt = 0;
 
	while(i < numOfInts){
		counts_transversions = x128;
		counts_gaps= x128;
		counts_transitions = x128;

		//TODO: review number of max iterations
		for(int j = 0;i<numOfInts && j < 31;i += 4){ //a maximum of 32 vectors allowed not to overflow things

				seq1 = *(__m128i*)&Achar[i];
				seq2 = *(__m128i*)&Bchar[i];

				gap1 = *(__m128i*)&Agap[i];
				gap2 = *(__m128i*)&Bgap[i];

                                // RMH: Brute-force gap initiation counting.  This is not
                                //      as elegant as an SSE implementation but it does
                                //      serve as a benchmark for improvement should one
                                //      be implemented.  This is only necessary ( for now )
                                //      if we are calculating the "MismatchesOneGap" 
                                //      distance metric.
                                if (this->corr_type == MismatchesOneGap)
                                  for ( int k = 0; k < 4; k++ ) // 32 bit words
                                    for ( int l = 3; l >= 0; l-- ) // 4 bytes each
                                      for ( int m = 3; m >= 0; m-- ){ // Only 4 bits of data used
                                        unsigned int bitmask = 1<<((l*8)+m);
                                        unsigned int Abit = Agap[i+k] & bitmask;
                                        unsigned int Bbit = Bgap[i+k] & bitmask;
                                        if ( Abit && Bbit ) 
                                        {
                                           inGap = 0;
                                        // Gap in B
                                        }else if ( Abit && !Bbit )
                                        {
                                          if ( inGap != 1 ) {
                                            gap1Cnt++;
                                          }  
                                          inGap = 1;
                                        // Gap in A
                                        }else if ( !Abit && Bbit )
                                        {
                                          if ( inGap != 2 ) {
                                            gap2Cnt++;
                                          }  
                                          inGap = 2;
                                        }
                                      }

			      count128(seq1,seq2,gap1, gap2, tmp,tmp2,tmp3,counts_transversions,counts_transitions, counts_gaps);
		 	      j+=4;
		}

		/*gather transversion counts*/

		counts_transversions = _mm_xor_si128(counts_transversions, x128);

		counts_transversions = _mm_sad_epu8 (counts_transversions, zero);
		tmp = _mm_shuffle_epi32(counts_transversions, _MM_SHUFFLE(1, 1, 1, 2));
		counts_transversions = _mm_add_epi16(counts_transversions, tmp);

		num_transversions +=  _mm_extract_epi16(counts_transversions, 0);


		/*gather transition counts*/

		counts_transitions = _mm_xor_si128(counts_transitions, x128);

		counts_transitions = _mm_sad_epu8 (counts_transitions, zero);
		tmp = _mm_shuffle_epi32(counts_transitions, _MM_SHUFFLE(1, 1, 1, 2));
		counts_transitions = _mm_add_epi16(counts_transitions, tmp);

		num_transitions += _mm_extract_epi16(counts_transitions, 0);

		/*gather gaps counts*/

		counts_gaps = _mm_xor_si128(counts_gaps, x128);

		counts_gaps = _mm_sad_epu8 (counts_gaps, zero);
		tmp = _mm_shuffle_epi32(counts_gaps, _MM_SHUFFLE(1, 1, 1, 2));
		counts_gaps = _mm_add_epi16(counts_gaps, tmp);

		gaps += _mm_extract_epi16(counts_gaps, 0);
               
	}

	length -= gaps;

	float dist = 0.0f;
	float maxscore =  (this->corr_type == none ? 1.0f : 3.0f);

	if(length == 0){
		dist = maxscore;
	}else{
		float p_f = (float)((float)num_transitions / (float)length);
		float q_f = (float)((float)num_transversions / (float)length);

		if (this->corr_type == JukesCantor)
			dist = (float)(-(0.75)*log((double)(1.0-(4.0/3.0)*(p_f+q_f))));
		else if (this->corr_type == Kimura2)
			dist = (float)(-0.5 * log(1.0 - 2*p_f - q_f) - 0.25 * log( 1.0-2*q_f ));
                else if (this->corr_type == MismatchesOneGap)// RMH 
                        dist = (float)(num_transitions + num_transversions + gap1Cnt + gap2Cnt) / (float)( length + gap1Cnt + gap2Cnt ); 
                else if (this->corr_type == none)
			dist = p_f + q_f;
	}

        //  NOTE: Here gaplen includes spacer gaps   
        //printf("RMH [ a=%d vs b=%d ]: transI = %d, transV = %d, len = %d (minus gaps), gaplen = %d, AGapCnt = %d, BGapCnt = %d, dist = %lf\n", a, b, num_transitions, num_transversions, length, gaps, gap2Cnt, gap1Cnt, dist); 

	double dist_d = (dist < maxscore ? dist : maxscore);

    return dist_d;

}


inline void DistanceCalculator::count128P(register __m128i &seq1, register __m128i &seq2,  register __m128i &gap1,  register __m128i &gap2, register __m128i &VALUES_0,  register __m128i &VALUES_1,  register __m128i &VALUES_2,  register __m128i &VALUES_3,  register __m128i &VALUES_4,  register __m128i &VALUES_5,  register __m128i &VALUES_6,  register __m128i &VALUES_7, register __m128i &sum, register __m128i &gap_count, register __m128i &tmp1, register __m128i &tmp2, int a, int b){

	tmp1 = _mm_min_epu8(seq1,seq2);
	tmp2 = _mm_max_epu8(seq1,seq2);

	gap1 = _mm_and_si128(gap1, gap2);

	// use tmp1, tmp2, gap1
	//tmp1 = min
	//tmp2 = max
	//gap1 = gaps

	seq2 = _mm_set1_epi8((char)29);

	seq2 = _mm_sub_epi8(seq2, tmp1); //29-min = index

	// use tmp1, tmp2, gap1, seq2
	// use tmp1, tmp2, gap1
	//tmp1 = min
	//tmp2 = max
	//gap1 = gaps
	//gap2 = index

	seq1 = _mm_slli_epi16(seq2, 8); //go left by 8 to erase higher 8

	//index1
	seq1 = _mm_srli_epi16(seq1, 8); //go back by 8

	gap2 = _mm_slli_epi16(tmp1, 8); //go left by 8 to erase higher 8

	//min1
	gap2 = _mm_srli_epi16(gap2, 8); //go back by 8

	seq1 = _mm_mullo_epi16(seq1, gap2); //multiply index1 and min1

	seq1 = _mm_srli_epi16(seq1, 1); //divide by 2

	seq1 = _mm_add_epi8(seq1, tmp2); //add seq1 to max

	seq1 = _mm_slli_epi16(seq1, 8); // go left by 8 to erase higher 8

	seq1 = _mm_srli_epi16(seq1, 8); // go back

	//used tmp1, tmp2, gap1, seq1
	//tmp1 = min
	//tmp2 = max
	//gap1 = gaps
	//seq1 = index1 final


	//min2
	gap2 = _mm_srli_epi16(tmp1, 8); //go right by 8

	__m128i min = tmp1; //new vector, store min

	tmp1 = _mm_cmpeq_epi8(tmp1, tmp2); //store equal proteins
	//after this I should not need max anymore

	//used tmp1, tmp2, gap1, gap2, seq1
	//tmp1 = equal
	//tmp2 = max
	//gap1 = gaps
	//gap2 = min2
	//seq1 = index1 final

	//index2
	seq2 = _mm_srli_epi16(seq2, 8); //go right by 8

	seq2 = _mm_mullo_epi16(gap2, seq2);

	seq2 = _mm_srli_epi16(seq2, 1); //divide by 2

	seq2 = _mm_slli_epi16(seq2, 8); //go left by 8

	seq2 = _mm_add_epi8(seq2, tmp2); // add max to it

	seq2 = _mm_srli_epi16(seq2, 8); //go right by 8 to clear

	seq2 = _mm_slli_epi16(seq2, 8); //go back

	//final index
	seq1 = _mm_or_si128(seq2, seq1);

	//used tmp1, tmp2, gap1, seq1
	//tmp1 = equal
	//tmp2 = max
	//gap1 = gaps
	//seq1 = index final
	//empty seq2, gap2

	//let`s make a srli_si128 that shifts bytes while all the other srli shift bits? sure, why not?
	gap2 = _mm_srli_epi16(seq1, 4);

	seq2 = _mm_set1_epi8(15);

	gap2 = _mm_and_si128(gap2, seq2);


	seq2 = _mm_and_si128(seq1, seq2);

	//use tmp1, tmp2, gap1, gap2, seq1, seq2
	//tmp1 = equal
	//tmp2 = max
	//gap1 = gaps
	//seq2 = used for next step, lower 4 bits
	//gap2 = used for next step, higher 4 bits
	//seq1 is free, I should not need the original index anymore

	__m128i sum_aux = _mm_setzero_si128(); //create new temp variable and set it to 0

	//for all 8 vectors

	seq1 = _mm_set1_epi8(0);

	seq1 = _mm_cmpeq_epi8(gap2, seq1);

	__m128i tmp3 = _mm_shuffle_epi8(VALUES_0, seq2);

	seq1 = _mm_and_si128(seq1, tmp3);

	sum_aux = _mm_add_epi8(sum_aux, seq1);


	seq1 = _mm_set1_epi8(1);

	seq1 = _mm_cmpeq_epi8(gap2, seq1);

	tmp3 = _mm_shuffle_epi8(VALUES_1, seq2);

	seq1 = _mm_and_si128(seq1, tmp3);

	sum_aux = _mm_add_epi8(sum_aux, seq1);


	seq1 = _mm_set1_epi8(2);

	seq1 = _mm_cmpeq_epi8(gap2, seq1);

	tmp3 = _mm_shuffle_epi8(VALUES_2, seq2);

	seq1 = _mm_and_si128(seq1, tmp3);

	sum_aux = _mm_add_epi8(sum_aux, seq1);


	seq1 = _mm_set1_epi8(3);

	seq1 = _mm_cmpeq_epi8(gap2, seq1);

	tmp3 = _mm_shuffle_epi8(VALUES_3, seq2);

	seq1 = _mm_and_si128(seq1, tmp3);

	sum_aux = _mm_add_epi8(sum_aux, seq1);


	seq1 = _mm_set1_epi8(4);

	seq1 = _mm_cmpeq_epi8(gap2, seq1);

	tmp3 = _mm_shuffle_epi8(VALUES_4, seq2);

	seq1 = _mm_and_si128(seq1, tmp3);

	sum_aux = _mm_add_epi8(sum_aux, seq1);



	seq1 = _mm_set1_epi8(5);

	seq1 = _mm_cmpeq_epi8(gap2, seq1);

	tmp3 = _mm_shuffle_epi8(VALUES_5, seq2);

	seq1 = _mm_and_si128(seq1, tmp3);

	sum_aux = _mm_add_epi8(sum_aux, seq1);



	seq1 = _mm_set1_epi8(6);

	seq1 = _mm_cmpeq_epi8(gap2, seq1);

	tmp3 = _mm_shuffle_epi8(VALUES_6, seq2);

	seq1 = _mm_and_si128(seq1, tmp3);

	sum_aux = _mm_add_epi8(sum_aux, seq1);



	seq1 = _mm_set1_epi8(7);

	seq1 = _mm_cmpeq_epi8(gap2, seq1);

	tmp3 = _mm_shuffle_epi8(VALUES_7, seq2);

	seq1 = _mm_and_si128(seq1, tmp3);

	sum_aux = _mm_add_epi8(sum_aux, seq1);

	//use tmp1, gap1, sum_aux, min
	//tmp1 = equal proteins
	//gap1 = gaps



	seq2 = _mm_set_epi8(22, 30, 18, 16, 22, 20, 16, 15, 24, 20, 15, 26, 20, 20, 18, 16); //EQUAL_VALUES


	seq1 = _mm_shuffle_epi8(seq2, min);
	seq1 = _mm_and_si128(tmp1, seq1);
	seq2 = _mm_set1_epi8(255);
	//this did not work like I envisioned with a nand on itself
	tmp2 = _mm_xor_si128(seq2, tmp1); // not on equal proteins
	sum_aux = _mm_and_si128(sum_aux, tmp2);
	sum_aux = _mm_or_si128(sum_aux, seq1); //add equal values
	sum_aux = _mm_and_si128(sum_aux, gap1);
	sum = _mm_add_epi8(sum_aux, sum);

/* for (int i=0;i<16;i++){
			char min2 = *((char*)(&min)+i);
			char value = *((char*)(&sum_aux)+i);
			printf("%d\n",(int) value);
	}*/

	sum_aux = _mm_setzero_si128();
	gap2 = _mm_set1_epi8(1);

	seq1 = _mm_cmpeq_epi8(gap1, sum_aux);
	seq1 = _mm_and_si128(seq1,gap2);
	gap_count = _mm_add_epi8(seq1, gap_count);

}
double DistanceCalculator::newCalcProtein(int a, int b){

	register __m128i seq1;
	register __m128i seq2;
	register __m128i gap1;
	register __m128i gap2;
	register __m128i tmp1;
	register __m128i tmp2;
	//register __m128i tmp3;
	register __m128i distance;
	register __m128i counts_gaps;



	int numOfInts = ceil((float)this->lengthOfSequences/4.0);
	if(numOfInts % 4 != 0)
		numOfInts += 4 - (numOfInts % 4);

	int length = this->lengthOfSequences;

	const unsigned int* Achar = this->convertedSequences[a];
	const unsigned int* Bchar = this->convertedSequences[b];

	const unsigned int* Agap = this->gapInTheSequences[a];
	const unsigned int* Bgap = this->gapInTheSequences[b];

	int sum = 0;
	int gaps = this->additionalGaps;

	int i = 0;

	while(i < numOfInts){ //a maximum of 8 vectors allowed not to overflow things
		distance = x128;
		counts_gaps= x128;

		for(int j = 0;i<numOfInts && j < 8; j++){

				seq1 = *(__m128i*)&Achar[i];
				seq2 = *(__m128i*)&Bchar[i];

				gap1 = *(__m128i*)&Agap[i];
				gap2 = *(__m128i*)&Bgap[i];

				count128P(seq1, seq2, gap1, gap2, this->VALUES_0, this->VALUES_1, this->VALUES_2, this->VALUES_3, this->VALUES_4, this->VALUES_5, this->VALUES_6, this->VALUES_7, distance, counts_gaps, tmp1, tmp2, a, b);

				i+=4;
		}

		/*gather distance*/

		distance = _mm_xor_si128(distance, x128);

		distance = _mm_sad_epu8 (distance, zero);
		tmp1 = _mm_shuffle_epi32(distance, _MM_SHUFFLE(1, 1, 1, 2));
		distance = _mm_add_epi16(distance, tmp1);

		sum +=  _mm_extract_epi16(distance, 0);


		/*gather gaps counts*/

		counts_gaps = _mm_xor_si128(counts_gaps, x128);

		counts_gaps = _mm_sad_epu8 (counts_gaps, zero);
		tmp1 = _mm_shuffle_epi32(counts_gaps, _MM_SHUFFLE(1, 1, 1, 2));
		counts_gaps = _mm_add_epi16(counts_gaps, tmp1);

		gaps += _mm_extract_epi16(counts_gaps, 0);
	}

	int relevant = length-gaps;

	if (relevant == 0)
		return -1; //max value

    return (sum/(relevant));
}
double DistanceCalculator::testDifferenceCluster(int a, int b){ //for test purposes only

	static const int bl62[20][20] = {
			{16, 6, 4, 4, 8, 6, 6, 8, 4, 6, 6, 6, 6, 4, 6, 10, 8, 2, 4, 8},
			 {6, 18, 8, 4, 2, 10, 8, 4, 8, 2, 4, 12, 6, 2, 4, 6, 6, 2, 4, 2},
			 {4, 8, 20, 10, 2, 8, 8, 8, 10, 2, 2, 8, 4, 2, 4, 10, 8, 0, 4, 2},
			 {4, 4, 10, 20, 2, 8, 12, 6, 6, 2, 0, 6, 2, 2, 6, 8, 6, 0, 2, 2},
			 {8, 2, 2, 2, 26, 2, 0, 2, 2, 6, 6, 2, 6, 4, 2, 6, 6, 4, 4, 6},
			 {6, 10, 8, 8, 2, 18, 12, 4, 8, 2, 4, 10, 8, 2, 6, 8, 6, 4, 6, 4},
			 {6, 8, 8, 12, 0, 12, 18, 4, 8, 2, 2, 10, 4, 2, 6, 8, 6, 2, 4, 4},
			 {8, 4, 8, 6, 2, 4, 4, 20, 4, 0, 0, 4, 2, 2, 4, 8, 4, 4, 2, 2},
			 {4, 8, 10, 6, 2, 8, 8, 4, 24, 2, 2, 6, 4, 6, 4, 6, 4, 4, 12, 2},
			 {6, 2, 2, 2, 6, 2, 2, 0, 2, 16, 12, 2, 10, 8, 2, 4, 6, 2, 6, 14},
			 {6, 4, 2, 0, 6, 4, 2, 0, 2, 12, 16, 4, 12, 8, 2, 4, 6, 4, 6, 10},
			 {6, 12, 8, 6, 2, 10, 10, 4, 6, 2, 4, 18, 6, 2, 6, 8, 6, 2, 4, 4},
			 {6, 6, 4, 2, 6, 8, 4, 2, 4, 10, 12, 6, 18, 8, 4, 6, 6, 6, 6, 10},
			 {4, 2, 2, 2, 4, 2, 2, 2, 6, 8, 8, 2, 8, 20, 0, 4, 4, 10, 14, 6},
			 {6, 4, 4, 6, 2, 6, 6, 4, 4, 2, 2, 6, 4, 0, 22, 6, 6, 0, 2, 4},
			 {10, 6, 10, 8, 6, 8, 8, 8, 6, 4, 4, 8, 6, 4, 6, 16, 10, 2, 4, 4},
			 {8, 6, 8, 6, 6, 6, 6, 4, 4, 6, 6, 6, 6, 4, 6, 10, 18, 4, 4, 8},
			 {2, 2, 0, 0, 4, 4, 2, 4, 4, 2, 4, 2, 6, 10, 0, 2, 4, 30, 12, 2},
			 {4, 4, 4, 2, 4, 6, 4, 2, 12, 6, 6, 4, 6, 14, 2, 4, 4, 12, 22, 6},
			 {8, 2, 2, 2, 6, 4, 4, 2, 2, 14, 10, 4, 10, 6, 4, 4, 8, 2, 6, 16}
	};
/*	static const int bl62_clusterized[16][16] = {
			{16, 6, 4, 4, 8, 6, 8, 4, 7, 6, 4, 6, 10, 8, 2, 4},
			 {6, 16, 8, 5, 2, 10, 4, 7, 3, 6, 2, 5, 7, 6, 2, 4},
			 {4, 8, 20, 10, 2, 8, 8, 10, 2, 4, 2, 4, 10, 8, 0, 4},
			 {4, 5, 10, 20, 2, 11, 6, 6, 1, 2, 2, 6, 8, 6, 0, 2},
			 {8, 2, 2, 2, 26, 1, 2, 2, 6, 6, 4, 2, 6, 6, 4, 4},
			 {6, 10, 8, 11, 1, 16, 4, 8, 3, 6, 2, 6, 8, 6, 3, 5},
			 {8, 4, 8, 6, 2, 4, 20, 4, 1, 2, 2, 4, 8, 4, 4, 2},
			 {4, 7, 10, 6, 2, 8, 4, 24, 2, 4, 6, 4, 6, 4, 4, 12},
			 {7, 3, 2, 1, 6, 3, 1, 2, 14, 11, 8, 3, 4, 7, 3, 6},
			 {6, 6, 4, 2, 6, 6, 2, 4, 11, 18, 8, 4, 6, 6, 6, 6},
			 {4, 2, 2, 2, 4, 2, 2, 6, 8, 8, 20, 0, 4, 4, 10, 14},
			 {6, 5, 4, 6, 2, 6, 4, 4, 3, 4, 0, 22, 6, 6, 0, 2},
			 {10, 7, 10, 8, 6, 8, 8, 6, 4, 6, 4, 6, 16, 10, 2, 4},
			 {8, 6, 8, 6, 6, 6, 4, 4, 7, 6, 4, 6, 10, 18, 4, 4},
			 {2, 2, 0, 0, 4, 3, 4, 4, 3, 6, 10, 0, 2, 4, 30, 12},
			 {4, 4, 4, 2, 4, 5, 2, 12, 6, 6, 14, 2, 4, 4, 12, 22}
					};*/

	//considering frequencies on distant calculation for the clusters
	static const int bl62_clusterized[16][16] = {
			{16, 6, 4, 4, 8, 6, 8, 4, 7, 6, 4, 6, 10, 8, 2, 4},
			 {6, 18, 8, 4, 2, 11, 4, 8, 2, 4, 2, 4, 6, 6, 2, 4},
			 {4, 8, 20, 10, 2, 8, 8, 10, 2, 2, 2, 4, 10, 8, 0, 4},
			 {4, 4, 10, 20, 2, 10, 6, 6, 2, 0, 2, 6, 8, 6, 0, 2},
			 {8, 2, 2, 2, 26, 1, 2, 2, 6, 6, 4, 2, 6, 6, 4, 4},
			 {6, 11, 8, 10, 1, 15, 4, 7, 3, 3, 2, 6, 8, 6, 2, 5},
			 {8, 4, 8, 6, 2, 4, 20, 4, 1, 0, 2, 4, 8, 4, 4, 2},
			 {4, 8, 10, 6, 2, 7, 4, 24, 2, 2, 6, 4, 6, 4, 4, 12},
			 {7, 2, 2, 2, 6, 3, 1, 2, 15, 11, 7, 3, 4, 7, 2, 6},
			 {6, 4, 2, 0, 6, 3, 0, 2, 11, 16, 8, 2, 4, 6, 4, 6},
			 {4, 2, 2, 2, 4, 2, 2, 6, 7, 8, 20, 0, 4, 4, 10, 14},//
			 {6, 4, 4, 6, 2, 6, 4, 4, 3, 2, 0, 22, 6, 6, 0, 2},
			 {10, 6, 10, 8, 6, 8, 8, 6, 4, 4, 4, 6, 16, 10, 2, 4},
			 {8, 6, 8, 6, 6, 6, 4, 4, 7, 6, 4, 6, 10, 18, 4, 4},
			 {2, 2, 0, 0, 4, 2, 4, 4, 2, 4, 10, 0, 2, 4, 30, 12},
			 {4, 4, 4, 2, 4, 5, 2, 12, 6, 6, 14, 2, 4, 4, 12, 22}


	};

	float dist = 0.0f;
	int dist_2 = 0;

	if (this->corr_type != none && this->corr_type != FastTree){
		fprintf(stderr, "illegal choice of correction method; must be 'n' or 's'");
		Exception::critical();
	}

	int count = 0;
	int length = (int)this->A[0]->length();
	const char* Achar = this->A[a]->c_str();
	const char* Bchar = this->A[b]->c_str();

	for (int i=0; i<length; i++){
		if (this->inv_alph[(int)Achar[i]] >= 0 && this->inv_alph[(int)Bchar[i]] >= 0) { // both are characters in the core alphabet
			//dist += bl62[this->protein_dict_original[(int)Achar[i]]][this->protein_dict_original[(int)Bchar[i]]];
			//printf("%d\n", this->protein_dict[(int)Bchar[i]]);
			dist_2 += bl62_clusterized[ this->protein_dict[(int)Achar[i]]][this->protein_dict[(int)Bchar[i]]];
			//printf("%c:%d %c:%d-- ", Achar[i], this->protein_dict_original[(int)Achar[i]], Bchar[i], this->protein_dict[(int)Bchar[i]]);
			//printf("%d %d || ", bl62[this->protein_dict_original[(int)Achar[i]]][this->protein_dict_original[(int)Bchar[i]]], bl62_clusterized[ this->protein_dict[(int)Achar[i]]][this->protein_dict[(int)Bchar[i]]]);
			count++;
		}
	}

	if (count != 0){
		//printf("%f %f %f %f %f\n", dist, dist_2, (dist - dist_2), ((dist - dist_2)/dist), (dist - dist_2)/(float)count);
		printf("%d\n", dist_2);
	}
	return dist - dist_2;
}
double DistanceCalculator::getMaxScore(){
	return (this->corr_type == none ? 1.0f : 3.0f);
}

double DistanceCalculator::calc (int a, int b){
	float dist = 0.0f;
	float maxscore =  (this->corr_type == none ? 1.0f : 3.0f);

	if(this->newCalculation && this->alph_type == this->dna){
		return newCalcDNA(a,b);
	}else if(this->newCalculation && this->alph_type == this->amino){
		//TODO: I should use a dissimilarity matrix for this
		dist =  newCalcProtein(a,b);
		if (dist == -1)
			dist = maxscore;
		dist = dist < 0.91 ? (float)(-1.3*log((double)(1.0 - dist))) : maxscore;
		return 	(double)-(dist < maxscore ? dist : maxscore);		
	}



	if (this->alph_type == amino) {
		if (this->corr_type != none && this->corr_type != FastTree){
			fprintf(stderr, "illegal choice of correction method; must be 'n' or 's'");
			Exception::critical();
		}

		int count = 0;
		const char* Achar = this->A[a]->c_str();
		const char* Bchar = this->A[b]->c_str();

		float localDist;

		for (int i=0; i<this->lengthOfSequences; i++){
//			if (this->inv_alph[(int)Achar[i]] >= 0 && this->inv_alph[(int)Bchar[i]] >= 0) { // both are characters in the core alphabet
//				dist += this->bl45[ this->inv_alph[(int)Achar[i]]  ]  [ this->inv_alph[(int)Bchar[i]]  ] ;
//				count++;
//			}
			localDist = this->bl45[(int)Achar[i]][(int)Bchar[i]];
			if (localDist != 0){
				dist += localDist;
				count++;
			}
		}

		if (count==0) {
			dist = maxscore;
		} else {
			dist /= count;
			if (this->corr_type == FastTree)
				dist = dist < 0.91 ? (float)(-1.3*log((double)(1.0 - dist))) : maxscore;
		}
	} else {
		if (this->corr_type == FastTree){
			fprintf(stderr, "illegal choice of correction method; must be 'n', 'j', or 'k'");
			Exception::critical();
		}


		int p=0; //transitions
		int q=0; //transversion
		int count = 0;
		const char* Achar = this->A[a]->c_str();
		const char* Bchar = this->A[b]->c_str();
		for (int i=0; i<this->lengthOfSequences; i++) {
			if (this->inv_alph[(int)Achar[i]] >= 0 && this->inv_alph[(int)Bchar[i]] >= 0) { // both are characters in the core alphabet
				count++;
				if (this->inv_alph[(int)Achar[i]] != this->inv_alph[(int)Bchar[i]]) {
					//count transitions (A-G, C-T) and transversions (others)
					if ( (this->inv_alph[(int)Achar[i]] <2 && this->inv_alph[(int)Bchar[i]] < 2) || // both A and G, and not equal
							(this->inv_alph[(int)Achar[i]] >1 && this->inv_alph[(int)Bchar[i]] > 1) ) // both C and T, not equal
						p++;
					else
						q++;
				}
			}
		}

		if (count == 0) {
				dist = maxscore;
		} else {
			float p_f = (float)p / count;
			float q_f = (float)q / count;

			if ( p_f+q_f == 0)
				dist = 0;
			else if (this->corr_type == JukesCantor)
				dist = (float)(-(0.75)*log((double)(1.0-(4.0/3.0)*(p_f+q_f))));
			else if (this->corr_type == Kimura2){
				dist = (float)(-0.5 * log(1.0 - 2*p_f - q_f) - 0.25 * sqrt( 1.0-2*q_f ));
			}else if (this->corr_type == none)
				dist = p_f + q_f;
		}

	}

	double dist_d = (dist < maxscore ? dist : maxscore);
	return dist_d;
}

DistanceCalculator::~DistanceCalculator(){
	delete[](this->inv_alph);
	if (this->gapInTheSequences != NULL){
		for(int i=0;i<this->numberOfSequences;i++){
			if(this->gapInTheSequences[i] != NULL)
				delete(this->gapInTheSequences[i]);
		}
	}
	if (this->convertedSequences != NULL){
		for(int i=0;i<this->numberOfSequences;i++){
			if(this->convertedSequences[i] != NULL)
				delete(this->convertedSequences[i]);
		}
	}
}
void DistanceCalculator::getBitsDNA(char* seq, int* size, unsigned int *seqOut, unsigned int *gapOut){
	/*
	 * transform the sequence of files into 2-bit-packed characters
	 *
	 * A = 00
	 * C = 01
	 * G = 10
	 * T = 11
	 *
	 * A is ignored in the loop, because it would add 0 to the integer.
	 * If the size left is less than 16, the other bits are all set to 0 by default,
	 * therefore it will not be counted as transitions or transversions. It will not
	 * affect the length of the sequence as well because we assume all sequences have
	 * the same length and this value is stored before the conversion to bits.
	 *
	 * The gaps (represented in ASCII as '-' ) are stored in gapOut, in the following manner:
	 *
	 * In each byte, the higher four bits are 0, and the lower 4 bits are either 1 or 0 if there is a gap or not, respectively.
	 *
	 */
	*seqOut = 0x0;
	*gapOut = 0xF0F0F0F; // initialize higher 4 bits as 0, and lower 4 bits as one

	if(*size<=0){
		return;
	}

	const static int numCharPerElement = 16;

	const static int whereToGo[] = {27,19,11,3}; //lookup table to figure where should one add the gap value

	const unsigned int powersOfTwo[] = {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,
			65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432,
			67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648};

	int i;

	for(i=0;i<*size && i<numCharPerElement;i++){ //goes until the sequence ends or 16 characters are converted
		if(seq[i] == 'C'){ //01
			*seqOut += powersOfTwo[((numCharPerElement-i)*2)-2];
		}else if(seq[i] == 'G'){ //10
			*seqOut += powersOfTwo[((numCharPerElement-i)*2)-1];
		}else if(seq[i] == 'T'){ //11
			*seqOut += (powersOfTwo[((numCharPerElement-i)*2)-2] + powersOfTwo[((numCharPerElement-i)*2)-1]);
		}else if(seq[i] == '-'){ //gap, 1 in one of the 4 lower bits of a byte
			*gapOut -= powersOfTwo[whereToGo[i/4]-(i%4)]; //get the right place to add and and the powers of 2 correspondent
		}
	}

	*size -= i;

}
void DistanceCalculator::generateProteinClusterDict(int* protein_dictionary){

	//new frequencies clusters
	//[['A'], ['R'], ['N'], ['D'], ['C'], ['Q', 'E', 'K'], ['G'], ['H'], ['I', 'V', 'M'], ['L'], ['F'], ['P'], ['S'], ['T'], ['W'], ['Y']]
	char proteins[20] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'K', 'G', 'H', 'I', 'V', 'M', 'L', 'F', 'P', 'S', 'T', 'W', 'Y'};
	protein_dictionary[(int)proteins[0]] = 0;
	protein_dictionary[(int)proteins[1]] = 1;
	protein_dictionary[(int)proteins[2]] = 2;
	protein_dictionary[(int)proteins[3]] = 3;
	protein_dictionary[(int)proteins[4]] = 4;
	protein_dictionary[(int)proteins[5]] = 5;
	protein_dictionary[(int)proteins[6]] = 5;
	protein_dictionary[(int)proteins[7]] = 5;
	protein_dictionary[(int)proteins[8]] = 6;
	protein_dictionary[(int)proteins[9]] = 7;
	protein_dictionary[(int)proteins[10]] = 8;
	protein_dictionary[(int)proteins[11]] = 8;
	protein_dictionary[(int)proteins[12]] = 8;
	protein_dictionary[(int)proteins[13]] = 9;
	protein_dictionary[(int)proteins[14]] = 10;
	protein_dictionary[(int)proteins[15]] = 11;
	protein_dictionary[(int)proteins[16]] = 12;
	protein_dictionary[(int)proteins[17]] = 13;
	protein_dictionary[(int)proteins[18]] = 14;
	protein_dictionary[(int)proteins[19]] = 15;
}
void DistanceCalculator::generateProteinOriginalDict(int* protein_dictionary){
	/*
	 * Clusters: {A} {R, N} {D} {C} {Q, E} {G} {H} {I, L, K} {M} {F} {S} {T} {W} {Y} {V}
	 * A R N D C Q E G H I L K M F P   S T  W Y V
	 * [['ala'], ['arg', 'lys'], ['asn'], ['asp'], ['cys'], ['gln', 'glu'], ['gly'], ['his'], ['ile', 'val', 'leu'], ['met'], ['phe'], ['pro'], ['ser'], ['thr'], ['trp'], ['tyr']]
	 */
	static const char proteins[20] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
	protein_dictionary[(int)proteins[0]] = 0;
	protein_dictionary[(int)proteins[1]] = 1;
	protein_dictionary[(int)proteins[2]] = 2;
	protein_dictionary[(int)proteins[3]] = 3;
	protein_dictionary[(int)proteins[4]] = 4;
	protein_dictionary[(int)proteins[5]] = 5;
	protein_dictionary[(int)proteins[6]] = 6;
	protein_dictionary[(int)proteins[7]] = 7;
	protein_dictionary[(int)proteins[8]] = 8;
	protein_dictionary[(int)proteins[9]] = 9;
	protein_dictionary[(int)proteins[10]] = 10;
	protein_dictionary[(int)proteins[11]] = 11;
	protein_dictionary[(int)proteins[12]] = 12;
	protein_dictionary[(int)proteins[13]] = 13;
	protein_dictionary[(int)proteins[14]] = 14;
	protein_dictionary[(int)proteins[15]] = 15;
	protein_dictionary[(int)proteins[16]] = 16;
	protein_dictionary[(int)proteins[17]] = 17;
	protein_dictionary[(int)proteins[18]] = 18;
	protein_dictionary[(int)proteins[19]] = 19;
}
void DistanceCalculator::getBitsProteinClustered(char* seq, int* size, unsigned int *seqOut, unsigned int *gapOut){
	/*
	 * The sequence and gap is initialized to zero. If it is not a gap or an undertermined protein, I set the
	 * gap to 1`s and the sequence according to generateProteinClusterDict.
	 */

	*seqOut = 0x0;
	*gapOut = 0x0;

	if(*size<=0){
		return;
	}

	const static int numCharPerElement = 4;

	const unsigned int gapValues[] = {255, 65280, 16711680, 4278190080};

	int i;

	for(i=0;i<*size && i<numCharPerElement;i++){
		if(seq[i] != '-' && seq[i] != 'X'){
			*seqOut += (this->protein_dict[(int)seq[i]]) << (i*8);
			*gapOut += gapValues[i];
		}
	}
	*size -= i;
}
unsigned int* DistanceCalculator::getProteinDic (std::string alph, int length) {
	unsigned int* inv_alph = new unsigned int[256];
	for (int i=0; i<256; i++)
		inv_alph[i]= 0xFFFFFFFF;
	for (int i=0; i<length; i++)
		inv_alph[(int)alph.at(i)] = i;
	return inv_alph;
}

void DistanceCalculator::convertAllProtein(){

	if (this->newCalculation == false){
		//convert this to be used with the bl45 matrix
		static const char proteins[22] = {'-', 'X', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
		this->protein_dict[(int)proteins[0]] = (char)0;
		this->protein_dict[(int)proteins[1]] = (char)0;
		this->protein_dict[(int)proteins[2]] = (char)1;
		this->protein_dict[(int)proteins[3]] = (char)2;
		this->protein_dict[(int)proteins[4]] = (char)3;
		this->protein_dict[(int)proteins[5]] = (char)4;
		this->protein_dict[(int)proteins[6]] = (char)5;
		this->protein_dict[(int)proteins[7]] = (char)6;
		this->protein_dict[(int)proteins[8]] = (char)7;
		this->protein_dict[(int)proteins[9]] = (char)8;
		this->protein_dict[(int)proteins[10]] = (char)9;
		this->protein_dict[(int)proteins[11]] = (char)10;
		this->protein_dict[(int)proteins[12]] = (char)11;
		this->protein_dict[(int)proteins[13]] = (char)12;
		this->protein_dict[(int)proteins[14]] = (char)13;
		this->protein_dict[(int)proteins[15]] = (char)14;
		this->protein_dict[(int)proteins[16]] = (char)15;
		this->protein_dict[(int)proteins[17]] = (char)16;
		this->protein_dict[(int)proteins[18]] = (char)17;
		this->protein_dict[(int)proteins[19]] = (char)18;
		this->protein_dict[(int)proteins[20]] = (char)19;
		this->protein_dict[(int)proteins[21]] = (char)20;
		char* seq;
		for(int i=0;i<this->numberOfSequences;i++){
			seq = (char*)this->A[i]->c_str();
			for(int j=0;j<this->lengthOfSequences;j++){
				seq[j] = (char)this->protein_dict[(int)seq[j]];
			}
		}
		return;
	}

	//calculate (k*(n-1)*n/2 k = length n = number of sequences
	long int k = this->lengthOfSequences;
	long int n = this->numberOfSequences;
	long int res = (k*(n-1)*n)/(long int)2; //number of pairs


	this->x128 = _mm_set1_epi8((int8_t) -128);
	this->zero = _mm_set1_epi8((int8_t) 0x00);

	generateProteinClusterDict(this->protein_dict);

	//the values are set in a peculiar way, it is actually the inverse of the expected, right to left
	this->VALUES_0 =_mm_set_epi8(4 ,2 ,8 ,10 ,6 ,4 ,6 ,7 ,4 ,8 ,6 ,8 ,4 ,4 ,6 ,0);
	this->VALUES_1 =_mm_set_epi8(2 ,10 ,4 ,2 ,6 ,6 ,4 ,2 ,4 ,2 ,8 ,4 ,11 ,2 ,4 ,8);
	this->VALUES_2 =_mm_set_epi8(2 ,6 ,6 ,10 ,2 ,4 ,0 ,8 ,10 ,4 ,2 ,2 ,2 ,10 ,8 ,8);
	this->VALUES_3 =_mm_set_epi8(6 ,6 ,2 ,4 ,6 ,6 ,2 ,2 ,1 ,2 ,0 ,6 ,8 ,6 ,2 ,0);
	this->VALUES_4 =_mm_set_epi8(2 ,0 ,1 ,4 ,5 ,2 ,6 ,8 ,6 ,2 ,3 ,3 ,7 ,4 ,4 ,4);
	this->VALUES_5 =_mm_set_epi8(3 ,7 ,11 ,12 ,4 ,4 ,6 ,4 ,6 ,2 ,2 ,2 ,4 ,4 ,8 ,4);
	this->VALUES_6 =_mm_set_epi8(6 ,14 ,10 ,4 ,4 ,0 ,6 ,4 ,6 ,4 ,2 ,8 ,6 ,2 ,7 ,4);
	this->VALUES_7 =_mm_set_epi8(0 ,0 ,0 ,0 ,0 ,0 ,0 ,12 ,4 ,4 ,4 ,2 ,10 ,2 ,0 ,6);


	//TODO: make sure all of these are aligned, so I can load them faster
	this->convertedSequences = new unsigned int*[this->numberOfSequences];
	this->gapInTheSequences = new unsigned int*[this->numberOfSequences];

	int allocSize = ceil((float)this->lengthOfSequences/4.0);
	if(allocSize % 4 != 0) //min size of 128bits
		allocSize += 4 - (allocSize % 4);
	int sizeLeft;
	for(int i=0;i<this->numberOfSequences;i++){
		this->convertedSequences[i] = new unsigned int[allocSize]; //min of 128bits, no need to change
		this->gapInTheSequences[i] = new unsigned int[allocSize]; //min of 128bits, no need to change
		sizeLeft = this->lengthOfSequences;
		for(int j=0;j<allocSize;j++){
			getBitsProteinClustered((char*)(this->A[i]->c_str()+(this->lengthOfSequences-sizeLeft)),&sizeLeft, &(this->convertedSequences[i][j]), &(this->gapInTheSequences[i][j]));
		}

	}

	this->additionalGaps = this->lengthOfSequences - (allocSize*4);

//	if (this->A != NULL){
//		for(int i=0;i<this->numberOfSequences;i++){
//			if(this->A[i] != NULL)
//				delete(this->A[i]);
//		}
//	}

}
void DistanceCalculator::convertAllDNA(){
	
	long int k = this->lengthOfSequences;
	long int n = this->numberOfSequences;
	long int res = (k*(n-1)*n)/(long int)2; //number of pairs

	//fprintf(stdout, "%ld",res);

	this->x128 = _mm_set1_epi8((int8_t) -128);
	this->zero = _mm_set1_epi8((int8_t) 0x00);
	this->COUNTS_MASK = _mm_set1_epi8((int8_t) 0xF);

	this->GAPS_COUNT_MASK = _mm_set_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);

	this->DECOMPRESSED_GAPS = _mm_set_epi8(255, 252, 243, 240, 207, 204, 195, 192, 63, 60, 51, 48, 15, 12, 3, 0);

	this->TRANSITIONS_MASK = _mm_set_epi8(0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0);

	this->TRANSVERSIONS_MASK = _mm_set_epi8(2, 1, 2, 1, 1, 0, 1, 0, 2, 1, 2, 1, 1, 0, 1, 0);

	//TODO: make sure all of these are aligned, so I can load them faster
	this->convertedSequences = new unsigned int*[this->numberOfSequences];
	this->gapInTheSequences = new unsigned int*[this->numberOfSequences];

	int allocSize = ceil((float)this->lengthOfSequences/16.0);
	if(allocSize % 4 != 0) //min size of 128bits
		allocSize += 4 - (allocSize % 4);
	int sizeLeft;
	for(int i=0;i<this->numberOfSequences;i++){
		this->convertedSequences[i] = new unsigned int[allocSize];
		this->gapInTheSequences[i] = new unsigned int[allocSize];

		sizeLeft = this->lengthOfSequences;
		for(int j=0;j<allocSize;j++){
			getBitsDNA((char*)(this->A[i]->c_str()+(this->lengthOfSequences-sizeLeft)),&sizeLeft, 
                                    &(this->convertedSequences[i][j]), &(this->gapInTheSequences[i][j]));
		}

	}

//	if (this->A != NULL){
//		for(int i=0;i<this->numberOfSequences;i++){
//			if(this->A[i] != NULL)
//				delete(this->A[i]);
//		}
//	}
}
