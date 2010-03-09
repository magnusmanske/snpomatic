#ifndef __SNPOMATIC_H__
#define __SNPOMATIC_H__

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h> 
#include <time.h>
#include <getopt.h>
#include <string.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std ;

// __________________________________________________________________________________________________________________________________________________________________ DEFINES

// DEBUGGING: set to 1 to get status information over STDOUT
#define DEBUGGING 1

// READ_CHAR_BUFFER: Buffer size for reading a FASTA/FASTQ line
#define READ_CHAR_BUFFER 1024

// MINQUAL: Make bases in output lowercase if their quality falls below this threshold
#define MINQUAL 20

// MAXLINES: Max number of lines in alignment
#define MAXLINES 70

// Misc stuff
#define IUPAC_A 1
#define IUPAC_C 2
#define IUPAC_G 4
#define IUPAC_T 8
#define IUPAC_ALL 15

// Macros (to save runtime)
#define MERGE_IUPAC(_a,_b) IUPAC2char[char2IUPAC[_a]|char2IUPAC[_b]]
#define twotoone(__lane,__pos) ((__pos)*MAXLINES+(__lane))

#define REVERSE_COMPLEMENT_CHARP(__start) { \
		char *__c , *__d , __e ; \
		for ( __c = __start ; *__c ; __c++ ) ; \
		for ( __c-- , __d = __start ; __d < __c ; __c-- , __d++ ) { \
			__e = base_rc[*__d] ; \
			*__d = base_rc[*__c] ; \
			*__c = __e ; \
		} \
		if ( __c == __d ) *__c = base_rc[*__c] ; \
}


// __________________________________________________________________________________________________________________________________________________________________ TYPEDEFS

class TPositionWithIndex ;
class TChromosome ;
class TCoverage ;

typedef unsigned char uchar ;
typedef unsigned int index2type ;
typedef vector <uint> vuint ;
typedef vector <TPositionWithIndex> pwi_vector ;
typedef vector <TChromosome> chrvec ;
typedef vector <TPositionWithIndex>::iterator TPWIi ;
typedef vector <string>::iterator TSVI ;
typedef vector <TCoverage> TVC ;

// __________________________________________________________________________________________________________________________________________________________________ GLOBAL VARIABLES


extern uchar char2num[256] , num2char[256] , to_lc[256] , to_uc[256] , char2IUPAC[256] , IUPAC2char[256] , isIUPAC[256] , base_rc[256] , IUPAC_variance[256] ;
extern int globally_allowed_mismatches ;

// __________________________________________________________________________________________________________________________________________________________________ GLOBAL FUNCTIONS

void Tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ") ;
void init_iupac () ;
bool matchIUPAC ( const char *c1 , const char *c2 ) ;
void load_all_chromosomes ( string filename , vector <TChromosome> &chrs , bool keep_original_snps = false ) ;
void incorporate_all_gff ( string filename , vector <TChromosome> &chrs ) ;
void incorporate_simple_snps ( string filename , vector <TChromosome> &chrs ) ;

// __________________________________________________________________________________________________________________________________________________________________ CLASSES

class TCoverage {
	public :
	TCoverage () { A = C = G = T = rc = 0 ; }
	uint A , C , G , T , rc ;
} ;

/// Stores and outputs the pileup.
class TAlignmentOutput {
	public :
	char *align ; ///< Two-dimensional alignment array for the bases.
	char *qalign ; ///< Same as *align, but for quality.
	vector <bool> align_full ; ///< Shortcut mechanism to speed up alignment of extreme depths.
	
	void init ( TChromosome *_chr ) ; ///< Initialization of align and qualign memory.

	void show_pileup ( FILE *pileup , bool snps_only ) ; ///< Outputs the alignment to the given file.
	void add_align ( const string &seq , const string &quality , uint pos , int chromosome ) ;
	TChromosome *chr ; ///< Pointer back to the chromosome.
} ;

/// Stores secondary index, position, chromosome, and reverse/complement information. Heavily used by TChromosomalIndices.
class TPositionWithIndex {
	public :
	index2type index2 ; ///< The secondary index.
	uint position ; ///< The position in the chromosome.
	bool reverse_complement ; ///< Is this reverse/complement?
	uint chromosome ; ///< Number of chromosome.
} ;


/// Calculates and stores the indices on the chromosomes.
class TChromosomalIndices {
	public:
	TChromosomalIndices ( int _index_length , int _index_length2 ) ;
	void run_reads ( string filename ) ;
	void run ( chrvec *_chrs ) ;
	void run_ends ( chrvec *_chrs , int min_length = 0 ) ;
	void run_regions ( chrvec *_chrs , string filename ) ;
	uint calculate_index ( const char *sequence_kmer ) ;
	index2type calculate_index2 ( const char *sequence_kmer ) ;
	void save_index ( string filename ) ;
	void load_index ( string filename ) ;
	void uniqueness ( string filename ) ;
	string decompile ( uint index , uint bases ) ;
	
	vector <pwi_vector> position ; ///< The complete index.
	int index_length , index_length2 , index_length_double , noshortcuts ;
	int max_snps_per_index ; ///< Maximal number of known SNPs per index (1+2); prevent excessive memory usage.
	chrvec *chrs ; ///< Pointer to vector with all the chromosomes (TChromosome).
	string index_file ;
	uint index_from , index_to ;
	uint index_steps ;

	private:
	int add_sequence ( char *begin , char *cur , int from , int chr , bool rc , int nos = 0 ) ;
	void add_index ( char *start , bool rc , uint pos , uchar chromosome ) ;
	void reassemble_reads ( FILE *file , uint i1 , index2type i2a , index2type i2b ) ;
	bool align_read2contigs ( char *read , vector <string> &vs ) ;
	
	char tmp[READ_CHAR_BUFFER] ;
} ;

/// Holds an individual chromosome.
class TChromosome {
	public:
	string name , sequence , original_sequence ;
	TAlignmentOutput ao ;
	TVC coverage ;
	vector <uint> uniqueness ;

	void read_chromosome ( string filename , string chromosome_name ) ;
	void incorporate_known_snps ( string filename ) ;
	void incorporate_known_snps_gff ( string filename ) ;
	void add_coverage ( const string &seq , uint pos , bool rc ) ;
} ;

/// Reads Solexa data from FASTA/FASTQ and finds perfect matches in the reference.
class TChromosomeAlign {
	public:
	TChromosomeAlign ( chrvec &_chrs , TChromosomalIndices &_index ) ;
	void init () ;
	void align_solexa_fasta ( string filename ) ;
	void align_solexa_fastq ( string filename ) ;
	void align_solexa_paired ( string filename , string filename2 , uint read_length_1 , uint fragment_length , uint range ) ;
	void align_solexa_fastq_variety ( string filename ) ;
	void show_pileup ( bool snps_only = false ) ;
	void align_contigs ( chrvec &contigs , TChromosomalIndices &contig_index , uint maxerr , string output_filename ) ;
	void add_nono_list ( string filename ) ;
	bool is_nono ( const string &s ) ;
	void dump_coverage () ;
	void finish_sqlite ( int fragment , int variance , int read_length ) ;
	void get_align_solexa_read ( const char *_seq , pwi_vector &results ) ;
	
	FILE *pileup , *binfile_no_match , *binfile_single_match , *binfile_multi_match , *snpsinreads ;
	FILE *binfile_iupac , *cigar , *wobble , *fragmentplot , *featuretable , *gffout , *coverage ;
	FILE *indelplot , *inversions , *sqlite , *sam , *spancontigs , *faceaway , *rpa ;
	bool using_bins , multimatch , singlematch , force_one_unique_match ;
	uint chop , wobblemax , single_read_length ;
	
	protected:
	
	int align_solexa_read ( char *_seq , const string &_quality , bool inv_compl = false ) ;
	void add_cigar_sqlite ( const string &sequence , const string &name , uint position , uint chromosome , char orientation ) ;
	virtual uint align_solexa_paired_read ( const string &name , const string &sequence , string quality , uint read_length_1 , uint fragment_length , uint range ) ;
	int paired_read_combine ( const pwi_vector &v1 , const pwi_vector &v2 , uint fragment_length , uint range , uint read_length_1 , const string &seq1 , const string &qual1 , const string &seq2 , const string &qual2 , bool wobbling = false , bool inverting = false , int order = 0 ) ;
	void write2bin ( const string &name , const string &sequence , const string &quality , int matches ) ;
	void run_wobble ( const string &seq1 , const string &seq2 , const string &qual1 , const string &qual2 , const pwi_vector &res1 , const pwi_vector &res2 , uint fragment_length ) ;
	void wobble4snps ( const pwi_vector &res , const string &seq , uint fragment_size , uint range , const string &anchor_seq , bool rc ) ;
	void add_coverage ( const string &seq , uint pos ) ;
	int try_matching_read ( const string & seq , string quality ) ;
	string get_ref_substr ( const TPositionWithIndex &p , int length ) ;
	void add_snpsinreads ( const string &seq , const string &quality , int chr , int pos , int readpart ) ;
	int count_snpsinreads ( const string &seq , int chr , int pos ) ;
	void add_sam ( const string &seq , const string &quality , int chr , int pos , int readpart , int matepos , int ins_size , bool rc , bool as_pair , bool mate_rc , bool unique_match ) ;
	void wobble_single_read ( char *seq , const string &_quality ) ;
	string matchIUPAC_tolerant ( const char *c1 , const char *c2 , int tolerance ) ;
	string generate_sqlite_paired_command ( int chr , string seq1 , int pos1 , string seq2 , int pos2 , char mode = ' ' , string add_keys = "" , string add_values = "" ) ;
	void add2sqlite_cache ( string s ) ;
	void rpa_out ( const TPositionWithIndex &a , const TPositionWithIndex &b , uint read_length_1 , const string &seq1 , const string &qual1 , const string &seq2 , const string &qual2 ) ;
	
	uint max_align ;
	string last_solexa_name ;
	vector <string> sqlite_out ;
	FILE *sqlite_cache ;
	string sqlite_cache_name ;
	chrvec *chrs ;
	TChromosomalIndices *index ;
	bool fasta , use_nono , sam_wrote ;
	pwi_vector wobble_results ;
	vector <string> nono ;
	uint fragment_range ;
	vector <int> fragment_stats ;
	bool snpsinreads_first ;
	bool sqlite_prefix_first ;
	bool rpa_header_written ;
	string sqlite_prefix ;
	string rgid ;
} ;

class prc_cache {
	public :
	TPositionWithIndex a , b ;
} ;


bool operator < ( const TPositionWithIndex &a , const TPositionWithIndex &b ) ;
bool tpi_full_compare ( const TPositionWithIndex a , const TPositionWithIndex b ) ;
string dor ( string s ) ;
string doc ( string s ) ;
string dorc ( string s ) ;

#endif
