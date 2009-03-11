#include "snpomatic.h"
#include <fstream>

int maxfrag = 0 ;

class TLeftovers {
	public :
	string seq ;
	int pos ; // Not used ATM
	bool used ;
	vector <uint> keys ;
	vector <uint> related ;
} ;

class TRegionScan : public TChromosomeAlign {
	public:
	TRegionScan ( chrvec &_chrs , TChromosomalIndices &_index ) : TChromosomeAlign ( _chrs , _index ) {} ;
	void assemble_leftovers ( int cnt ) ;
	vector <TLeftovers> leftovers ;
	
	protected:
	virtual uint align_solexa_paired_read ( const string &name , const string &sequence , string quality , uint read_length_1 , uint fragment_length , uint range ) ;
	int get_overlap ( TLeftovers &l1 , TLeftovers &l2 ) ;
	TLeftovers merge_reads ( TLeftovers &l1 , TLeftovers &l2 , int overlap ) ;
} ;

bool operator < ( const TLeftovers &a , const TLeftovers &b ) {
	return a.seq < b.seq ;
}

TLeftovers TRegionScan::merge_reads ( TLeftovers &l1 , TLeftovers &l2 , int overlap ) {
	TLeftovers l ;
	l.seq = l2.seq.substr ( 0 , overlap ) + l1.seq ;
	for ( int n = 0 ; n + 16 < l.seq.length() ; n++ ) {
		l.keys.push_back ( index->calculate_index2 ( l.seq.c_str() + n ) ) ;
	}
	l.used = false ;
	return l ;
}

int TRegionScan::get_overlap ( TLeftovers &l1 , TLeftovers &l2 ) {
	int a, b ;
	int o = -1 ;
	for ( a = 0 ; a < l1.keys.size() ; a++ ) {
		for ( b = a+1 ; b < l2.keys.size() ; b++ ) {
			if ( l1.keys[a] != l2.keys[b] ) continue ;
			int diff = b - a ;
			if ( o == -1 ) o = diff ;
			else if ( o != diff ) return -1 ; // More than one overlap
		}
	}
/*
	string s1 = l1.seq.substr ( o , l1.seq.length() - o ) ;
	string s2 = l2.seq ;
	
	if ( s1.length() < s2.length() ) s2 = s2.substr ( 0 , s1.length() ) ;
	else if ( s2.length() < s1.length() ) s1 = s1.substr ( 0 , s2.length() ) ;
	if ( s1 != s2 ) return -1 ; // No match
*/
	return o ;
}

void TRegionScan::assemble_leftovers ( int cnt ) {
	uint a , b ;
	vector <TLeftovers> lo ;
	
	cout << "Assembling, depth " << cnt << endl ;
	
	for ( a = 0 ; a+1 < leftovers.size() ; a++ ) {
		for ( b = a+1 ; b < leftovers.size() ; b++ ) {
			int o = get_overlap ( leftovers[a] , leftovers[b] ) ;
			if ( o < 0 ) continue ; // No overlap, or more than one
			leftovers[a].used = true ;
			leftovers[b].used = true ;
			lo.push_back ( merge_reads ( leftovers[a] , leftovers[b] , o ) ) ;
//			cout << leftovers[a].seq << " <=> " << leftovers[b].seq << " : " << o << endl ;
		}
	}
	if ( cnt <= 0 ) return ;
	
	for ( a = 0 ; a < leftovers.size() ; a++ ) {
		if ( leftovers[a].used ) continue ;
		lo.push_back ( leftovers[a] ) ;
	}
	
	// Add unique sequences as new source set
	sort ( lo.begin() , lo.end() ) ;
	leftovers.clear() ;
	for ( a = 0 ; a < lo.size() ; a++ ) {
		if ( a == 0 || lo[a-1].seq != lo[a].seq ) leftovers.push_back ( lo[a] ) ;
	}
	
	assemble_leftovers ( cnt - 1 ) ;
}

uint TRegionScan::align_solexa_paired_read ( const string &name , const string &sequence , string quality , uint read_length_1 , uint fragment_length , uint range ) {
	string seq1 = sequence.substr ( 0 , read_length_1 ) ;
	string seq2 = sequence.substr ( read_length_1 , sequence.length() - read_length_1 ) ;
	
	while ( quality.length() < sequence.length() ) quality += (char) 33+30 ;
	string qual1 = quality.substr ( 0 , read_length_1 ) ;
	string qual2 = quality.substr ( read_length_1 , sequence.length() - read_length_1 ) ;
	
	pwi_vector res1 , res2 ;
	get_align_solexa_read ( seq1.c_str() , res1 ) ;
	get_align_solexa_read ( seq2.c_str() , res2 ) ;
	if ( res1.size() + res2.size() == 0 ) return 0 ;

	sort ( res1.begin() , res1.end() , tpi_full_compare ) ;
	sort ( res2.begin() , res2.end() , tpi_full_compare ) ;

	int ret = 0 ;
	ret += paired_read_combine ( res1 , res2 , maxfrag/2 , maxfrag/2 , read_length_1 , seq1 , qual1 , seq2 , qual2 ) ;
	ret += paired_read_combine ( res2 , res1 , maxfrag/2 , maxfrag/2 , read_length_1 , seq2 , qual2 , seq1 , qual1 ) ;
	
	if ( ret == 0 && res1.size() + res2.size() == 1 ) {
		string seq ;
		int newpos , chr ;
		if ( res1.size() ) {
			chr = res1[0].chromosome ;
			if ( res1[0].reverse_complement ) { seq = seq2 ; newpos = res1[0].position - fragment_length ; }
			else { seq = dorc ( seq2 ) ; newpos = (int) res1[0].position + fragment_length ; }
		} else {
			chr = res2[0].chromosome ;
			if ( res2[0].reverse_complement ) { seq = seq1 ; newpos = (int) res2[0].position - fragment_length ; }
			else { seq = dorc ( seq1 ) ; newpos = res2[0].position + fragment_length ; }
		}
		if ( newpos > 0 || newpos + seq.length() < (*chrs)[chr].sequence.length() ) {
			TLeftovers l ;
			l.seq = seq ;
			l.pos = newpos ;
			l.used = false ;
			for ( int n = 0 ; n + 16 < seq.length() ; n++ ) {
				l.keys.push_back ( index->calculate_index2 ( seq.c_str() + n ) ) ;
			}
			leftovers.push_back ( l ) ;
		}
	}
	
	return ret ;
}

// __________________________________________________________________________________________________________________________________________________________________ MAIN

/// Shows possible parameters, and an (optional) error message.
void die_with_description ( const char *msg = NULL ) {
	if ( msg ) cerr << msg << endl ;
	cerr << "Parameters : \n" ;
	exit ( 1 ) ;
}

/// Well, guess what.
int main ( int argc , char **argv ) {
	string genome_file , gff_file , simple_snp_file , fastq_file , chromosome , name , coverage_file ;
	string bin_prefix , binmask ;
	uint mspi = 8 ;
	uint pair_length = 0 ;
	uint fragment_length = 0 ;
	uint variance = 0 ;
	uint index_length_1 = 10 ;
	uint index_length_2 = 16 ;
	int noshortcuts = false ;
	uint from = 0 , to = 0 ;

	// Parameters
	static struct option long_options[] = {
		{ "genome" , required_argument , 0 , 'g' } ,
		{ "fastq" , optional_argument , 0 , 'q' } ,
		{ "snps" , optional_argument , 0 , 's' } ,
		{ "gff" , optional_argument , 0 , 'f' } ,
		{ "pair" , optional_argument , 0 , 'i' } ,
		{ "fragment" , optional_argument , 0 , 't' } ,
		{ "variance" , optional_argument , 0 , 'v' } ,
		{ "coverage" , optional_argument , 0 , 'e' } ,
		{ "mspi" , optional_argument , 0 , 'm' } ,
		{ "chr" , optional_argument , 0 , 'c' } ,
		{ "from" , optional_argument , 0 , 'r' } ,
		{ "to" , optional_argument , 0 , 'o' } ,
		{ "name" , optional_argument , 0 , 'n' } ,
		{ "bins" , optional_argument , 0 , 'b' } ,
		{ "binmask" , optional_argument , 0 , 'z' } ,
		{ 0 , 0 , 0 , 0 }
	} ;

	int c ;
	int option_index = 0;
	while ( -1 != ( c = getopt_long (argc, argv, "g:a::q::s::f::p::b::r::w::n::",long_options, NULL) ) ) {
		
		switch ( c ) {
			case 0 : break ;
			case 'g' : genome_file = optarg ; break ;
			case 'q' : fastq_file = optarg ; break ;
			case 's' : simple_snp_file = optarg ; break ;
			case 'f' : gff_file = optarg ; break ;
			case 'c' : chromosome = optarg ; break ;
			case 'n' : name = optarg ; break ;
			case 'e' : coverage_file = optarg ; break ;
			case 'i' : pair_length = atoi ( optarg ) ; break ;
			case 't' : fragment_length = atoi ( optarg ) ; break ;
			case 'v' : variance = atoi ( optarg ) ; break ;
			case 'm' : mspi = atoi ( optarg ) ; break ;
			case 'r' : from = atoi ( optarg ) ; break ;
			case 'o' : to = atoi ( optarg ) ; break ;
			case 'b' : bin_prefix = optarg ; break ;
			case 'z' : binmask = optarg ; break ;
			default : die_with_description("Unknown option") ;
		}
	}
	binmask += "0000" ; // Ignore bins if shorter string was supplied

	if ( genome_file.empty() ) die_with_description("No genome file given!") ; // Need this
	if ( fastq_file.empty() ) die_with_description("No Solexa file given!") ; // Need one
	if ( chromosome.empty() ) die_with_description("No chromosome given!") ; // Need one
	if ( name.empty() ) die_with_description("No name given!") ; // Need one
	if ( from == 0 ) die_with_description("Needs --from!") ; // Need one
	if ( to == 0 ) die_with_description("Needs --to!") ; // Need one
	

	// Incorporating known SNPs and indexing chromosomes
	init_iupac () ;
	TChromosomalIndices index ( index_length_1 , index_length_2 ) ;
	index.noshortcuts = noshortcuts ;
	index.max_snps_per_index = mspi ;
	vector <TChromosome> chrs ;
	load_all_chromosomes ( genome_file , chrs ) ;
	
	if ( !gff_file.empty() ) incorporate_all_gff ( gff_file , chrs ) ;
	if ( !simple_snp_file.empty() ) incorporate_simple_snps ( simple_snp_file , chrs ) ;

	// Remove all but one chromosome
	uint a ;
	for ( a = 0 ; a < chrs.size() ; a++ ) {
		if ( chrs[a].name == chromosome ) continue ;
		for ( int b = a ; b+1 < chrs.size() ; b++ ) chrs[b] = chrs[b+1] ;
		chrs.pop_back() ;
		a-- ;
	}
	chrs[0].sequence = chrs[0].sequence.substr ( from-1 , to-from ) ;
	maxfrag = chrs[0].sequence.length() ;
	cout << "Keeping " << chrs[0].sequence.length() << " bases from " << chromosome << endl ;

	index.run ( &chrs ) ;
	
	// Now look at the SOlexa reads and do the thing
	TRegionScan ca ( chrs , index ) ;
	if ( !coverage_file.empty() ) ca.coverage = fopen ( coverage_file.c_str() , "w" ) ;

	ca.init () ;
/*
	if ( !bin_prefix.empty() ) {
		string ext = ".fasta" ;
		if ( !fastq_file.empty() ) ext = ".fastq" ;
		cout << "Using binmask " << binmask << endl ;
		if ( binmask[0] == '1' ) ca.binfile_no_match = fopen ( string ( bin_prefix + "_no_match" + ext ).c_str() , "w" ) ;
		if ( binmask[1] == '1' ) ca.binfile_single_match = fopen ( string ( bin_prefix + "_single_match" + ext ).c_str() , "w" ) ;
		if ( binmask[2] == '1' ) ca.binfile_multi_match = fopen ( string ( bin_prefix + "_multi_match" + ext ).c_str() , "w" ) ;
		if ( binmask[3] == '1' ) ca.binfile_iupac = fopen ( string ( bin_prefix + "_iupac" + ext ).c_str() , "w" ) ;
		ca.using_bins = true ;
	}
*/
	ca.align_solexa_paired ( fastq_file , pair_length , fragment_length , variance ) ;
//	ca.assemble_leftovers ( 2 ) ;
	
	string fasta_output = name + ".fasta" ;
	fstream filestr;
	filestr.open ( fasta_output.c_str() , fstream::out ); //| fstream::app);
	for ( a = 0 ; a < ca.leftovers.size() ; a++ ) {
		filestr << ">" << a << endl << ca.leftovers[a].seq << endl ;
	}
	filestr.close();

	ca.dump_coverage() ;
//	cout << chrs[0].sequence << endl ;
}
