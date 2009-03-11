#include "snpomatic.h"

// __________________________________________________________________________________________________________________________________________________________________ GLOBAL VARIABLES


uchar char2num[256] , num2char[256] , to_lc[256] , to_uc[256] , char2IUPAC[256] , IUPAC2char[256] , isIUPAC[256] , base_rc[256] , IUPAC_variance[256] ;
int globally_allowed_mismatches = 0 ;

// __________________________________________________________________________________________________________________________________________________________________ GLOBAL FUNCTIONS


/// Breaks up a string in a vector of strings at given sepatator
void Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters)
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

/// Initializes IUPAC lookup tables
void init_iupac () {
	int a ;
	for ( a = 0 ; a < 256 ; a++ ) char2IUPAC[a] = 0 ;
	char2IUPAC['A'] = IUPAC_A ;
	char2IUPAC['C'] = IUPAC_C ;
	char2IUPAC['G'] = IUPAC_G ;
	char2IUPAC['T'] = IUPAC_T ;
//	char2IUPAC['U'] = IUPAC_T ; // U <=> T
	char2IUPAC['R'] = IUPAC_A|IUPAC_G ; // R	Purine (A or G)
	char2IUPAC['Y'] = IUPAC_C|IUPAC_T ; // Y	Pyrimidine (C, T)
	char2IUPAC['M'] = IUPAC_C|IUPAC_A ; // M	C or A
	char2IUPAC['K'] = IUPAC_T|IUPAC_G ; // K	T, or G
	char2IUPAC['W'] = IUPAC_T|IUPAC_A ; // W	T, or A
	char2IUPAC['S'] = IUPAC_C|IUPAC_G ; // S	C or G
	char2IUPAC['B'] = IUPAC_C|IUPAC_T|IUPAC_G ; // B	C, T, or G (not A)
	char2IUPAC['D'] = IUPAC_T|IUPAC_A|IUPAC_G ; // D	A, T, or G (not C)
	char2IUPAC['H'] = IUPAC_T|IUPAC_A|IUPAC_C ; // H	A, T, or C (not G)
	char2IUPAC['V'] = IUPAC_A|IUPAC_C|IUPAC_G ; // V	A, C, or G (not T)
	char2IUPAC['N'] = IUPAC_A|IUPAC_C|IUPAC_G|IUPAC_T ; // N	Any base (A, C, G, T)
	for ( a = 0 ; a < 256 ; a++ ) IUPAC2char[char2IUPAC[a]] = a ;
	char2IUPAC['X'] = IUPAC_A|IUPAC_C|IUPAC_G|IUPAC_T ; // X	Any base (A, C, G, T)

	for ( a = 0 ; a < 256 ; a++ ) isIUPAC[a] = 1 ;
	isIUPAC['A'] = isIUPAC['C'] = isIUPAC['G'] = isIUPAC['T'] = 0 ;
	isIUPAC[IUPAC_A] = isIUPAC[IUPAC_C] = isIUPAC[IUPAC_G] = isIUPAC[IUPAC_T] = 0 ; // 1,2,4,8 - isIUPAC is safe for double use...
	for ( a = 0 ; a < 256 ; a++ ) base_rc[a] = ' ' ;
	base_rc['A'] = 'T' ;
	base_rc['T'] = 'A' ;
	base_rc['C'] = 'G' ;
	base_rc['G'] = 'C' ;
	base_rc['N'] = 'N' ;
	base_rc['X'] = 'X' ;
	for ( a = 0 ; a < 256 ; a++ ) IUPAC_variance[a] = 0 ;
	for ( a = 0 ; a < 256 ; a++ ) {
		if ( ( char2IUPAC[a] & IUPAC_A ) > 0 ) IUPAC_variance[a]++ ;
		if ( ( char2IUPAC[a] & IUPAC_C ) > 0 ) IUPAC_variance[a]++ ;
		if ( ( char2IUPAC[a] & IUPAC_G ) > 0 ) IUPAC_variance[a]++ ;
		if ( ( char2IUPAC[a] & IUPAC_T ) > 0 ) IUPAC_variance[a]++ ;
	}

	for ( a = 0 ; a < 256 ; a++ ) to_lc[a] = a ;
	for ( a = 'A' ; a <= 'Z' ; a++ ) to_lc[a] = a - 'A' + 'a' ;

	for ( a = 0 ; a < 256 ; a++ ) to_uc[a] = a ;
	for ( a = 'a' ; a <= 'z' ; a++ ) to_uc[a] = a - 'a' + 'A' ;

	for ( a = 0 ; a < 256 ; a++ ) char2num[a] = 0 ;
	char2num[(uchar)'A'] = 1 ;
	char2num[(uchar)'C'] = 2 ;
	char2num[(uchar)'G'] = 3 ;
	char2num[(uchar)'T'] = 4 ;
	char2num[(uchar)'N'] = 5 ;
	char2num[(uchar)'X'] = 5 ;
	for ( a = 0 ; a < 256 ; a++ ) num2char[char2num[a]] = a ;
}

/// IUPAC-tolerant basematch for two char strings
bool matchIUPAC ( const char *c1 , const char *c2 ) {
	int left = globally_allowed_mismatches ;
	while ( *c1 && *c2 ) {
		if ( char2IUPAC[*c1++] & char2IUPAC[*c2++] ) continue ;
		if ( --left >= 0 ) continue ;
		return false ;
	}
	return true ;
}




/// Loads an entire FASTA file into a vector of TChromosome.
void load_all_chromosomes ( string filename , vector <TChromosome> &chrs , bool keep_original_snps ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Reading all chromosomes from %s ... " , filename.c_str() ) ; fflush(stdout); }
	string *seq = NULL ;
	FILE * file = fopen ( filename.c_str() , "r" ) ;
	if ( NULL == file ) {
		cerr << "file could not be accessed - aborting.\n" ;
		exit ( 1 ) ;
	}
	char dummy[READ_CHAR_BUFFER] , *c1 , *c2 ;
	while ( !feof(file) ) {
		*dummy = 0 ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		if ( *dummy == '>' ) {
			for ( c1 = dummy+1 ; *c1 == ' ' ; c1++ ) ;
			for ( c2 = c1 ; *c2 && *c2 != 10 && *c2 != 13 && *c2 != ' ' && *c2 != '|' ; c2++ ) ;
			*c2 = 0 ;
			chrs.push_back ( TChromosome() ) ;
			chrs[chrs.size()-1].name = c1 ;
			seq = &( chrs[chrs.size()-1].sequence ) ;
		} else {
			for ( c1 = dummy ; *c1 ; c1++ ) {
				switch ( *c1 ) {
					case 'a' : *c1 = 'A' ; break ;
					case 'c' : *c1 = 'C' ; break ;
					case 'g' : *c1 = 'G' ; break ;
					case 't' : *c1 = 'T' ; break ;
					case 'A' : break ;
					case 'C' : break ;
					case 'G' : break ;
					case 'T' : break ;
					case 10 : *c1 = 0 ; break ; // Line end
					case 13 : *c1 = 0 ; break ; // Line end
					default : if ( !keep_original_snps ) *c1 = 'X' ; // Setting everything else to 'X'
				}
//				if ( !*c1 ) break ;
			}
			*seq += dummy ;
		}
	}
	fclose ( file ) ;

	if ( DEBUGGING ) { fprintf ( stdout , "read %d chromosomes.\n" , chrs.size() ) ; fflush(stdout); }
}

/// Includes SNPs from a .gff file into chromosomes.
void incorporate_all_gff ( string filename , vector <TChromosome> &chrs ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Reading SNPS from %s ... ", filename.c_str() ) ; fflush(stdout); }

	// Backup
	for ( uint a = 0 ; a < chrs.size() ; a++ ) {
		if ( chrs[a].original_sequence.empty() )
			chrs[a].original_sequence = chrs[a].sequence ;
	}

	FILE * file = fopen ( filename.c_str() , "r" ) ;
	int snpcnt = 0 ;
	char dummy[READ_CHAR_BUFFER] , *c1 , *c2 ;
	while ( !feof(file) ) {
		vector <string> parts , parts2 ;
		*dummy = 0 ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		Tokenize ( dummy , parts , "\t" ) ;
		if ( parts.size() < 7 ) continue ; // These is not the data you are looking for

		int i ;
		for ( i = 0 ; i < chrs.size() && chrs[i].name != parts[0] ; i++ ) ;
		if ( i >= chrs.size() ) continue ; // Chromosome unknown, ignore this one

		if ( parts[2] != "SNP" ) continue ; // Not a SNP
		if ( parts[3] != parts[4] ) continue ; // Not a single point
		if ( parts[5] != "+" ) continue ; // Wrong strand (?)
		int pos = atoi ( parts[3].c_str() ) ;
		
		char merge = chrs[i].sequence[pos-1] ;
		Tokenize ( parts[6].c_str() , parts2 , "\"\" \"\"" ) ;
		bool found_reference = false ;
		for ( int b = 0 ; b < parts2.size() ; b++ ) {
			int l = parts2[b].length() ;
			if ( l < 2 ) continue ;
			if ( parts2[b][l-2] != ':' ) continue ;
			char c = parts2[b][l-1] ;
			if ( c == chrs[i].sequence[pos-1] ) found_reference = true ;
			merge = MERGE_IUPAC(merge,c) ;
		}
		if ( !found_reference ) { // Strange; maybe original sequence contains "N"
			// We should count this!
			continue ;  
		}
		chrs[i].sequence[pos-1] = merge ;
		snpcnt++ ;
	}
	fclose ( file ) ;
	if ( DEBUGGING ) { fprintf ( stdout , "incorporated %d SNPs.\n" , snpcnt ) ; fflush(stdout); }
}

/// Include SNPs from a simple SNP list into chromosome.
void incorporate_simple_snps ( string filename , vector <TChromosome> &chrs ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Reading SNPS from %s ... ", filename.c_str() ) ; fflush(stdout); }

	// Backup
	for ( uint a = 0 ; a < chrs.size() ; a++ ) {
		if ( chrs[a].original_sequence.empty() )
			chrs[a].original_sequence = chrs[a].sequence ;
	}

	FILE * file = fopen ( filename.c_str() , "r" ) ;
	int snpcnt = 0 , already_snp = 0 ;
	char dummy[READ_CHAR_BUFFER] , *c1 , *c2 ;
	uint chr_unknown = 0 ;
	uint reference_mismatch = 0 ;

	while ( !feof(file) ) {
		bool found_reference = false ;
		vector <string> parts ;
		*dummy = 0 ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		if ( !*dummy ) break ;
//		cout << dummy << endl ;

		Tokenize ( dummy , parts , "\t" ) ;
		if ( parts.size() < 4 ) {
			parts.clear() ;
			Tokenize ( dummy , parts , " " ) ;
			if ( parts.size() < 4 ) continue ; // These is not the data you are looking for
		}
		
		int i ;
		for ( i = 0 ; i < chrs.size() && chrs[i].name != parts[0] ; i++ ) ;
		if ( i >= chrs.size() ) { chr_unknown++ ; continue ; } // Chromosome unknown, ignore this one
		
		int pos = atoi ( parts[1].c_str() ) ;
		if ( pos < 1 || pos >= chrs[i].sequence.length() ) continue ; // Out-of-range
		char c = parts[2][0] ;
		char d = chrs[i].original_sequence[pos-1] ;
		if ( d == c ) found_reference = true ;
		else if ( d != 'A' && d != 'C' && d != 'G' && d != 'T' ) {
			found_reference = true ;
			already_snp++ ;
		} else {
			reference_mismatch++ ;
			//cerr << "Called as " << c << ", is really " << d << endl ;
		}
		c = parts[3][0] ;
		char merge = MERGE_IUPAC ( chrs[i].sequence[pos-1] , c ) ;
		chrs[i].sequence[pos-1] = merge ;
		snpcnt++ ;
	}

	fclose ( file ) ;
	if ( DEBUGGING ) { fprintf ( stdout , "incorporated %d SNPs (%d unknown chromosomes, %d reference mismatches, %d already marked as SNPs in genome).\n" , snpcnt , chr_unknown , reference_mismatch , already_snp ) ; fflush(stdout); }
}



// _________________________________________________________________________________________________________________________________________________________________ OPERATORS

/// Compares two indices; important for finding second-level indices only, so chromosome is not important.
bool operator < ( const TPositionWithIndex &a , const TPositionWithIndex &b ) {
	return a.index2 < b.index2 ;
}
