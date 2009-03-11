#include "snpomatic.h"

/// Reads a single chromosome from a FASTA file; superceded by load_all_chromosomes().
void TChromosome::read_chromosome ( string filename , string chromosome_name ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Reading %s from %s ... " , chromosome_name.c_str() , filename.c_str() ) ; fflush(stdout); }
	
	name = chromosome_name ;
	FILE * file = fopen ( filename.c_str() , "r" ) ;
	char dummy[READ_CHAR_BUFFER] , *c1 , *c2 ;
	uchar store = 0 ;
	while ( !feof(file) ) {
		*dummy = 0 ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		if ( *dummy == '>' ) {
			if ( store ) break ; // Read the chromosome
			for ( c1 = dummy+1 ; *c1 == ' ' ; c1++ ) ;
			for ( c2 = c1 ; *c2 && *c2 != 10 && *c2 != 13 && *c2 != ' ' && *c2 != '|' ; c2++ ) ;
			*c2 = 0 ;
			if ( c1 == chromosome_name ) store = 1 ;
		} else if ( store ) {
			for ( c1 = dummy ; *c1 ; c1++ ) {
				switch ( *c1 ) {
					case 'a' : *c1 = 'A' ; break ;
					case 'c' : *c1 = 'C' ; break ;
					case 'g' : *c1 = 'G' ; break ;
					case 't' : *c1 = 'T' ; break ;
				}

//				if ( *c1 == 'A' || *c1 == 'C' || *c1 == 'G' || *c1 == 'T' ) {
					sequence += *c1 ;
//				} else if ( *c1 != 13 && *c1 != 10 ) {
//					fprintf ( stdout , "Genome data contains incompatible character '%c' - aborting.\n" , *c1 ) ;
//					exit ( 1 ) ;
//				}

			}
		}
	}
	fclose ( file ) ;
	
	if ( DEBUGGING ) { fprintf ( stdout , "read %d nucleotides.\n" , sequence.length() ) ; fflush(stdout); }
}

/// Reads a SNP list in .gff and puts the SNPs as IUPAC bases into the chromosome; superceded by incorporate_all_gff().
void TChromosome::incorporate_known_snps_gff ( string filename ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Reading SNPS into %s from %s ... " , name.c_str() , filename.c_str() ) ; fflush(stdout); }

	FILE * file = fopen ( filename.c_str() , "r" ) ;
	int snpcnt = 0 ;
	char dummy[READ_CHAR_BUFFER] , *c1 , *c2 ;
	while ( !feof(file) ) {
		vector <string> parts , parts2 ;
		*dummy = 0 ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		Tokenize ( dummy , parts , "\t" ) ;
		if ( parts.size() < 7 ) continue ; // These is not the data you are looking for
		if ( parts[0] != name ) continue ; // Wrong chromosome
		if ( parts[2] != "SNP" ) continue ; // Not a SNP
		if ( parts[3] != parts[4] ) continue ; // Not a single point
		if ( parts[5] != "+" ) continue ; // Wrong strand (?)
		int pos = atoi ( parts[3].c_str() ) ;
		
		char merge = sequence[pos-1] ;
		Tokenize ( parts[6].c_str() , parts2 , "\"\" \"\"" ) ;
		bool found_reference = false ;
		for ( int b = 0 ; b < parts2.size() ; b++ ) {
			int l = parts2[b].length() ;
			if ( l < 2 ) continue ;
			if ( parts2[b][l-2] != ':' ) continue ;
			char c = parts2[b][l-1] ;
			if ( c == sequence[pos-1] ) found_reference = true ;
			merge = MERGE_IUPAC(merge,c) ;
		}
		if ( !found_reference ) { // Strange; maybe original sequence contains "N"
			// We should count this!
			continue ;  
		}
		sequence[pos-1] = merge ;
		snpcnt++ ;
	}
	fclose ( file ) ;
	if ( DEBUGGING ) { fprintf ( stdout , "incorporated %d SNPs.\n" , snpcnt ) ; fflush(stdout); }
}

/// Reads a SNP list in simple format and puts the SNPs as IUPAC bases into the chromosome.
void TChromosome::incorporate_known_snps ( string filename ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Reading from %s ... " , filename.c_str() ) ; fflush(stdout); }

	FILE * file = fopen ( filename.c_str() , "r" ) ;
	char dummy[READ_CHAR_BUFFER] , *c1 , *c2 ;
	int snpcnt = 0 ;
	while ( !feof(file) ) {
		vector <string> parts ;
		*dummy = 0 ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		Tokenize ( dummy , parts , "\t" ) ;
		if ( parts.size() != 4 ) continue ;
		if ( name != parts[0] ) continue ;
		int pos = atoi ( parts[1].c_str() ) ;
		uchar orig = parts[2][0] ;
		uchar snp = parts[3][0] ;
		sequence[pos-1] = MERGE_IUPAC(orig,snp) ;
		snpcnt++ ;
	}
	fclose ( file ) ;

	if ( DEBUGGING ) { fprintf ( stdout , "incorporated %d SNPs.\n" , snpcnt ) ; fflush(stdout); }
}

/// Adds a (single) read to the coverage statistics
void TChromosome::add_coverage ( const string &seq , uint pos , bool rc ) {
	if ( coverage.size() == 0 ) { // Set up vector
		TCoverage dc ;
		dc.A = dc.C = dc.T = dc.G = dc.rc = 0 ; // Paranoia
		coverage.resize ( sequence.length() , dc ) ;
	}

	for ( uint p = 0 ; p < seq.length() ; p++ ) {
		if ( p + pos >= sequence.length() ) break ; // End of the line, my friend...
		if ( rc ) coverage[p+pos].rc++ ;
		switch ( seq[p] ) {
			case 'A' : coverage[p+pos].A++ ; break ;
			case 'C' : coverage[p+pos].C++ ; break ;
			case 'G' : coverage[p+pos].G++ ; break ;
			case 'T' : coverage[p+pos].T++ ; break ;
			default : cerr << seq[p] << endl ;
		} ;
	}
}
