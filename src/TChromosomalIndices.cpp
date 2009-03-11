#include "snpomatic.h"

/// Constructor.
TChromosomalIndices::TChromosomalIndices ( int _index_length , int _index_length2 ) {
	index_length = _index_length ;
	index_length2 = _index_length2 ;
	index_length_double = index_length + index_length2 ;
	uint size = 1 ;
	for ( int a = 0 ; a < index_length ; a++ ) size *= 4 ;
	position.assign ( size , pwi_vector() ) ;
	noshortcuts = 0 ; // Take shortcuts by default
	max_snps_per_index = 8 ;
	index_from = 0 ;
	index_to = 0 ;
	index_steps = 1 ;
}

/// Indexes reads for reassembly.
void TChromosomalIndices::run_reads ( string filename ) {
	FILE * file = fopen ( filename.c_str() , "r" ) ;
	char dummy[READ_CHAR_BUFFER] , *c ;
	ulong count = 0 , a , from , to , l ;
	ulong filepos = 0 ;
	
	for ( a = 0 ; a < position.size() ; a++ ) position[a].reserve ( 95 ) ; // Speedup cheat
	
	while ( !feof(file) ) {
		count++ ;
		filepos = ftell ( file ) ;
		if ( count % 1000000 == 0 ) cerr << count << endl ;

		*dummy = 0 ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		if ( !*dummy ) break ;
		if ( *dummy != '>' ) break ; // Not a FASTA file, apparently
		*dummy = 0 ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		if ( !*dummy ) break ;
		
		l = 0 ;
		for ( c = dummy ; *c > 13 ; c++ ) {
			if ( isIUPAC[*c] ) { // No IUPAC reads, please...
				*dummy = 0 ;
				break ; 
			}
			l++ ;
		}
		if ( !*dummy ) continue ;
		*c = 0 ;
		
		from = 3 ; // 0 = no offset
		to = l - from - index_length_double ;
		
		for ( a = from ; a <= to ; a++ ) add_index ( dummy + from , false , filepos , from ) ;
		REVERSE_COMPLEMENT_CHARP ( dummy ) ;
		for ( a = from ; a <= to ; a++ ) add_index ( dummy + from , true , filepos , from ) ;
	}
	
	// Sort
	for ( a = 0 ; a < position.size() ; a++ ) sort ( position[a].begin() , position[a].end() ) ;
	
	// Analyze
	uint min = 5 ;
	uint max = 100 ;
	for ( a = 0 ; a < position.size() ; a++ ) {
		uint b , last = 0 , cnt = 0 ;
		for ( b = 0 ; b < position[a].size() ; b++ ) {
			if ( position[a][b].index2 == position[a][last].index2 ) {
				cnt++ ;
			} else {
				if ( cnt >= min && cnt <= max ) reassemble_reads ( file , a , last , b-1 ) ;
				last = b ;
				cnt = 1 ;
			}
		}
		if ( cnt >= min && cnt <= max ) reassemble_reads ( file , a , last , b-1 ) ;
	}
}

bool TChromosomalIndices::align_read2contigs ( char *read , vector <string> &vs ) {
	uint a , vsmax = vs.size() ;
	int l , b , c ;
	bool ret = false ;
	for ( l = 0 ; *(read+l) ; l++ ) ;
	int from = -10 ;
	int to = 20 ;
	for ( a = 0 ; a < vsmax ; a++ ) {
		if ( vs[a] == read ) return ret ;
		for ( b = from ; b < to ; b++ ) {
			bool match = true ;
			for ( c = 0 ; match && c < l ; c++ ) {
				if ( b + c < 0 ) continue ;
				if ( b + c >= vs[a].length() ) continue ;
				if ( *(read+c) != vs[a][b+c] ) match = false ;
			}
			if ( match ) {
				string m ;
				string n ;
				for ( c = b ; c < 0 ; c++ ) n += ' ' ;
				for ( c = 0 ; c < b ; c++ ) m += ' ' ;
				n += vs[a] ;
				m += read ;
				if ( n.length() == vs[a].length() ) continue ;
				
				while ( n.length() < m.length() ) n += ' ' ;
				while ( m.length() < n.length() ) m += ' ' ;
				for ( c = 0 ; c < n.length() ; c++ ) {
					if ( n[c] == ' ' ) n[c] = m[c] ;
				}
				vs.push_back ( n ) ;
				
				break ; // One match will do just fine...
			}
		}
	}
	if ( !ret ) vs.push_back ( read ) ; // New contig since it doesn't match old ones
	return ret ;
}

void TChromosomalIndices::reassemble_reads ( FILE *file , uint i1 , index2type i2a , index2type i2b ) {
	char dummy[READ_CHAR_BUFFER] , *c ;
	vector <string> contigs ;
	uint maxlen = 0 ;
	for ( index2type i2 = i2a ; i2 <= i2b ; i2++ ) {
		fseek ( file , position[i1][i2].position , SEEK_SET ) ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		uint l = 0 ;
		for ( c = dummy ; *c > 13 ; l++ , c++ ) ; // All possible hits are already IUPAC-filtered, no need to do it again
		*c = 0 ;
		if ( l > maxlen ) maxlen = l ;
		if ( position[i1][i2].reverse_complement ) REVERSE_COMPLEMENT_CHARP ( dummy ) ;
		
		align_read2contigs ( dummy , contigs ) ;
	}
	for ( uint a = 0 ; a < contigs.size() ; a++ ) {
		if ( contigs[a].length() <= maxlen ) continue ;
		cout << contigs[a] << endl ;
	}
}

void TChromosomalIndices::add_index ( char *start , bool rc , uint pos , uchar chromosome ) {
	uint key = calculate_index ( start ) ;
	TPositionWithIndex pwi ;
	pwi.position = pos ;
	pwi.index2 = calculate_index2 ( start + index_length ) ;
	pwi.reverse_complement = rc ;
	pwi.chromosome = chromosome ;
	position[key].push_back ( pwi ) ;
}


/// Creates index for all regions, mark regions N
void TChromosomalIndices::run_regions ( chrvec *_chrs , string filename ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Indexing regions ... " ) ; fflush(stdout); }
	chrs = _chrs ;
	int a , chr ;

//	for ( a = 0 ; a < position.size() ; a++ ) position[a].reserve ( 100 ) ; // Pre-reserve memory, rule-of-thumb
	
//	cout << endl ;
	// Go through region file
	FILE * file = fopen ( filename.c_str() , "r" ) ;
	char dummy[READ_CHAR_BUFFER] ;
	while ( !feof(file) ) {
		*dummy = 0 ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		if ( !*dummy ) break ;
		vector <string> parts ;
		Tokenize ( dummy , parts , "\t" ) ;
		if ( parts.size() < 3 ) continue ;
		
		for ( chr = 0 ; chr < chrs->size() && (*chrs)[chr].name != parts[0] ; a++ ) ;
		if ( chr == chrs->size() ) {
			cout << "Unknown chromosome " << parts[0] << " in region list!\n" ;
			continue ;
		}

//		cout << dummy ;
		uint from = atoi ( parts[1].c_str() ) - 1 ;
		uint to = atoi ( parts[2].c_str() ) - 1 ;
		
		// Mark region as "N"
		for ( a = from ; a <= to ; a++ ) (*chrs)[chr].sequence[a] = 'N' ;
		
		uint offset = 20 ;
		string s ;
		uint pl = (*chrs)[chr].sequence.length() - index_length_double ;
		string *p = & ( (*chrs)[chr].sequence ) ;
		
		for ( a = from - index_length_double - offset ; a <= to ; a++ ) {
			if ( a < 0 ) continue ;
			if ( a >= pl ) continue ;
			s = p->substr ( a , index_length_double ) ;
			
			char *sc = (char*) s.c_str() ;
			add_sequence ( sc , sc , a , chr , false ) ;
			add_sequence ( sc , sc , a+index_length_double , chr , true ) ;
		}
		
	}

	// Sort all sub-lists
	for ( a = 0 ; a < position.size() ; a++ ) {
		sort ( position[a].begin() , position[a].end() ) ;
	}

	if ( DEBUGGING ) { fprintf ( stdout , "done.\n" ) ; fflush(stdout); }
}

/// Create index for beginning and end of all contigs.
void TChromosomalIndices::run_ends ( chrvec *_chrs , int min_length ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Indexing contigs ... " ) ; fflush(stdout); }
	chrs = _chrs ;
	uint a , chr ;

//	for ( a = 0 ; a < position.size() ; a++ ) position[a].reserve ( 100 ) ; // Pre-reserve memory, rule-of-thumb
	
	for ( chr = 0 ; chr < chrs->size() ; chr++ ) {
		if ( (*chrs)[chr].sequence.length() < min_length ) continue ;
		string s , *p ;
		uint pl = (*chrs)[chr].sequence.length() - index_length_double ;
		char *sc = (char*) (*chrs)[chr].sequence.c_str() ;

		for ( uint a = 0 ; a < pl ; a++ ) {
			add_sequence ( sc , sc , a , chr , false ) ;
//			add_sequence ( sc , sc , a+index_length_double , chr , true ) ;
		}
	}
	
/*	
	if ( DEBUGGING ) { fprintf ( stdout , "sorting ... " ) ; fflush(stdout); }
	// Sort all sub-lists
	for ( a = 0 ; a < position.size() ; a++ ) {
		sort ( position[a].begin() , position[a].end() ) ;
	}
*/
	if ( DEBUGGING ) { fprintf ( stdout , "done.\n" ) ; fflush(stdout); }
}

string TChromosomalIndices::decompile ( uint index , uint bases ) {
	string ret ;
	while ( bases > 0 ) {
		ret += num2char[(index&3)+1] ;
		index >>= 2 ;
		bases-- ;
	}
	reverse ( ret.begin() , ret.end() ) ;
	return ret ;
}


/// Create index for all chromosomes.
void TChromosomalIndices::run ( chrvec *_chrs ) {
	chrs = _chrs ;
	if ( !index_file.empty() && fopen ( index_file.c_str() , "r" ) ) {
		load_index ( index_file ) ;
		return ;
	}
	
	if ( DEBUGGING ) {
		if ( index_from == 0 && index_to == 0 ) fprintf ( stdout , "Indexing %d chromosomes ... " , chrs->size() ) ;
		else fprintf ( stdout , "Indexing %d chromosomes (%d-%d) ... " , chrs->size() , index_from , index_to ) ;
		fflush(stdout);
	}
	
	uint a , chr ;

	for ( chr = 0 ; chr < chrs->size() ; chr++ ) {
		string s , *p ;
		char *sc = (char*) (*chrs)[chr].sequence.c_str() ;
		uint pl = (*chrs)[chr].sequence.length() ;
		if ( index_to > 0 && index_to < pl ) pl = index_to ;
		sc += index_from ;
		for ( a = index_from ; a < pl ; a += index_steps ) {
			int r = 0 ;
			r += add_sequence ( sc , sc , a , chr , false ) ;
			r += add_sequence ( sc , sc , a+index_length_double , chr , true ) ;
			sc += index_steps ;
		}
	}
	
	if ( DEBUGGING ) { fprintf ( stdout , "sorting ... " ) ; fflush(stdout); }
	
	// Nuke huge indices; speedup: 50%
	if ( !noshortcuts ) {
		pwi_vector largest ;
		for ( a = 0 ; a < position.size() ; a++ ) {
			if ( largest.size() < 20 ) {
				TPositionWithIndex tmp ;
				tmp.index2 = position[a].size() ;
				tmp.position = a ;
				largest.push_back ( tmp ) ;
				sort ( largest.begin() , largest.end() ) ;
				continue ;
			}
			if ( largest[0].index2 >= position[a].size() ) continue ;

			largest[0].index2 = position[a].size() ;
			largest[0].position = a ;
			sort ( largest.begin() , largest.end() ) ;
		}

		for ( a = 0 ; a < largest.size() ; a++ ) {
			if ( largest[a].index2 < largest.back().index2 / 4 ) continue ;
			position[largest[a].position].clear() ;
		}
	}

	// Sort all sub-lists
	for ( a = 0 ; a < position.size() ; a++ ) {
		sort ( position[a].begin() , position[a].end() ) ;
	}

	if ( DEBUGGING ) { fprintf ( stdout , "done.\n" ) ; fflush(stdout); }

	if ( !index_file.empty() ) {
		save_index ( index_file ) ;
		return ;
	}
}


/// Recurses/iterates through sequence with possible IUPAC codes, creates all possible combinations, and adds them to the index list.
int TChromosomalIndices::add_sequence ( char *begin , char *cur , int from , int chr , bool rc , int nos ) {
	if ( nos > max_snps_per_index ) return 0 ;
	int ret = 0 ;
	while ( index_length_double > cur - begin ) {
		if ( isIUPAC[*cur] ) {
			if ( *cur == 'X' ) return ret ; // Contains X, aborting
			uchar c = char2IUPAC[*cur] ;
			uchar oc = *cur ;
			if ( c & IUPAC_A ) { *cur = 'A' ; ret += add_sequence ( begin , cur + 1 , from , chr , rc , nos+1 ) ; }
			if ( c & IUPAC_C ) { *cur = 'C' ; ret += add_sequence ( begin , cur + 1 , from , chr , rc , nos+1 ) ; }
			if ( c & IUPAC_G ) { *cur = 'G' ; ret += add_sequence ( begin , cur + 1 , from , chr , rc , nos+1 ) ; }
			if ( c & IUPAC_T ) { *cur = 'T' ; ret += add_sequence ( begin , cur + 1 , from , chr , rc , nos+1 ) ; }
			*cur = oc ;
			return ret ;
		}
		cur++ ;
	}
	
	char *sc = begin ;
	if ( rc ) {
		sc = tmp ;
		memcpy ( tmp , begin , index_length_double ) ;
		tmp[index_length_double] = 0 ;
		REVERSE_COMPLEMENT_CHARP ( sc ) ;
	}

	uint key = calculate_index ( sc ) ;
	TPositionWithIndex pwi ;
	pwi.position = from ;
	pwi.index2 = calculate_index2 ( sc + index_length ) ;
	pwi.reverse_complement = rc ;
	pwi.chromosome = chr ;
	position[key].push_back ( pwi ) ;
	return ret + 1 ;
}

/// sequence_kmer MUST only contain As, Cs, Gs, or Ts
uint TChromosomalIndices::calculate_index ( const char *sequence_kmer ) {
	uint ret = 0 ;
	for ( int a = 0 ; a < index_length ; a++ ) {
		ret = ( ret << 2 ) | ( ((uint)char2num[(uint)sequence_kmer[a]])-1 ) ;
	}
	return ret ;
}

/// sequence_kmer MUST only contain As, Cs, Gs, or Ts
index2type TChromosomalIndices::calculate_index2 ( const char *sequence_kmer ) {
	uint ret = 0 ;
	for ( int a = 0 ; a < index_length2 ; a++ ) {
		ret = ( ret << 2 ) | ( ((uint)char2num[(uint)sequence_kmer[a]])-1 ) ;
	}
	return ret ;
}

void TChromosomalIndices::save_index ( string filename ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Saving index file ... " ) ; fflush(stdout); }
	FILE *file = fopen ( string ( filename.c_str() ).c_str() , "w" ) ;
	uint a , b ;
	
	uint temp = 1 ; // Version
	fwrite ( &temp , sizeof ( temp ) , 1 , file ) ;
	temp = position.size() ;
//	cout << endl << position.size() << endl << temp << endl ;
	fwrite ( &temp , sizeof ( temp ) , 1 , file ) ;
	for ( a = 0 ; a < position.size() ; a++ ) {
		temp = position[a].size() ;
		fwrite ( &temp , sizeof ( temp ) , 1 , file ) ;
		for ( b = 0 ; b < position[a].size() ; b++ ) {
			fwrite ( &position[a][b].index2 , sizeof ( position[a][b].index2 ) , 1 , file ) ;
			fwrite ( &position[a][b].position , sizeof ( position[a][b].position ) , 1 , file ) ;
			fwrite ( &position[a][b].chromosome , sizeof ( position[a][b].chromosome ) , 1 , file ) ;
			fwrite ( &position[a][b].reverse_complement , sizeof ( position[a][b].reverse_complement ) , 1 , file ) ;
		}
	}
	
	fclose ( file ) ;
	if ( DEBUGGING ) { fprintf ( stdout , "done.\n" ) ; fflush(stdout); }
}

void TChromosomalIndices::load_index ( string filename ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Loading index file ... " ) ; fflush(stdout); }
	uint temp ;
	FILE *file = fopen ( string ( filename.c_str() ).c_str() , "r" ) ;
	if ( !file ) exit ( 1 ) ;
	fread ( &temp , sizeof ( temp ) , 1 , file ) ; // Version
	fread ( &temp , sizeof ( temp ) , 1 , file ) ; // position size
	if ( position.size() != temp ) exit ( 1 ) ;

	uint a , b ;
	for ( a = 0 ; a < position.size() ; a++ ) {
		fread ( &temp , sizeof ( temp ) , 1 , file ) ; // Sub-index size
		position[a].reserve ( temp ) ;
		for ( b = 0 ; b < temp ; b++ ) {
			TPositionWithIndex i ;
			fread ( &i.index2 , sizeof ( i.index2 ) , 1 , file ) ;
			fread ( &i.position , sizeof ( i.position ) , 1 , file ) ;
			fread ( &i.chromosome , sizeof ( i.chromosome ) , 1 , file ) ;
			fread ( &i.reverse_complement , sizeof ( i.reverse_complement ) , 1 , file ) ;
			position[a].push_back ( i ) ;
//			if ( i.chromosome > 100 ) cerr << "!" << i.chromosome << endl ;
		}
	}

	fclose ( file ) ;
	if ( DEBUGGING ) { fprintf ( stdout , "done.\n" ) ; fflush(stdout); }
}

void TChromosomalIndices::uniqueness ( string filename ) {
	uint a , b , c ;
	
	for ( a = 0 ; a < chrs->size() ; a++ ) {
		(*chrs)[a].uniqueness.resize ( (*chrs)[a].sequence.length() , 0 ) ;
	}
	
	for ( a = 0 ; a < position.size() ; a++ ) {
		for ( b = 0 ; b < position[a].size() ; b++ ) {
			for ( c = b+1 ; c < position[a].size() && position[a][b].index2 == position[a][c].index2 ; c++ ) ;
			if ( c == b+1 ) continue ;
			uint count = c - b ;
			cout << count << endl ;
			for ( c = b ; c < position[a].size() && position[a][b].index2 == position[a][c].index2 ; c++ ) {
				uint p = position[a][c].position ;
				uint chr = position[a][c].chromosome ;
				if ( chr > chrs->size() ) {
					cerr << position[a][c].chromosome << "/" << c << " : " << chr << endl ;
					exit ( 1 ) ;
				}
				uint si = (*chrs)[chr].uniqueness.size() ;
				if ( si == 0 ) {
					(*chrs)[chr].uniqueness.resize ( (*chrs)[chr].sequence.length() , 0 ) ;
					si = (*chrs)[chr].uniqueness.size() ;
				}
				for ( uint u = 0 ; u < index_length_double ; u++ ) {
					if ( u + p >= si ) break ;
					(*chrs)[chr].uniqueness[u+p] += count ;
				}
			}
			b = c - 1 ;
		}
	}

	for ( a = 0 ; a < chrs->size() ; a++ ) {
		for ( b = 0 ; b < (*chrs)[a].uniqueness.size() ; b++ ) {
			if ( (*chrs)[a].uniqueness[b] == 0 )
				(*chrs)[a].uniqueness[b] = 1 ;
		}
	}
	
	FILE *u = fopen ( filename.c_str() , "w" ) ;
	for ( a = 0 ; a < chrs->size() ; a++ ) {
		for ( b = 0 ; b < (*chrs)[a].uniqueness.size() ; b++ ) {
			fprintf ( u , "%s\t%d\t%d\n" , (*chrs)[a].name.c_str() , b + 1 , (*chrs)[a].uniqueness[b] ) ;
		}
	}
	fclose ( u ) ;
}
