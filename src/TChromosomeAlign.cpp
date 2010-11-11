#include "snpomatic.h"

string dor ( string s ) { // Do Reverse
	reverse ( s.begin() , s.end() ) ;
	return s ;
}

string doc ( string s ) { // Do Complementary
	for ( int a = 0 ; a < s.length() ; a++ ) s[a] = base_rc[s[a]] ;
	return s ;
}

string dorc ( string s ) { // Do Reverse/Complementary
	for ( int a = 0 ; a < s.length() ; a++ ) s[a] = base_rc[s[a]] ;
	reverse ( s.begin() , s.end() ) ;
	return s ;
}

bool tpi_full_compare ( const TPositionWithIndex a , const TPositionWithIndex b ) {
	return ( a.chromosome != b.chromosome ) ? a.chromosome < b.chromosome : a.position < b.position ;
}




/// Constructor
TChromosomeAlign::TChromosomeAlign ( chrvec &_chrs , TChromosomalIndices &_index ) {
	pileup = binfile_no_match = binfile_single_match = binfile_multi_match = binfile_iupac = cigar = NULL ;
	wobble = fragmentplot = featuretable = gffout = coverage = snpsinreads = indelplot = inversions = NULL ;
	sqlite = sam = spancontigs = faceaway = rpa = NULL ;
	chrs = &_chrs ;
	index = &_index ;
	chop = 0 ;
	wobblemax = 2 ;
	max_align = (uint) 2 << 30 ; // Just a very large number
	fasta = false ;
	use_nono = false ;
	using_bins = false ;
	snpsinreads_first = true ;
	multimatch = false ;
	singlematch = false ;
	force_one_unique_match = false ;
	sqlite_prefix_first = true ;
	rpa_header_written = false ;
	single_read_length = 0 ;
}

void TChromosomeAlign::dump_coverage () {
	if ( !coverage ) return ;
	fprintf ( coverage , "Chromosome\tPosition\tReference\tA\tC\tG\tT\tReverse-complement\n" ) ;
	for ( uint chr = 0 ; chr < chrs->size() ; chr++ ) {
		if ( (*chrs)[chr].coverage.size() == 0 ) continue ; // Paranoia
		for ( uint pos = 0 ; pos < (*chrs)[chr].sequence.length() ; pos++ ) {
			if ( pos >= (*chrs)[chr].coverage.size() ) {
			fprintf ( coverage , "%s\t%d\t%c\t%d\t%d\t%d\t%d\t%d\n" , 
				(*chrs)[chr].name.c_str() ,
				pos + 1 ,
				(*chrs)[chr].original_sequence.length() > 0 ? (*chrs)[chr].original_sequence[pos] : (*chrs)[chr].sequence[pos] ,
				0 , 0 , 0 , 0 , 0 ) ;
			}
			fprintf ( coverage , "%s\t%d\t%c\t%d\t%d\t%d\t%d\t%d\n" , 
				(*chrs)[chr].name.c_str() ,
				pos + 1 ,
				(*chrs)[chr].original_sequence.length() > 0 ? (*chrs)[chr].original_sequence[pos] : (*chrs)[chr].sequence[pos] ,
				(*chrs)[chr].coverage[pos].A ,
				(*chrs)[chr].coverage[pos].C ,
				(*chrs)[chr].coverage[pos].G ,
				(*chrs)[chr].coverage[pos].T ,
				(*chrs)[chr].coverage[pos].rc ) ;
		}
	}
	fclose ( coverage ) ;
	coverage = NULL ; // Paranoia
}

/// Initialization (after constructor)
void TChromosomeAlign::init () {
	int i ;
	if ( pileup ) for ( i = 0 ; i < chrs->size() ; i++ ) (*chrs)[i].ao.init(&((*chrs)[i])) ;
	if ( sam ) {
		fprintf ( sam , "@HD\tVN:1.0\n" ) ;
		for ( i = 0 ; i < chrs->size() ; i++ ) {
			fprintf ( sam , "@SQ\tSN:%s\tLN:%d\n" , (*chrs)[i].name.c_str() , (*chrs)[i].sequence.length() ) ;
		}
	}
}

void TChromosomeAlign::add_nono_list ( string filename ) {
	FILE * file = fopen ( filename.c_str() , "r" ) ;
	char dummy[READ_CHAR_BUFFER] , *c1 , *c2 ;
	while ( !feof(file) ) {
		*dummy = 0 ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		if ( !*dummy ) break ;
		
		c1 = dummy ;
		if ( *c1 == '>' || *c1 == '@' ) c1++ ;
		while ( *c1 == ' ' ) c1++ ;
		for ( c2 = c1 ; *c2 > 13 ; c2++ ) ;
		*c2 = 0 ;
		nono.push_back ( c1 ) ;
	}
	sort ( nono.begin() , nono.end() ) ;
	use_nono = true ;
}

bool TChromosomeAlign::is_nono ( const string &s ) {
	pair <TSVI , TSVI> ss = equal_range ( nono.begin() , nono.end() , s ) ;
	return ss.first != ss.second ;
}

void TChromosomeAlign::write2bin ( const string &name , const string &sequence , const string &quality , int matches ) {
	string out ;
	if ( fasta ) {
		out = ">" + name + "\n" + sequence + "\n" ;
	} else {
		out = "@" + name + "\n" + sequence + "\n+\n" + quality + "\n" ;
	}

	if ( matches < 0 ) { if ( binfile_iupac ) fprintf ( binfile_iupac , "%s" , out.c_str() ) ; }
	else if ( matches == 0 ) { if ( binfile_no_match ) fprintf ( binfile_no_match , "%s" , out.c_str() ) ; }
	else if ( matches == 1 ) { if ( binfile_single_match ) fprintf ( binfile_single_match , "%s" , out.c_str() ) ; }
	else { if ( binfile_multi_match ) fprintf ( binfile_multi_match , "%s" , out.c_str() ) ; }
}

uint TChromosomeAlign::align_solexa_paired_read ( const string &name , const string &sequence , string quality , uint read_length_1 , uint fragment_length , uint range ) {
	string seq1 = sequence.substr ( 0 , read_length_1 ) ;
	string seq2 = sequence.substr ( read_length_1 , sequence.length() - read_length_1 ) ;
	
	while ( quality.length() < sequence.length() ) quality += (char) 33+30 ;
	string qual1 = quality.substr ( 0 , read_length_1 ) ;
	string qual2 = quality.substr ( read_length_1 , sequence.length() - read_length_1 ) ;
	
	pwi_vector res1 , res2 ;
	get_align_solexa_read ( seq1.c_str() , res1 ) ;
	if ( !sam && !wobble && !chop && !sqlite && res1.size() == 0 ) return 0 ;
	get_align_solexa_read ( seq2.c_str() , res2 ) ;
	if ( !sam && !wobble && !chop && !sqlite && res2.size() == 0 ) return 0 ;

	sort ( res1.begin() , res1.end() , tpi_full_compare ) ;
	sort ( res2.begin() , res2.end() , tpi_full_compare ) ;

	int ret = 0 , old_fragment , old_variance ;

	// Spanning contigs
	if ( ( rpa || spancontigs ) && res1.size() == 1 && res2.size() == 1 && res1[0].chromosome != res2[0].chromosome ) {
		if ( spancontigs ) {
			fprintf ( spancontigs , "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\n" , 
				last_solexa_name.c_str() ,
				(*chrs)[res1[0].chromosome].name.c_str() , res1[0].position , res1[0].reverse_complement , seq1.c_str() , qual1.c_str() ,
				(*chrs)[res2[0].chromosome].name.c_str() , res2[0].position , res2[0].reverse_complement , seq2.c_str() , qual2.c_str()
				) ;
		}
		if ( rpa ) rpa_out ( res1[0] , res2[0] , read_length_1 , seq1 , qual1 , seq2 , qual2 ) ;
		return 0 ;
	}
	
	if ( ( sqlite || sam ) && res1.size() + res2.size() == 1 ) { // Unique single read
		TPositionWithIndex r = res1.size() ? res1[0] : res2[0] ;
		string seq = res1.size() ? seq1 : seq2 , s ;
		if ( sqlite ) {
			if ( r.reverse_complement ) s = generate_sqlite_paired_command ( r.chromosome , dorc ( seq ) , r.position - seq.length() , "" , 0 , '1' ) ;
			else s = generate_sqlite_paired_command ( r.chromosome , seq , r.position , "" , 0 , '1' ) ;
			add2sqlite_cache ( s ) ;
		}
		if ( sam ) {
			string qual = res1.size() ? qual1 : qual2 , s ;
			int readpart = res1.size() ? 1 : 2 ;
			if ( r.reverse_complement ) add_sam ( dorc ( seq ) , dor ( qual ) , r.chromosome , r.position - seq.length() , readpart , 0 , 0 , r.reverse_complement , false , false , true ) ;
			else add_sam ( seq , qual , r.chromosome , r.position , readpart , 0 , 0 , r.reverse_complement , false , false , true ) ;
		}
		return 0 ;
	}
	
	
	if ( sqlite ) { // Clear temporary cache, increase fragment/variance temporarily
		sqlite_out.clear () ;
	}
	sam_wrote = false ;
	ret += paired_read_combine ( res1 , res2 , fragment_length , range , read_length_1 , seq1 , qual1 , seq2 , qual2 ) ;
	ret += paired_read_combine ( res2 , res1 , fragment_length , range , read_length_1 , seq2 , qual2 , seq1 , qual1 ) ;
	
	if ( sqlite && sqlite_out.size() == 1 ) {
		add2sqlite_cache ( sqlite_out[0] ) ;
		return ret ;
	}
	
	// Try all fragment sizes if sqlite and no in-range find
	if ( ( sam && !sam_wrote ) || ( sqlite && sqlite_out.size() == 0 ) ) {
		old_fragment = fragment_length ;
		old_variance = range ;
		fragment_length = 90000000 ; // Long...
		range = 90000000 ; // Looooooong...
		bool done = false ;

		paired_read_combine ( res1 , res2 , fragment_length , range , read_length_1 , seq1 , qual1 , seq2 , qual2 ) ;
		paired_read_combine ( res2 , res1 , fragment_length , range , read_length_1 , seq2 , qual2 , seq1 , qual1 ) ;
		if ( sqlite_out.size() == 1 ) {
			add2sqlite_cache ( sqlite_out[0] ) ;
			done = true ;
		}

		fragment_length = old_fragment ; // Ah. Just kidding.
		range = old_variance ;
		if ( done ) return 0 ;
	}

	// Try inversions
	if ( ret == 0 && ( rpa || inversions || sqlite || ( sam && !sam_wrote ) ) ) {
		bool fake_long_region = true ; // sqlite ;

		if ( fake_long_region ) {
			old_fragment = fragment_length ;
			old_variance = range ;
			fragment_length = 90000000 ; // Long...
			range = 90000000 ; // Looooooong...
		}
		
		paired_read_combine ( res1 , res2 , fragment_length , range , read_length_1 , seq1 , qual1 , seq2 , qual2 , false , true , 1 ) ;
		paired_read_combine ( res2 , res1 , fragment_length , range , read_length_1 , seq2 , qual2 , seq1 , qual1 , false , true , 2 ) ;

		if ( fake_long_region ) {
			fragment_length = old_fragment ; // Ah. Just kidding.
			range = old_variance ;
		}

		if ( sqlite_out.size() == 1 ) {
			add2sqlite_cache ( sqlite_out[0] ) ;
			return 0 ;
		}
	}
	
	if ( ret > 0 ) return ret ;
	
	// Chop off read ends and try again (optional!)
	if ( chop > 0 ) {
		uint temp_chop = chop ;

		if ( 1 ) { // Chop off read end
			string nseq = seq1.substr ( 0 , seq1.length() - chop ) + seq2.substr ( chop , seq2.length() - chop ) ;
			string nqual = qual1.substr ( 0 , qual1.length() - chop ) + qual2.substr ( chop , qual2.length() - chop ) ;
			chop = 0 ;
			ret = align_solexa_paired_read ( name , nseq , nqual , read_length_1 - temp_chop , fragment_length , range ) ;
			chop = temp_chop ;
		}
		
		if ( ret == 0 ) { // Still no luck, try to chop off the beginning
			string nseq = seq1.substr ( chop , seq1.length() - chop ) + seq2.substr ( 0 , seq2.length() - chop ) ;
			string nqual = qual1.substr ( chop , qual1.length() - chop ) + qual2.substr ( 0 , qual2.length() - chop ) ;
			chop = 0 ;
			ret = align_solexa_paired_read ( name , nseq , nqual , read_length_1 - temp_chop , fragment_length , range ) ;
			chop = temp_chop ;
		}

		if ( ret == 0 ) { // Still no luck, try to chop off some of both beginning and end
			uint chop1 = chop / 2 ;
			uint chop2 = chop - chop1 ;
			string nseq = seq1.substr ( chop1 , seq1.length() - chop ) + seq2.substr ( chop2 , seq2.length() - chop ) ;
			string nqual = qual1.substr ( chop1 , qual1.length() - chop ) + qual2.substr ( chop2 , qual2.length() - chop ) ;
			chop = 0 ;
			ret = align_solexa_paired_read ( name , nseq , nqual , read_length_1 - temp_chop , fragment_length , range ) ;
			chop = temp_chop ;
		}

	}
	else if ( wobble ) run_wobble ( seq1 , seq2 , qual1 , qual2 , res1 , res2 , fragment_length ) ; // Don't do in chopped and unchopped state

	return ret ;
}

void TChromosomeAlign::get_align_solexa_read ( const char *_seq , pwi_vector &results ) {
	string rc , rcq ;
	uint l = strlen ( _seq ) ;
//	results.reserve ( 50 ) ;

	for ( uint offset = 0 ; offset < index->index_steps ; offset++ ) {
		uint i1 = index->calculate_index ( _seq + offset ) ;
		if ( !index->position[i1].size() ) continue ;
		
		// Binary search
		TPositionWithIndex findme ;
		findme.index2 = index->calculate_index2 ( _seq + offset + index->index_length ) ;
		pair <TPWIi , TPWIi> pp = equal_range ( index->position[i1].begin() , index->position[i1].end() , findme ) ;

		for ( TPWIi it = pp.first ; it < pp.second ; it++ ) {
			int p = it->position ;
			if ( it->reverse_complement ) {
				if ( p+offset < l ) continue ;
				if ( rc.empty() ) { // Creating the reverse-complement
					rc = string ( _seq ) ;
					REVERSE_COMPLEMENT_CHARP ( (char*)rc.c_str() ) ;
				}
				char *chrseq = (char*) (*chrs)[it->chromosome].sequence.c_str() ;
				if ( !matchIUPAC ( rc.c_str() , chrseq+p-l+offset ) ) continue ;
				p += offset ;
			} else {
				if ( p + l > (*chrs)[it->chromosome].sequence.length() ) continue ;
				if ( p < offset ) continue ;
				char *chrseq = (char*) (*chrs)[it->chromosome].sequence.c_str() ;
				if ( !matchIUPAC ( _seq , chrseq+p-offset ) ) continue ;
				p -= offset ;
			}
			results.push_back ( *it ) ;
			if ( offset ) results.back().position = p ;
		}
	}
}

void TChromosomeAlign::run_wobble ( const string &seq1 , const string &seq2 , const string &qual1 , const string &qual2 , const pwi_vector &res1 , const pwi_vector &res2 , uint fragment_length ) {
	uint cnt , frag , range = fragment_range ;
	uint read_length_1 = seq1.length() ;
	uint insertion_length = 10000 ;

	if ( res1.size() + res2.size() > 10 ) return ; // Too many (number's just a whim)
	
	if ( res1.size() * res2.size() == 0 ) { // One read matches perfectly but the other does not
	
		if ( res1.size() + res2.size() == 1 ) { // Single match for one side, none for the other
			int the_pos , the_chr , match_pos ;
			string out ;
			bool rc ;
			if ( res1.size() == 1 ) {
				rc = res1[0].reverse_complement ;
				match_pos = res1[0].position ;
				the_pos = res1[0].position + ( fragment_length - seq1.length() * 2 ) * ( rc ? -1 : 1 ) ;
				the_chr = res1[0].chromosome ;
				out = !rc ? dorc ( seq2 ) : seq2 ;
			} else {
				rc = res2[0].reverse_complement ;
				match_pos = res2[0].position ;
				the_pos = res2[0].position - ( fragment_length - seq1.length() * 2 ) * ( rc ? -1 : 1 ) ;
				the_chr = res2[0].chromosome ;
				out = rc ? dorc ( seq1 ) : seq1 ;
			}
			fprintf ( wobble , "MIS\t%s\t%d\t%d\t%s\n" , (*chrs)[the_chr].name.c_str() , match_pos , the_pos , out.c_str() ) ;
		} else if ( res1.size() == 1 || res2.size() == 1 ) {
			pwi_vector single = res1.size() == 1 ? res1 : res2 ;
			pwi_vector multi = res1.size() == 1 ? res2 : res1 ;
			bool good = true ;
			int avg_pos = 0 ;
			int avg_cnt = 0 ;
			for ( int a = 0 ; a < multi.size() ; a++ ) {
				if ( single[0].chromosome != multi[a].chromosome ) {
					good = false ;
					break ;
				}
				if ( abs ( ((int)single[0].position) - ((int)multi[a].position) ) > fragment_length * 2 ) {
					good = false ;
					break ;
				}
				if ( single[0].reverse_complement == multi[a].reverse_complement ) {
					good = false ;
					break ;
				}
				avg_pos += multi[a].position ;
				avg_cnt++ ;
			}
			if ( good ) {
				int the_chr = single[0].chromosome ;
				int match_pos = single[0].position ;
				int the_pos = avg_pos / avg_cnt ;
				string out = res1.size() == 1 ? seq2 : seq1 ;
				out = multi[0].reverse_complement ? dorc ( out ) : out ;
				fprintf ( wobble , "MUL\t%s\t%d\t%d\t%s\n" , (*chrs)[the_chr].name.c_str() , match_pos , the_pos , out.c_str() ) ;
			}
		}
	
		if ( sqlite ) sqlite_out.clear() ;
	
		range = fragment_range ;
//		if ( res1.size() + res2.size() == 1 ) {
			if ( res1.size() > 0 ) wobble4snps ( res1 , seq2 , fragment_length , range , seq1 , false ) ;
			if ( res2.size() > 0 ) wobble4snps ( res2 , seq1 , fragment_length , range , seq2 , true ) ;
//		}

		if ( sqlite_out.size() == 1 ) {
			add2sqlite_cache ( sqlite_out[0] ) ;
		}

		return ;
	}
	return ;
	
	
	uint a ;
	cout << endl << seq1 << "\t" ;
	for ( a = 0 ; a < res1.size() ; a++ ) cout << res1[a].chromosome << "@" << res1[a].position << ":" << res1[a].reverse_complement << ", " ;
	cout << endl ;
	cout << seq2 << "\t" ;
	for ( a = 0 ; a < res2.size() ; a++ ) cout << res2[a].chromosome << "@" << res2[a].position << ":" << res2[a].reverse_complement << ", " ;
	cout << endl ;

	uint b ;
	for ( a = 0 ; a < res1.size() ; a++ ) {
		for ( b = 0 ; b < res2.size() ; b++ ) {
			if ( res1[a].chromosome != res2[b].chromosome ) continue ;
			
			TPositionWithIndex pa , pb ;
			if ( res1[a].position > res2[b].position ) {
				pa = res2[b] ;
				pb = res1[a] ;
			} else {
				pa = res1[a] ;
				pb = res2[b] ;
			}
			
			pa.position += (int) ( pa.reverse_complement ? -1 : 1 ) * seq1.length() ;
			pb.position += (int) ( pb.reverse_complement ? -1 : 1 ) * seq2.length() ;
			
//			uint diff = pb.position - pa.position ;
			fprintf ( wobble , "DEL\t%s\t%d\t%d\t%d\n" , (*chrs)[pa.chromosome].name.c_str() , pa.position , pb.position , pb.position - pa.position ) ;
//			cout << res1[a].position << ":" << res1[a].reverse_complement << " <=> " ;
//			cout << res2[b].position << ":" << res2[b].reverse_complement << "\t" << diff << endl ;
		}
	}


return;

/*

	// Deletions
	frag = fragment_length / 2 ;
	range = frag - read_length_1 ;
	wobble_results.clear() ;
	cnt = paired_read_combine ( res1 , res2 , frag , range , read_length_1 , seq1 , qual1 , seq2 , qual2 , true ) ;
	while ( wobble_results.size() ) {
		TPositionWithIndex p1 , p2 ;
		p2 = wobble_results.back() ; // Reverse order
		wobble_results.pop_back() ;
		p1 = wobble_results.back() ; // than pushing
		wobble_results.pop_back() ;
		if ( p1.position > p2.position ) { TPositionWithIndex p3 = p1 ; p1 = p2 ; p2 = p1 ; }
		fprintf ( wobble , "DEL\t%s\t%d\t%d\t%d\n" , (*chrs)[p1.chromosome].name.c_str() , p1.position , p1.position + fragment_length , p2.position - p1.position ) ;
	}
	
	
	// Insertions
	frag = fragment_length + insertion_length / 2 ;
	range = insertion_length / 2 ;
	wobble_results.clear() ;
	cnt = paired_read_combine ( res1 , res2 , frag , range , read_length_1 , seq1 , qual1 , seq2 , qual2 , true ) ;
	while ( wobble_results.size() ) {
		TPositionWithIndex p1 , p2 ;
		p2 = wobble_results.back() ; // Reverse order
		wobble_results.pop_back() ;
		p1 = wobble_results.back() ; // than pushing
		wobble_results.pop_back() ;
		if ( p1.position > p2.position ) { TPositionWithIndex p3 = p1 ; p1 = p2 ; p2 = p1 ; }
		fprintf ( wobble , "INS\t%s\t%d\t%d\t%d\n" , (*chrs)[p1.chromosome].name.c_str() , p1.position , p1.position + fragment_length , p2.position - p1.position ) ;
	}
	
	*/
}

void TChromosomeAlign::wobble4snps ( const pwi_vector &res , const string &seq , uint fragment_size , uint range , const string &anchor_seq , bool rc ) {
	if ( res.empty() ) return ;
	
	string seq_rc = dorc ( seq ) ;
	int from , to , c ;
	vector <string> out ;
	char dummy[READ_CHAR_BUFFER] ;
	int total_matches = 0 ;
	
//	string rc = dorc ( seq ) ;
	int sl = seq.length() ;
	for ( uint a = 0 ; a < res.size() ; a++ ) {
		string *cseq = &(*chrs)[res[a].chromosome].sequence ;
		if ( !(*chrs)[res[a].chromosome].original_sequence.empty() ) cseq = &(*chrs)[res[a].chromosome].original_sequence ;
//		if ( !cseq->empty() ) cseq = &(*chrs)[res[a].chromosome].sequence ;
		if ( res[a].reverse_complement ) from = res[a].position - fragment_size - range ;
		else from = res[a].position + fragment_size - range ;
		if ( from < 0 ) from = 0 ;
		to = from + range * 2 ;
		if ( to >= cseq->length() ) to = cseq->length() - 1 ;
		to -= sl ;

		// Add general variation range
//		fprintf ( wobble , "VAR\t%s\t%d\t%d\n" , (*chrs)[res[a].chromosome].name.c_str() , (int) ((from+to)/2-seq.length()/2) , (int) ((from+to)/2+seq.length()/2) ) ;
//		out.push_back ( dummy ) ;
		
		for ( int b = from ; b <= to ; b++ ) {
			vector <int> pos ;
			vector <char> o , n ;
			pos.reserve ( wobblemax ) ;
			o.reserve ( wobblemax ) ;
			n.reserve ( wobblemax ) ;
			const char *t = res[a].reverse_complement ? seq.c_str() : seq_rc.c_str() ;
			int mismatches = 0 ;
			const char *orig = cseq->c_str() + b ;
			for ( c = 0 ; c < sl ; c++ , t++ , orig++ ) {
				if ( char2IUPAC[*t] & char2IUPAC[*orig] ) continue ;
				if ( ++mismatches > wobblemax ) break ;
				pos.push_back ( b + c + 1 ) ;
				o.push_back ( *orig ) ;
				n.push_back ( *t ) ;
			}
			if ( mismatches == 0 ) continue ; // Weird...
			if ( mismatches > wobblemax ) continue ;
			
//			total_matches++ ;
//			if ( total_matches > 1 ) break ;
			if ( ++total_matches > 1 ) return ;
			
			for ( c = 0 ; c < pos.size() ; c++ ) {
				sprintf ( dummy , "SNP\t%s\t%d\t%c\t%c\n" , (*chrs)[res[a].chromosome].name.c_str() , pos[c] , o[c] , n[c] ) ;
				out.push_back ( dummy ) ;
			}
			
			if ( sqlite ) {
				uint opos = res[a].position ;
				if ( res[a].reverse_complement ) opos -= anchor_seq.length() ;
				string s = generate_sqlite_paired_command ( res[a].chromosome , 
								res[a].reverse_complement ? dorc ( anchor_seq ) : anchor_seq , 
								opos , 
								res[a].reverse_complement ? seq : seq_rc , 
								b , 'S' ) ;
/*				cout << s << endl ;
				cout << ( res[a].reverse_complement ? "RC" : "NORM" ) << endl ;
				cout << anchor_seq << endl ;
				cout << dorc ( anchor_seq ) << endl ;
				cout << dor ( anchor_seq ) << endl ;
				cout << doc ( anchor_seq ) << endl ;
				cout << cseq->substr ( opos , anchor_seq.length() ) << endl ;
				cout << endl ;*/
				
				sqlite_out.push_back ( s ) ;
			}
		}
	}
	
	if ( total_matches > 1 ) return ;
	for ( c = 0 ; c < out.size() ; c++ ) {
		fprintf ( wobble , "%s" , out[c].c_str() ) ;
	}
}

string TChromosomeAlign::generate_sqlite_paired_command ( int chr , string seq1 , int pos1 , string seq2 , int pos2 , char mode , string add_keys , string add_values ) {
	string ret ;
	string *refp = (*chrs)[chr].original_sequence.empty() ? &((*chrs)[chr].sequence) : &((*chrs)[chr].original_sequence) ;
	string ref1 = refp->substr ( pos1 , seq1.length() ) ;
	string ref2 = refp->substr ( pos2 , seq2.length() ) ;
	if ( seq1 == ref1 && seq2 == ref2 && mode == ' ' ) { // No SNPs
		char tmp[1000] , *chrn = (char*) (*chrs)[chr].name.c_str() ;
		sprintf ( tmp , "/*P|%s|%9d*/ INSERT INTO %s_perfect_match (read_name,pos1,pos2) VALUES (\"%s\",%d,%d);" , chrn , pos1 , chrn , last_solexa_name.c_str() , pos1+1 , pos2+1 ) ;
		ret = tmp ;
	} else { // SNP(s)
		char tmp[1000] , *chrn = (char*) (*chrs)[chr].name.c_str() ;
		int a ;
		string out1 , out2 ;
		string table = "_snp_match" ;
		if ( mode == ' ' ) mode = 'S' ;
		else if ( mode == 'I' ) table = "_inversions" ;
		else if ( mode == '1' ) table = "_single" ;
		if ( seq1 != ref1 ) {
			out1 = seq1 ;
			for ( a = 0 ; a < out1.length() ; a++ ) out1[a] = seq1[a] != ref1[a] ? seq1[a] : seq1[a] - 'A' + 'a' ;
		}
		if ( seq2 != ref2 && mode != '1' ) {
			out2 = seq2 ;
			for ( a = 0 ; a < out2.length() ; a++ ) out2[a] = seq2[a] != ref2[a] ? seq2[a] : seq2[a] - 'A' + 'a' ;
		}
		
		if ( mode == '1' ) {
			sprintf ( tmp , "/*%c|%s|%9d*/ INSERT INTO %s%s (read_name,seq1,pos1) VALUES (\"%s\",\"%s\",%d);" , 
								mode , chrn , pos1 , chrn , table.c_str() , last_solexa_name.c_str() , 
								out1.c_str() , 
								pos1+1 ) ;
		} else {
			sprintf ( tmp , "/*%c|%s|%9d*/ INSERT INTO %s%s (read_name,seq1,pos1,seq2,pos2%s) VALUES (\"%s\",\"%s\",%d,\"%s\",%d%s);" , 
								mode , chrn , pos1 , chrn , table.c_str() , 
								add_keys.c_str(),
								last_solexa_name.c_str() , 
								out1.c_str() , 
								pos1+1 , 
								out2.c_str() , 
								pos2+1 ,
								add_values.c_str()
								) ;
		}
		ret = tmp ;
	}
	return ret ;
}

int TChromosomeAlign::paired_read_combine ( const pwi_vector &v1 , const pwi_vector &v2 , uint fragment_length , uint range , 
											uint read_length_1 , const string &seq1 , const string &qual1 , const string &seq2 , const string &qual2 , 
											bool wobbling , bool inverting , int order ) {
	int ret = 0 ;
	vector <prc_cache> cache ;
	uint a , b , c ;
	uint seq1_length = seq1.length() ;
	uint seq2_length = seq2.length() ;
	
	if ( force_one_unique_match && v1.size() > 1 && v2.size() > 1 ) return ret ;

	a = 0 ;
	b = 0 ;
	while ( a < v1.size() && b < v2.size() ) {
		int chr = v1[a].chromosome ;
		int l1 = v1[a].position + fragment_length - range - read_length_1 ;
		
		if ( l1 < v1[a].position + read_length_1 ) l1 = v1[a].position + read_length_1 ;

		while ( b < v2.size() && ( v2[b].chromosome < chr || ( v2[b].chromosome == chr && v2[b].position < l1 ) ) ) b++ ;
		if ( b == v2.size() ) break ;

		int l2 = v1[a].position + fragment_length + range ;


		for ( c = b ; c < v2.size() && chr == v2[c].chromosome && v2[c].position <= l2 ; c++ ) {
			if ( v2[c].position < l1 ) continue ;
			if ( inverting && v1[a].reverse_complement == v2[c].reverse_complement ) { // INVERSION! FROM MARS!!!
				if ( inversions ) fprintf ( inversions , "%s\t%d\t%d\n" , (*chrs)[v1[a].chromosome].name.c_str() , v1[a].position , v2[c].position ) ;
				if ( sqlite ) {
					string s , t ;
					if ( v1[a].reverse_complement ) {
						t = ",\"--\"" ;
						t[2] = '0' + order ;
						s = generate_sqlite_paired_command ( chr , dorc ( seq1 ) , v1[a].position - seq1_length , 
																	dorc ( seq2 ) , v2[c].position - seq2_length , 
																	'I' , ",dir" , t ) ;
					} else {
						t = ",\"++\"" ;
						t[2] = '0' + order ;
						s = generate_sqlite_paired_command ( chr , seq1 , v1[a].position , seq2 , v2[c].position , 'I' , ",dir" , t ) ;
					}
					sqlite_out.push_back ( s ) ;
				}

				if ( sam ) {
					int position_a = v1[a].reverse_complement ? v1[a].position - seq1_length : v1[a].position ;
					int position_b = v2[c].reverse_complement ? v2[c].position - seq2_length : v2[c].position ;
					int dist = v2[c].position > v1[a].position ? v2[c].position - v1[a].position : v1[a].position - v2[c].position ;
					bool unique = true ; // Not true, but too complicated to fix...
					if ( v1[a].reverse_complement ) {
						add_sam ( dorc ( seq1 ) , dor ( qual1 ) , chr , position_a , 1 , position_b , dist , true , true , true , unique ) ;
						add_sam ( dorc ( seq2 ) , dor ( qual2 ) , chr , position_b , 2 , position_a , dist , true , true , true , unique ) ;
					} else {
						add_sam ( seq1 , qual1 , chr , position_a , 1 , position_b , dist , false , true , false , unique ) ;
						add_sam ( seq2 , qual2 , chr , position_b , 2 , position_a , dist , false , true , false , unique ) ;
					}
					sam_wrote = true ;
				}
				continue ;
			}
			if ( v1[a].reverse_complement == v2[c].reverse_complement ) continue ;

			ret++ ;
			
			// Add to output cache
			prc_cache ca ;
			ca.a = v1[a] ;
			ca.b = v2[c] ;
			cache.push_back ( ca ) ;

			if ( gffout ) {
				fprintf ( gffout , "\"%s\"\t\"SNP-O-MATIC\"\t\"misc_feature\"\t%d\t%d\t\"-\"\t\"%c\"\t\"-\"\t\"ID=%s_a;Gap=M%d\"\n" , 
															(*chrs)[v1[a].chromosome].name.c_str() , 
															v1[a].position , 
															v1[a].position + read_length_1 , 
															v1[a].reverse_complement ? '-' : '+' ,
															last_solexa_name.c_str() ,
															read_length_1
															) ;
				fprintf ( gffout , "\"%s\"\t\"SNP-O-MATIC\"\t\"misc_feature\"\t%d\t%d\t\"-\"\t\"%c\"\t\"-\"\t\"ID=%s_b;Gap=M%d\"\n" , 
															(*chrs)[v2[c].chromosome].name.c_str() , 
															v2[c].position , 
															v2[c].position + read_length_1 , 
															v2[c].reverse_complement ? '-' : '+' ,
															last_solexa_name.c_str() ,
															read_length_1
															) ;
			}


			if ( wobbling ) {
				TPositionWithIndex p1 , p2 ;
				p1.chromosome = chr ;
				p2.chromosome = chr ;
				if ( v1[a].reverse_complement ) p1.position = v1[a].position - seq1_length ;
				else p1.position = v1[a].position ;
				if ( v2[c].reverse_complement ) p2.position = v2[c].position - seq2_length ;
				else p2.position = v2[c].position ;
				wobble_results.push_back ( p1 ) ;
				wobble_results.push_back ( p2 ) ;
			}

		}
		a++ ;
	}
	
	if ( multimatch && cache.size() > 1 ) { // Force to map the read in one of its possible locations
		cache[0] = cache[rand()%cache.size()] ;
		cache.resize ( 1 ) ;
	}
	
	if ( singlematch && cache.size() > 1 ) return ret ;
	
	for ( uint n = 0 ; n < cache.size() ; n++ ) {
		prc_cache *p = &cache[n] ;
		int chr = p->a.chromosome ; // Same as p->b.chromosome
		
		int position_a = p->a.reverse_complement ? p->a.position - seq1_length : p->a.position ;
		int position_b = p->b.reverse_complement ? p->b.position - seq2_length : p->b.position ;
		int dist = p->b.position > p->a.position ? p->b.position - p->a.position : p->a.position - p->b.position ;

		if ( faceaway ) {
			if ( p->a.reverse_complement && !p->b.reverse_complement && p->a.position < p->b.position ) {
				fprintf ( faceaway , "%s\t%s" , last_solexa_name.c_str() , (*chrs)[p->a.chromosome].name.c_str() ) ;
				fprintf ( faceaway , "\t%d" , position_a ) ;
				fprintf ( faceaway , "\t%s" , p->a.reverse_complement ? dorc(seq1).c_str() : seq1.c_str() ) ;
				fprintf ( faceaway , "\t%d" , position_b ) ;
				fprintf ( faceaway , "\t%s" , p->b.reverse_complement ? dorc(seq2).c_str() : seq2.c_str() ) ;
				fprintf ( faceaway , "\n" ) ;
			}
		}
		
		if ( rpa ) {
			rpa_out ( p->a , p->b , read_length_1 , seq1 , qual1 , seq2 , qual2 ) ;
		}

		if ( sqlite ) {
			string s ;
			if ( p->a.reverse_complement ) s = generate_sqlite_paired_command ( chr , dorc ( seq1 ) , position_a , seq2 , position_b ) ;
			else s = generate_sqlite_paired_command ( chr , seq1 , position_a , dorc ( seq2 ) , position_b ) ;
			sqlite_out.push_back ( s ) ;
		}
		
		if ( fragmentplot ) {
			while ( fragment_stats.size() <= dist ) fragment_stats.push_back ( 0 ) ;
			fragment_stats[dist]++ ;
		}
		if ( indelplot ) {
			fprintf ( indelplot , "%s\t%d\t%d\t%d\n" , (*chrs)[p->a.chromosome].name.c_str() , p->a.position , p->b.position , dist ) ;
		}
		if ( coverage ) {
			if ( p->a.reverse_complement ) (*chrs)[chr].add_coverage ( dorc ( seq1 ) , position_a , true ) ;
			else (*chrs)[chr].add_coverage ( seq1 , position_a , false ) ;
			if ( p->b.reverse_complement ) (*chrs)[chr].add_coverage ( dorc ( seq2 ) , position_b , true ) ;
			else (*chrs)[chr].add_coverage ( seq2 , position_b , false ) ;
		}
		if ( sam ) {
			bool unique = cache.size() == 1 ;
			if ( p->a.reverse_complement ) add_sam ( dorc ( seq1 ) , dor ( qual1 ) , chr , position_a , 1 , position_b , dist , p->a.reverse_complement , true , p->b.reverse_complement , unique ) ;
			else add_sam ( seq1 , qual1 , chr , position_a , 1 , position_b , dist , p->a.reverse_complement , true , p->b.reverse_complement , unique ) ;
			if ( p->b.reverse_complement ) add_sam ( dorc ( seq2 ) , dor ( qual2 ) , chr , position_b , 2 , position_a , dist , p->b.reverse_complement , true , p->a.reverse_complement , unique ) ;
			else add_sam ( seq2 , qual2 , chr , position_b , 2 , position_a , dist , p->b.reverse_complement , true , p->a.reverse_complement , unique ) ;
			sam_wrote = true ;
		}
		if ( snpsinreads ) {
			if ( p->a.reverse_complement ) add_snpsinreads ( dorc ( seq1 ) , dor ( qual1 ) , chr , position_a , 1 ) ;
			else add_snpsinreads ( seq1 , qual1 , chr , position_a , 1 ) ;
			if ( p->b.reverse_complement ) add_snpsinreads ( dorc ( seq2 ) , dor ( qual2 ) , chr , position_b , 2 ) ;
			else add_snpsinreads ( seq2 , qual2 , chr , position_b , 2 ) ;
		}
		if ( pileup ) {
			if ( p->a.reverse_complement ) (*chrs)[chr].ao.add_align ( dorc ( seq1 ) , dor ( qual1 ) , position_a , chr ) ;
			else (*chrs)[chr].ao.add_align ( seq1 , qual1 , position_a , chr ) ;
			if ( p->b.reverse_complement ) (*chrs)[chr].ao.add_align ( dorc ( seq2 ) , dor ( qual2 ) , position_b , chr ) ;
			else (*chrs)[chr].ao.add_align ( seq2 , qual2 , position_b , chr ) ;
		}
		if ( featuretable ) {
			fprintf ( featuretable , "gene\t%d-%d\n\t/product=\"%s\"\n\t/chromosome=\"%s\"\n" , p->a.position , p->a.position + read_length_1 , last_solexa_name.c_str() , (*chrs)[p->a.chromosome].name.c_str() ) ;
			fprintf ( featuretable , "gene\t%d-%d\n\t/product=\"%s\"\n\t/chromosome=\"%s\"\n" , p->b.position , p->b.position + read_length_1 , last_solexa_name.c_str() , (*chrs)[p->b.chromosome].name.c_str() ) ;
		}
	}
	
	return ret ;
}

void TChromosomeAlign::rpa_out ( const TPositionWithIndex &a , const TPositionWithIndex &b , uint read_length_1 , const string &seq1 , const string &qual1 , const string &seq2 , const string &qual2 ) {
	if ( !rpa_header_written ) {
		fprintf ( rpa , "#Readname\tr1_seq\tr1_qual\tr1_chr\tr1_pos\tr1_rc\tr2_seq\tr2_qual\tr2_chr\tr2_pos\tr2_rc\n" ) ;
		rpa_header_written = true ;
	}
	fprintf ( rpa , "%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\n" ,
		last_solexa_name.c_str() ,
		seq1.c_str() ,
		qual1.c_str() ,
		(*chrs)[a.chromosome].name.c_str() ,
		a.position ,
		a.reverse_complement ? 1 : 0 ,
		seq2.c_str() ,
		qual2.c_str() ,
		(*chrs)[b.chromosome].name.c_str() ,
		b.position ,
		b.reverse_complement ? 1 : 0
	) ;
}

void TChromosomeAlign::add_sam ( const string &seq , const string &quality , int chr , int pos , int readpart , int matepos , int ins_size , bool rc , bool as_pair , bool mate_rc , bool unique_match ) {
	char cigar[100] ;
	sprintf ( cigar , "%dM" , (int) seq.length() ) ;
	
	unsigned int flag = 0 ;
	flag |= 0x0001 ; // Paired in sequencing
	if ( as_pair ) flag |= 0x0002 ; // Mapped as pair
	else flag |= 0x0008 ; // Mate unmapped
	if ( rc ) flag |= 0x0010 ; // On reverse strand
	if ( as_pair ) {
		if ( mate_rc ) flag |= 0x0020 ; // Mate RC
		if ( 1 == readpart ) flag |= 0x0040 ; // First read in a pair
		if ( 2 == readpart ) flag |= 0x0080 ; // Second read in a pair
	}
	if ( !unique_match ) flag |= 0x0100 ; // Not uniquely mapped
	
	char mismatches[1000] ;
	sprintf ( mismatches , "\tNM:i:%d" , count_snpsinreads ( seq , chr , pos ) ) ;
	string tags ; // FILL ME
//	tags += "\tRG:Z:" + rgid ;
//	tags += "\tPU:Z:" + rgid ;
	tags += "\tPG:Z:SNP-o-matic" ;
	tags += mismatches ;
	
	if ( pos > matepos ) ins_size = -ins_size ;
	
	fprintf ( sam , "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s%s\n" ,
				last_solexa_name.c_str() , // Read name
				flag , // flag
				(*chrs)[chr].name.c_str() , // Chromosome name
				pos + 1 , // Alignment start position
				255 , // PHRED alignment score dummy
				cigar , // CIGAR
				"=" , // Mate
				matepos , // Mate position
				ins_size , // ISIZE (actually, fragment size)
				seq.c_str() , // Sequence
				quality.c_str() , // Quality
				tags.c_str() // Tag
				) ;
}

int TChromosomeAlign::count_snpsinreads ( const string &seq , int chr , int pos ) {
	uint sl = (*chrs)[chr].original_sequence.length() ;
	if ( sl == 0 ) return 0 ; // No SNPs here
	
	int ret = 0 ;
	int seql = seq.length() , n , cnt ;
	const char *os = (*chrs)[chr].original_sequence.c_str() ;
	const char *ss = (*chrs)[chr].sequence.c_str() ;
	for ( int p = 0 ; p < seql && p + pos < sl ; p++ ) {
		if ( !isIUPAC[ss[p+pos]] || os[p+pos] == 'X' ) continue ;
		if ( seq[p] != os[p+pos] ) ret++ ;
	}
	return ret ;
}

void TChromosomeAlign::add_snpsinreads ( const string &seq , const string &quality , int chr , int pos , int readpart ) {
	uint sl = (*chrs)[chr].original_sequence.length() ;
	if ( sl == 0 ) return ; // No SNPs here
	
	if ( snpsinreads_first ) { // Output header
		snpsinreads_first = false ;
		fprintf ( snpsinreads , "#read_name\tread_part\tsnp_name\tchromosome\tposition_on_chromosome\tposition_on_read\tobserved_base\treference_base\tbase_quality\tmin_total\tavg_total\tmin_5\tavg_5\n" ) ;
	}
	
	int min_total = -1 , avg_total , min_5 , avg_5 ;
	int seql = seq.length() , n , cnt ;
	
	const char *os = (*chrs)[chr].original_sequence.c_str() ;
	const char *ss = (*chrs)[chr].sequence.c_str() ;
	for ( int p = 0 ; p < seql && p + pos < sl ; p++ ) {
//		if ( os[p+pos] == seq[p] || os[p+pos] == 'X' ) continue ;
		if ( !isIUPAC[ss[p+pos]] || os[p+pos] == 'X' ) continue ;
		
		// Preparing total quality scores
		if ( min_total == -1 ) {
			min_total = 1000 ; // Some high number
			avg_total = 0 ;
			for ( n = 0 ; n < seql ; n++ ) {
				if ( quality[n] < min_total ) min_total = quality[n] ;
				avg_total += quality[n] ;
			}
			min_total -= 33 ;
			avg_total = avg_total / seql - 33 ;
		}
		
		// Prepare surrounding quality scores
		min_5 = 1000 ; // Some high number
		avg_5 = 0 ;
		cnt = 0 ;
		for ( n = ( p >= 5 ) ? p - 5 : 0 ; n < p + 5 && n < seql ; n++ ) {
			if ( n == p ) continue ; // Do not use actual SNP base for quality
			cnt++ ;
			if ( quality[n] < min_5 ) min_5 = quality[n] ;
			avg_5 += quality[n] ;
		}
		if ( cnt > 0 ) {
			min_5 -= 33 ;
			avg_5 = avg_5 / cnt - 33 ;
		} else {
			avg_5 = 0 ;
			cerr << "TChromosomeAlign::add_snpsinreads : cnt == 0\n" ;
		}
		
		fprintf ( snpsinreads ,
							"%s\t%d\t%s:%d\t%s\t%d\t%d\t%c\t%c\t%d\t%d\t%d\t%d\t%d\n" , 
							last_solexa_name.c_str() ,
							readpart ,
							(*chrs)[chr].name.c_str() ,
							pos+p+1 ,
							(*chrs)[chr].name.c_str() ,
							pos+p+1 ,
							p+1 ,
							seq[p] ,
							os[p+pos] ,
							(unsigned char) quality[p] - 33 ,
							min_total ,
							avg_total ,
							min_5 ,
							avg_5
							) ;
	}
}

void TChromosomeAlign::add2sqlite_cache ( string s ) {
	if ( sqlite_cache_name.empty() ) {
		char *c = new char[10000] ;
		strcpy ( c , "sqlite_cache_XXXXXXXX" ) ;
		close ( mkstemp ( c ) ) ;
		//tmpnam ( c ) ;
		sqlite_cache_name = c ;
		sqlite_cache = fopen ( c , "w" ) ;
		delete c ;
	}
	fprintf ( sqlite_cache , "%s\n" , s.c_str() ) ;
}

void TChromosomeAlign::finish_sqlite ( int fragment , int variance , int read_length ) {
	if ( !sqlite ) return ;
	int a , b ;
	fprintf ( sqlite , "BEGIN EXCLUSIVE TRANSACTION;\n" ) ;
	
	// Meta-data
	fprintf ( sqlite , "CREATE TABLE meta ( key VARCHAR[64] , value VARCHAR[256] ) ;\n" ) ;
	fprintf ( sqlite , "INSERT INTO meta ( key , value ) VALUES ( \"dbversion\" , \"3\" ) ;\n" ) ;
	fprintf ( sqlite , "INSERT INTO meta ( key , value ) VALUES ( \"%s\" , \"%d\" ) ;\n" , "fragment" , fragment ) ;
	fprintf ( sqlite , "INSERT INTO meta ( key , value ) VALUES ( \"%s\" , \"%d\" ) ;\n" , "variance" , variance ) ;
	fprintf ( sqlite , "INSERT INTO meta ( key , value ) VALUES ( \"%s\" , \"%d\" ) ;\n" , "read_length" , read_length ) ;
	fprintf ( sqlite , "INSERT INTO meta ( key , value ) VALUES ( \"%s\" , \"%s\" ) ;\n" , "name_prefix" , sqlite_prefix.c_str() ) ;

	// Create tables
	fprintf ( sqlite , "CREATE TABLE chromosomes ( name VARCHAR[256] , size INTEGER ) ;\n" ) ;
	for ( a = 0 ; a < chrs->size() ; a++ ) {
		fprintf ( sqlite , "CREATE TABLE %s_perfect_match ( read_name VARCHAR[32], pos1 INTEGER, pos2 INTEGER ) ;\n" , (*chrs)[a].name.c_str() ) ;
		fprintf ( sqlite , "CREATE TABLE %s_snp_match ( read_name VARCHAR[32], pos1 INTEGER, seq1 VARCHAR[64] , pos2 INTEGER , seq2 VARCHAR[64] ) ;\n" , (*chrs)[a].name.c_str() ) ;
		fprintf ( sqlite , "CREATE TABLE %s_inversions ( read_name VARCHAR[32], pos1 INTEGER, seq1 VARCHAR[64] , pos2 INTEGER , seq2 VARCHAR[64] , dir VARCHAR[4] ) ;\n" , (*chrs)[a].name.c_str() ) ;
		fprintf ( sqlite , "CREATE TABLE %s_single ( read_name VARCHAR[32], pos1 INTEGER, seq1 VARCHAR[64] ) ;\n" , (*chrs)[a].name.c_str() ) ;
		fprintf ( sqlite , "INSERT INTO chromosomes ( name , size ) VALUES ( \"%s\" , \"%d\") ;\n" , (*chrs)[a].name.c_str() , (*chrs)[a].sequence.length() ) ;
	}
	
	// Write contents
//	if ( DEBUGGING ) { fprintf ( stdout , "Sorting sqlite output ... " ) ; fflush(stdout); }
//	sort ( sqlite_cache.begin() , sqlite_cache.end() ) ;
//	if ( DEBUGGING ) { fprintf ( stdout , "writing ... " ) ; fflush(stdout); }

	fclose ( sqlite_cache ) ;
	sqlite_cache = fopen ( sqlite_cache_name.c_str() , "r" ) ;
	char dummy[READ_CHAR_BUFFER] , *c1 , *c2 ;

	string pref = "\"" + sqlite_prefix ;
//	for ( a = 0 ; a < sqlite_cache.size() ; a++ ) {
//		string s = sqlite_cache[a] ;
	while ( !feof ( sqlite_cache ) ) {
		fgets ( dummy , READ_CHAR_BUFFER , sqlite_cache ) ;
		for ( c1 = dummy ; *c1 > 13 ; c1++ ) ;
		*c1 = 0 ;
		string s ( dummy ) ;
		if ( s.empty() ) continue ;
		
		// Removing initial SQL comment to save disk space; we could leave it in, and sqlite would ignore it
		b = s.find ( "*/" ) ;
		if ( b > 0 && b < s.length() ) {
			b += 2 ;
			while ( s[b] == ' ' ) b++ ;
			s = s.substr ( b ) ;
		}
		
		// Remove common read name prefix
		int position = s.find ( pref ) ; // find first space
		while ( position != string::npos ) {
			s.replace( position, pref.length(), "\"" );
			position = s.find( pref, position + 1 );
		} 
		
		fprintf ( sqlite , "%s\n" , s.c_str() ) ;
	}

	// Create indices
	fprintf ( sqlite , "COMMIT;\n" ) ;
	fprintf ( sqlite , "BEGIN EXCLUSIVE TRANSACTION;\n" ) ;
	for ( a = 0 ; a < chrs->size() ; a++ ) {
		char *c = (char*) (*chrs)[a].name.c_str() ;
		fprintf ( sqlite , "CREATE INDEX %s_perfect_index ON %s_perfect_match ( pos1 , pos2 ) ;\n" , c , c ) ;
		fprintf ( sqlite , "CREATE INDEX %s_snp_index ON %s_snp_match ( pos1 , pos2 ) ;\n" , c , c ) ;
		fprintf ( sqlite , "CREATE INDEX %s_inv_index ON %s_inversions ( pos1 , pos2 ) ;\n" , c , c ) ;
		fprintf ( sqlite , "CREATE INDEX %s_sin_index ON %s_single ( pos1 ) ;\n" , c , c ) ;
	}
	
	// Close up
	fprintf ( sqlite , "COMMIT;\n" ) ;
	fclose ( sqlite_cache ) ;
	remove ( sqlite_cache_name.c_str() ) ;
	if ( DEBUGGING ) { fprintf ( stdout , "done.\n" ) ; fflush(stdout); }
}

void TChromosomeAlign::align_solexa_paired ( string filename , string filename2 , uint read_length_1 , uint fragment_length , uint range ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Reading solexa pair data from %s ... " , filename.c_str() ) ; fflush(stdout); }

	FILE * file = fopen ( filename.c_str() , "r" ) ;
	if ( !file ) {
		cerr << "Could not open FASTQ file \"" << filename << "\" - aborting." << endl ;
		exit ( 1 ) ;
	}

	FILE *file2 = NULL ;
	if ( !filename2.empty() ) {
		file2 = fopen ( filename2.c_str() , "r" ) ;
		if ( !file2 ) {
			cerr << "Could not open FASTQ file \"" << filename2 << "\" - aborting." << endl ;
			exit ( 1 ) ;
		}
	}

	if ( sam ) {
		for ( int a = 0 ; a < filename.length() ; a++ ) {
			if ( filename[a] == '/' ) rgid.clear() ;
			else rgid += filename[a] ;
		}
		
		char predicted_insert_size[100] ;
		sprintf ( predicted_insert_size , "\tPI:%d" , fragment_length - read_length_1 * 2 ) ;
		
/*		string rg ( "@RG" ) ;
		rg += "\tID:" + rgid ;
		rg += "\tPU:" + rgid ;
		rg += predicted_insert_size ;
		rg += "\tPL:SOLEXA" ;
		fprintf ( sam , "%s\n" , rg.c_str() ) ;
*/
		string pg ( "@PG\tID:SNP-o-matic" ) ;
		fprintf ( sam , "%s\n" , pg.c_str() ) ;
	}

	char dummy[READ_CHAR_BUFFER] , *c1 , *c2 ;
	string lastseq , lastqual ;
	fasta = true ;
	bool last_was_quality = false , lastone = false , eof = false , skip = false ;
	uint readcount = 0 , matched_reads = 0 , total_matches = 0 , skipped = 0 ;

	if ( range == 0 ) range = fragment_length / 4 ;
	fragment_range = range ;
	
	int linecounter = 0 ;
	
	while ( !eof || lastone ) {
//		if ( readcount > 10000 ) break ; // TESTING
		if ( lastone ) {
			strcpy ( dummy , "@dummy" ) ;
			if ( fasta ) *dummy = '>' ;
		} else {
			*dummy = 0 ;
			fgets ( dummy , READ_CHAR_BUFFER , file ) ;
			if ( file2 ) {
				for ( c1 = dummy ; *c1 > 13 ; c1++ ) ;
				fgets ( c1 , READ_CHAR_BUFFER , file2 ) ;
			}
		}
		
		// Remove EOL
		for ( c1 = dummy ; *c1 > 13 ; c1++ ) ;
		*c1 = 0 ;
	
		if ( linecounter == 0 && ( ( fasta && *dummy == '>' ) || *dummy == '@' ) ) {
			if ( *dummy == '@' ) fasta = false ;
			
			if ( !skip ) {

				if ( sqlite && !last_solexa_name.empty() ) {
					if ( sqlite_prefix_first ) {
						sqlite_prefix_first = false ;
						sqlite_prefix = last_solexa_name ;
					} else {
						int a ;
						for ( a = 1 ; a < sqlite_prefix.length() && a < last_solexa_name.length() && last_solexa_name[a] == sqlite_prefix[a] ; a++ ) ;
						if ( sqlite_prefix.length() != a ) sqlite_prefix = sqlite_prefix.substr ( 0 , a ) ;
					}
				}

				for ( c1 = (char*) lastseq.c_str() ; *c1 && !isIUPAC[*c1] ; c1++ ) ;
				if ( *c1 ) { // Contains IUPAC
					if ( using_bins ) write2bin ( last_solexa_name , lastseq , lastqual , -1 ) ;
				} else if ( !lastseq.empty() ) {
					uint matches = align_solexa_paired_read ( last_solexa_name , lastseq , lastqual , read_length_1 , fragment_length , range ) ;
					if ( matches > 0 ) matched_reads++ ;
					total_matches += matches ;
					if ( using_bins ) write2bin ( last_solexa_name , lastseq , lastqual , (int) matches ) ;
				}
			}

			if ( lastone ) break ;
			
			// Clear data
			c1 = dummy ;
			if ( *c1 == '>' || *c1 == '@' ) c1++ ;
			while ( *c1 == ' ' ) c1++ ;
			last_solexa_name = c1 ;
			readcount++ ;
			last_was_quality = false ;
			lastseq.clear() ;
			lastqual.clear() ;
			
			if ( use_nono ) {
				skip = is_nono ( last_solexa_name ) ;
				if ( skip ) skipped++ ;
			}
		} else if ( *dummy == '+' ) {
			last_was_quality = true ;
		} else if ( last_was_quality ) {
			lastqual += dummy ;
		} else {
			lastseq += dummy ;
		}
		
		linecounter++ ;
		linecounter %= 4 ;
		
		eof = feof ( file ) ;
		if ( eof ) lastone = true ; // Fake last read
	}

	if ( fragmentplot ) {
		for ( uint a = 0 ; a < fragment_stats.size() ; a++ ) {
			fprintf ( fragmentplot , "%d\t%d\n" , a , fragment_stats[a] ) ;
		}
	}

	if ( DEBUGGING ) { fprintf ( stdout , "scanned %d solexa reads, %d (%2.2f%%) matched , %d total matches, %d skipped.\n" , readcount , matched_reads , (float) matched_reads * 100 / readcount , total_matches , skipped ) ; fflush(stdout); }
}

void TChromosomeAlign::align_contigs ( chrvec &contigs , TChromosomalIndices &contig_index , uint maxerr , string output_filename ) {
	uint a , b , c ;
	vector <pwi_vector> vp ;
	vector <vuint> vi ;
	for ( a = 0 ; a < contig_index.position.size() ; a++ ) {
		for ( b = 0 ; b < contig_index.position[a].size() ; b++ ) {
			uint chr = contig_index.position[a][b].chromosome ;
			while ( vp.size() <= chr ) vp.push_back ( pwi_vector () ) ;
			while ( vi.size() <= chr ) vi.push_back ( vuint () ) ;
			vp[chr].push_back ( contig_index.position[a][b] ) ;
			vi[chr].push_back ( a ) ;
		}
	}

	uint mismatch = 0 ;
	vector <string> out ;
	char c2[50000] ;
	
	for ( a = 0 ; a < vp.size() ; a++ ) {
		string refseq = contigs[a].sequence ;
		for ( b = 0 ; b < vp[a].size() ; b++ ) {
			vector <string> matches ;
			uint i1 = vi[a][b] ;
			pair <TPWIi , TPWIi> pp = equal_range ( index->position[i1].begin() , index->position[i1].end() , vp[a][b] ) ;

			for ( TPWIi it = pp.first ; it < pp.second ; it++ ) {
				uint chr = it->chromosome ;
				string s ;
				uint from , len = contigs[a].sequence.length() ;
				if ( vp[a][b].reverse_complement ) {
					if ( it->reverse_complement ) {
						from = it->position - index->index_length_double ;
						from -= vp[a][b].position - index->index_length_double ;
					} else {
						from = it->position + index->index_length_double - contigs[a].sequence.length() ;
						from += vp[a][b].position - index->index_length_double ;
					}
				} else {
					if ( it->reverse_complement ) {
						from = it->position - contigs[a].sequence.length() ;
						from += vp[a][b].position ;
					} else {
						from = it->position ;
						from -= vp[a][b].position ;
					}
				}
				
				bool rc = false ;
				s = (*chrs)[chr].sequence.substr ( from , len ) ;
				if ( s.length() != refseq.length() ) continue ; // Sanity check
				
				if ( it->reverse_complement ) {
					REVERSE_COMPLEMENT_CHARP ( (char*)s.c_str() ) ;
					rc = !rc ;
				}
				if ( vp[a][b].reverse_complement ) {
					REVERSE_COMPLEMENT_CHARP ( (char*)s.c_str() ) ;
					rc = !rc ;
				}
				if ( s == refseq ) continue ; // No need to bother with a perfect match - or is there? Where did it come from?

				// Check for errors
				uint errcnt = 0 ;
				for ( c = 0 ; errcnt <= maxerr && c < s.length() ; c++ ) {
					if ( s[c] != refseq[c] ) {
						errcnt++ ;
					} else { // Lower case for all matching chars
						s[c] = s[c] - 'A' + 'a' ;
					}
				}
				if ( errcnt > maxerr ) continue ;

				sprintf ( c2 , "%s\t%d\t%d\t%s\t%c" , s.c_str() , chr , from , refseq.c_str() , rc ? 'X' : ' ' ) ;
				string s2 ( c2 ) ;
				for ( c = 0 ; c < matches.size() && matches[c] != s2 ; c++ ) ;
				if ( c < matches.size() ) continue ; // We have that one already
				matches.push_back ( s2 ) ;
			}
		
			int lastpos , lastchr = -1 ;
			for ( c = 0 ; c < matches.size() ; c++ ) {
				if ( c > 0 && matches[c-1] == matches[c] ) continue ;

				int d ;
				vector <string> vs ;
				Tokenize ( matches[c] , vs , "\t" ) ;
				int _chr = atoi ( vs[1].c_str() ) ;
				int _pos = atoi ( vs[2].c_str() ) ;
				if ( lastchr == _chr && lastpos == _pos ) continue ;
				
				lastchr = _chr ;
				lastpos = _pos ;
				
				bool rc = vs[4] == "X" ;
				
				int thepos ;
				for ( d = 0 ; d < vs[0].length() ; d++ ) {
					if ( vs[0][d] < 'A' || vs[0][d] > 'Z' ) continue ; // Perfect match
					if ( rc ) thepos = _pos + ( vs[0].length() - d ) ;
					else thepos = _pos + 1 + d ;
					char the_c = rc ? base_rc[vs[0][d]] : vs[0][d] ;
					if ( (*chrs)[_chr].sequence[thepos-1] != the_c ) mismatch++ ;
					if ( the_c == vs[3][d] ) continue ; // Pa-ra-no-ia
					sprintf ( c2 , "%s\t%d\t%c\t%c\n" , (*chrs)[_chr].name.c_str() , thepos , the_c , vs[3][d] ) ;
					out.push_back ( c2 ) ;
				}
			}
			
		}
	}
	
	sort ( out.begin() , out.end() ) ;
	FILE * file = fopen ( output_filename.c_str() , "w" ) ;
	for ( a = 0 ; a < out.size() ; a++ ) {
		if ( a > 0 && out[a] == out[a-1] ) continue ;
		fwrite ( out[a].c_str() , 1 , out[a].length() , file ) ;
	}
	fclose ( file ) ;
	
//	if ( DEBUGGING ) cerr << mismatch << " mismatches\n" ;
}


/// Gets Solexa reads from a FASTQ file and matches them against the sequence.
void TChromosomeAlign::align_solexa_fastq ( string filename ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Reading solexa data from %s ... " , filename.c_str() ) ; fflush(stdout); }

	FILE * file = fopen ( filename.c_str() , "r" ) ;
	if ( !file ) {
		cerr << "Could not open FASTQ file \"" << filename << "\" - aborting." << endl ;
		exit ( 1 ) ;
	}
	
	char dummy[READ_CHAR_BUFFER] , *c1 , *c2 ;
	int solexacount = 0 , total_matches = 0 , matched_reads = 0 ;
	string seq , q ;
	char *c ;
	bool invalid ;
	fasta = false ;

	if ( !pileup ) max_align = 2 ; // If we don't do the pileup but just bin, 2 matches are all we need (=> multiple_bin)

	// Run through file
	while ( !feof(file) ) {
		invalid = false ;
		*dummy = 0 ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		if ( !*dummy ) break ;
		
		// Read name
		if ( *dummy != '@' ) { cerr << "Expected '@'\n" << dummy ; exit ( 1 ) ; }
		for ( c = dummy ; *c > 13 ; c++ ) ;
		*c = 0 ;
		last_solexa_name = dummy + 1 ;
		
		// Sequence
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		for ( c = dummy ; !isIUPAC[*c] ; c++ ) ;
		if ( *c == 10 || *c == 13 ) *c = 0 ;
		if ( *c ) {
			invalid = true ;
			while ( *c > 13 ) c++ ;
			*c = 0 ;
		}
		seq = dummy ;

		fgets ( dummy , READ_CHAR_BUFFER , file ) ;

		if ( *dummy != '+' ) { cerr << "Expected '+'\n" ; exit ( 1 ) ; }
		for ( c = dummy ; *c > 13 ; c++ ) ;
		*c = 0 ;
		string qname = dummy + 1 ;
		if ( qname.empty() ) qname = last_solexa_name ;
		if ( last_solexa_name != qname ) { cerr << "Name mismatch\n" << last_solexa_name << " : " << qname ; exit ( 1 ) ; }
		
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		for ( c = dummy ; *c > 13 ; c++ ) ;
		*c = 0 ;
		q = dummy ;
		
		if ( !single_read_length ) single_read_length = seq.length() ;
		int m = invalid ? -1 : align_solexa_read ( (char*) seq.c_str() , q ) ;
		if ( m > 0 ) matched_reads++ ;
		if ( m > 0 ) total_matches += m ;
		if ( using_bins ) { // Write to bins
			if ( last_solexa_name.empty() ) {} // Paranoia
			else write2bin ( last_solexa_name , seq , q , m ) ;
/*			else if ( invalid ) fprintf ( binfile_iupac , "@%s\n%s\n+%s\n%s\n" , last_solexa_name.c_str() , seq.c_str() , last_solexa_name.c_str() , q.c_str() ) ;
			else if ( m == 0 ) fprintf ( binfile_no_match , "@%s\n%s\n+%s\n%s\n" , last_solexa_name.c_str() , seq.c_str() , last_solexa_name.c_str() , q.c_str() ) ;
			else if ( m == 1 ) fprintf ( binfile_single_match , "@%s\n%s\n+%s\n%s\n" , last_solexa_name.c_str() , seq.c_str() , last_solexa_name.c_str() , q.c_str() ) ;
			else fprintf ( binfile_multi_match , "@%s\n%s\n+%s\n%s\n" , last_solexa_name.c_str() , seq.c_str() , last_solexa_name.c_str() , q.c_str() ) ;*/
			last_solexa_name.clear() ;
		}
		if ( !invalid ) solexacount++ ;
	}
	fclose ( file ) ;
	
	if ( DEBUGGING ) { fprintf ( stdout , "scanned %d solexa reads, %d (%2.2f%%) matched , %d total matches.\n" , solexacount , matched_reads , (float) matched_reads * 100 / solexacount , total_matches ) ; fflush(stdout); }
}


/// Gets Solexa reads from a FASTA file and matches them against the sequence. Should not be used anymore.
void TChromosomeAlign::align_solexa_fasta ( string filename ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Reading solexa data from %s ... " , filename.c_str() ) ; fflush(stdout); }

	FILE * file = fopen ( filename.c_str() , "r" ) ;
	if ( !file ) {
		cerr << "Could not open FASTA file \"" << filename << "\" - aborting." << endl ;
		exit ( 1 ) ;
	}
	char dummy[READ_CHAR_BUFFER] , *c ;
	int solexacount = 0 , total_matches = 0 , matched_reads = 0 ;
	bool invalid ;
	fasta = true ;
	
	// Artificial quality strings
	vector <string> q ;
	for ( int a = 0 ; a < 100 ; a++ ) q.push_back ( string ( a , (char) 30+33 ) ) ;
	
	if ( !pileup ) max_align = 2 ; // If we don't do the pileup but just bin, 2 matches are all we need (=> multiple_bin)
	
	// Run through file
	while ( !feof(file) ) {
		*dummy = 0 ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		if ( !*dummy ) break ;
		if ( *dummy == '>' ) {
			if ( using_bins || cigar ) { // Writing to binfiles, storing solexa name
				for ( c = dummy ; *c > 13 ; c++ ) ;
				*c = 0 ;
				last_solexa_name = dummy + 1 ;
			}
			continue ;
		}
		
		invalid = false ;
		for ( c = dummy ; !isIUPAC[*c] ; c++ ) ;
		if ( *c == 10 || *c == 13 ) *c = 0 ;

		if ( *c ) {
			invalid = true ;
			while ( *c > 13 ) c++ ;
			*c = 0 ;
		}

		if ( !invalid && !single_read_length ) single_read_length = strlen ( dummy ) ;
		int m = invalid ? -1 : align_solexa_read ( dummy , q[strlen(dummy)] ) ;
		if ( m > 0 ) matched_reads++ ;
		if ( m > 0 ) total_matches += m ;
		if ( using_bins ) { // Write to bins
			if ( last_solexa_name.empty() ) {} // Paranoia
			else write2bin ( last_solexa_name , string ( dummy ) , "" , m ) ;
/*			else if ( invalid ) fprintf ( binfile_iupac , ">%s\n%s\n" , last_solexa_name.c_str() , dummy ) ;
			else if ( m == 0 ) fprintf ( binfile_no_match , ">%s\n%s\n" , last_solexa_name.c_str() , dummy ) ;
			else if ( m == 1 ) fprintf ( binfile_single_match , ">%s\n%s\n" , last_solexa_name.c_str() , dummy ) ;
			else fprintf ( binfile_multi_match , ">%s\n%s\n" , last_solexa_name.c_str() , dummy ) ;*/
			last_solexa_name.clear() ;
		}
		solexacount++ ;
	}
	fclose ( file ) ;
	
	if ( DEBUGGING ) { fprintf ( stdout , "scanned %d solexa reads, %d (%2.2f%%) matched , %d total matches.\n" , solexacount , matched_reads , (float) matched_reads * 100 / solexacount , total_matches ) ; fflush(stdout); }
}


int TChromosomeAlign::align_solexa_read ( char *_seq , const string &_quality , bool inv_compl ) {
	uint ret = 0 ;
	uint i1 = index->calculate_index ( _seq ) ;
	if ( !index->position[i1].size() ) return 0 ;
	string rc , rcq ;
	uint l = strlen ( _seq ) ;

	// Binary search
	TPositionWithIndex findme ;
	findme.index2 = index->calculate_index2 ( _seq + index->index_length ) ;
	pair <TPWIi , TPWIi> pp = equal_range ( index->position[i1].begin() , index->position[i1].end() , findme ) ;
	if ( sqlite ) sqlite_out.clear () ;

	for ( TPWIi it = pp.first ; it < pp.second ; it++ ) {
		uint chr = it->chromosome ;
		char *chrseq = (char*) (*chrs)[chr].sequence.c_str() ;
		int p = it->position ;
		if ( it->reverse_complement ) {
			if ( p < l ) continue ;
			if ( rc.empty() ) { // Creating the reverse-complement
				rc = string ( _seq ) ;
				REVERSE_COMPLEMENT_CHARP ( (char*)rc.c_str() ) ;
			}
			if ( !matchIUPAC ( rc.substr(0,l-index->index_length_double).c_str() , chrseq+p-l ) ) continue ;
			if ( cigar || sqlite ) add_cigar_sqlite ( rc , last_solexa_name , p-l , chr , '-' ) ;
			if ( coverage ) (*chrs)[chr].add_coverage ( rc , p-l , true ) ;
			if ( pileup ) {
				if ( rcq.empty() ) {
					rcq = _quality ;
					reverse ( rcq.begin() , rcq.end() ) ;
				}
				(*chrs)[chr].ao.add_align ( rc , rcq , p-l , chr ) ;
			}
			if ( snpsinreads ) add_snpsinreads ( rc , rcq , chr , p-l , 1 ) ;
			if ( sqlite ) sqlite_out.push_back ( generate_sqlite_paired_command ( chr , rc , p-l , "" , 0 , '1' ) ) ;
		} else {
			uint chrl = (*chrs)[chr].sequence.length() ;
			if ( p + l > chrl ) continue ;
			if ( !matchIUPAC ( _seq+index->index_length_double , chrseq+p+index->index_length_double ) ) continue ;
			if ( pileup ) (*chrs)[chr].ao.add_align ( _seq , _quality , p , chr ) ;
			if ( cigar || sqlite ) add_cigar_sqlite ( _seq , last_solexa_name , p , chr , '+' ) ;
			if ( coverage ) (*chrs)[chr].add_coverage ( _seq , p , false ) ;
			if ( snpsinreads ) add_snpsinreads ( _seq , _quality , chr , p , 1 ) ;
			if ( sqlite ) sqlite_out.push_back ( generate_sqlite_paired_command ( chr , _seq , p , "" , 0 , '1' ) ) ;
		}
		ret++ ;
		if ( ret >= max_align && !sqlite ) return ret ; // Shortcut for binning-only
	}
	
	if ( wobble && ret == 0 ) wobble_single_read ( _seq , _quality ) ;

	if ( sqlite && sqlite_out.size() == 1 ) add2sqlite_cache ( sqlite_out[0] ) ;

	return ret ;
}

string TChromosomeAlign::matchIUPAC_tolerant ( const char *c1 , const char *c2 , int tolerance ) {
	string ret ;
	while ( *c1 && *c2 ) {
		if ( char2IUPAC[*c1++] & char2IUPAC[*c2++] ) {
			ret += *(c2-1) ;
			continue ;
		}
		if ( tolerance-- == 0 ) return "" ;
		ret += *(c2-1) ;
	}
	return ret ;
}



void TChromosomeAlign::wobble_single_read ( char *seq , const string &_quality ) {
	uint i1 = index->calculate_index ( seq ) ;
	uint l = strlen ( seq ) ;
	string rc , rcq , query , realseq ;
	int start , dir ;

	// Binary search
	TPositionWithIndex findme ;
	findme.index2 = index->calculate_index2 ( seq + index->index_length ) ;
	pair <TPWIi , TPWIi> pp = equal_range ( index->position[i1].begin() , index->position[i1].end() , findme ) ;

	for ( TPWIi it = pp.first ; it < pp.second ; it++ ) {
		uint chr = it->chromosome ;
		char *chrseq = (char*) (*chrs)[chr].sequence.c_str() ;
		int p = it->position ;
		string found ;
		if ( it->reverse_complement ) {
			if ( p < l ) continue ;
			if ( rc.empty() ) { // Creating the reverse-complement
				rc = string ( seq ) ;
				REVERSE_COMPLEMENT_CHARP ( (char*)rc.c_str() ) ;
			}
			query = rc.substr(0,l-index->index_length_double) ;
			start = p-l ;
			found = matchIUPAC_tolerant ( query.c_str() , chrseq+start , wobblemax ) ;
			if ( found.empty() ) continue ;
			realseq = dorc(seq).substr(0,l-index->index_length_double) ;
		} else {
			uint chrl = (*chrs)[chr].sequence.length() ;
			if ( p + l > chrl ) continue ;
			query = seq+index->index_length_double ;
			start = p+index->index_length_double ;
			found = matchIUPAC_tolerant ( query.c_str() , chrseq+start , wobblemax ) ;
			if ( found.empty() ) continue ;
			realseq = query ;
		}
		dir = 1 ;
		for ( int a = 0 ; a < query.length() && a < found.length() ; a++ ) {
			char c = (*chrs)[chr].original_sequence.empty() ? (*chrs)[chr].sequence[start + a * dir] : (*chrs)[chr].original_sequence[start + a * dir] ; 
//			cout << start + a * dir + 1 << "\t" << c << "\t" << query[a] << "\t" << found[a] << endl ;
			if ( char2IUPAC[c] & char2IUPAC[query[a]] ) continue ;
//			if ( c == found[a] ) continue ;
			fprintf ( wobble , "SNP\t%s\t%d\t%c\t%c\n" , (*chrs)[chr].name.c_str() , start + a * dir + 1 , c , query[a] ) ;
		}
//		cout << query << " O" << endl << found << " N" << endl << endl ;
	}
}

void TChromosomeAlign::add_cigar_sqlite ( const string &sequence , const string &name , uint position , uint chromosome , char orientation ) {
	char cig[READ_CHAR_BUFFER] ;
	string cigar_text ;
	
	int a ;
	bool iu = false ;
	for ( a = 0 ; a < sequence.length() && !iu ; a++ ) {
		if ( isIUPAC[(*chrs)[chromosome].sequence[a+position]] ) iu = true ;
	}
	if ( iu ) {
		string s ;
		for ( a = 0 ; a < sequence.length() ; a++ ) {
			if ( (*chrs)[chromosome].sequence[a+position] != sequence[a] ) s += "D" ;
			else s += "M" ;
		}
		
		// Compressing CIGAR
		char last = ' ' ;
		uint cnt = 0 ;
		for ( a = 0 ; a < s.length() ; a++ ) {
			if ( s[a] != last ) {
				if ( last != ' ' ) {
					if ( !cigar_text.empty() ) cigar_text += " " ;
					sprintf ( cig , "%c %d" , last , cnt ) ;
					cigar_text += cig ;
					if ( last == 'D' ) {
						sprintf ( cig , " I %d" , cnt ) ;
						cigar_text += cig ;
					}
				}
				last = s[a] ;
				cnt = 1 ;
			} else cnt++ ;
		}
		sprintf ( cig , " %c %d" , last , cnt ) ;
		cigar_text += cig ;
		if ( last == 'D' ) {
			sprintf ( cig , " I %d" , cnt ) ;
			cigar_text += cig ;
		}
		
//		cout << cigar_text << endl ;
	} else {
		sprintf ( cig , "M %d" , (int)sequence.length() ) ;
		cigar_text = cig ;
	}

	if ( cigar ) {
		fprintf ( cigar , "%s %d %d %c %s %d %d + %d %s\n" ,
								name.c_str() , // Read name
								1 , // Read begin
								sequence.length() , // Read end
								orientation , // Read orientation
								(*chrs)[chromosome].name.c_str(), // Chromosome name
								position + 1 , // Begin on chromosome
								position + sequence.length() , // End on chromosome
								99 , // Fake Smith-Waterman score
								cigar_text.c_str() // CIGAR commands
								) ;
	}
}

void TChromosomeAlign::show_pileup ( bool snps_only ) {
	if ( !pileup ) return ;
	for ( int chromosome = 0 ; chromosome < chrs->size() ; chromosome++ ) (*chrs)[chromosome].ao.show_pileup ( pileup , snps_only ) ;
//	cerr << "SNPS searched : " << search_snps << ", found : " << found_snps << " (" << (found_snps*100/search_snps) << "%)" << endl ;
//	if ( bogus ) cerr << "Bogus : " << bogus << endl ;
}



void TChromosomeAlign::align_solexa_fastq_variety ( string filename ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Reading solexa data from %s ... " , filename.c_str() ) ; fflush(stdout); }

	FILE * file = fopen ( filename.c_str() , "r" ) ;
	if ( !file ) {
		cerr << "Could not open FASTQ file \"" << filename << "\" - aborting." << endl ;
		exit ( 1 ) ;
	}
	char dummy[READ_CHAR_BUFFER] , *c1 , *c2 ;
	int solexacount = 0 , total_matches = 0 , matched_reads = 0 ;

	string seq , q ;
	char *c ;
	bool invalid ;
	fasta = false ;

	if ( !pileup ) max_align = 2 ; // If we don't do the pileup but just bin, 2 matches are all we need (=> multiple_bin)

	// Run through file
	while ( !feof(file) ) {
		invalid = false ;
		*dummy = 0 ;
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		if ( !*dummy ) break ;
		// Read name
		if ( *dummy != '@' ) { cerr << "Expected '@'\n" << dummy ; exit ( 1 ) ; }
		for ( c = dummy ; *c > 13 ; c++ ) ;
		*c = 0 ;
		last_solexa_name = dummy + 1 ;

		// Sequence
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		for ( c = dummy ; !isIUPAC[*c] ; c++ ) ;
		if ( *c == 10 || *c == 13 ) *c = 0 ;
		if ( *c ) {
			invalid = true ;
			while ( *c > 13 ) c++ ;
			*c = 0 ;
		}
		seq = dummy ;

		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		if ( *dummy != '+' ) { cerr << "Expected '+'\n" ; exit ( 1 ) ; }
		for ( c = dummy ; *c > 13 ; c++ ) ;
		*c = 0 ;
		string qname = dummy + 1 ;
		if ( qname.empty() ) qname = last_solexa_name ;
		if ( last_solexa_name != qname ) { cerr << "Name mismatch\n" << last_solexa_name << " : " << qname ; exit ( 1 ) ; }
		
		fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		for ( c = dummy ; *c > 13 ; *c++ ) ;
		*c = 0 ;
		q = dummy ;
		
		try_matching_read ( seq , q ) ;
		
		int m = invalid ? -1 : try_matching_read ( seq , q ) ;
		if ( m > 0 ) matched_reads++ ;
		if ( m > 0 ) total_matches += m ;
		if ( using_bins ) { // Write to bins
			if ( last_solexa_name.empty() ) {} // Paranoia
			else write2bin ( last_solexa_name , seq , q , m ) ;
			last_solexa_name.clear() ;
		}
		if ( !invalid ) solexacount++ ;
	}
	fclose ( file ) ;
	
	if ( DEBUGGING ) { fprintf ( stdout , "scanned %d solexa reads, %d (%2.2f%%) matched , %d total matches.\n" , solexacount , matched_reads , (float) matched_reads * 100 / solexacount , total_matches ) ; fflush(stdout); }
}

int TChromosomeAlign::try_matching_read ( const string & seq , string quality ) {
	vector <uint> indices ;
	uint idd = index->index_length_double ;
	
	// Beginning
	indices.push_back ( index->calculate_index ( seq.c_str() ) ) ;
	indices.push_back ( index->calculate_index2 ( seq.c_str() + index->index_length ) ) ;
	
	// Middle
//	indices.push_back ( index->calculate_index ( seq.c_str() + seq.length() / 2 + idd / 2 ) ) ;
//	indices.push_back ( index->calculate_index2 ( seq.c_str() + seq.length() / 2 + idd / 2 + index->index_length ) ) ;
	
	// End
	indices.push_back ( index->calculate_index ( seq.c_str() + seq.length() - 1 - idd ) ) ;
	indices.push_back ( index->calculate_index2 ( seq.c_str() + seq.length() - 1 - idd + index->index_length ) ) ;
	
	uint m , n ;
	pwi_vector out ;
	
	for ( n = 0 ; n < indices.size() ; n += 2 ) {
		uint i1 = indices[n] ;
		TPositionWithIndex findme ;
		findme.index2 = indices[n+1] ;
		pair <TPWIi , TPWIi> pp = equal_range ( index->position[i1].begin() , index->position[i1].end() , findme ) ;

		for ( TPWIi it = pp.first ; it < pp.second ; it++ ) {
			out.push_back ( *it ) ;
		}
	}
	
	if ( out.size() < 2 ) return 0 ;
	
	sort ( out.begin() , out.end() , tpi_full_compare ) ;
	
	
	int cnt = 0 , start = -1 ;
	for ( n = 0 ; n < out.size() ; n++ ) {
		bool found = false ;
		for ( m = n + 1 ; m < out.size() ; m++ ) {
			if ( out[m].chromosome != out[n].chromosome ) break ;
			if ( out[m].reverse_complement != out[n].reverse_complement ) break ;
			
			int diff = out[m].position - out[n].position ;
			if ( diff < 0 ) diff = -diff ;
			
			if ( diff > seq.length() ) continue ;
			found = true ;
			start = out[m].position < out[n].position ? m : n ;
			break ;
		}
		if ( !found ) continue ;
		cnt++ ;
//		cout << out[n].chromosome << "\t" << out[n].position << "\t" << out[n].reverse_complement << endl ;
	}
//	cout << endl ;

	if ( cnt == 1 ) {
		string ref = get_ref_substr ( out[start] , seq.length() ) ;
		string s = out[start].reverse_complement ? dorc ( seq ) : seq ; // ????
		for ( int a = 0 ; a < s.length() ; a++ ) if ( s[a] != ref[a] ) s[a] = s[a] - 'A' + 'a' ;
		cout << out[start].chromosome << "\t" << out[start].position << "\t" << ref << "\t" << s << "\t" << out.size() << "\t" << cnt << "\t" << out[start].reverse_complement << endl ;
	}

	return cnt ;
}

string TChromosomeAlign::get_ref_substr ( const TPositionWithIndex &p , int length ) {
	int start = p.position ;
	if ( p.reverse_complement ) start -= index->index_length_double + 1 ;
	if ( start < 0 ) return "" ;
	string ret = (*chrs)[p.chromosome].sequence.substr ( start , length ) ;
//	if ( p.reverse_complement ) return doc ( ret ) ;
	return ret ;
}
