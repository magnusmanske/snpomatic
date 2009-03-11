/*! \mainpage SNP-O-MATIC
 *
 * \section authors Authors
 * 
 *  - Magnus Manske (mm6@sanger.ac.uk)
 * 
 * \section intro_sec Introduction
 *
 * SNP-O-MATIC can align Solexa reads against a reference sequence.
 * Known SNPs can be added to the reference sequence as IUPAC bases.
 * Only perfectly matching Solexa reads, with respect to the IUPAC bases, will be aligned.
 * Reads can be binned into those that do not match the reference, those that match uniquely, and those that match two or more times.
 *
 * \section install_sec Installation
 * 
 * SNP-O-MATIC comes with a Makefile. The Intel compiler is used by default for performance reasons, but the g++ variant is in the Makefile as well in the form of comments.
 */

#include "snpomatic.h"
#include <algorithm>
#include <map>
#include <fstream>

#define MAX_NEW 1

TChromosomalIndices *idx ;

class TRead ;

class TPart {
	public :
	void do_index ( uint id , bool rc = false ) ;
	string seq ;
} ;

class TRead {
	public :
	TRead () {} ;
	TRead ( const string &seq ) ;
	TRead ( const string &seq1 , const string &seq2 ) ;
	inline void reset ( uint id  ) ;
	bool used ;
	vector <uint> related ;
	TPart p1 , p2 ;
} ;

TRead::TRead ( const string &seq ) {
	int l = seq.length() / 2 ;
	p1.seq = seq.substr ( 0 , l ) ;
	p2.seq = seq.substr ( l , l ) ;
}

TRead::TRead ( const string &seq1 , const string &seq2 ) {
	p1.seq = seq1 ;
	p2.seq = seq2 ;
}

void TRead::reset ( uint id  ) {
	p1.do_index ( id , false ) ;
	p2.do_index ( id , true ) ;
	used = false ;
	related.clear() ;
}

void TPart::do_index ( uint id , bool rc ) {
	for ( int a = 0 ; a + idx->index_length_double < seq.length() ; a++ ) {
		uint i1 = idx->calculate_index ( (char*) seq.c_str() + a ) ;
		TPositionWithIndex p ;
		p.index2 = idx->calculate_index2 ( (char*) seq.c_str() + a + idx->index_length ) ;
		p.position = a ;
		p.reverse_complement = rc ;
		p.chromosome = id ;
		idx->position[i1].push_back ( p ) ;
		break ;
	}
}



string pad ( int a ) {
	string ret ;
	while ( ret.length() < a ) ret += " " ;
	return ret ;
}



/*
int remove_full_sequences ( vector <TRead> &vr , vector <string> &vs ) {
	int min = 30 ;
	int cnt = 0 ;
	for ( int a = 0 ; a < vr.size() ; a++ ) {
		if ( vr[a].p1.seq.length() < min ) continue ;
		string m = vr[a].p1.seq.substr ( vr[a].p1.seq.length() - min , min ) ;
		int n = vr[a].p2.seq.find ( m , 0 ) ;
		if ( n == std::string::npos ) continue ;
		m = vr[a].p1.seq + vr[a].p2.seq.substr ( n ) ;
		vs.push_back ( m ) ;
		cnt++ ;
		vr[a].p1.seq = "" ;
		vr[a].p2.seq = "" ;
	}
	return cnt ;
}

void create_contigs ( vector <TRead> &vr , vector <TRead> &vn ) {
	for ( uint r1 = 0 ; r1 < vr.size() ; r1++ ) {
		if ( vr[r1].used ) continue ;
		for ( uint r2x = 0 ; r2x < vr[r1].related.size() ; r2x++ ) {
			uint r2 = vr[r1].related[r2x] ;
			if ( vr[r1].p1.seq == vr[r2].p1.seq ) continue ;
			if ( vr[r1].p2.seq == vr[r2].p2.seq ) continue ;
			string t1 = combine_strings ( vr[r1].p1.seq , vr[r2].p1.seq ) ;
			if ( t1.empty() ) continue ;
			string t2 = combine_strings ( vr[r1].p2.seq , vr[r2].p2.seq ) ;
			if ( t2.empty() ) continue ;
			vr[r1].used = true ;
			vr[r2].used = true ;
			vn.push_back ( TRead ( t1 , t2 ) ) ;
//			cout << t1 << "\t" << t2 << endl ;
			break ;
		}
	}
}

void join_single_sequences ( vector <string> &vs ) {
	vector <string> vs2 ;
	vs.swap ( vs2 ) ;
	do {
		for ( int a = 0 ; a < vs.size() ; a++ ) {
			if ( vs[a].empty() ) continue ;
			vs2.push_back ( vs[a] ) ;
		}
		vs.swap ( vs2 ) ;
		vs2.clear() ;
		for ( int a = 0 ; a < vs.size() ; a++ ) {
			if ( vs[a].length() < 100 ) continue ;
			for ( int b = a+1 ; b < vs.size() ; b++ ) {
				if ( vs[b].empty() ) continue ;
				string m = combine_strings ( vs[a] , vs[b] , 100 ) ;
				if ( m.empty() ) continue ;
				vs[a].clear() ;
				vs[b].clear() ;
				vs2.push_back ( m ) ;
			}
		}
	} while ( vs2.size() > 0 ) ;
}

void interlink_reads ( vector <TRead> &vr ) {
	uint p1s = idx->position.size() ;
	for ( uint i1 = 0 ; i1 < p1s ; i1++ ) {
		uint last = 0 ;
		uint p2s = idx->position[i1].size() ;
		for ( uint a = 1 ; a < p2s ; a++ ) {
			if ( idx->position[i1][last].index2 == idx->position[i1][a].index2 ) {
				uint n2 = idx->position[i1][a].chromosome ;
				for ( uint b = last ; b < a ; b++ ) {
					uint n1 = idx->position[i1][b].chromosome ;
					vr[n2].related.push_back ( n1 ) ;
					vr[n1].related.push_back ( n2 ) ;
				}
			} else {
				last = a ;
			}
		}
	}
	
	uint cnt = 0 ;
	for ( uint a = 0 ; a < vr.size() ; a++ ) {
		sort ( vr[a].related.begin() , vr[a].related.end() ) ;
		vector <uint> vi ;
		vi.reserve ( vr[a].related.size() ) ;
		for ( uint b = 0 ; b < vr[a].related.size() ; b++ ) {
			if ( b > 0 && vr[a].related[b-1] == vr[a].related[b] ) continue ;
			vi.push_back ( vr[a].related[b] ) ;
		}
		cnt += vi.size() ;
		vr[a].related.swap ( vi ) ;
	}
	cerr << cnt << " crosslinks found\n" ;
}


void align_perfect ( vector <TRead> &reads , int offset ) {
	int a , b ;
	if ( DEBUGGING ) { fprintf ( stdout , "Indexing ... " ) ; fflush(stdout); }
	for ( a = 0 ; a < reads.size() ; a++ ) reads[a].reset ( a ) ;
	if ( DEBUGGING ) { fprintf ( stdout , "done.\n" ) ; fflush(stdout); }
	
	if ( DEBUGGING ) { fprintf ( stdout , "Scanning with offset %d ... " , offset ) ; fflush(stdout); }
	
	for ( uint i1 = 0 ; i1 < idx->position.size() ; i1++ ) {
		uint i1s = idx->position[i1].size() ;
		for ( uint a = 0 ; a < i1s ; a++ ) {
			if ( idx->position[i1][a].position > offset ) continue ;
			uint i2 = idx->position[i1][a].index2 ;
			for ( uint b = a+1 ; b < i1s && idx->position[i1][b].index2 == i2 ; b++ ) {
				if ( idx->position[i1][b].position > offset ) continue ;
				uint r1 = idx->position[i1][a].chromosome ;
				uint r2 = idx->position[i1][b].chromosome ;
				if ( reads[r2].p1.seq != reads[r1].p1.seq && reads[r2].p1.seq != reads[r1].p2.seq && 
					reads[r2].p2.seq != reads[r1].p1.seq && reads[r2].p2.seq != reads[r1].p2.seq ) continue ;
			}
		}
	}
	
	if ( DEBUGGING ) { fprintf ( stdout , "done.\n" ) ; fflush(stdout); }
}

string combine_strings ( const string &s1 , const string &s2 , int max = -1 ) {
	if ( max == -1 ) max = idx->index_length_double ;
	string ret ;
	int s1l = s1.length() ;
	int s2l = s2.length() ;
	int maxa = s2l - max ;
	int a , b ;
	for (  a = - ( s1l - max ) ; a < maxa ; a++ ) {
		int maxb = s2l < s1l-a ? s2l : s1l-a ;
		for ( b = a<0?-a:0 ; b < maxb && s1[a+b] == s2[b] ; b++ ) ;
		if ( b < maxb ) continue ;
		
		if ( !ret.empty() ) return string () ; // Do not allow multiple hits
		
		if ( a < 0 ) {
			ret = s2.substr ( 0 , -a ) + s1 ;
			if ( s1l-a < s2l ) ret += s2.substr ( s1l-a ) ;
		} else {
			ret = s1.substr ( 0 , a ) + s2 ;
			if ( s2l+a < s1l ) ret += s1.substr ( s2l+a ) ;
		}
	}
	return ret ;
}

*/

string merge_strings ( string s1 , string s2 ) {
	string ret ;
	if ( s1.length() > s2.length() ) {
		string tmp = s1 ;
		s1 = s2 ;
		s2 = tmp ;
	}
	int sl1 = s1.length() ;
	int sl2 = s2.length() ;
	
	int best = sl1 + sl2 ;
	int bestpos = 0 ;
	
	int a , b ;
	int min_offset = 5 ;
	for ( a = -sl1 + min_offset ; a < sl2 - min_offset ; a++ ) {
		int match = 0 ;
		int mismatch = 0 ;
		for ( b = 0 ; b < sl1 ; b++ ) {
			int p2 = a + b ;
			if ( p2 < 0 || p2 >= sl2 ) continue ;
			if ( s1[b] == 'N' || s2[p2] == 'N' ) continue ;
			if ( s1[b] == s2[p2] ) match++ ;
			else mismatch++ ;
		}
		if ( mismatch > match / 3 ) continue ;
		int total = match + mismatch ;
		if ( mismatch < best ) {
			best = mismatch ;
			bestpos = a ;
		}
	}
	
	if ( best > 5 ) return ret ; // Hard cutoff
	
	ret = s2 + pad ( sl1 ) ;
	if ( bestpos < 0 ) {
		ret = pad ( -bestpos ) + ret ;
		bestpos = 0 ;
	}
	
	// Merge
	for ( a = 0 ; a < sl1 ; a++ ) {
		b = bestpos + a ;
		if ( s1[a] == ret[b] ) continue ;
		if ( ret[b] != 'N' && s1[a] != ' ' ) {
			if ( ret[b] == ' ' ) ret[b] = s1[a] ;
			else ret[b] = 'N' ;
		}
	}
	
	// Trim blanks
	for ( a = 0 ; ret[a] == ' ' ; a++ ) ;
	if ( a > 0 ) ret = ret.substr ( a , ret.length()-a ) ;
	for ( a = 0 ; a < ret.length() && ret[a] != ' ' ; a++ ) ;
	ret = ret.substr ( 0 , a ) ;
	
	return ret ;
}


void read_fastq ( string filename , vector <TRead> &vr , int min_qual = 0 ) {
	if ( DEBUGGING ) { fprintf ( stdout , "Reading solexa pair data from %s ... " , filename.c_str() ) ; fflush(stdout); }

	FILE * file = fopen ( filename.c_str() , "r" ) ;
	if ( !file ) {
		cerr << "Could not open FASTQ file \"" << filename << "\" - aborting." << endl ;
		exit ( 1 ) ;
	}
	char dummy[READ_CHAR_BUFFER] , *c1 , *c2 ;
	string lastseq , lastqual , last_solexa_name ;
	bool fasta = true ;
	bool last_was_quality = false , lastone = false , eof = false , skip = false , goodqual = true ;
	uint readcount = 0 , matched_reads = 0 , total_matches = 0 , skipped = 0 , badqualcnt = 0 ;

	while ( !eof || lastone ) {
		if ( lastone ) {
			strcpy ( dummy , "@dummy" ) ;
			if ( fasta ) *dummy = '>' ;
		} else {
			*dummy = 0 ;
			fgets ( dummy , READ_CHAR_BUFFER , file ) ;
		}
		
		// Remove EOL
		for ( c1 = dummy ; *c1 > 13 ; c1++ ) ;
		*c1 = 0 ;
		
		if ( ( fasta && *dummy == '>' ) || *dummy == '@' ) {
			if ( *dummy == '@' ) fasta = false ;
			
			if ( !skip ) {
				for ( c1 = (char*) lastseq.c_str() ; *c1 && !isIUPAC[*c1] ; c1++ ) ;
				if ( *c1 ) { // Contains IUPAC
				} else if ( !goodqual ) {
					badqualcnt++ ;
				} else if ( !lastseq.empty() ) {
					vr.push_back ( TRead ( lastseq ) ) ;
				}
			}

			goodqual = true ;
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
			
		} else if ( *dummy == '+' ) {
			last_was_quality = true ;
		} else if ( last_was_quality ) {
			if ( min_qual > 0 ) {
				goodqual = true ;
				for ( c1 = dummy ; goodqual && *c1 > 20 ; c1++ ) {
//					cout << (int) *c1-33 << ", " ;
					goodqual = ( *c1-33 ) >= min_qual ;
				}
//				cout << endl ;
			}
			last_was_quality = false ;
//			lastqual += dummy ;
		} else {
			lastseq += dummy ;
		}
		
		eof = feof ( file ) ;
//		if ( vr.size() >= 100000 ) eof = 1 ; // TESTING!!!!
		if ( eof ) lastone = true ; // Fake last read
	}
	
	if ( DEBUGGING ) { fprintf ( stdout , "scanned %d solexa reads, %d (%2.2f%%) matched , %d total matches, %d skipped, %d bad quality.\n" , 
		readcount , matched_reads , (float) matched_reads * 100 / readcount , total_matches , skipped , badqualcnt ) ; fflush(stdout); }
}


// __________________________________________________________________________________________________________________________________________________________________ MAIN

/// Shows possible parameters, and an (optional) error message.
void die_with_description ( const char *msg = NULL ) {
	if ( msg ) cerr << msg << endl ;
	cerr << "Parameters : \n" ;
	cerr << "--genome=GENOME_FILE     FASTA file with chromosomes (mandatory)\n" ;
	cerr << "--fasta=FASTA_FILE       FASTA file with Solexa reads (mandatory, except when --fastq or --index is used)\n" ;
	cerr << "--fastq=FASTQ_FILE       FASTQ file with Solexa reads (mandatory, except when --fasta or --index is used)\n" ;
	cerr << "--nono=FILENAME          File with list of read names (!) to ignore (optional)\n" ;
	cerr << "--regions=REGION_FILE    Region file for finding new SNPs (optional) [DEPRECATED]\n" ;
	cerr << "--snps=SNP_FILE          Simple SNP file (optional)\n" ;
	cerr << "--gff=GFF_FILE           GFF file with SNPs (optional)\n" ;
	cerr << "--uniqueness=FILE        Output a uniqueness data file for the reference; no Solexa reads needed; implies --noshortcuts (optional)\n" ;
	cerr << "--pileup=FILE            Outputs complete pileup into that file (optional)\n" ;
	cerr << "--cigar=FILE             Outputs alignment in CIGAR format (optional)\n" ;
	cerr << "--gffout=FILE            Outputs reads in GFF format (optional)\n" ;
	cerr << "--coverage=FILENAME      Outputs (high depth) coverage data (optional)\n" ;
	cerr << "--wobble=FILENAME        Outputs a list of possible variations (optional; paired reads only) [UNDER CONSTRUCTION]\n" ;
	cerr << "--fragmentplot=FILENAME  Outputs a plot of fragment size distribution to a file (optional)\n" ;
	cerr << "--snpsinreads=FILENAME   Outputs a list of reads containing known SNPs to a file (optional)\n" ;
	cerr << "--indelplot=FILENAME     Outputs indel data to a file (optional)\n" ;
	cerr << "--inversions=FILENAME    For paired reads, writes read matches indicating  inversions into a file (optional)\n" ;
	cerr << "--bins=FILE_PREFIX       Outputs no-match, single-match and multi-match Solexa reads into prefixed files (optional)\n" ;
	cerr << "--binmask=MASK           Mask of 1s and 0s to turn off individual bins. Order: No match, single match, multi-match, IUPAC.\n" ;
	cout << "                         Example: 0100 creates only single-match bin. (optional; default:1111)\n" ;
	cerr << "--pair=NUMBER            For paired reads, the length of the first part of the read (mandatory for paired reads)\n" ;
	cerr << "--fragment=NUMBER        For paired reads, the average fragment length (mandatory for paired reads)\n" ;
	cerr << "--variance=NUMBER        For paired reads, the variance of the fragment length to either side (optional; default: 1/4 of fragment size)\n" ;
	cerr << "--wobblemax=NUMBER       Maximum number of mismatches for wobble (optional; default 2; see --wobble)\n" ;
	cerr << "--mspi=NUMBER            Maximum number of SNPs per chromosomal index (optional; default:8)\n" ;
	cerr << "--index=FILENAME         Index filename (index will be created if it does not exist; optional)\n" ;
	cerr << "--noshortcuts            Will process all chrososomal regions, even those with lots'o'repeats (optional; no value)\n" ;
	cerr << "--snpsonly               Only lists found SNPs in the pileup (optional; no value)\n" ;
	cerr << "--chromosome=NAME        Discards all chromosomes but NAME prior to run (optional)\n" ;
	cerr << "--chop=NUMBER            For paired reads, if one but not the other matches, shorten the other by NUMBER bases (optional)\n" ;
	cerr << "--index1=NUMBER          Length of internal index 1 (optional; default:10)\n" ;
	cerr << "--index2=NUMBER          Length of internal index 2 (optional; default:16)\n" ;
	cerr << "--multimatch             Puts a multiple-matching read to a random position (optional) [currently paired reads only]\n" ;
	cerr << "--foum                   For paired reads, at least one read has to match uniquely in the genome (force one unique match) (optional)\n" ;
//	cerr << "--featuretable=OUTPUT_FILE\tOutputs a DDBJ/EMBL/GenBank feature table for matched reads (optional; paired reads only for now)\n" ;
	exit ( 0 ) ;
}


class TContig {
	public :
	TContig () { remove = false ; }
	string sequence ;
//	vector <unsigned int> reads ;
	bool remove ;
} ;

vector <TContig> contigs ;
map <string,bool> cache ;
vector <TPositionWithIndex> contig_idx ;

void initialize_contigs ( vector <TRead> & vr , string start ) {
	unsigned int a , b ;
	for ( a = 0 ; a < vr.size() ; a++ ) {
//		char *cc = (char*) vr[a].p1.seq.c_str() ;
//		if ( *(cc) != 'A'  || *(cc+1) != 'T' || *(cc+2) != 'G' ) continue ;

	
		if ( vr[a].p1.seq.substr(0,3) == start && cache.end() == cache.find ( vr[a].p1.seq ) ) {
			cache[vr[a].p1.seq] = true ;
			TContig c ;
			c.sequence = vr[a].p1.seq ;
			contigs.push_back ( c ) ;
			vr[a].used = true ;
		}

		if ( vr[a].p2.seq.substr(0,3) == start && cache.end() == cache.find ( vr[a].p2.seq ) ) {
			cache[vr[a].p2.seq] = true ;
			TContig c ;
			c.sequence = vr[a].p2.seq ;
			contigs.push_back ( c ) ;
			vr[a].used = true ;
		}

	}
	cout << contigs.size() << " contigs created (from " << vr.size() << " reads).\n" ;
}

void add_contig ( string seq , int readnum ) {
	if ( cache.end() == cache.find ( seq ) ) {
		TContig c ;
		c.sequence = seq ;
		contigs.push_back ( c ) ;
		cache[seq] = 1 ;
	} else {
	}
}


int grow_contigs ( vector <TRead> & vr , int pos , string start ) {
	cout << "Growing at position " << pos << " ... " ;
	int read_length = vr[0].p1.seq.length() ;
	int ret = 0 ;
	
	contig_idx.clear() ;
	for ( unsigned int a = 0 ; a < contigs.size() ; a++ ) {
		TPositionWithIndex p ;
		p.index2 = idx->calculate_index2 ( (char*) contigs[a].sequence.c_str() ) ;
		p.chromosome = a ;
		contig_idx.push_back ( p ) ;
	}
	sort ( contig_idx.begin() , contig_idx.end() ) ;
	
	for ( unsigned int a = 0 ; a < contigs.size() ; a++ ) {
		if ( contigs[a].remove ) continue ; // Marked for deletion
		string cseq = contigs[a].sequence ;
		if ( read_length + pos > cseq.length() + MAX_NEW ) continue ;
		
		int len = cseq.length() - pos ;
		if ( len > read_length ) len = read_length ;
		string csub = cseq.substr ( pos , len ) ;

		contigs[a].remove = true ;
		bool altered = false ;
		
		if ( csub.substr ( 0 , 3 ) == start ) { // Elongate with another "contig"
			csub = cseq.substr ( pos , cseq.length()-pos ) ;
			unsigned int i2 = idx->calculate_index2 ( (char*) csub.c_str() ) ;
			unsigned int b ;
			for ( unsigned int b2 = 0 ; b2 < contig_idx.size() ; b2++ ) {
				if ( contig_idx[b2].index2 > i2 ) break ;
				if ( contig_idx[b2].index2 < i2 ) continue ;
				b = contig_idx[b2].chromosome ;
				if ( a == b ) continue ; // No self-elongation
				// if ( contigs[b].remove ) continue ;
				if ( contigs[b].sequence.length() <= csub.length() ) continue ; // Contig smaller than fragment
				if ( contigs[b].sequence.substr ( 0 , csub.length() ) == csub ) {
					string s = cseq.substr ( 0 , pos ) + contigs[b].sequence ;
					add_contig ( s , -1 ) ;
					altered = true ;
					contigs[b].remove = true ;
					break ;
				}
			}
		} else { // Try to elongate with reads
			uint i1 = idx->calculate_index ( (char*) csub.c_str() ) ;
			uint i2 = idx->calculate_index2 ( (char*) csub.c_str() + idx->index_length ) ;
			
			
			for ( unsigned d = 0 ; d < idx->position[i1].size() ; d++ ) {
				if ( idx->position[i1][d].index2 > i2 ) break ;
				if ( idx->position[i1][d].index2 < i2 ) continue ;
				unsigned int b = idx->position[i1][d].chromosome ;
				if ( vr[b].used ) continue ;
				
				string *sp = idx->position[i1][d].reverse_complement ? &(vr[b].p2.seq) : &(vr[b].p1.seq) ;
				
				if ( sp->substr ( 0 , len ) != csub ) continue ;
				string s = cseq.substr ( 0 , pos ) + *sp ;
				if ( s.length() <= cseq.length() ) continue ; // No new length increase
				add_contig ( s , b ) ;
				vr[b].used = true ;
				altered = true ;
			}
		}

		if ( altered ) ret++ ;
		else contigs[a].remove = false ;
	}
	
	// Cleanup
	for ( int a = 0 ; a < contigs.size() ; a++ ) {
		if ( !contigs[a].remove ) continue ;
		cache.erase ( contigs[a].sequence ) ;
		contigs[a] = contigs[contigs.size()-1] ;
		contigs.pop_back() ;
		a-- ;
	}
	
	cout << ret << " contigs grew (" << contigs.size() << " contigs total)." << endl ;
	return ret ;
}

void assemble_contigs ( vector <TRead> &vr , vector <TChromosome> &chrs , string start ) {
	initialize_contigs ( vr , start ) ;
	int a = 1 , growing , growth ;
	do {
		growing = grow_contigs ( vr , a++ , start ) ;
		if ( growing > 0 ) growth = 10 ;
		else growth-- ;
	} while ( growth > 0 ) ;
	
	
	cout << "Writing contigs to contigs.fasta\n" ;
	ofstream fout( "contigs.fasta", ios::out|ios::trunc );
	for ( a = 0 ; a < contigs.size() ; a++ ) {
		if ( contigs[a].sequence.length() == vr[0].p1.seq.length() ) continue ;
		
		string s = contigs[a].sequence ;
		//new_contigs.push_back ( s ) ;
		
		// FIXME : Add to chrs
		
		fout << ">contig_" << a << "_" << s.length() << endl ;
		while ( !s.empty() ) {
			if ( s.length() >= 60 ) {
				fout << s.substr ( 0 , 60 ) << endl ;
				s = s.substr ( 60 ) ;
			} else {
				fout << s << endl ;
				break ;
			}
		}
	}
	fout.close() ;
	
	int used = 0 ;
	for ( a = 0 ; a < vr.size() ; a++ ) {
		if ( vr[a].used ) used++ ;
		vr[a].used = false ; // Resetting
	}
	cout << used << " of " << vr.size() << " reads (" << 100 * used / vr.size() << "%) used in contigs\n" ;
	
	// Cleanup
	contigs.clear() ;
	cache.clear() ;
	contig_idx.clear() ;
}

void write_chromosomes_to_fasta ( vector <TChromosome> &chrs ) {
	cout << "Writing contigs to chrs.fasta\n" ;
	ofstream fout( "chrs.fasta", ios::out|ios::trunc );
	map <string,bool> seen ;
	for ( int a = 0 ; a < chrs.size() ; a++ ) {
		string s = chrs[a].sequence ;
		if ( seen.end() != seen.find ( s ) ) continue ; // Wrote that one already
		seen[s] = true ;
		fout << ">" << chrs[a].name << endl ;
		while ( !s.empty() ) {
			if ( s.length() >= 60 ) {
				fout << s.substr ( 0 , 60 ) << endl ;
				s = s.substr ( 60 ) ;
			} else {
				fout << s << endl ;
				break ;
			}
		}
	}
	fout.close() ;
}

void make_read_indices ( vector <TRead> &vr ) {
	// Create indices
	for ( uint a = 0 ; a < vr.size() ; a++ ) {
		vr[a].reset ( a ) ;
	}
	
	// Sort indices
	for ( uint a = 0 ; a < idx->position.size() ; a++ ) {
		sort ( idx->position[a].begin() , idx->position[a].end() ) ;
	}
}



void join_chromosomes ( vector <TRead> &vr , uint index_length_1 , uint index_length_2 , string genome_file ) {

	for ( int cnt = 0 ; cnt < 100 ; cnt++ ) {
		cout << "Chromosome joining, round " << cnt << endl ;

		vector <TChromosome> chrs ;
		TChromosomalIndices contig_index ( index_length_1 , index_length_2 ) ;
//		contig_index.noshortcuts = noshortcuts ;
		idx = &contig_index ;
		
		load_all_chromosomes ( genome_file , chrs , true ) ;
		contig_index.run ( &chrs ) ;

		TChromosomeAlign ca ( chrs , contig_index ) ;
		ca.init () ;
		for ( uint a = 0 ; a < vr.size() ; a++ ) {
			pwi_vector res1 , res2 ;
			ca.get_align_solexa_read ( vr[a].p1.seq.c_str() , res1 ) ;
			if ( res1.size() != 1 ) continue ;
			ca.get_align_solexa_read ( vr[a].p2.seq.c_str() , res2 ) ;
			if ( res2.size() != 1 ) continue ;
			if ( res1[0].chromosome == res2[0].chromosome ) continue ;
			if ( res1[0].reverse_complement == res2[0].reverse_complement ) continue ;
			
			// Abusing uniqueness vector to store connections
			chrs[res1[0].chromosome].uniqueness.push_back ( res2[0].chromosome ) ;
			chrs[res2[0].chromosome].uniqueness.push_back ( res1[0].chromosome ) ;
		}
		
		for ( uint a = 0 ; a < chrs.size() ; a++ ) {
			sort ( chrs[a].uniqueness.begin() , chrs[a].uniqueness.end() ) ;
			chrs[a].uniqueness.push_back ( chrs.size()+1 ) ; // Dummy, will be ignored
			uint last = chrs.size()+1 ;
			int cnt = -1 ;
			for ( uint b = 0 ; b < chrs[a].uniqueness.size() ; b++ ) {
				if ( chrs[a].uniqueness[b] == last ) {
					cnt++ ;
				} else {
					if ( cnt >= 0 ) {
						if ( cnt > 5 ) {
							string joined = merge_strings ( chrs[a].sequence , chrs[last].sequence ) ;
							if ( !joined.empty() ) {
								chrs[a].sequence = joined ;
								chrs[last].sequence = joined ;
							}
//						cout << a << " => " << last << " (" << cnt << "x)\n" ;
//						cout << chrs[a].sequence << endl << chrs[last].sequence << endl ;
//						cout << merge_strings ( chrs[a].sequence , chrs[last].sequence ) << endl ;
						}
					}
					cnt = 1 ;
					last = chrs[a].uniqueness[b] ;
				}
			}
		}
		
		write_chromosomes_to_fasta ( chrs ) ;
	}
}




/// Well, guess what.
int main ( int argc , char **argv ) {
	string genome_file , gff_file , simple_snp_file , fastq_file , fasta_file , pileup_file , bin_prefix , regions_file ;
	string cigar_file , wobble_file , nono , fragmentplot , featuretable , gffout , coverage_file , uniqueness , snpsinreads ;
	string index_file , indelplot_file , inversions_file , chromosome ;
	string binmask ( "1111" ) ;
	uint wobblemax = 2 ;
	uint mspi = 8 ;
	uint pair_length = 0 ;
	uint fragment_length = 0 ;
	uint index_length_1 = 10 ;
	uint index_length_2 = 16 ;
	int noshortcuts = false ;
	int snpsonly = false ;
	int chop = 0 ;
	uint variance = 0 ;
	int multimatch = false ;
	int force_one_unique_match = false ;

	// Parameters
	static struct option long_options[] = {
		{ "genome" , required_argument , 0 , 'g' } ,
		{ "fasta" , optional_argument , 0 , 'a' } ,
		{ "fastq" , optional_argument , 0 , 'q' } ,
		{ "snps" , optional_argument , 0 , 's' } ,
		{ "coverage" , optional_argument , 0 , 'e' } ,
		{ "gff" , optional_argument , 0 , 'f' } ,
		{ "pileup" , optional_argument , 0 , 'p' } ,
		{ "bins" , optional_argument , 0 , 'b' } ,
		{ "binmask" , optional_argument , 0 , 'z' } ,
		{ "index" , optional_argument , 0 , 'x' } ,
		{ "nono" , optional_argument , 0 , 'n' } ,
		{ "uniqueness" , optional_argument , 0 , '0' } ,
		{ "wobble" , optional_argument , 0 , 'w' } ,
		{ "pair" , optional_argument , 0 , 'i' } ,
		{ "fragment" , optional_argument , 0 , 't' } ,
		{ "featuretable" , optional_argument , 0 , 'u' } ,
		{ "gffout" , optional_argument , 0 , 'y' } ,
		{ "variance" , optional_argument , 0 , 'v' } ,
		{ "mspi" , optional_argument , 0 , 'm' } ,
		{ "cigar" , optional_argument , 0 , 'c' } ,
		{ "fragmentplot" , optional_argument , 0 , 'o' } ,
		{ "snpsinreads" , optional_argument , 0 , 'd' } ,
		{ "regions" , optional_argument , 0 , 'r' } ,
		{ "index1" , optional_argument , 0 , '1' } ,
		{ "index2" , optional_argument , 0 , '2' } ,
		{ "chop" , optional_argument , 0 , '3' } ,
		{ "indelplot" , optional_argument , 0 , '4' } ,
		{ "wobblemax" , optional_argument , 0 , '5' } ,
		{ "inversions" , optional_argument , 0 , '6' } ,
		{ "chromosome" , optional_argument , 0 , '7' } ,
		{ "multimatch" , optional_argument , &multimatch , true } ,
		{ "noshortcuts" , optional_argument , &noshortcuts , true } ,
		{ "foum" , optional_argument , &force_one_unique_match , true } ,
		{ "snpsonly" , optional_argument , &snpsonly , true } ,
		{ 0 , 0 , 0 , 0 }
	} ;
	
	
	int c ;
	int option_index = 0;
	while ( -1 != ( c = getopt_long (argc, argv, "g:a::q::s::f::p::b::r::w::n::",long_options, NULL) ) ) {
		
		switch ( c ) {
			case 0 : break ;
			case 'g' : genome_file = optarg ; break ;
			case 'a' : fasta_file = optarg ; break ;
			case 'q' : fastq_file = optarg ; break ;
			case 's' : simple_snp_file = optarg ; break ;
			case 'f' : gff_file = optarg ; break ;
			case 'e' : coverage_file = optarg ; break ;
			case 'p' : pileup_file = optarg ; break ;
			case 'b' : bin_prefix = optarg ; break ;
			case 'z' : binmask = optarg ; break ;
			case '0' : uniqueness = optarg ; break ;
			case 'w' : wobble_file = optarg ; break ;
			case 'n' : nono = optarg ; break ;
			case 'x' : index_file = optarg ; break ;
			case 'i' : pair_length = atoi ( optarg ) ; break ;
			case 'o' : fragmentplot = optarg ; break ;
			case 'u' : featuretable = optarg ; break ;
			case 'y' : gffout = optarg ; break ;
			case 'd' : snpsinreads = optarg ; break ;
			case 'v' : variance = atoi ( optarg ) ; break ;
			case 't' : fragment_length = atoi ( optarg ) ; break ;
			case 'm' : mspi = atoi ( optarg ) ; break ;
			case '1' : index_length_1 = atoi ( optarg ) ; break ;
			case '2' : index_length_2 = atoi ( optarg ) ; break ;
			case '3' : chop = atoi ( optarg ) ; break ;
			case '4' : indelplot_file = optarg ; break ;
			case '5' : wobblemax = atoi ( optarg ) ; break ;
			case '6' : inversions_file = optarg ; break ;
			case '7' : chromosome = optarg ; break ;
			case 'c' : cigar_file = optarg ; break ;
			case 'r' : regions_file = optarg ; break ;
			default : die_with_description("Unknown option") ;
		}
	}
	binmask += "0000" ; // Ignore bins if shorter string was supplied
	
	init_iupac () ;
	TChromosomalIndices main_index ( index_length_1 , index_length_2 ) ;
	main_index.noshortcuts = noshortcuts ;
	idx = &main_index ;
	
	// FASTA sequence? Just try to join them up and ignore everything smaller than 100 bases.
/*	if ( !fasta_file.empty() ) {
		vector <TChromosome> chrs ;
		load_all_chromosomes ( fasta_file , chrs ) ;
		vector <string> vs ;
		for ( int a = 0 ; a < chrs.size() ; a++ ) vs.push_back ( chrs[a].sequence ) ;
		join_single_sequences ( vs ) ;
		for ( int a = 0 ; a < vs.size() ; a++ ) {
			if ( vs[a].length() < 100 ) continue ; // Hard cutoff
			cout << ">CONTIG" << a << "_(" << vs[a].length() << ")" << endl << vs[a] << endl ;
		}
		return 0 ;
	}*/
	
	vector <TRead> vr ;
	

	read_fastq ( fastq_file , vr , 5 ) ;
	make_read_indices ( vr ) ;
	vector <TChromosome> chrs ;
	assemble_contigs ( vr , chrs , "ATG" ) ;

//	join_chromosomes ( vr , index_length_1 , index_length_2 , genome_file ) ;
	





	
	
//	align_perfect ( vr , 0 ) ;
	

	
/*
	vector <string> vs ;
	vector <TRead> vr , vn ;
	read_fastq ( fastq_file , vr ) ;
	vn.swap ( vr ) ;
	int iteration = 0 ;
	do {
		cerr << "Iteration " << ++iteration << endl ;
		cerr << vn.size() << " items left" << endl ;
		vn.swap ( vr ) ;
		vn.clear() ;
		vn.reserve ( vr.size() ) ;
		interlink_reads ( vr ) ;
		for ( uint a = 0 ; a < idx->position.size() ; a++ ) {
			idx->position[a].clear() ;
		}
		create_contigs ( vr , vn ) ;

		// for ( int a = 0 ; a < vr.size() ; a++ ) {
			// if ( vr[a].used ) continue ;
			// vn.push_back ( vr[a] ) ;
		// }

		for ( int a = 0 ; a < vn.size() ; a++ ) {
			vn[a].reset ( a ) ;
		}
		
		if ( remove_full_sequences ( vn , vs ) == 0 && iteration > 5 ) break ;
		if ( iteration > 15 ) break ; // Hard cutoff
	} while ( 1 ) ;
	
//	join_single_sequences ( vs ) ;
	
	for ( int a = 0 ; a < vs.size() ; a++ ) {
		if ( vs[a].length() < 100 ) continue ; // Hard cutoff
		cout << ">CONTIG" << a << "_(" << vs[a].length() << ")" << endl << vs[a] << endl ;
	}*/
	return 0 ;
}


// rm -f reassemble ; make clean ; make reassemble ; ./reassemble --fastq=/nfs/repository/d0071/SLX_ST_1451_BK1/1547_s_2.fastq --pair=76 --genome=contigs.fasta
