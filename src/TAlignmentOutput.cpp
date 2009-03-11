#include "snpomatic.h"

#define MINOCCUR 1

void TAlignmentOutput::init ( TChromosome *_chr ) {
	chr = _chr ;
	align_full.assign ( chr->sequence.length() , false ) ;
	
	uint size = MAXLINES * chr->sequence.length() ;
	align = new char[size] ;
	qalign = new char[size] ;
	
	for ( uint i = 0 ; i < size ; i++ ) align[i] = qalign[i] = ' ' ;
}

void TAlignmentOutput::add_align ( const string &seq , const string &quality , uint pos , int chromosome ) {
	if ( align_full[pos] ) return ;
	int lane , seql = seq.length() ;
	for ( lane = 0 ; lane < MAXLINES ; lane++ ) {
		if ( align[twotoone(lane,pos)] + align[twotoone(lane,pos+seql/2-1)] + align[twotoone(lane,pos+seql-1)] == 96 ) break ;
	}
	
	if ( lane >= MAXLINES ) {
		align_full[pos] = true ;
		return ; // HARD CUT OFF
	}
	
	for ( int i = 0 ; i < seql ; i++ ) {
		uint p = twotoone(lane,pos+i) ;
		align[p] = seq[i] ;
		qalign[p] = quality[i] ;
	}
}

void TAlignmentOutput::show_pileup ( FILE *pileup , bool snps_only ) {
	int search_snps = 0 , found_snps = 0 , bogus = 0 ;
	uint base , l , cnt , p , lastns ;
	uchar common ;
	char *s2 = new char[MAXLINES+5] ;
	char *s3 = new char[MAXLINES+5] ;
	char *t , *u ;
	
	uint count[256] ;

	for ( base = 0 ; base < chr->sequence.length() ; base++ ) {
		if ( snps_only && !isIUPAC[chr->sequence[base]] ) continue ;
		if ( isIUPAC[chr->sequence[base]] ) search_snps++ ;
		count['A']=count['C']=count['G']=count['T']=0 ;
		cnt = 0 ;
		common = ' ' ;
		lastns = 0 ;
		p = twotoone(0,base) ;
		t = s2 ;
		u = s3 ;
		char *c = align + p ;
		char *q = qalign + p ;
		for ( l = 0 ; l < MAXLINES ; l++ , c++ , q++ ) {
			if ( *c && *c != ' ' ) {
				lastns = l ;
				cnt++ ;
				count[*c]++ ;
				common = MERGE_IUPAC ( *c , common ) ;
				if ( *q < MINQUAL + 33 ) *t++ = to_lc[*c] ;
				else *t++ = *c ;
				*u++ = *q ;
			} else {
				*t++ = ' ' ;
				*u++ = ' ' ;
			}
		}
		s2[lastns+1] = 0 ;
		s3[lastns+1] = 0 ;
		
		bool confirmed_snp = false ;
		uint n = 0 ;
		if ( count['A'] >= MINOCCUR ) n++ ;
		if ( count['C'] >= MINOCCUR ) n++ ;
		if ( count['G'] >= MINOCCUR ) n++ ;
		if ( count['T'] >= MINOCCUR ) n++ ;
		if ( n > 1 ) {
			confirmed_snp = true ;
			found_snps++ ;
		}
		
/*		
		if ( common != ' ' && isIUPAC[common] ) {
			confirmed_snp = true ;
			found_snps++ ;
		}
*/

		if ( common != ' ' && isIUPAC[common] && !isIUPAC[chr->sequence[base]] ) bogus++ ;

		if ( snps_only && common == ' ' ) continue ;
		
		char mark = confirmed_snp ? '*' : ' ' ;
		
		string ref ;
		if ( !chr->original_sequence.empty() ) {
			ref += chr->original_sequence[base] ;
			ref += "\t" ;
		}
		ref += chr->sequence[base] ;

		fprintf ( pileup , "%s\t%s\t%d\t%c%c\t%d\t" , chr->name.c_str() , ref.c_str() , base+1 , common , mark , cnt ) ;
		fprintf ( pileup , "%s\t%s\n" , s2 , s3 ) ;
	}
}
