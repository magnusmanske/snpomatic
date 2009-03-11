#include "snpomatic.h"

/// Shows possible parameters, and an (optional) error message.
void die_with_description ( const char *msg = NULL ) {
	if ( msg ) cerr << msg << endl ;
	cerr << "Parameters : \n" ;
	cerr << "--genome=GENOME_FILE\tFASTA file with chromosomes (mandatory)\n" ;
	cerr << "--fasta=FASTA_FILE\tFASTA file with Solexa reads (mandatory, except when --fastq is used)\n" ;
	cerr << "--output=SNP_FILE\tWhere the output goes (mandatory)\n" ;
	cerr << "--mincontigsize=NUMBER\tMinimal size of contigs (optional)\n" ;
	exit ( 0 ) ;
}

void analyze_contigs ( vector <TChromosome> &contigs )  {
	vector <int> sizes ;
	int a ;
	for ( a = 0 ; a < contigs.size() ; a++ ) {
		int l = contigs[a].sequence.length() ;
		if ( l > 99 ) l = 100 ;
		while ( sizes.size() <= l ) sizes.push_back ( 0 ) ;
		sizes[l]++ ;
	}
	
	for ( a = 0 ; a < sizes.size() ; a++ ) {
		if ( sizes[a] == 0 ) continue ;
		cout << a << "\t" << sizes[a] << endl ;
	}
	exit ( 0 ) ;
}

int main ( int argc , char **argv ) {
	string genome_file , fasta_file , output_file ;
	uint min_contig_size = 0 ;
	uint index_length_1 = 8 ;
	uint index_length_2 = 16 ;

	// Parameters
	static struct option long_options[] = {
		{ "fasta" , optional_argument , 0 , 'a' } ,
		{ "genome" , optional_argument , 0 , 'g' } ,
		{ "output" , optional_argument , 0 , 'o' } ,
		{ "mincontigsize" , optional_argument , 0 , 'm' } ,
		{ 0 , 0 , 0 , 0 }
	} ;
	
		int c ;
	int option_index = 0;
	while ( -1 != ( c = getopt_long (argc, argv, "g:a::o::s::f::p::b::m::",long_options, NULL) ) ) {
		
		switch ( c ) {
			case 0 : break ;
			case 'a' : fasta_file = optarg ; break ;
			case 'g' : genome_file = optarg ; break ;
			case 'o' : output_file = optarg ; break ;
			case 'm' : min_contig_size = atoi ( optarg ) ; break ;
			default : die_with_description("Unknown option") ;
		}
	}

	if ( genome_file.empty() ) die_with_description("No genome file given!") ; // Need one
	if ( fasta_file.empty() ) die_with_description("No Solexa file given!") ; // Need one
	if ( output_file.empty() ) die_with_description("No output file given!") ; // Need one

	init_iupac () ;

	vector <TChromosome> contigs ;
	load_all_chromosomes ( fasta_file , contigs ) ;
	TChromosomalIndices contig_index ( index_length_1 , index_length_2 ) ;
	contig_index.noshortcuts = true ;
	contig_index.run_ends ( &contigs , min_contig_size ) ;
	
//	analyze_contigs ( contigs ) ;
	
	vector <TChromosome> chrs ;
	load_all_chromosomes ( genome_file , chrs ) ;
	TChromosomalIndices genome_index ( index_length_1 , index_length_2 ) ;
	genome_index.noshortcuts = true ;
	genome_index.run ( &chrs ) ;

	TChromosomeAlign ca ( chrs , genome_index ) ;
	ca.align_contigs ( contigs , contig_index , 3 , output_file ) ;

}

/*
./velveth temp1 21 -shortPaired ~/som++/q1_no_match.fasta
./velvetg temp1

rm mapcontigs ; make mapcontigs ; time ./mapcontigs --genome=/nfs/malaria/data2/3D7/Pf.3D7.v2.1.chromosomes.dna.fa --fasta=../velvet/velvet_0.5.05/temp1/contigs.fa --out=new_snps_60 --mincontigsize=60
*/