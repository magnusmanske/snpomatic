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




// __________________________________________________________________________________________________________________________________________________________________ MAIN

/// Shows possible parameters, and an (optional) error message.
void die_with_description ( const char *msg = NULL ) {
	if ( msg ) cerr << msg << endl ;
	cerr << "Parameters : \n" ;
	cerr << "--genome=GENOME_FILE\tFASTA file with chromosomes (mandatory)\n" ;
	cerr << "--fasta=FASTA_FILE\tFASTA file with Solexa reads (mandatory, except when --fastq or --index is used)\n" ;
	cerr << "--fastq=FASTQ_FILE\tFASTQ file with Solexa reads (mandatory, except when --fasta or --index is used)\n" ;
	cerr << "--nono=FILENAME\tFile with list of read names (!) to ignore (optional)\n" ;
	cerr << "--regions=REGION_FILE\tRegion file for finding new SNPs (optional)\n" ;
	cerr << "--snps=SNP_FILE\t\tSimple SNP file (optional)\n" ;
	cerr << "--gff=GFF_FILE\t\tGFF file with SNPs (optional)\n" ;
	cerr << "--uniqueness=OUTPUT_FILE\tOutput a uniqueness data file for the reference; no Solexa reads needed; implies --noshortcuts (optional)\n" ;
	cerr << "--pileup=OUTPUT_FILE\tOutputs complete pileup into that file (optional)\n" ;
	cerr << "--cigar=OUTPUT_FILE\tOutputs alignment in CIGAR format (optional)\n" ;
	cerr << "--gffout=OUTPUT_FILE\tOutputs reads in GFF format (optional)\n" ;
//	cerr << "--featuretable=OUTPUT_FILE\tOutputs a DDBJ/EMBL/GenBank feature table for matched reads (optional; paired reads only for now)\n" ;
	cerr << "--bins=FILE_PREFIX\tOutputs no-match, single-match and multi-match Solexa reads into prefixed files (optional)\n" ;
	cerr << "--binmask=MASK\tMask of 1s and 0s to turn off individual bins. Order: No match, single match, multi-match, IUPAC. Example: 0100 creates only single-match bin. (optional; default:1111)\n" ;
	cerr << "--pair=NUMBER\t\tFor paired reads, the length of the first part of the read (mandatory for paired reads)\n" ;
	cerr << "--fragment=NUMBER\tFor paired reads, the average fragment length (mandatory for paired reads)\n" ;
	cerr << "--variance=NUMBER\tFor paired reads, the variance of the fragment length to either side (optional; default: 1/4 of fragment size)\n" ;
	cerr << "--wobble=FILENAME\tOutputs a list of possible variations (optional; paired reads only) [UNDER CONSTRUCTION]\n" ;
	cerr << "--coverage=FILENAME\tOutputs (high depth) coverage data (optional)\n" ;
	cerr << "--index=FILENAME\tIndex filename (index will be created if it does not exist; optional)\n" ;
	cerr << "--mspi=NUMBER\t\tMaximum number of SNPs per chromosomal index (optional; default:8)\n" ;
	cerr << "--noshortcuts\t\tWill process all chrososomal regions, even those with lots'o'repeats (optional; no value)\n" ;
	cerr << "--snpsonly\t\tOnly lists found SNPs in the pileup (optional; no value)\n" ;
	cerr << "--fragmentplot=FILENAME\tWrites a plot of fragment size distribution to a file (optional)\n" ;
	exit ( 0 ) ;
}

/// Well, guess what.
int main ( int argc , char **argv ) {
	string genome_file , fastq_file ;
	uint index_length_1 = 4 ;
	uint index_length_2 = 10 ;
	int noshortcuts = true ;
	uint mspi = 8 ;

	static struct option long_options[] = {
		{ "genome" , required_argument , 0 , 'g' } ,
		{ "fastq" , optional_argument , 0 , 'q' } ,
		{ 0 , 0 , 0 , 0 }
	} ;

	int c ;
	int option_index = 0;
	while ( -1 != ( c = getopt_long (argc, argv, "g:a::q::s::f::p::b::r::w::n::",long_options, NULL) ) ) {
		
		switch ( c ) {
			case 0 : break ;
			case 'g' : genome_file = optarg ; break ;
			case 'q' : fastq_file = optarg ; break ;
		}
	}
	
	init_iupac () ;
	TChromosomalIndices index ( index_length_1 , index_length_2 ) ;
	index.noshortcuts = noshortcuts ;
	index.max_snps_per_index = mspi ;
	vector <TChromosome> chrs ;
	load_all_chromosomes ( genome_file , chrs ) ;
	index.run ( &chrs ) ;

	TChromosomeAlign ca ( chrs , index ) ;
	ca.align_solexa_fastq_variety ( fastq_file ) ;


/*
	string genome_file , gff_file , simple_snp_file , fastq_file , fasta_file , pileup_file , bin_prefix , regions_file ;
	string cigar_file , wobble_file , nono , fragmentplot , featuretable , gffout , coverage_file , uniqueness ;
	string index_file ;
	string binmask ( "1111" ) ;
	uint mspi = 8 ;
	uint pair_length = 0 ;
	uint fragment_length = 0 ;
	uint index_length_1 = 10 ;
	uint index_length_2 = 16 ;
	int noshortcuts = false ;
	int snpsonly = false ;
	uint variance = 0 ;

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
		{ "regions" , optional_argument , 0 , 'r' } ,
		{ "noshortcuts" , optional_argument , &noshortcuts , true } ,
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
			case 'v' : variance = atoi ( optarg ) ; break ;
			case 't' : fragment_length = atoi ( optarg ) ; break ;
			case 'm' : mspi = atoi ( optarg ) ; break ;
			case 'c' : cigar_file = optarg ; break ;
			case 'r' : regions_file = optarg ; break ;
			default : die_with_description("Unknown option") ;
		}
	}
	binmask += "0000" ; // Ignore bins if shorter string was supplied
	

	// Sanity checks
	if ( genome_file.empty() ) die_with_description("No genome file given!") ; // Need this
	if ( index_file.empty() ) {
		if ( fastq_file.empty() && fasta_file.empty() && uniqueness.empty() ) die_with_description("No Solexa file given!") ; // Need one
		if ( !fastq_file.empty() && !fasta_file.empty() ) die_with_description("Two Solexa files given!") ; // Can't have both
	}
	if ( !regions_file.empty() || !uniqueness.empty() ) noshortcuts = true ;

	// Incorporating known SNPs and indexing chromosomes
	init_iupac () ;
	TChromosomalIndices index ( index_length_1 , index_length_2 ) ;
	index.noshortcuts = noshortcuts ;
	index.max_snps_per_index = mspi ;
	index.index_file = index_file ;
	vector <TChromosome> chrs ;
	load_all_chromosomes ( genome_file , chrs ) ;

	if ( !gff_file.empty() ) incorporate_all_gff ( gff_file , chrs ) ;
	if ( !simple_snp_file.empty() ) incorporate_simple_snps ( simple_snp_file , chrs ) ;
	
	if ( !regions_file.empty() ) index.run_regions ( &chrs , regions_file ) ;
	else index.run ( &chrs ) ;
	
	if ( !uniqueness.empty() ) {
		index.uniqueness ( uniqueness ) ;
		if ( fastq_file.empty() && fasta_file.empty() ) return 0 ; // My work here is done
	}
	

	// Now look at the SOlexa reads and do the thing
	TChromosomeAlign ca ( chrs , index ) ;
	if ( !pileup_file.empty() ) ca.pileup = fopen ( pileup_file.c_str() , "w" ) ;
	if ( !cigar_file.empty() ) ca.cigar = fopen ( cigar_file.c_str() , "w" ) ;
	if ( !fragmentplot.empty() ) ca.fragmentplot = fopen ( fragmentplot.c_str() , "w" ) ;
	if ( !gffout.empty() ) ca.gffout = fopen ( gffout.c_str() , "w" ) ;
	if ( !coverage_file.empty() ) ca.coverage = fopen ( coverage_file.c_str() , "w" ) ;
	if ( !featuretable.empty() ) {
		ca.featuretable = fopen ( featuretable.c_str() , "w" ) ;
		fprintf ( ca.featuretable , "Key\tLocation/Qualifiers\n" ) ;
	}
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
	ca.init () ;
	if ( !nono.empty() ) ca.add_nono_list ( nono ) ;
	if ( pair_length == 0 ) {
		if ( !fasta_file.empty() ) ca.align_solexa_fasta ( fasta_file ) ;
		if ( !fastq_file.empty() ) ca.align_solexa_fastq ( fastq_file ) ;
	} else {
		ca.wobble = fopen ( string ( wobble_file ).c_str() , "w" ) ;
		if ( !fasta_file.empty() ) ca.align_solexa_paired ( fasta_file , pair_length , fragment_length , variance ) ;
		if ( !fastq_file.empty() ) ca.align_solexa_paired ( fastq_file , pair_length , fragment_length , variance ) ;
	}
	ca.show_pileup ( !regions_file.empty() || snpsonly ) ;
	ca.dump_coverage () ;
*/
	return 0 ;
}

// Example runs:
// rm test ; make clean ; make ; /usr/bin/time ./findknownsnps --genome=/nfs/malaria/data2/3D7/Pf.3D7.v2.1.chromosomes.dna.fa --gff=/nfs/malaria/data2/3D7/known_snps/PfalciparumCombinedSNPs.gff --fasta=/nfs/malaria/data2/3D7/all.solexa.fasta --noshortcuts --bins=test
// make clean ; make ; time ./findknownsnps --genome=/nfs/malaria/data2/human_sequence/Homo_sapiens.NCBI36.48.dna.chromosome.22.fa --fasta=/nfs/malaria/data2/3D7/all.solexa.fasta --noshortcuts --bins=human
// make clean ; make ; /usr/bin/time ./findknownsnps --genome=/nfs/malaria/data2/3D7/Pf.3D7.v2.1.chromosomes.dna.fa --fasta=/nfs/malaria/data2/3D7/all.solexa.fasta --bins=q2
// make clean ; make ; /usr/bin/time ./findknownsnps --genome=/nfs/malaria/data2/3D7/Pf.3D7.v2.1.chromosomes.dna.fa --fasta=/nfs/malaria/data2/3D7/all.solexa.fasta --gff=/nfs/malaria/data2/3D7/known_snps/PfalciparumCombinedSNPs.gff --cigar=test.cigar
// make clean ; make ; /usr/bin/time ./findknownsnps --genome=/nfs/malaria/data2/3D7/Pf.3D7.v2.1.chromosomes.dna.fa --fasta=q1_no_match.fasta --snps=new_snps_60_uniq --bins=q3
// make clean ; make ; /usr/bin/time ./findknownsnps --genome=/nfs/faculty/mm6/3D7_pm.fa --fastq=$MAL_SOLEXA_HOME/run368/368_s_3.fastq --pair=35 --fragment=300 --bins=q3


/* TEST:
	make clean ; make variety ; ./variety --genome=tools/MAL3.925000-927000.fasta --fastq=tools/585_s_5.bin2_no_match.fastq > v.out
*/
