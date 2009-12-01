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
#include <fstream>

using namespace std ;



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
	cerr << "--inversions=FILENAME    For paired reads, writes read matches indicating inversions into a file (optional)\n" ;
	cerr << "--faceaway=FILENAME      For paired reads, writes read matches that \"face away\" from each other into a file (optional)\n" ;
	cerr << "--sqlite=FILENAME        Creates a sqlite text file with alignment data [EXPERIMENTAL] (optional)\n" ;
	cerr << "--sam=FILENAME           Creates a SAM alignment file (optional)\n" ;
	cerr << "--spancontigs=FILENAME   Outputs read pairs where \"half\" reads map uniquely to different contigs (optional)\n" ;
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
	cerr << "--index_from=NUMBER      Starts indexing at this position on all chromosomes (optional)\n" ;
	cerr << "--index_to=NUMBER        Stops indexing at this position on all chromosomes (optional)\n" ;
	cerr << "--chop=NUMBER            For paired reads, if one but not the other matches, shorten the other by NUMBER bases (optional)\n" ;
	cerr << "--index1=NUMBER          Length of internal index 1 (optional; default:10)\n" ;
	cerr << "--index2=NUMBER          Length of internal index 2 (optional; default:16)\n" ;
	cerr << "--memory_save=NUMBER     Indexes the genome every NUMBER of positions; saves memory and runtime, but can have strange side effects (optional)\n" ;
	cerr << "--multimatch             Puts a multiple-matching read to a random position (optional) [currently paired reads only]\n" ;
	cerr << "--singlematch            Only performs additional output functions for single matches (optional) [currently paired reads only]\n" ;
	cerr << "--foum                   For paired reads, at least one read has to match uniquely in the genome (force one unique match) (optional)\n" ;
	cerr << "--mismatch               The number of mismatches allowed outside the index (index1+index2) (optional)\n" ;
	cerr << "--rpa=FILENAME           Writes all read pair alignments into a file (optional)\n" ;
//	cerr << "--featuretable=OUTPUT_FILE\tOutputs a DDBJ/EMBL/GenBank feature table for matched reads (optional; paired reads only for now)\n" ;
	exit ( 0 ) ;
}

/// Well, guess what.
int main ( int argc , char **argv ) {
	string genome_file , gff_file , simple_snp_file , fastq_file , fastq2_file , fasta_file , pileup_file , bin_prefix , regions_file ;
	string cigar_file , wobble_file , nono , fragmentplot , featuretable , gffout , coverage_file , uniqueness , snpsinreads ;
	string index_file , indelplot_file , inversions_file , chromosome , sqlite , sam_file , spancontigs , face_away , rpa ;
	string binmask ( "1111" ) ;
	uint wobblemax = 2 ;
	uint mspi = 8 ;
	uint index_from = 0 ;
	uint index_to = 0 ;
	uint pair_length = 0 ;
	uint fragment_length = 0 ;
	uint index_length_1 = 10 ;
	uint index_length_2 = 16 ;
	int noshortcuts = false ;
	int snpsonly = false ;
	int chop = 0 ;
	uint variance = 0 ;
	int multimatch = false ;
	int singlematch = false ;
	int force_one_unique_match = false ;
	int memory_save = 1 ;

	// Parameters
	static struct option long_options[] = {
		{ "genome" , required_argument , 0 , 'g' } ,
		{ "fasta" , optional_argument , 0 , 'a' } ,
		{ "fastq" , optional_argument , 0 , 'q' } ,
		{ "fastq2" , optional_argument , 0 , 'Q' } ,
		{ "snps" , optional_argument , 0 , 's' } ,
		{ "coverage" , optional_argument , 0 , 'e' } ,
		{ "sam" , optional_argument , 0 , 'S' } ,
		{ "gff" , optional_argument , 0 , 'f' } ,
		{ "rpa" , optional_argument , 0 , 'R' } ,
		{ "faceaway" , optional_argument , 0 , 'F' } ,
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
		{ "sqlite" , optional_argument , 0 , 'P' } ,
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
		{ "index_from" , optional_argument , 0 , '8' } ,
		{ "index_to" , optional_argument , 0 , '9' } ,
		{ "mismatch" , optional_argument , 0 , 'M' } ,
		{ "spancontigs" , optional_argument , 0 , 'X' } ,
		{ "memory_save" , optional_argument , 0 , 'Y' } ,
		{ "multimatch" , optional_argument , &multimatch , true } ,
		{ "singlematch" , optional_argument , &singlematch , true } ,
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
			case 'S' : sam_file = optarg ; break ;
			case 'Q' : fastq2_file = optarg ; break ;
			case 's' : simple_snp_file = optarg ; break ;
			case 'f' : gff_file = optarg ; break ;
			case 'e' : coverage_file = optarg ; break ;
			case 'p' : pileup_file = optarg ; break ;
			case 'b' : bin_prefix = optarg ; break ;
			case 'z' : binmask = optarg ; break ;
			case 'F' : face_away = optarg ; break ;
			case '0' : uniqueness = optarg ; break ;
			case 'w' : wobble_file = optarg ; break ;
			case 'n' : nono = optarg ; break ;
			case 'x' : index_file = optarg ; break ;
			case 'i' : pair_length = atoi ( optarg ) ; break ;
			case 'o' : fragmentplot = optarg ; break ;
			case 'u' : featuretable = optarg ; break ;
			case 'y' : gffout = optarg ; break ;
			case 'X' : spancontigs = optarg ; break ;
			case 'd' : snpsinreads = optarg ; break ;
			case 'v' : variance = atoi ( optarg ) ; break ;
			case 'Y' : memory_save = atoi ( optarg ) ; break ;
			case 't' : fragment_length = atoi ( optarg ) ; break ;
			case 'm' : mspi = atoi ( optarg ) ; break ;
			case '1' : index_length_1 = atoi ( optarg ) ; break ;
			case '2' : index_length_2 = atoi ( optarg ) ; break ;
			case '3' : chop = atoi ( optarg ) ; break ;
			case '4' : indelplot_file = optarg ; break ;
			case '5' : wobblemax = atoi ( optarg ) ; break ;
			case '6' : inversions_file = optarg ; break ;
			case '7' : chromosome = optarg ; break ;
			case '8' : index_from = atoi ( optarg ) ; break ;
			case '9' : index_to = atoi ( optarg ) ; break ;
			case 'M' : globally_allowed_mismatches = atoi ( optarg ) ; break ;
			case 'c' : cigar_file = optarg ; break ;
			case 'r' : regions_file = optarg ; break ;
			case 'P' : sqlite = optarg ; break ;
			case 'R' : rpa = optarg ; break ;
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
	index.index_from = index_from ;
	index.index_to = index_to ;
	vector <TChromosome> chrs ;
	load_all_chromosomes ( genome_file , chrs ) ;

	if ( !gff_file.empty() ) incorporate_all_gff ( gff_file , chrs ) ;
	if ( !simple_snp_file.empty() ) incorporate_simple_snps ( simple_snp_file , chrs ) ;
	
	if ( !chromosome.empty() ) {
		for ( int a = 0 ; a < chrs.size() ; a++ ) {
			if ( chrs[a].name == chromosome ) continue ;
			chrs.erase ( chrs.begin() + a ) ;
			a-- ;
		}
		if ( chrs.size() == 0 ) {
			cout << "No chromosome with that name (" << chromosome << ") in file -- aborting.\n" ;
			exit ( 1 ) ;
		}
		cout << "Keeping " << chrs.size() << " chromosomes" << endl ;
	}
	
	index.index_steps = memory_save ;
	
	if ( !regions_file.empty() ) index.run_regions ( &chrs , regions_file ) ;
	else index.run ( &chrs ) ;
	
	if ( !uniqueness.empty() ) {
		index.uniqueness ( uniqueness ) ;
		if ( fastq_file.empty() && fasta_file.empty() ) return 0 ; // My work here is done
	}
	
	if ( !fastq_file.empty() && !fastq2_file.empty() && pair_length == 0 ) { // Auto-detect read length
		ifstream in ( fastq_file.c_str() ) ;
		char tmp[10000] ;
		in.getline ( tmp , 10000 ) ;
		in.getline ( tmp , 10000 ) ;
		in.close() ;
		pair_length = strlen ( tmp ) ;
		cout << "Setting read length to " << pair_length << endl ;
	}
	

	// Now look at the Solexa reads and do the thing
	TChromosomeAlign ca ( chrs , index ) ;
	ca.wobblemax = wobblemax ;
	ca.multimatch = multimatch ;
	ca.singlematch = singlematch ;
	ca.force_one_unique_match = force_one_unique_match ;
	if ( !pileup_file.empty() ) ca.pileup = fopen ( pileup_file.c_str() , "w" ) ;
	if ( !cigar_file.empty() ) ca.cigar = fopen ( cigar_file.c_str() , "w" ) ;
	if ( !fragmentplot.empty() ) ca.fragmentplot = fopen ( fragmentplot.c_str() , "w" ) ;
	if ( !indelplot_file.empty() ) ca.indelplot = fopen ( indelplot_file.c_str() , "w" ) ;
	if ( !gffout.empty() ) ca.gffout = fopen ( gffout.c_str() , "w" ) ;
	if ( !snpsinreads.empty() ) ca.snpsinreads = fopen ( snpsinreads.c_str() , "w" ) ;
	if ( !face_away.empty() ) ca.faceaway = fopen ( face_away.c_str() , "w" ) ;
	if ( !coverage_file.empty() ) ca.coverage = fopen ( coverage_file.c_str() , "w" ) ;
	if ( !inversions_file.empty() ) ca.inversions = fopen ( inversions_file.c_str() , "w" ) ;
	if ( !sqlite.empty() ) ca.sqlite = fopen ( sqlite.c_str() , "w" ) ;
	if ( !sam_file.empty() ) ca.sam = fopen ( sam_file.c_str() , "w" ) ;
	if ( !spancontigs.empty() ) ca.spancontigs = fopen ( spancontigs.c_str() , "w" ) ;
	if ( !rpa.empty() ) ca.rpa = fopen ( rpa.c_str() , "w" ) ;
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
	ca.chop = chop ;
	ca.init () ;
	if ( !nono.empty() ) ca.add_nono_list ( nono ) ;
	if ( pair_length == 0 ) {
		if ( !wobble_file.empty() ) ca.wobble = fopen ( string ( wobble_file ).c_str() , "w" ) ;
		if ( !fasta_file.empty() ) ca.align_solexa_fasta ( fasta_file ) ;
		if ( !fastq_file.empty() ) ca.align_solexa_fastq ( fastq_file ) ;
	} else {
		if ( fragment_length == 0 && variance == 0 ) { // If no fragment size is given, do entire chromosome
			for ( uint a = 0 ; a < chrs.size() ; a++ ) {
				if ( fragment_length < chrs[a].sequence.length() ) fragment_length = chrs[a].sequence.length() ;
			}
			fragment_length = fragment_length / 2 ;
			variance = fragment_length ;
		}
		if ( !wobble_file.empty() ) ca.wobble = fopen ( string ( wobble_file ).c_str() , "w" ) ;
		if ( !fasta_file.empty() ) ca.align_solexa_paired ( fasta_file , fastq2_file , pair_length , fragment_length , variance ) ;
		if ( !fastq_file.empty() ) ca.align_solexa_paired ( fastq_file , fastq2_file , pair_length , fragment_length , variance ) ;
	}
	ca.show_pileup ( !regions_file.empty() || snpsonly ) ;
	ca.dump_coverage () ;
	ca.finish_sqlite ( fragment_length , variance , pair_length != 0 ? pair_length : ca.single_read_length ) ;
	
	return 0 ;
}

// Example runs:
// rm test ; make clean ; make ; /usr/bin/time ./findknownsnps --genome=/nfs/malaria/data2/3D7/Pf.3D7.v2.1.chromosomes.dna.fa --gff=/nfs/malaria/data2/3D7/known_snps/PfalciparumCombinedSNPs.gff --fasta=/nfs/malaria/data2/3D7/all.solexa.fasta --noshortcuts --bins=test
// make clean ; make ; time ./findknownsnps --genome=/nfs/malaria/data2/human_sequence/Homo_sapiens.NCBI36.48.dna.chromosome.22.fa --fasta=/nfs/malaria/data2/3D7/all.solexa.fasta --noshortcuts --bins=human
// make clean ; make ; /usr/bin/time ./findknownsnps --genome=/nfs/malaria/data2/3D7/Pf.3D7.v2.1.chromosomes.dna.fa --fasta=/nfs/malaria/data2/3D7/all.solexa.fasta --bins=q2
// make clean ; make ; /usr/bin/time ./findknownsnps --genome=/nfs/malaria/data2/3D7/Pf.3D7.v2.1.chromosomes.dna.fa --fasta=/nfs/malaria/data2/3D7/all.solexa.fasta --gff=/nfs/malaria/data2/3D7/known_snps/PfalciparumCombinedSNPs.gff --cigar=test.cigar
// make clean ; make ; /usr/bin/time ./findknownsnps --genome=/nfs/malaria/data2/3D7/Pf.3D7.v2.1.chromosomes.dna.fa --fasta=q1_no_match.fasta --snps=new_snps_60_uniq --bins=q3
// make clean ; make ; /usr/bin/time ./findknownsnps --genome=/nfs/faculty/mm6/3D7_pm.fa --fastq=$MAL_SOLEXA_HOME/run368/368_s_3.fastq --pair=35 --fragment=300 --bins=q3
// make clean ; make ; /usr/bin/time ./findknownsnps --genome=../3D7_pm.fa --snps=MERGED_CAPILLARY_MAQ.tab.subset  --fastq=$MAL_SOLEXA_DATA/986_s_8/986_s_8.fastq --pair=37 --fragment=200 --variance=200 --sqlite=test.sqlite.text

/*
	// Test params, obsolete
	genome_file = "/nfs/malaria/data2/3D7/Pf.3D7.v2.1.chromosomes.dna.fa" ;
	gff_file = "/nfs/malaria/data2/3D7/known_snps/PfalciparumCombinedSNPs.gff" ;
	fasta_file = "/nfs/malaria/data2/3D7/all.solexa.fasta" ;
	//simple_snp_file = "/nfs/malaria/data2/3D7/known_snps/MAL1.snps" ;
	//fastq_file = "/nfs/malaria/data2/3D7/fastq/PFCLIN_slx.fastq" ;
*/
