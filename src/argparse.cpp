
/* argument parsing ------------------------------------------------------------------------------------------------- */

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <fstream>

using namespace std;

#include "argparse.hpp"
#include "stringBuffer.h"
#include "configparser.hpp"


// some crazy shit for stringifying preprocessor directives
#define STRIFY(x) #x
#define TO_STR(x) STRIFY(x)

namespace argparse {
const char usage[] =
"usage: " PROGNAME " [-h] "
"[-v] "
"[-o OUTPUT] "
"[-r REFERENCE] "
"[-s SCORE] "
"[-t DATATYPE] "
"[-l LOCAL_ALIGNMENT] "
"[-f FORMAT] "
"[-S SPACE] "
"[-a] "
"[-q] "
"[-I] "
"[-R] "
"[FASTA]\n";

const char help_msg[] =
"perform a pairwise alignment between a reference sequence and a set of other sequences\n"
"\n"
"optional arguments:\n"
"  -h, --help               show this help message and exit\n"
"  -v, --version            show " TO_STR (PROGNAME) " version \n"
"  -o OUTPUT                direct the output to a file named OUTPUT (default=stdout)\n"
"  -r REFERENCE             read the reference sequence from this file (default=" TO_STR (DEFAULT_REFERENCE)")\n"
"                           first checks to see if the filepath exists, if not looks inside the res/references directory\n"
"                           relative to the install path (/usr/local/share/cawlign by default)\n"
"  -s SCORE                 read the scoring matrices and options from this file (default=" TO_STR (DEFAULT_SCORING)")\n"
"                           first checks to see if the filepath exists, if not looks inside the res/scoring directory\n"
"                           relative to the install path (/usr/local/share/cawlign by default)\n"
"                           see INSERT URL later for how to read and create alignment scoring files\n"
"  -t DATATYPE              datatype (default=" TO_STR( DEFAULT_DATA_TYPE ) ")\n"
"                           nucleotide : align sequences in the nucleotide space;\n"
"                           protein    : align sequences in the protein space;\n"
"                           codon: align sequences in the codon space (reference must be in frame; stop codons are defined in the scoring file);\n"
"  -R REVERSE_COMPLEMENT    options of reverse complementation [rc] (default=" TO_STR( DEFAULT_RC_TYPE ) ")\n"
"                           none       : do not consider reverse complements of sequences;\n"
"                           silent     : align both the sequence and its rc to the reference, select the one with the highest score and report it;\n"
"                           annotated  : align both the sequence and its rc to the reference, select the one with the highest score and report it\n"
"                                        annotate sequences whose reverse complements were reported in the FASTA by appending '|RC' to the sequence name;\n"
"  -l LOCAL_ALIGNMENT       global/local alignment (default=" TO_STR (DEFAULT_LOCAL_TYPE)")\n"
"                           global : full string alignment; all gaps in the alignments are scored the same\n"
"                           local  : partial string local (smith-waterman type) alignment which maximizes the alignment score\n"
"  -f FORMAT                controls the format of the output (default=" TO_STR( DEFAULT_OUTPUT_FORMAT ) ")\n"
"                           refmap   : aligns query sequences to the reference and does NOT retain instertions relative to the reference;\n"
"                           pairwise : aligns query sequences to the reference and DOES retain instertions relative to the reference;\n"
"                                      no MSA is generated, but rather pair-wise alignments are all reported (2x the number of sequences);\n"
"  -S SPACE                 which version of the algorithm to use (an integer >0, default=" TO_STR( DEFAULT_SPACE ) "):\n"
"                           quadratic : build the entire dynamic programming matrix (NxM);\n"
"                           linear    : use the divide and conquer recursion to keep only 6 columns in memory (~ max (N,M));\n"
"                                       NOT IMPLEMENTED FOR CODON DATA\n"
"  -a                       do NOT use affine gap scoring (use by default)\n"
"  -I                       write out the reference sequence for refmap and refalign output options (default = no) \n"
"  FASTA                    read sequences to compare from this file (default=stdin)\n";

    inline
    void help()
    {
        fprintf( stderr, "%s\n%s", usage, help_msg );
        exit( 1 );
    }

    inline
    void version()
    {
        fprintf( stderr, "%s\n", VERSION_NUMBER);
        exit( 0 );
    }

    /**
     * Prints an error message to stderr and terminates the program.
     *
     * @param msg The format string for the error message, similar to printf formatting.
     * @param ... The additional arguments corresponding to the format specifiers in the `msg` string.
     */
    inline
    void ERROR( const char * msg, ... )
    {
        va_list args;
        fprintf( stderr, "%s" PROGNAME ": error: ", usage );
        va_start( args, msg );
        vfprintf( stderr, msg, args );
        va_end( args );
        fprintf( stderr, "\n" );
        exit( 1 );
    }

    /**
     * Prints an error message to stderr and terminates the program without printing usage information.
     *
     * @param msg The format string for the error message, similar to printf formatting.
     * @param ... The additional arguments corresponding to the format specifiers in the `msg` string.
     */
    void ERROR_NO_USAGE ( const char * msg, ... )
    {
        va_list args;
        fprintf( stderr, "" PROGNAME ": error: ");
        va_start( args, msg );
        vfprintf( stderr, msg, args );
        va_end( args );
        fprintf( stderr, "\n" );
        exit( 1 );
    }

    /**
     * Checks if a file exists at the given path, and if not, attempts to find it in a subpath of the library.
     *
     * @param path The path to the file.
     * @param subpath The subpath within the library where the file may be located.
     * @return A pointer to a FILE if the file is found, or NULL if not.
     */
    FILE* check_file_path (const char* path, const char * subpath) {
        FILE * test = fopen (path, "rb");
        if (!test) {
            StringBuffer sb;
            sb.appendBuffer (LIBRARY_PATH);
            if (sb.getString()[sb.length()-1] != '/') {
                sb.appendChar('/');
            }
            sb.appendBuffer (subpath);
            sb.appendChar('/');
            sb.appendBuffer(path);
            test = fopen (sb.getString(), "rb");
        }
        return test;
    }

    /**
     * Checks if a file exists at the given path (using C++ streams), and if not, attempts to find it in a subpath of the library.
     *
     * @param path The path to the file.
     * @param subpath The subpath within the library where the file may be located.
     * @return An ifstream object with the file opened if found, or not opened if the file does not exist.
     */
    ifstream check_file_path_stream (const char* path, const char * subpath) {
        ifstream test;
        test.open (path);
        if (!test.is_open()) {
            StringBuffer sb;
            sb.appendBuffer (LIBRARY_PATH);
            if (sb.getString()[sb.length()-1] != '/') {
                sb.appendChar('/');
            }
            sb.appendBuffer (subpath);
            sb.appendChar('/');
            sb.appendBuffer(path);
            test.open (sb.getString());
        }
        return test;
    }

    /**
     * Retrieves the next argument from the command-line arguments.
     *
     * @param i The current index in the argument list, which will be incremented.
     * @param argc The total number of command-line arguments.
     * @param argv The array of command-line arguments.
     * @return The next command-line argument.
     */
    const char * next_arg (int& i, const int argc, const char * argv[]) {
        i++;
        if (i == argc)
            ERROR ("ran out of command line arguments");
        
        return argv[i];
        
    }

    /**
     * Constructor for args_t, which parses command-line arguments and sets up configuration options.
     *
     * This constructor processes the command-line arguments, setting up the program's input, output,
     * reference files, and various other configuration options such as data type, space type, and output format.
     * If necessary, default values are assigned to some options.
     *
     * @param argc The number of command-line arguments.
     * @param argv The array of command-line arguments.
     */
    args_t::args_t( int argc, const char * argv[] ) :
    output (stdout),
    reference (nullptr),
    input (stdin),
    scores (nullptr),
    data_type (DEFAULT_DATA_TYPE),
    local_option (DEFAULT_LOCAL_TYPE),
    space_type (DEFAULT_SPACE),
    out_format (DEFAULT_OUTPUT_FORMAT),
    reverse_complement (DEFAULT_RC_TYPE),
    quiet (false),
    affine (true),
    include_reference (false) {
        // skip arg[0], it's just the program name
        for (int i = 1; i < argc; ++i ) {
            const char * arg = argv[i];
            
            if ( arg[0] == '-' && arg[1] == '-' ) {
                if ( !strcmp( &arg[2], "help" ) ) help();
                else if ( !strcmp( &arg[2], "version" ) ) version();
                else
                    ERROR( "unknown argument: %s", arg );
            }
            else if ( arg[0] == '-' ) {
                if ( !strcmp( &arg[1], "h" ) ) help();
                else if (  arg[1] == 'v' ) version();
                else if (  arg[1] == 'o' ) parse_output ( next_arg (i, argc, argv) );
                else if (  arg[1] == 'r' ) parse_reference ( next_arg (i, argc, argv) );
                else if (  arg[1] == 's')  parse_scores( next_arg (i, argc, argv) );
                else if (  arg[1] == 't')  parse_data_t( next_arg (i, argc, argv) );
                else if (  arg[1] == 'f')  parse_out_format_t( next_arg (i, argc, argv) );
                else if (  arg[1] == 'S')  parse_space_t( next_arg (i, argc, argv) );
                else if (  arg[1] == 'l')  parse_local_t( next_arg (i, argc, argv) );
                else if (  arg[1] == 'R')  parse_rc ( next_arg (i, argc, argv) );
                else if (  arg[1] == 'a')  parse_affine ();
                else if (  arg[1] == 'I')  parse_include_ref ();
                else if (  arg[1] == 'q')  parse_quiet ();
                else
                    ERROR( "unknown argument: %s", arg );
            }
            else
                if (i == argc-1) {
                    parse_input (arg);
                } else {
                    ERROR( "unknown argument: %s", arg );
                }
        }
        if ( !reference ) {
            parse_reference ( DEFAULT_REFERENCE );
        }
    }

    /**
     * Destructor for args_t, responsible for cleaning up any resources used (files or memory).
     * It closes the input/output/reference files and deletes the scores object, if applicable.
     */
    args_t::~args_t() {
        if ( output && output != stdout )
            fclose( output );
        
        if ( input && input != stdin)
            fclose (input);
        
        if ( reference )
            fclose (reference);
        
        if ( scores ) {
            delete scores;
        }
    }

    /**
     * Parses the output file path from a command-line argument.
     * Opens the output file for writing. If the argument is "-", stdout is used.
     *
     * @param str The path to the output file.
     */
    void args_t::parse_output( const char * str )
    {
        if ( str && strcmp( str, "-" ) )
            output = fopen( str, "wb" );
        else
            output = stdout;
        
        if ( !output )
            ERROR( "failed to open the OUTPUT file %s", str );
    }

    /**
     * Parses the input file path from a command-line argument.
     * Opens the input file for reading. If the argument is "-", stdin is used.
     *
     * @param str The path to the input file.
     */
    void args_t::parse_input( const char * str )
    {
        if ( str && strcmp( str, "-" ) )
            input = fopen( str, "rb" );
        else
            input = stdin;
        
        if ( !input )
            ERROR( "failed to open the INPUT file %s", str );
    }

    /**
     * Parses the reference file path from a command-line argument.
     * Attempts to open the reference file in the given subpath.
     *
     * @param str The path to the reference file.
     */
    void args_t::parse_reference ( const char * str ) {
        if ( str ) {
            reference = check_file_path (str, REF_SUBPATH);
            if ( ! reference )
                ERROR( "failed to open the REFERENCE file %s", str );
        }
    }

    /**
     * Parses the scores file path from a command-line argument.
     * Opens the scores file using an ifstream and initializes a ConfigParser.
     *
     * @param str The path to the scores file.
     */
    void args_t::parse_scores ( const char * str ) {
        if ( str ) {
            ifstream score_stream = check_file_path_stream(str, SCORES_SUBPATH);
            if ( ! score_stream.is_open() )
                ERROR( "failed to open the SCORES file %s", str );
            scores = new ConfigParser (score_stream);
            
        }
    }

    /**
     * Parses the space type from a command-line argument.
     * Valid options are "linear" or "quadratic".
     *
     * @param str The space type argument.
     */
    void args_t::parse_space_t( const char * str ) {
        if (!strcmp (str, "linear")) {
            space_type = linear;
        } else if (!strcmp (str, "quadratic")) {
            space_type = quadratic;
        } else  {
            ERROR( "invalid algorithm type: %s", str );
        }
    }

    /**
     * Parses the data type from a command-line argument.
     * Valid options are "nucleotide", "codon", or "protein".
     *
     * @param str The data type argument.
     */
    void args_t::parse_data_t( const char * str ) {
        if (!strcmp (str, "nucleotide")) {
            data_type = nucleotide;
        } else if (!strcmp (str, "codon")) {
            data_type = codon;
        } else if (!strcmp (str, "protein")) {
            data_type = protein;
        } else  {
            ERROR( "invalid data type: %s", str );
        }
    }

    /**
     * Parses the local alignment type from a command-line argument.
     * Valid options are "trim", "global", or "local".
     *
     * @param str The local type argument.
     */
    void args_t::parse_local_t( const char * str ) {
        if (!strcmp (str, "trim")) {
            local_option = trim;
        } else if (!strcmp (str, "global")) {
            local_option = global;
        } else if (!strcmp (str, "local")) {
            local_option = local;
        } else  {
            ERROR( "invalid local type: %s", str );
        }
    }

    /**
     * Parses the reverse complement options from a command-line argument.
     * Valid options are "none", "silent", or "annotated".
     *
     * @param str The rc option.
     */
    void args_t::parse_rc ( const char * str ) {
        if (!strcmp (str, "none")) {
            reverse_complement = none;
        } else if (!strcmp (str, "silent")) {
            reverse_complement = silent;
        } else if (!strcmp (str, "annotated")) {
            reverse_complement = annotated;
        } else  {
            ERROR( "invalid reverse complement option: %s", str );
        }
    }

    /**
     * Parses the output format type from a command-line argument.
     * Valid options are "refmap", "refalign", or "pairwise".
     *
     * @param str The output format argument.
     */
    void args_t::parse_out_format_t( const char * str ) {
        if (!strcmp (str, "refmap")) {
            out_format = refmap;
        } else if (!strcmp (str, "refalign")) {
            out_format = refalign;
        } else if (!strcmp (str, "pairwise")) {
            out_format = pairwise;
        } else  {
            ERROR( "invalid output format: %s", str );
        }
    }

    /**
     * Enables the inclusion of the reference file in the output.
     */
    void args_t::parse_include_ref() {
        include_reference = true;
    }

    /**
     * Disables the affine flag, indicating a non-affine alignment.
     */
    void args_t::parse_affine() {
        affine = false;
    }

    /**
     * Enables quiet mode, suppressing certain outputs.
     */
    void args_t::parse_quiet() {
        quiet = true;
    }

}
