
#ifndef ARGPARSE_H
#define ARGPARSE_H

#include <stdio.h>
#include <configparser.hpp>
// argument defaults

#define PROGNAME                 "cawlign"
#define DEFAULT_DATA_TYPE        nucleotide
#define DEFAULT_REFERENCE        "HXB2_pol"
#define DEFAULT_SCORING          "Nucleotide-BLAST"
#define DEFAULT_SPACE            quadratic
#define DEFAULT_LOCAL_TYPE       trim
#define DEFAULT_OUTPUT_FORMAT    refmap
#define DEFAULT_RC_TYPE          none

#include "stringBuffer.h"

#ifndef VERSION_NUMBER
    #define VERSION_NUMBER            "0.0.1"
#endif

#ifndef LIBRARY_PATH
    #define LIBRARY_PATH            "/usr/local/shares/cawlign/"
#endif

#define SCORES_SUBPATH "scoring"
#define REF_SUBPATH    "references"

namespace argparse
{
  
    enum data_t {
      nucleotide,
      codon,
      protein
    };

    enum local_t {
      trim,
      global,
      local
    };

    enum space_t {
        quadratic,
        linear
    };

    enum out_format_t {
        refmap,
        refalign,
        pairwise
    };

    enum rc_t {
        none,
        silent,
        annotated
    };

   class args_t {
    public:
 
        FILE            * output,
                        * reference,
                        * input;
       
        ConfigParser    * scores;
             
        data_t          data_type;
        local_t         local_option;
        space_t         space_type;
        out_format_t    out_format;
        rc_t            reverse_complement;
        
        bool            quiet;
        bool            affine;
        bool            include_reference;
       
       StringBuffer*   memory_ref;
        
      
        args_t( int, const char ** );
        ~args_t();
        
    private:
        void parse_input        ( const char * );
        void parse_reference    ( const char * );
        void parse_output       ( const char * );
        void parse_scores       ( const char * );
        void parse_quiet        ( void );
        void parse_affine       ( void );
        void parse_include_ref  ( void );
        void parse_rc           ( const char * );
        void parse_space_t      ( const char * );
        void parse_data_t       ( const char * );
        void parse_local_t      ( const char * );
        void parse_out_format_t ( const char * );

    };

    void ERROR_NO_USAGE ( const char * msg, ... );
}

#endif // ARGPARSE_H
