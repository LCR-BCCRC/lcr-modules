/**************************************************************************
#
#  This software/database is "United States Government Work" under the terms of
#  the United States Copyright Act.  It was written as part of the authors'
#  official duties for the United States Government and thus cannot be
#  copyrighted.  This software/database is freely available to the public for
#  use without a copyright notice.  Restrictions cannot be placed on its present
#  or future use.
# 
#  Although all reasonable efforts have been taken to ensure the accuracy and
#  reliability of the software and data, the National Human Genome Research
#  Institute (NHGRI) and the U.S. Government does not and cannot warrant the
#  performance or results that may be obtained by using this software or data.
#  NHGRI and the U.S.  Government disclaims all warranties as to performance,
#  merchantability or fitness for any particular purpose.
# 
#  In any work or product derived from this material, proper attribution of the
#  authors as the source of the software or data should be made, using "NHGRI
#  Genome Technology Branch" as the citation.
#
**************************************************************************/

#include "printcomp.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

Params *parameters;

get_params(argc, argv)
    int argc;
    char **argv;
{
    int i;
    char *emptyptr;

    parameters = (Params *) malloc(sizeof(Params));

    emptyptr = (char *) malloc(sizeof(char));
    strcpy(emptyptr, "");
    parameters->region = emptyptr;
    parameters->bedfile = emptyptr;
    parameters->minqual = 0;
    parameters->mapqual = 0;

    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-fasta")) {
            parameters->fasta = (char *) malloc((strlen(argv[++i]) + 1) * sizeof(char)); 
            strcpy(parameters->fasta, argv[i]);
        }
        else if (!strcmp(argv[i], "-bam1")) {
            parameters->bam1 = (char *) malloc((strlen(argv[++i]) + 1) * sizeof(char)); 
            strcpy(parameters->bam1, argv[i]);
        }
        else if (!strcmp(argv[i], "-bam2")) {
            parameters->bam2 = (char *) malloc((strlen(argv[++i]) + 1) * sizeof(char)); 
            strcpy(parameters->bam2, argv[i]);
        }
        else if (!strcmp(argv[i], "-region")) {
            parameters->region = (char *)malloc((strlen(argv[++i]) + 1) * sizeof(char)); 
            strcpy(parameters->region, argv[i]);
        }
        else if (!strcmp(argv[i], "-minqual")) {
            parameters->minqual = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-mapqual")) {
            parameters->mapqual = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "-bedfile")) {
            parameters->bedfile = (char *)malloc((strlen(argv[++i]) + 1) * sizeof(char));
            strcpy(parameters->bedfile, argv[i]);
        }
    }
}

void free_params() {
    free(parameters->fasta);
    free(parameters->bam1);
    free(parameters->bam2);
    free(parameters->region);
    free(parameters->bedfile);
}
