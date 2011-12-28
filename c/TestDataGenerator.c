//
//  TestDataGenerator.c
//  
//
//  Created by Lifeng Tian on 12/26/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "common.h"
#include "options.h"

/* command line option specifications 
 * copied from Jim Kent
 */
static struct optionSpec optionSpecs[] = {
    {"AlleleBalance", OPTION_INT},
    {"Depth", OPTION_INT},
    {"MAFinPercent_ALL", OPTION_INT},
    {"ExonicFunc",OPTION_STRING},
    {"dbSNP", OPTION_DOUBLE},
    {NULL, 0}
};


/* variables ...
 */
int allelebalance = 0;
int depth = 0;
int MAFinPercent_ALL = 0;
char *exonicFunction = '';
double dbsnp = 0.1;


int RandDbl (void)
{
    return rand() / ((double)RAND_MAX + 1.0);
}

void usage()
{
    printf("blah
           blah
           blah\n");
}

int main(int argc, char ** argv )
{
    optionInit(&argc, argv, optionSpecs);

    if ( optionExists("dbSNP") )
    {
        dbSNP = optionDouble("dbSNP",dbSNP);
        printf("dbSNP=%f\n",dbSNP);
    }

    
    return 1;
}