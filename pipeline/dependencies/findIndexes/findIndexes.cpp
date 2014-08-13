/*
 * Copyright (C) 2012  <Jose Fernandez Navarro>
 * <jose.fernandez.navarro@scilifelab.com>
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * 
 * Author : Jose Fernandez Navarro     jose.fernandez.navarro@scilifelab.com
 * 
 * Spatial Transcriptomics group.
 * 
 * Spatial Transcriptomics Viewer (stVi).
 * 
 */

/* This tool is based on a tool written by Paul Costea and Pelin Akan : 
 * https://github.com/pelinakan/UBD
 */
  
#include <stdio.h>  
#include "kmerFinder.h"

using namespace std;

static int print_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: findIndexes \n");
	fprintf(stderr, "Contact: Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>\n\n");
	fprintf(stderr, "Usage:   findIndexes [options] <ids.txt> <transcripts.fastq>\n\n");
	fprintf(stderr, "Options: \n");
	fprintf(stderr, "         -o STR     output name->id hash for found reads to STR file\n");
	fprintf(stderr, "         -m INT     allowed mismatches [3]\n");
	fprintf(stderr, "         -k INT     kMer length [7]\n");
	fprintf(stderr, "         -s INT     start position of ID [0]\n");
	fprintf(stderr, "         -l INT     length of ID [0]\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{

	FILE *out = NULL;
    int Nmismatch = 3;
    int kLen = 7;
    int probeStartPos = 0;
    int probeLength = 18;
    
	int arg;
	//Get args
	while ((arg = getopt(argc, argv, "o:m:k:s:l")) >= 0) 
    {
		switch (arg) 
        {
            case 'o':
                if ((out = fopen(optarg, "w")) == 0) 
                {
                    fprintf(stderr, "Failed to open for writing: %s\n", optarg);
                    out = NULL;
                }
                fprintf(stdout,"Printing map to %s\n",optarg);
                break;
            case 'm': Nmismatch = atoi(optarg); break;
            case 'k': kLen = atoi(optarg); break;
            case 's': probeStartPos = atoi(optarg); break;
            case 'l': probeLength = atoi(optarg); break;
            default:
                fprintf(stderr,"Read wrong arguments! \n");
                break;
		}
	}

	if (argc-optind != 2) 
    {
		//Not enough paramters!
		print_usage();
		return 1;
	}
	if (probeLength <= 0) 
    {
		fprintf(stderr,"No support for non-positional ID's");
		return 1;
	}
	//Seed rand
	srand((unsigned)time(0));
	//Open files!
	FILE *ids,*inF;
	if ((ids = fopen(argv[optind], "r")) == 0) 
    {
		fprintf(stderr, "Failed to open file %s\n", argv[optind]);
		return 1;
	}
	if ((inF = fopen(argv[optind+1], "r")) == 0) 
    {
		fprintf(stderr, "Failed to open file %s\n", argv[optind+1]);
		return 1;
	}

	kmerFinder kmerfinder;
    kmerfinder.kLen = kLen;
    kmerfinder.Nmismatch = Nmismatch;
    kmerfinder.probeLength = probeLength;
    kmerfinder.probeStartPos = probeStartPos;
    kmerfinder.ids = ids;
    kmerfinder.inF = inF;
    kmerfinder.out = out;
    
    pthread_t threadList[8];
//     threadArgs args;
//     args.inFile = inF;
//     args.outFile = out;
    
    /*DO THREAD STUFF*/
    for (int i=0; i<8; ++i) 
    {
        pthread_create(&threadList[i], NULL, kmerfinder::idSearch, NULL);
    }
    //Now join them...
    for (int i=0; i<8; ++i) 
    {
        pthread_join(threadList[i],NULL);
    }
    
    
    fprintf(stdout, "Found: %ld\n",kmerfinder.counter);
    fprintf(stdout, "Perfect match: %ld\n",kmerfinder.perfectMatch);
    fprintf(stdout, "Ambiguous ID's: %ld\n",kmerfinder.ambiguous);
    fprintf(stdout, "Edit distance exceeding limit: %ld\n",kmerfinder.editDistanceTooBig);
    fprintf(stdout, "Containing polyT after ID: %ld\n",kmerfinder.hasPolyTAferID);
    fprintf(stdout, "Only contains polyT further down: %ld\n",kmerfinder.hasPolyTFurtherDown);
    fprintf(stdout, "Containing polyA after ID: %ld\n",kmerfinder.hasPolyAAfterID);
    fprintf(stdout, "Containing polyA after ID_polyT: %ld\n",kmerfinder.hasManyPolyA);
    fprintf(stdout, "Has many T's even after the needed polyT: %ld\n",kmerfinder.hasManyPolyT);
    fprintf(stdout, "Deletions in ID's : %d\nInsertions in ID's: %d\n",kmerfinder.probeDeletions,probeInsertions);
    fprintf(stdout, "Different Ids: %d\n",kmerfinder.diffIds);
    fprintf(stdout, "Different Ids with perfect matches: %d\n",kmerfinder.diffIdsOnPerfectMatching);
    fprintf(stdout, "Total reads: %ld\n",kmerfinder.totalReads);
    
	fclose(ids);
	fclose(inF);
	if (out != NULL) fclose(out);
	return 0;
}


