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

#ifndef KMERFINDER_H
#define KMERFINDER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include <time.h>

#ifdef __APPLE__
#include "boost/unordered_map.hpp"
#else
#include <boost/unordered_map.hpp>
#endif

#include <pthread.h>
#include "pol_fastq.h"

pthread_mutex_t readMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t writeMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t countersMutex = PTHREAD_MUTEX_INITIALIZER;

#define GUARDED_INC(x) \
pthread_mutex_lock(&countersMutex); \
++x; \
pthread_mutex_unlock(&countersMutex);

typedef struct {
    int x,y;
} Location;

typedef struct {
    FILE *inFile,*outFile;
}threadArgs;

class IdStruct {
public:
    IdStruct() {
    };
    
    IdStruct(const IdStruct &other) {
        id = other.id;
        loc = other.loc;
    };
    
    Location loc;
    std::string id;
    
};

class kmerFinder
{

	public:
	
        kmerFinder(int Nmismatch = 3, int klen = 7, int probeStartPos = 0, int probeLength = 18);
        kmerFinder(const kmerFinder & other);
        virtual ~ kmerFinder();
        void run(FILE *ids, FILE *inF, FILE *out);
       
    private:
    
        bool loadIds(FILE* snpFile, boost::unordered_map<std::string,Location>& map, int mat[][MATRIX_Y]);
        
        int editDistance(const std::string &s1,const std::string &s2,
                         const std::string &qual, int maxDist=1000, int* alignScore=0);
        
        void createKmerMap(int kLen);
        void registerHit(pol_util::FastqEntry* e, const std::string &id, int x, int y, bool perfect = false);
        void buildKMerMap(const std::string &id,  
                        boost::unordered_map<std::string,int> &kMerHitMap,std::list<std::string>* orderedIds);
        
        bool searchBestHit(const std::list<std::string> &orderedIds, const std::string &id,
                           const std::string &qual, int &min_ed, std::string &found_id, Location &l);

        void* idSearch();
        
        //globals..really?
        int qualMin; 
        int qualMax;
        
        //Parameters
        int Nmismatch;
        int kLen;
        int probeStartPos;
        int probeLength;
        FILE *ids; 
        FILE *inF; 
        FILE *out;
        
        //for debuging
        long int counter;
        long int perfectMatch;
        long int totalReads;
        long int ambiguous;
        long int editDistanceTooBig;

        //internal parameters
        int pos;
        int leftShift;
        int rightShift;
        int maxKmerHits;

        boost::unordered_map<std::string,Location> mainIds_map;
        boost::unordered_map<std::string,std::list<IdStruct> > kIds_map;
    
    //int scoreDistr[1000];
};

#endif	/* // KMERFINDER_H */
