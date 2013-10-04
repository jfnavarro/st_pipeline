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

#include "kmerFinder.h"


static const int NUMBER_TRIES = 5;

kmerFinder::kmerFinder()
{
    //logging
    counter = 0;
    perfectMatch = 0;
    totalReads = 0;
    ambiguous = 0;
    editDistanceTooBig = 0;
    
    //globals..really? NOTE try to put them into editDistance
    qualMin = 1000;
    qualMax = 0;
    
    Nmismatch = 3;
    kLen = 7;
    probeStartPos = 0;
    probeLength = 18;
}

kmerFinder::kmerFinder(const kmerFinder & other)
{

}

kmerFinder::~kmerFinder()
{
    //TODO should destroy matrices
}

//assumes input files are correct
void kmerFinder::run(FILE *ids, FILE *inF, FILE *out) 
{
    this->ids = ids;
    this->inF = inF;
    this->out = out;
    //internal parameters
    pos = probeStartPos+probeLength;
    leftShift = (probeStartPos-2 >= 0)? 2 : 0;
    rightShift = 2;
    maxKmerHits = (probeLength+rightShift) - (leftShift);
    
    //Load IDS!
    loadIds(ids,mainIds_map);
    
    //Create kMer map
    createKmerMap(kLen); 

}
//input = ids file, output = hash map of id -> (x,y)
bool kmerFinder::loadIds(FILE* snpFile, boost::unordered_map<std::string,Location>& map)
{
    std::ifstream file(snpFile);
    std::string line;
    while(std::getline(file, line))
    {
        std::stringstream linestream(line);
        std::string id;
        int x;
        int y;
        //std::getline(linestream, data, '\t');
        linestream >> id >> x >> y;
        map[id] = Location(x,y);
    }
    return true;
}

//compute penalty distance scores using dynamic programming.
//input : two ids strings, the quality string
//output : the computed score and minimum distance
int kmerFinder::editDistance(const std::string &s1, const std::string &s2, 
                             const std::string &qual, int maxDist, int &alignScore)
{
    int xLen = s1.size();
    int yLen = s2.size();
    
    unsigned int **d = new unsigned int*[xLen+1];
    
    for (int i=0;i<=xLen;++i) 
    {
        d[i] = new unsigned int[yLen+1];
    }
    
    d[0][0] = 0;
    for(int x = 1; x <= xLen; ++x) d[x][0] = 0;
    for(int x = 1; x <= yLen; ++x) d[0][x] = x;
    
    for(int x = 1; x <= xLen; ++x)
    {
        for(int y = 1; y <= yLen; ++y) 
        {
            if (qual[x-1] > qualMax) 
            {
                qualMax = qual[x-1];
            } 
            else if (qual[x-1] < qualMin) 
            {
                qualMin = qual[x-1];
            }
            d[x][y] = std::min( std::min(d[x - 1][y] + 1,d[x][y - 1] + 1),
                                d[x - 1][y - 1] + (s1[x - 1] == s2[y - 1] ? 0 : 1) );
        }
    }
    //Find min for sub-global alignment
    int min = 100;
    int iPos = 0;
    for (int i = xLen; i>0/*xLen-extraBases-2*/;--i) 
    {
        if (d[i][yLen] < min) 
        {
            min = d[i][yLen];
            iPos = i;
        }
    }
    
    if (min <= maxDist) // && alignScore != 0
    {   //Do traceback
        //Compute score!
        alignScore = 0;
        int i = iPos;
        int j = yLen;
        
        while ((i > 0) && (j > 0))
        {
            int Score = d[i][j];
            int ScoreDiag = d[i - 1][j - 1];
            int ScoreUp = d[i][j - 1];
            int ScoreLeft = d[i - 1][j];
            if (Score == ScoreDiag + (s1[i-1] == s2[j-1] ? 0 : 1))
            {
                if (s1[i-1] != s2[j-1])
                {
                    alignScore -= qual[i-1]-33;
                }
                else
                {
                    alignScore += qual[i-1]-33;
                }
                i = i - 1;
                j = j - 1;
            }
            else if (Score == ScoreLeft + 1)
            {
                //Add deletion score
                alignScore -= 100;
                i = i - 1;
            }
            else if (Score == ScoreUp + 1)
            {
                //Add insertions score
                alignScore -= 100;
                j = j - 1;
            }
        }
    }
    for (int i=0;i<=xLen;++i) 
    {
        delete[] d[i];
    }
    delete[] d;
    
    //FAAAKE!
    //alignScore = 0;
    return min;
}

void kmerFinder::createKmerMap(int kLen)
{
    boost::unordered_map<std::string,Location>::const_iterator it;
    for (it = mainIds_map.begin(); it != mainIds_map.end(); ++it) 
    {//For each, introduce only part of it into the new map!
        for (unsigned int i=0; i<=(it->first.size()-kLen); ++i) 
        {
            std::string newKey = it->first.substr(i,kLen);
            IdStruct idS;
            idS.id = it->first;
            idS.loc = it->second;
            kIds_map[newKey].push_back(idS);
        }
    }
}

void kmerFinder::registerHit(pol_util::FastqEntry* e, const std::string &id, int x, int y, bool perfect)
{
    GUARDED_INC(counter)
    
    if(perfect)
    {
        GUARDED_INC(perfectMatch)
    }
    
    //NOTE add option to save fastq
    
    //Print this to output map!
    pthread_mutex_lock(&writeMutex);
    fprintf(out,"%s\t%s\t%d\t%d\t%s\t%s\n",e->getName().c_str(),
            id.c_str(),x,y,e->getQuality().c_str(),e->getSequence().c_str());
    pthread_mutex_unlock(&writeMutex);
}

void kmerFinder::buildKMerMap(const std::string &id,  
                              std::list<std::string>* orderedIds)
{           
    
    boost::unordered_map<std::string,std::list<IdStruct> >::iterator splitIt;
    std::list<IdStruct>::const_iterator it;
    boost::unordered_map<std::string,int>::iterator kIt;
    boost::unordered_map<std::string,int> kMerHitMap;
    
    int x = 0;
    bool forced = false;
    while ( (x <= (id.size()-kLen) ) && (!forced) ) 
    {
        std::string id1 = id.substr(x,kLen);
        splitIt = kIds_map.find(id1);
        if (splitIt != kIds_map.end()) 
        {
            //Add to counted map
            std::list<IdStruct>* l = &(splitIt->second);
            for (it = l->begin(); it != l->end(); ++it) 
            {
                std::string locId = it->id;
                kIt = kMerHitMap.find(locId);
                if (kIt != kMerHitMap.end()) 
                {
                    kMerHitMap[locId]+=1;
                } 
                else 
                {   //Add
                    kMerHitMap[locId]=1;
                }
            }
        }
        //Make sure we get the last overlapping k-mer
        x += kLen;
        if ( x > (id.size()-kLen) )
        {
            x = id.size()-kLen;
            forced = true;
        }
    }
    
    //Go through map and find maximum.
    //Order hits
    for (kIt = kMerHitMap.begin(); kIt != kMerHitMap.end(); ++kIt) 
    {
        orderedIds[kIt->second].push_back(kIt->first);
    }
}

//input = list of ids, id and qual to search. Output = min distance, id and location of found match
bool kmerFinder::searchBestHit(const std::list<std::string> &orderedIds, const std::string &id,
                               const std::string &qual, int &min_ed, std::string &found_id, Location &l)
{
    bool goodHit = false;
    int max_score = 0;
    int goOn = 0;
    int score = 0;
    
    //Take only first x most kMer containing Id's
    for (int i = maxKmerHits-1; i >= 0; i--) 
    {
        //goOn = 0;//TOTRY
        
        //This is the list of most hits
        for (std::list<std::string>::const_iterator it = orderedIds[i].begin(); 
             it != orderedIds[i].end() && (goOn <= NUMBER_TRIES); ++it) 
            {
                double ed = editDistance(id,(*it),qual,Nmismatch,&score);

                if (ed < min_ed) 
                {   //New min
                    min_ed = ed;
                    max_score = score;
                    goodHit = true;
                    found_id = (*it);
                    l = mainIds_map[found_id];
                } 
                else if (ed == min_ed) 
                {   //Two with same min :(
                    //Discriminate by score!
                    if (score == max_score) 
                    {
                        goodHit = false;
                        goOn = 0;
                    } 
                    else if (score > max_score) 
                    {
                        goodHit = true;
                        max_score = score;
                        found_id = (*it);
                        l = mainIds_map[found_id];
                    }
                }
            }
            
         ++goOn;
            
        //orderedIds[i].clear(); ??
    }
    
    return goodHit;
}

/**
 * Thread entry
 */
void* kmerFinder::idSearch()
{
    //boost::unordered_map<std::string,int> kMerHitMap;
    pol_util::FastqEntry* e;
    
    while (true) 
    {
        //Lock file mutex
        pthread_mutex_lock(&readMutex);
        e = pol_util::FastqEntry::readEntry(inF);
        pthread_mutex_unlock(&readMutex);
        
        if (e == NULL) break; //NOTE use another exit method
        
        GUARDED_INC(totalReads)
        //ID should be left of LP sequence
        std::string id = e->getSequence(pos-probeLength,pos);
        std::string qual = e->getQuality(pos-probeLength,pos);
        //Search ID in map!
        boost::unordered_map<std::string,Location>::iterator it = mainIds_map.find(id);
        if (it != mainIds_map.end()) 
        {    //Found!
            registerHit(e,id,it->second.x,it->second.y,true);
        } 
        else 
        {   //Probably with mismatch!
            //Do kmer search
            //Get longer sequence!!!!
            id = e->getSequence(probeStartPos-leftShift,probeStartPos+probeLength+rightShift);
            qual = e->getQuality(probeStartPos-leftShift,probeStartPos+probeLength+rightShift);
            // build kmer map and order it
            std::list<std::string>* orderedIds = new std::list<std::string>[maxKmerHits];
            buildKMerMap(id,orderedIds); //kMerHitMap
            
            int min_ed = 100;
            std::string found_id("");
            Location l;
            l.x = 0; l.y = 0;
            
            bool goodHit = searchBestHit(orderedIds,id,qual,min_ed,found_id,l);
            
            if (goodHit && min_ed <= Nmismatch) 
            {
                //We have a "best" hit!
                registerHit(e,found_id,l.x,l.y,false);
            } 
            else 
            {
                if (!goodHit) 
                {
                    GUARDED_INC(ambiguous)
                } 
                else 
                {
                    GUARDED_INC(editDistanceTooBig)
                }
            }
            
            delete[] orderedIds;
            //kMerHitMap.clear();
        }
    }
    
    delete e;
    return NULL;
}