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


using namespace std;

static int print_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: findIndexes \n");
    fprintf(stderr, "Contact: Paul Costea <paul.igor.costea@scilifelab.se>\n\n");
    fprintf(stderr, "Usage:   findIndexes [options] <ids.txt> <fw_1.fastq>\n\n");
    fprintf(stderr, "Options: \n");
    fprintf(stderr, "         -o STR     output name->id hash for found reads to STR file\n");
    fprintf(stderr, "         -d STR     output fastq records for not found reads to STR file\n");
    fprintf(stderr, "         -m INT     allowed mismatches [2]\n");
    fprintf(stderr, "         -k INT     kMer length [8]\n");
    fprintf(stderr, "         -s INT     start position of ID [0]\n");
    fprintf(stderr, "         -l INT     length of ID [0]\n");
    fprintf(stderr, "         -p         print species\n");
    fprintf(stderr, "\n");
    return 1;
}

typedef struct {
    int x,y;
} Location;

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

int probeDeletions;
int probeInsertions;
long int counter = 0;
long int perfectMatch = 0;
long int hasPolyTAferID = 0;
long int hasPolyTFurtherDown = 0;
long int hasManyPolyT = 0;
long int hasManyPolyA = 0;
long int hasPolyAAfterID = 0;
long int totalReads = 0;

long int ambiguous = 0;
long int editDistanceTooBig = 0;

//Parameter
int Nmismatch = 2;
int kLen = 8;
int probeStartPos = 0;
int probeLength = 18;
bool outputDiscarded = false;

int qualMin = 1000;
int qualMax = 0;

int badMapping = 0;
int recoverable = 0;

pthread_mutex_t readMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t writeMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t countersMutex = PTHREAD_MUTEX_INITIALIZER;

#define GUARDED_INC(x) \
pthread_mutex_lock(&countersMutex); \
++x; \
pthread_mutex_unlock(&countersMutex);

/////////////////Hash declaration//////////////////////
boost::unordered_map<std::string,Location> mainIds_map;
boost::unordered_map<std::string,std::list<IdStruct> > kIds_map;
///////////////////////////////////////////////////////

#define MATRIX_X 700
#define MATRIX_Y 700

int posMatrix[MATRIX_X][MATRIX_Y];
int perfectMatchMatrix[MATRIX_X][MATRIX_Y];
int uniqueIdMatrix[MATRIX_X][MATRIX_Y];

int scoreDistr[1000];

bool loadIds(FILE* snpFile, boost::unordered_map<std::string,Location>& map, int mat[][MATRIX_Y])
{
    char line[100000];
    char *tok = new char[10000];
    for (int i=0; i<10000; ++i)
        tok[i]='\0';
    int unique = 0;
    while (fgets(line,100000,snpFile)) {
        ++unique;
        const char *t = line;
        int pos = 0;
        std::string id = "";
        Location p;
        p.x=-1;p.y=-1;
        while (*t) {
            t = pol_util::toksplit(t, '\t', tok, 10000);
            if (pos == 0) {//ID
                id = tok;
            } else if (pos == 1) {
                p.x = atoi(tok);
            } else if (pos == 2) {
                p.y = atoi(tok);
            }
            ++pos;
        }
        mat[p.x][p.y]=unique;
        map[id] = p;
    }
    delete[] tok;
    
    return true;
}

void editDistanceWithTraceback(std::string s1, std::string s2)
{
    
    int xLen = s1.size();
    int yLen = s2.size();
    int extraBases = xLen-yLen;
    
    int x,y;
    unsigned int **d = new unsigned int*[xLen+1];
    for (int i=0;i<=xLen;++i) {
        d[i] = new unsigned int[yLen+1];
    }
    d[0][0] = 0;
    for(int x = 1; x <= xLen; ++x) d[x][0] = 0;
    for(int x = 1; x <= yLen; ++x) d[0][x] = x;
    
    for(x = 1; x <= xLen; ++x) {
        for(y = 1; y <= yLen; ++y) {
            d[x][y] = std::min( std::min(d[x - 1][y] + 1,d[x][y - 1] + 1),d[x - 1][y - 1] + (s1[x - 1] == s2[y - 1] ? 0 : 1) );
            //fprintf(stdout,"%d\t",d[x][y]);
        }
        //fprintf(stdout,"\n");
    }
    //fprintf(stdout,"\n");
    //Find min for sub-global alignment
    unsigned int min = 100;
    int iPos = 0;
    for (int i=xLen;i>=1;--i) {
        if (d[i][yLen] < min) {
            min = d[i][yLen];
            iPos = i;
        }
    }
    
    std::string aS1 = "";
    std::string aS2 = "";
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
            aS1 = s1[i-1] + aS1;
            aS2 = s2[j-1] + aS2;
            i = i - 1;
            j = j - 1;
        }
        else if (Score == ScoreLeft + 1)
        {
            aS1 = s1[i-1] + aS1;
            aS2 = "-" + aS2;
            i = i - 1;
            GUARDED_INC(probeInsertions);
        }
        else if (Score == ScoreUp + 1)
        {
            aS1 = "-" + aS1;
            aS2 = s2[j-1] + aS2;
            j = j - 1;
            GUARDED_INC(probeDeletions);
        }
    }
    /*while (i > 0)
     *    {
     *        aS1 = s1[i-1] + aS1;
     *        aS2 = "-" + aS2;
     *        i = i - 1;
     }*/
    while (j > 0)
    {
        aS1 = "-" + aS1;
        aS2 = s2[j-1] + aS2;
        j = j - 1;
    }
    
    fprintf(stdout,"%s\n%s\n\n",aS1.c_str(),aS2.c_str());
    for (int i=0;i<=xLen;++i) {
        delete[] d[i];
    }
    delete[] d;
    
     }
     
     int editDistance(std::string s1, std::string s2, std::string qual="", int* alignScore=NULL, int maxDist=1000)
     {
         int xLen = s1.size();
         int yLen = s2.size();
         
         unsigned int **d = new unsigned int*[xLen+1];
         for (int i=0;i<=xLen;++i) {
             d[i] = new unsigned int[yLen+1];
         }
         d[0][0] = 0;
         for(int x = 1; x <= xLen; ++x) d[x][0] = 0;
         for(int x = 1; x <= yLen; ++x) d[0][x] = x;
         
         for(int x = 1; x <= xLen; ++x)
             for(int y = 1; y <= yLen; ++y) {
                 if (qual[x-1] > qualMax) {
                     qualMax = qual[x-1];
                 } else if (qual[x-1] < qualMin) {
                     qualMin = qual[x-1];
                 }
                 d[x][y] = std::min( std::min(d[x - 1][y] + 1,d[x][y - 1] + 1),d[x - 1][y - 1] + (s1[x - 1] == s2[y - 1] ? 0 : 1) );
             }
             
             //Find min for sub-global alignment
             int min = 100;
             int iPos = 0;
             for (int i=xLen;i>0/*xLen-extraBases-2*/;--i) {
                 if (d[i][yLen] < min) {
                     min = d[i][yLen];
                     iPos = i;
                 }
             }
             
             if ((alignScore != NULL) && (min <= maxDist)){//Do traceback
                 //Compute score!
                 *alignScore = 0;
                 
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
                             *alignScore -= qual[i-1]-33;
                         else
                             *alignScore += qual[i-1]-33;
                         i = i - 1;
                         j = j - 1;
                     }
                     else if (Score == ScoreLeft + 1)
                     {
                         //Add deletion score
                         *alignScore -= 100;
                         i = i - 1;
                     }
                     else if (Score == ScoreUp + 1)
                     {
                         //Add insertions score
                         *alignScore -= 100;
                         j = j - 1;
                     }
                 }
     }
     for (int i=0;i<=xLen;++i) {
         delete[] d[i];
     }
     delete[] d;
     
     //FAAAKE!
     *alignScore = 0;
     
     return min;
}

void createKmerMap(int kLen)
{
    boost::unordered_map<std::string,Location>::const_iterator it;
    for (it = mainIds_map.begin(); it != mainIds_map.end(); ++it) {//For each, introduce only part of it into the new map!
        for (unsigned int i=0; i<=(it->first.size()-kLen); ++i) {
            std::string newKey = it->first.substr(i,kLen);
            IdStruct idS;
            idS.id = it->first;
            idS.loc = it->second;
            kIds_map[newKey].push_back(idS);
        }
    }
    
}

bool doFullSearch(std::string s, std::string qual, int __mismatch, std::string& found_id)
{
    boost::unordered_map<std::string,Location>::iterator mapIt;
    //Go through the entire list
    int score;
    int min_ed = 100;
    int max_score = 0;
    bool goodHit = false;
    for (mapIt = mainIds_map.begin(); mapIt != mainIds_map.end(); ++mapIt) {
        int ed = editDistance(s,mapIt->first,qual,&score,__mismatch);
        if (ed < min_ed) {//New min
            min_ed = ed;
            max_score = score;
            goodHit = true;
            found_id = mapIt->first;
        } else if (ed == min_ed) {//Two with same min :(
            //Discriminate by score!
            if (score == max_score) {
                goodHit = false;
            } else if (score > max_score) {
                goodHit = true;
                max_score = score;
                found_id = mapIt->first;
            }
        }
    }
    return goodHit;
}

void uniqueJoinList(std::list<IdStruct>& mainList,const std::list<IdStruct>* list)
{
    bool found = false;
    std::list<IdStruct>::iterator mainIt;
    std::list<IdStruct>::const_iterator it;
    for (it = list->begin(); it != list->end(); ++it) {
        found = false;
        for (mainIt = mainList.begin(); mainIt != mainList.end(); ++mainIt) {
            if (it->id.compare(mainIt->id) == 0) {
                found = true;
                break;
            }
        }
        if (!found) {
            mainList.push_back((*it));
        }
    }
}

typedef struct {
    FILE *inFile,*outFile,*outDiscarded;
}threadArgs;

/**
 * Thread entry
 */
void* idSearch(void* arguments)
{
    boost::unordered_map<std::string,Location>::iterator it;
    //std::list<IdStruct> all = std::list<IdStruct>();
    boost::unordered_map<std::string,int> kMerHitMap;
    pol_util::FastqEntry* e;
    threadArgs args = *(threadArgs*)arguments;
    while (true) {
        //Lock file mutex
        pthread_mutex_lock(&readMutex);
        e = pol_util::FastqEntry::readEntry(args.inFile);
        pthread_mutex_unlock(&readMutex);
        if (e == NULL)
            break;
        int pos = probeStartPos+probeLength;
        bool has_seq = true;
        if (has_seq) {
            GUARDED_INC(totalReads)
            //ID should be left of LP sequence
            std::string id = e->getSequence(pos-probeLength,pos);
            std::string qual = e->getQuality(pos-probeLength,pos);
            //Search ID in map!
            it = mainIds_map.find(id);
            if (it != mainIds_map.end()) {//Found!
                GUARDED_INC(counter)
                GUARDED_INC(perfectMatch)
                int tempP = 0;
                bool hasT = pol_util::find_substring(e->getSequence(probeLength,probeLength + 21),"TTTTTTTTTTTTTTTTTTTT",tempP,8);
                bool hasManyT = pol_util::find_substring(e->getSequence(probeLength + 21,probeLength + 21 + 20),"TTTTTTTTTTTTTTTTTTTT",tempP,5);
                bool hasManyA = pol_util::find_substring(e->getSequence(probeLength + 21,probeLength + 21 + 20),"AAAAAAAAAAAAAAAAAAAA",tempP,5);
                bool hasA = pol_util::find_substring(e->getSequence(probeLength,probeLength + 31),"AAAAAAAAAAAAAAAAAAAA",tempP,3);
                if (hasT && !hasManyT) {
                    GUARDED_INC(hasPolyTAferID)
                } else if (hasT && hasManyT) {
                    GUARDED_INC(hasManyPolyT)
                } else if (hasManyT) {
                    GUARDED_INC(hasPolyTFurtherDown)
                }
                if (hasA) {
                    GUARDED_INC(hasPolyAAfterID)
                }
                if (hasManyA) {
                    GUARDED_INC(hasManyPolyA)
                }
                GUARDED_INC(posMatrix[it->second.x][it->second.y]);
                GUARDED_INC(perfectMatchMatrix[it->second.x][it->second.y]);
                //Print this to output map!
                if (args.outFile != NULL) {
                    pthread_mutex_lock(&writeMutex);
                    fprintf(args.outFile,"%s\t%s\t%d\t%d\t%s\t%s\n",e->getName().c_str(),id.c_str(),it->second.x,it->second.y,e->getQuality().c_str(),e->getSequence().c_str());
                    pthread_mutex_unlock(&writeMutex);
                }
            } else {//Probably with mismatch!
                //////////////////////////////////////////////////////////////////////////
                ////////////////////////////Do kmer search////////////////////////////////
                //Get longer sequence!!!!
                int leftShift = (probeStartPos-2 >= 0)? 2 : 0;
                int rightShift = 2;
                id = e->getSequence(probeStartPos-leftShift,probeStartPos+probeLength+rightShift);
                qual = e->getQuality(probeStartPos-leftShift,probeStartPos+probeLength+rightShift);
                
                boost::unordered_map<std::string,std::list<IdStruct> >::iterator splitIt;
                std::list<IdStruct>::const_iterator it;
                boost::unordered_map<std::string,int>::iterator kIt;
                
                int x = 0;
                bool forced = false;
                while ( (x <= id.size()-kLen) && (!forced)) {
                    std::string id1 = id.substr(x,kLen);
                    splitIt = kIds_map.find(id1);
                    if (splitIt != kIds_map.end()) {
                        //Add to counted map
                        std::list<IdStruct>* l = &(splitIt->second);
                        for (it = l->begin(); it != l->end(); ++it) {
                            std::string locId = it->id;
                            kIt = kMerHitMap.find(locId);
                            if (kIt != kMerHitMap.end()) {
                                kMerHitMap[locId]+=1;
                            } else {//Add
                                kMerHitMap[locId]=1;
                            }
                        }
                    }
                    //Make sure we get the last overlapping k-mer
                    x+=kLen;
                    if (x > id.size()-kLen) {
                        x = id.size()-kLen;
                        forced = true;
                    }
                }
                
                
                ///////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////
                Location l;
                l.x = 0; l.y = 0;
                bool goodHit = false;
                int min_ed = 100;
                int max_score = 0;
                std::string found_id = "";
                //Go through map and find maximum.
                int maxKmerHits = (probeLength+rightShift) - (leftShift);
                std::list<std::string>* orderedIds = new std::list<std::string>[maxKmerHits];
                //Order hits
                for (kIt = kMerHitMap.begin(); kIt != kMerHitMap.end(); ++kIt) {
                    orderedIds[kIt->second].push_back(kIt->first);
                }
                int goOn = 0;
                int score = 0;
                int searching = 0;
                //Take only first x most kMer containing Id's
                for (int i=maxKmerHits-1;i>=0;i--) {
                    if (orderedIds[i].size() != 0 && ((goOn <= 5) )) {//This is the list of most hits
                        for (std::list<std::string>::const_iterator it=orderedIds[i].begin(); it != orderedIds[i].end(); ++it) {
                            double ed = editDistance(id,(*it),qual,&score,Nmismatch);
                            ++searching;
                            if (ed < min_ed) {//New min
                                min_ed = ed;
                                max_score = score;
                                goodHit = true;
                                found_id = (*it);
                                l = mainIds_map[found_id];
                            } else if (ed == min_ed) {//Two with same min :(
                                //Discriminate by score!
                                if (score == max_score) {
                                    goodHit = false;
                                    goOn = 0;
                                } else if (score > max_score) {
                                    goodHit = true;
                                    max_score = score;
                                    found_id = (*it);
                                    l = mainIds_map[found_id];
                                }
                            }
                        }
                        ++goOn;
                    }
                    orderedIds[i].clear();
                }
                
                if (goodHit && min_ed<= Nmismatch /*&& score >= 200*/) {
                    //We have a "best" hit!
                    //Add score to quality count
                    GUARDED_INC(scoreDistr[score])
                    GUARDED_INC(counter)
                    int tempP = 0;
                    bool hasT = pol_util::find_substring(e->getSequence(probeLength,probeLength + 21),"TTTTTTTTTTTTTTTTTTTT",tempP,8);
                    bool hasManyT = pol_util::find_substring(e->getSequence(probeLength + 21,probeLength + 21 + 20),"TTTTTTTTTTTTTTTTTTTT",tempP,5);
                    bool hasManyA = pol_util::find_substring(e->getSequence(probeLength + 21,probeLength + 21 + 20),"AAAAAAAAAAAAAAAAAAAA",tempP,5);
                    bool hasA = pol_util::find_substring(e->getSequence(probeLength,probeLength + 31),"AAAAAAAAAAAAAAAAAAAA",tempP,3);
                    if (hasT && !hasManyT) {
                        GUARDED_INC(hasPolyTAferID)
                    } else if (hasT && hasManyT) {
                        GUARDED_INC(hasManyPolyT)
                    } else if (hasManyT) {
                        GUARDED_INC(hasPolyTFurtherDown)
                    }
                    if (hasA) {
                        GUARDED_INC(hasPolyAAfterID)
                    }
                    
                    if (hasManyA) {
                        GUARDED_INC(hasManyPolyA)
                    }
                    
                    GUARDED_INC(posMatrix[l.x][l.y]);
                    if (args.outFile != NULL) {
                        pthread_mutex_lock(&writeMutex);
                        //Adding sequence and quality to the output
                        fprintf(args.outFile,"%s\t%s\t%d\t%d\t%s\t%s\n",e->getName().c_str(),found_id.c_str(),l.x,l.y,e->getQuality().c_str(),e->getSequence().c_str());
                        pthread_mutex_unlock(&writeMutex);
                    }
                } else {
                    //discarded
                    if (args.outDiscarded != NULL) {
                        fprintf(args.outDiscarded,"%s\n%s\n%s\n",e->getName().c_str(), e->getSequence().c_str(), e->getQuality().c_str());
                    }
                    
                    if (!goodHit) {
                        GUARDED_INC(ambiguous)
                    } else {
                        GUARDED_INC(editDistanceTooBig)
                    }
                }
                
                delete[] orderedIds;
                kMerHitMap.clear();
                
            }
        }
        delete e;
    }
    return NULL;
}

/**
 * Main of app
 */
int main(int argc, char *argv[])
{
    
    probeDeletions = 0;
    probeInsertions = 0;
    
    for (int i=0; i<MATRIX_X; ++i)
        for (int j=0; j<MATRIX_Y; ++j){
            posMatrix[i][j]=0;
            uniqueIdMatrix[i][j]=0;
            perfectMatchMatrix[i][j]=0;
        }
        
    FILE *out = NULL;
    FILE *outDiscarded = NULL;
    
    int arg;
    //Get args
    while ((arg = getopt(argc, argv, "o:d:m:k:s:l:p")) >= 0) {
        switch (arg) {
            case 'o':
                if ((out = fopen(optarg, "w")) == 0) {
                    fprintf(stderr, "Failed to open for writing: %s\n", optarg);
                    out = NULL;
                }
                fprintf(stdout,"Printing map to %s\n",optarg);
                break;
            case 'd':
                if ((outDiscarded = fopen(optarg, "w")) == 0) {
                    fprintf(stderr, "Failed to open for writing: %s\n", optarg);
                    outDiscarded = NULL;
                }
                fprintf(stdout,"Printing discarded reads to %s\n",optarg);
                break;
            case 'm': Nmismatch = atoi(optarg); break;
            case 'k': kLen = atoi(optarg);
                                                         if (kLen > 18) {//Ids are 18...
                                                             fprintf(stderr,"Id's are 18 in lenght. Can't have a k-mer larger than that!\nUsing default [11]\n");
                    kLen = 11;
                                                         }
                                                         break;
            case 's': probeStartPos = atoi(optarg); break;
            case 'l': probeLength = atoi(optarg); break;
            default:
                fprintf(stderr,"Read wrong arguments! \n");
                break;
        }
    }
    
    if (argc-optind != 2) {
        //Not enough paramters!
        print_usage();
        return 1;
    }
    if (probeLength <= 0) {
        fprintf(stderr,"No support for non-positional ID's");
        return 1;
    }
    //Seed rand
    srand((unsigned)time(0));
    //Open files!
    FILE *ids,*inF;
    if ((ids = fopen(argv[optind], "r")) == 0) {
        fprintf(stderr, "Failed to open file %s\n", argv[optind]);
        return 1;
    }
    if ((inF = fopen(argv[optind+1], "r")) == 0) {
        fprintf(stderr, "Failed to open file %s\n", argv[optind+1]);
        return 1;
    }
    
    //Load IDS!
    loadIds(ids,mainIds_map,uniqueIdMatrix);
    
    //Create kMer map
    createKmerMap(kLen);
    
    pthread_t threadList[8];
    threadArgs args;
    args.inFile = inF;
    args.outFile = out;
    args.outDiscarded = outDiscarded;
    
    /*DO THREAD STUFF*/
    for (int i=0; i<8; ++i) {
        pthread_create( &threadList[i], NULL, idSearch, (void*) &args);
    }
    //Now join them...
    for (int i=0; i<8; ++i) {
        pthread_join(threadList[i],NULL);
    }
    
    fprintf(stdout,"Found: %ld\n",counter);
    fprintf(stdout, "Perfect match: %ld\n",perfectMatch);
    fprintf(stdout, "Ambiguous ID's: %ld\n",ambiguous);
    fprintf(stdout, "Edit distance exceeding limit: %ld\n",editDistanceTooBig);
    fprintf(stdout, "Containing polyT after ID: %ld\n",hasPolyTAferID);
    fprintf(stdout, "Only contains polyT further down: %ld\n",hasPolyTFurtherDown);
    fprintf(stdout, "Containing polyA after ID: %ld\n",hasPolyAAfterID);
    fprintf(stdout, "Containing polyA after ID_polyT: %ld\n",hasManyPolyA);
    fprintf(stdout, "Has many T's even after the needed polyT: %ld\n",hasManyPolyT);
    fprintf(stdout, "Deletions in ID's : %d\nInsertions in ID's: %d\n",probeDeletions,probeInsertions);
    
    int diffIds = 0;
    int diffIdsOnPerfectMatching = 0;
    for (int i=0; i<MATRIX_X; ++i) {
        for (int j=0; j<MATRIX_Y; ++j) {
            if (posMatrix[i][j] != 0) {
                ++diffIds;
            }
            if (perfectMatchMatrix[i][j] != 0) {
                ++diffIdsOnPerfectMatching;
            }
        }
    }
    
    fprintf(stdout,"Different Ids: %d\n",diffIds);
    fprintf(stdout,"Different Ids with perfect matches: %d\n",diffIdsOnPerfectMatching);
    fprintf(stdout, "Total reads: %ld\n",totalReads);
    fclose(ids);
    fclose(inF);
    if (out != NULL) fclose(out);
    if (outDiscarded != NULL) fclose(outDiscarded);
    
    return 0;
}


