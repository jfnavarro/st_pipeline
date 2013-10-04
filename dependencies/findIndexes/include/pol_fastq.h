/* The MIT License                      
   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOS
*/

/*Contact <paul.igor.costea@scilifelab.se>*/

#include <string>
#include "util.h"
#include <sys/types.h>

namespace pol_util
{

  ///> Defines
  #define MAX_LINE 1000

  /**
   * @description Fastq entry reading/writing + manipulation class
   */
  class FastqEntry{
  public:
    /**
     * @brief Entry reader and generator. The only way to get a new instance of a FastqEntry.
     * @brief Call this function to read another entry into memory and get a pointer to it.
     * @return Pointer to newly created FastqEntry
     */
    static FastqEntry* readEntry(FILE *file) {
      char line[4][MAX_LINE];
      int x;
      for (x=0; x<4; ++x) {//Fastq entry has 4 lines.
	      if (NULL == fgets(line[x],MAX_LINE,file)) {
	      //EOF reached
	      return NULL;
	}
	//Remove end line character   
	line[x][strlen(line[x])-1] = '\0';

      }
      //DO some checking
      if (line[0][0] != '@') {
	fprintf(stderr,"[pol_util]:readEntry -> File seems to not contain fastq entries\nRead:\n%s\n%s\n%s\n",line[0],line[1],line[2]);
	return NULL;
      }
      //Entry read successfully. return FastqEntry
      return new FastqEntry(line);

    };

    /**
     * @brief Destructor
     */
    ~FastqEntry(){
      delete[] m_strName;
      delete[] m_strSeq;
      delete[] m_strPlus;
      delete[] m_strQual;
    };

    /**
     * @brief Return entry formatted as a 4 line string.
     */
    std::string toString() {
      std::string entry = "";
      entry += m_strName;
      entry += '\n';
      entry += m_strSeq;
      entry += '\n';
      entry += m_strPlus;
      entry += '\n';
      entry += m_strQual;
      entry += '\n';
      return entry;
    };

    /**
     * @brief toString wrapper for direct file writing.
     */
    void write(FILE* file) {
      fprintf(file,"%s",this->toString().c_str());
    };

    void writeFake(FILE* file) {
      fprintf(file,"%s\n%s\n",m_strName,"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n+\nBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB");
    }

    /**
     * @brief return part of sequnce. from start to end
     * @param end -> default return all: -1
     * @param start -> start of seq. default 0
     */
    std::string getSequence(int start = 0, int end = -1) {
      if (end == -1) end = strlen(m_strSeq);
      std::string seq = "";
      for (int i=start; i<end; ++i) seq += m_strSeq[i];
      return seq;
    };
    /*updated from git*/
    std::string getName() {
      std::string name = m_strName;
      return name;
    };
    
    uint getSize() {
      return strlen(m_strSeq);
    };

    std::string getQuality(int start = 0, int end = -1) {
      if (end == -1) end = strlen(m_strQual);
      std::string qual = "";
      for (int i=start; i<end; ++i) qual += m_strQual[i];
      return qual;
    }
    /*updated from git*/
    void setOptional(const std::string &info) {
      std::string s = m_strPlus;
      s += " " + info;
      delete[] m_strPlus;
      m_strPlus = new char[s.size()+1];
      strncpy(m_strPlus,s.c_str(),s.size()+1);
    }
    /*updated from git*/
    std::string getOptional() {
      return m_strPlus;
    }

    /**
     * @brief Trim read the specified quality and minimum lenght
     * @return True if read is still usable
     */
    bool trim(int qual, int min_length, bool isSolexa=false)
    {
      int charMove = 64;
      if (isSolexa) charMove = 33;
      int trim_qual = qual;
      int s = 0, l, max = 0, max_l = strlen(m_strQual) - 1;

      for (l = strlen(m_strQual) - 1; l >= 1; --l) {
	if (m_strQual[l] - charMove < 0) {
	  fprintf(stderr,"Base quality under 0! Make sure you are using the correct quality qualifications\n");
	}
	s += trim_qual - (m_strQual[l] - charMove);
	if (s < 0) break;
	if (s > max) {
	  max = s; max_l = l;
	}
      }

      //Ignore this if len < cutoff
      if (max_l < min_length) {
	return false;
      }
      //Trim seq and quality strings
      m_strQual[max_l] = '\0';
      m_strSeq[max_l] = '\0';
      return true;
    };

    /**
     * @breif Cut PolyX's from the beginning/end of read
     * @param nucl       - the nucleotide to cut poly regions of
     * @param min_length - minimum lenght accepted for output read. Won't be trimmed if under this size!
     * @param from_end   - cut from the end of read if True, otherwise from the beginning [Def: true]
     * @param count      - how many bases != from nucl to consider before stopping [Def: 2]
     * @return True is read after cutting N is longer than or equal to min_length
     */
    bool removePoly(char nucl, uint min_length, bool from_end = true, uint count = 2)
    {
      uint l;
      if (from_end) {//Removing from the end
	uint c = 0;
	for (l = strlen(m_strSeq) - 1; l >= min_length; --l) {
	  if (m_strSeq[l] != nucl) {
	    c++;
	    if (c >= count) {
	      //This is it! Recuperate pervious base.
	      l++;
	      break;
	    }
	  } else {
	    //Reset c just in case.
	    c = 0;
	  }
	}
	//Ignore this if len < cutoff
	if (l < min_length) {
	  return false;
	}
	m_strQual[l] = '\0';
	m_strSeq[l] = '\0';
	return true;

      } else {//Removing from the beginning.
	uint c = 0;
	for (l = 0; l < strlen(m_strSeq); ++l) {
          if (m_strSeq[l] != nucl) {
            c++;
            if (c >= count) {
              //This is it! Recuperate pervious base.
              l--;
              break;
            }
          } else {
            //Reset c just in case.
            c = 0;
          }
        }
	if (strlen(m_strQual)-l < min_length) {
	  return false;
	}
	//Move mem! Keep /0 termiator
	memmove(m_strQual, m_strQual+l, strlen(m_strQual)-l+1);
	memmove(m_strSeq, m_strSeq+l, strlen(m_strSeq)-l+1);
	//Check is give nucleodite is not overrepresented in the sequence.
	c=0;
	for (l=0; l < strlen(m_strSeq); ++l)
	  if (m_strSeq[l]==nucl) {
	    ++c;
	  }
	//Is it more than 70% of seqnece? then drop it!
	if (c >= l/2)
	  return false;
	return true;
      }
    };

    /**
     * @brief Check if entry contains given sequence.
     * @param seq -> Sequence to search.
     * @param mismatch -> maximum mismatches allowed.
     * @return True is seq found, false otherwise.
     */
    bool contains(std::string seq, int &pos, int mismatch = 0)
    {
      pos = 0;
      bool res = pol_util::find_substring(m_strSeq,seq,pos,mismatch);
      return res;
    };

  

  private:
    /**
     * @brief Constructor
     */
    FastqEntry(char entry[][MAX_LINE]){
      m_strName = new char[strlen(entry[0])+1];
      m_strSeq = new char[strlen(entry[1])+1];
      m_strPlus = new char[strlen(entry[2])+1];
      m_strQual = new char[strlen(entry[3])+1];

      strncpy(m_strName,entry[0],strlen(entry[0])+1);
      strncpy(m_strSeq,entry[1],strlen(entry[1])+1);
      strncpy(m_strPlus,entry[2],strlen(entry[2])+1);
      strncpy(m_strQual,entry[3],strlen(entry[3])+1);
    };

    ///> Sequence name
    char* m_strName;
    ///> Actual sequence of read
    char* m_strSeq;
    ///> Extra information line
    char* m_strPlus;
    ///> Quality string
    char* m_strQual;
  };

}
