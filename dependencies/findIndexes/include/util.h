#include <string>

namespace pol_util
{

  /**
   * @biref Match substring with k missmatches.
   * @return True is substring found!
   * @out pos = Alignment position. Can range from -k to len(ref)-len(needle)
   */
  bool find_substring(std::string ref, std::string needle,int& pos, int k = 0) {
    if (needle.size() > ref.size()+k)//If size of substring exceeds size of ref+number of missmatches
      return false;

    //Cover left overflow.
    /* For example:
       ref:     GATCGAAAG
       needle: AGATCG
     */
    for (int i=0; i<=k; ++i) {
      int mismatch = i;
      for (unsigned int j=0; j<needle.size()-i; ++j) {
	if (needle[i+j] != ref[j])
	  ++mismatch;
	if (mismatch > k)
	  break;
      }
      if (mismatch <= k) {
	pos = -i;
	return true;
      }
    }

    for (unsigned int i=1; i<ref.size()-needle.size()+k;++i) {//For all from second starting position in ref, since first was covered above!
      int mismatch = 0;
      for (unsigned int j=0; j<needle.size();++j) {
	if (needle[j] != ref[i+j])
	  ++mismatch;
	if (mismatch > k)
	  break;
      }
      if (mismatch <= k) {//Matched string with less than of equal K, return position!
	pos = (int)i;
	return true;
      }
    }
    return false;
  };

  /**
   * @brief Search for adapter in read and remove it if found.
   * @retrun True is adapter found and removed from read. False otherwise
   * @out read = clipped read
   */
  bool remove_adapter(char*& read, std::string adapter, int m = 0) {
    int pos = 0;
    //Get match position
    bool found = find_substring(read,adapter,pos,m);
    if (!found) //Not found, forget it!
      return false;

    if ((strlen(read)- adapter.size() ) > 0) {
      //Just remove the adapter from this sequence. NOT anything else.                                                                                                                                                          
      if ((int)strlen(read)-(int)adapter.size() <= pos) {//Adapter pases the end of the string. Simply clip at postion i.                                                                                           
	read[pos] = 0;
      } else {//Adapter is somewhere inside the string.                                                                                                                                                                         
	unsigned int aLen = adapter.size();
	if (pos < 0) {//Match oveshoots to the left. Cut less than length of adapter by i.                                                                                                                                        
	  aLen += pos;
	  pos = 0;
	}
	unsigned int limit = strlen(read)-aLen;
	for (unsigned int k=pos; k<limit; ++k) {
	  read[k] = read[k+aLen];
	}
	read[limit] = 0;
      }
    } else {
      if (pos < 0) pos=0;
      //Just trim the string after this position                                                                                                                                                                                
      read[pos] = 0;
    }
    return true;

  };

  /**
   * @brief Convert chromosome name from bwa standard to TopHat.
   * @brief Basically remove 'chr'. Except of chrM which should be translated to MT
   * @param chr -> Chromosome name
   * @param bwaToTopHat -> default true. Chop "chr". If false, make TopHat to bwa standard
   * @return Translated string
   */
  std::string translate_chr(const std::string chr, bool bwaToTopHat = true) {
    if (bwaToTopHat) {
      if (chr.compare("chrM") == 0)
	return "MT";
      else {
	std::string ret = chr;
	ret.replace(0,3,"");
	return ret;
      }
    } else {
      if (chr.compare("MT") == 0)
	return "chrM";
      else {
	std::string ret = "chr";
	ret += chr;
	return ret;
      }
    }
  };

#ifdef USE_SAM

  /**
   * Open a .sam/.bam file.
   * @returns NULL is open failed.
   */
  samfile_t * openAlignmentFile(std::string path)
  {
    samfile_t * fp = NULL;
    std::string flag = "r";
    if (path.substr(path.size()-3).compare("bam") == 0) {
      //BAM file!
      flag += "b";
    }
    if ((fp = samopen(path.c_str(), flag.c_str() , 0)) == 0) {
      fprintf(stderr, "Failed to open file %s\n", path.c_str());
    }
    return fp;
  };
#endif

  /**
   * @brief Proper tokanizer!
   */
  const char *toksplit(
		       const char *src, /* Source of tokens */
		       char tokchar, /* token delimiting char */
		       char *token, /* receiver of parsed token */
		       size_t lgh) /* length token can receive */
    /* not including final '\0' */
  {
    if (src) {
      while (' ' == *src) src++;
      while (*src && (tokchar != *src)) {
	if (lgh) {
	  *token++ = *src;
	  --lgh;
	}
	src++;
      }
      if (*src && (tokchar == *src)) src++;
    }
    *token = '\0';
    return src;
  }

}
