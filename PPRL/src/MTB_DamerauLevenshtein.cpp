/*
 * damerauLevenshtein.cpp
 *
 *  Created on: 21.02.2017
 *      Author: schnell-42
 */

# include "MTB_DamerauLevenshtein.h"

int mini(int a, int b, int c){
  return(min(a, min(b,c)));
}

// Damerau-Levenshtein Edit Distance
int damlevdist(string word1, string word2){
  /// Damerau-Levenshtein distance
  ///  Please use lower-case strings
  /// word1 : first word
  /// word2 : second word
  ///
  int size1 = word1.size(), size2 = word2.size();
  int suppr_dist, insert_dist, subs_dist, val;
  int* dist = new int[(size1+1)*(size2+1)];

  for(int i=0; i<size1+1; ++i)
    dist[(size2+1)*i] = i;
  for(int j=0; j<size2+1; ++j)
    dist[j] = j;
  for(int i=1; i<size1+1; ++i){
    for(int j=1; j<size2+1; ++j){
      suppr_dist = dist[(size2+1)*(i-1)+j] + 1;
      insert_dist = dist[(size2+1)*i+j-1] + 1;
      subs_dist = dist[(size2+1)*(i-1)+j-1];
      if(word1[i-1]!=word2[j-1])  // word indexes are implemented differently.
        subs_dist += 1;
      val = mini(suppr_dist, insert_dist, subs_dist);
      if(((i>=2) && (j>=2)) && ((word1[i-1]==word2[j-2]) && (word1[i-2]==word2[j-1])))
        val = min(dist[(size2+1)*(i-2)+j-2]+1, val);
      dist[(size2+1)*i+j] = val;
    }
  }

  int res = dist[(size1+1)*(size2+1) - 1];
  delete[] dist;
  return(res);
}

string MTB_DamerauLevenshteinAlgorithm::getName(){
  return "DamerauLevenshtein";
}

double MTB_DamerauLevenshteinAlgorithm::getRelativeValue(string o1, string o2){
  return calcFromDistanceToInterval(o1, o2, damlevdist(o1,o2));
}

double MTB_DamerauLevenshteinAlgorithm::getAbsoluteValue(string o1, string o2){
  return damlevdist(o1,o2);
}
