#include "dnastats.h"
#include <iostream>
#include <string>
#include <time.h>
using namespace std;


/*
Skip invalid lines
Bigrams are 2 only. (paired exclucively (eg 'ATCG' has AT and CG, but not TC))
assume nucleotides are even

*/
int main(int argc, char const *argv[]) {

string filename = "ZDNATEXT.txt";

srand(time(0));
DNAStats stats = DNAStats(filename);


cout << "3 sample DNA" << endl;
cout << "1: " << endl << stats.createDNA() << endl;
cout << "2: " << endl << stats.createDNA() << endl;
cout << "3: " << endl << stats.createDNA() << endl;
  return 0;
}
