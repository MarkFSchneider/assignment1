#include "DNAStats.h"
#include <iostream>
#include <string>
#include <time.h>
using namespace std;


/*
Skip invalid lines
Bigrams are 2 only. (paired exclucively (eg 'ATCG' has AT and CG, but not TC))
assume nucleotides are even

*/
string filename = "defaultDNAin.txt";

bool static programExit()
{
  char answer;
  cout << "Process another list? Y/N" << endl;
  cin >> answer;

  if(toupper(answer) == 'Y')
  {
    cout << "Enter the path of the new list" << endl;
    cin >> filename;
    return false;
  }

  if(toupper(answer) == 'N')
  {
    return true;
  }

  return programExit();


}


int main(int argc, char const *argv[]) {
if(argc > 1)
{
  filename = argv[1];
  ifstream infile2(filename);
  if(!infile2.good())
  {
    cout << "file not found. Program Terminated";
  }

}
else
{
  cout << "Invalid source file. Program Terminated";
}


srand(time(0));
bool exitProgram = false;

string outputFile = "schneider.out";


while(!exitProgram){

DNAStats stats = DNAStats(filename);


ofstream outFile;
outFile.open(outputFile);





//I was having trouble creating strings with variables
//I'm not proud of what follows
if(stats.caseError)
{
  outFile << "WARNING: A lowercase letter was found. Bigram counts may be inaccurate" << endl;
  cout << "WARNING: A lowercase letter was found. Bigram counts may be inaccurate" << endl;
}

if(stats.letterError)
{
  outFile << "WARNING: A non-ATCG character was found. Bigram counts may be inaccurate" << endl;
  cout << "WARNING: A non-ATCG character was found. Bigram counts may be inaccurate" << endl;
}

outFile << "sum DNA Length: ," << stats.sumDNALength << endl
<< "DNA Length Varience: " << stats.varianceDNALength << endl
<< "DNA Length Standard Deviation: " << sqrt(stats.varianceDNALength) << endl
<< endl
<< "DNA Neucleotide Probability:" << endl
<< "A: " << stats.aP * 100 << endl
<< "T: " << stats.tP * 100 << endl
<< "C: " << stats.cP * 100 << endl
<< "G: " << stats.gP * 100 << endl
<< endl
<< "DNA Bigram Probability:" << endl
<< "AA: " << stats.aaP * 100 << endl
<< "AT: " << stats.atP * 100 << endl
<< "AC: " << stats.acP * 100 << endl
<< "AG: " << stats.agP * 100 << endl
<< endl
<< "TA: " << stats.taP * 100 << endl
<< "TT: " << stats.ttP * 100 << endl
<< "TC: " << stats.tcP * 100 << endl
<< "TG: " << stats.tgP * 100 << endl
<< endl
<< "CA: " << stats.caP * 100 << endl
<< "CT: " << stats.ctP * 100 << endl
<< "CC: " << stats.ccP * 100 << endl
<< "CG: " << stats.cgP * 100 << endl
<< endl
<< "GA: " << stats.gaP * 100 << endl
<< "GT: " << stats.gtP * 100 << endl
<< "GC: " << stats.gcP * 100 << endl
<< "GG: " << stats.ggP * 100 << endl
<< endl;

for(int i = 0; i<1000; i++)
{
  outFile << stats.createDNA() << endl;
}

exitProgram = programExit();


ifstream infile(filename);
if(!infile.good())
{
  cout << "file not found. Program Exit";
}

}

  return 0;





}
