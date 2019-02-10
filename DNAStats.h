#include <string>
#include <iostream>
#include <random>
#include <fstream>
#include <math.h>
#include <locale>
using namespace std;

class DNAStats
{
public:
  DNAStats(string filename);
  ~DNAStats();

  int countLetter(string fullString, char chosenLetter);
  int countBigram(string fullString, string chosenBigram);
  void zeroVariables();
  void updateNumbers(string fullString);
  float divide(int numerator, int denomenator);
  void updatePercents();

  double calculateStandardGaussian();
  float calculateVariance(string filename);

  int newDNALength();
  string chooseNucleotide();
  string createDNA();



  int aN;
  int tN;
  int cN;
  int gN;

  int aaN;
  int atN;
  int acN;
  int agN;

  int taN;
  int ttN;
  int tcN;
  int tgN;

  int caN;
  int ctN;
  int ccN;
  int cgN;

  int gaN;
  int gtN;
  int gcN;
  int ggN;

  float aP;
  float tP;
  float cP;
  float gP;

  float aaP;
  float atP;
  float acP;
  float agP;

  float taP;
  float ttP;
  float tcP;
  float tgP;

  float caP;
  float ctP;
  float ccP;
  float cgP;

  float gaP;
  float gtP;
  float gcP;
  float ggP;

  int sumDNALength;
  int numberDNALines;
  float sumSquaredDNALengths;
  float meanDNALength;


  float varianceDNALength;
  //int stdDeviationDNALength;
  //Variance is The average of the squared differences from the Mean.
  //Take each difference, square it, and average the results
  //Standard deviation is the root of variance




};
