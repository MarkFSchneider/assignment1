#include "DNAStats.h"

DNAStats::DNAStats(string filename)
{
  //clear memory issues
  zeroVariables();

  //Read each line, and calculate it
  ifstream file(filename);
  string line;
  while(getline(file, line))
  {
    //checks if a line is empty before calculating
    if(line != "\n" || "")
    {
    updateNumbers(line);
    }
  }

  //updates the ratios.
  updatePercents();

  //calculate the variance. Can't be done during updateNumbers because it requires re-reading the file after knowing the mean
  calculateVariance(filename);

}

DNAStats::~DNAStats()
{

}

void DNAStats::zeroVariables()
{
  //helps clear memory errors
  aN = 0;
  tN  = 0;
  cN  = 0;
  gN  = 0;

  aaN  = 0;
  atN  = 0;
  acN  = 0;
  agN  = 0;

  taN  = 0;
  ttN  = 0;
  tcN  = 0;
  tgN  = 0;

  caN  = 0;
  ctN  = 0;
  ccN  = 0;
  cgN  = 0;

  gaN  = 0;
  gtN  = 0;
  gcN  = 0;
  ggN  = 0;

  aP = 0;
  tP  = 0;
  cP  = 0;
  gP  = 0;

  aaP  = 0;
  atP  = 0;
  acP  = 0;
  agP  = 0;

  taP  = 0;
  ttP  = 0;
  tcP  = 0;
  tgP  = 0;

  caP  = 0;
  ctP  = 0;
  ccP  = 0;
  cgP  = 0;

  gaP  = 0;
  gtP  = 0;
  gcP  = 0;
  ggP  = 0;

  sumDNALength  = 0;
  numberDNALines  = 0;
  sumSquaredDNALengths =0;
  meanDNALength = 0;

  varianceDNALength  = 0;

  accurateFile = false;
  bigramError = false;
  caseError = false;
  letterError = false;
}

int DNAStats::countLetter(string fullString, char chosenLetter)
{
  //counts the number of chosen letter chars in a string
  int count = 0;
  for(unsigned i=0; i<fullString.length(); ++i)
  {

    if(toupper(fullString.at(i))==chosenLetter)
    {
      count += 1;
    }

    if(!toupper(fullString.at(i)) == 'A' || 'T' || 'C' || 'G')
    {
      //checks if a letter is ATCG and flags an error if that doesn't happen (Unexpected letters mess with the bigram)
      letterError = true;
    }

    if(toupper(fullString.at(i))==chosenLetter && fullString.at(i)!= chosenLetter)
    {
      //checks for lowercase letters and flags them (lowercase letters mess with the bigram)
       caseError = true;
    }



  }
  return count;
}

int DNAStats::countBigram(string fullString, string chosenBigram)
{

  //reads a string and checks for a given two characters in a row, each two characters
  int count = 0;
  for(unsigned i=0; i<fullString.length(); i+=2)
  {
    if(fullString.substr(i,2)==chosenBigram){
    //if(("" + toupper(fullString.at(i) + "" + toupper(fullString.at(i+1)) == chosenBigram){
      count += 1;
    }
  }
  return count;
}

void DNAStats::updateNumbers(string fullString)
{
/*looks at a string and updates internal variables based on that string.
Ignores characters that arent ATCG, but doesn't break in case of an unexpected
character.

Variables include the count of each letter and bigram, the total DNA length,
number of lines, and the mean length of each line
*/
  sumDNALength += fullString.length();
  numberDNALines += 1;
  meanDNALength = (float)sumDNALength / numberDNALines;

  //sumSquaredDNALengths += thisLength * thisLength;

  //varianceDNALength
  //stdDeviationDNALength;

  aN += countLetter(fullString, 'A');
  tN += countLetter(fullString, 'T');
  cN += countLetter(fullString, 'C');
  gN += countLetter(fullString, 'G');

  aaN += countBigram(fullString, "AA");
  atN += countBigram(fullString, "AT");
  acN += countBigram(fullString, "AC");
  agN += countBigram(fullString, "AG");

  taN += countBigram(fullString, "TA");
  ttN += countBigram(fullString, "TT");
  tcN += countBigram(fullString, "TC");
  tgN += countBigram(fullString, "TG");

  caN += countBigram(fullString, "CA");
  ctN += countBigram(fullString, "CT");
  ccN += countBigram(fullString, "CC");
  cgN += countBigram(fullString, "CG");

  gaN += countBigram(fullString, "GA");
  gtN += countBigram(fullString, "GT");
  gcN += countBigram(fullString, "GC");
  ggN += countBigram(fullString, "GG");


}

float DNAStats::divide(int numerator, int denomenator)
{
  //divides 2 ints into a float. I'm not sure what I was thinking
  return (float)numerator / denomenator;
}

void DNAStats::updatePercents()
{
  aP = divide(aN, sumDNALength);
  tP = divide(tN, sumDNALength);
  cP = divide(cN, sumDNALength);
  gP = divide(gN, sumDNALength);

  aaP = divide(aaN, sumDNALength);
  atP = divide(atN, sumDNALength);
  acP = divide(acN, sumDNALength);
  agP = divide(agN, sumDNALength);

  taP = divide(taN, sumDNALength);
  ttP = divide(ttN, sumDNALength);
  tcP = divide(tcN, sumDNALength);
  tgP = divide(tgN, sumDNALength);

  caP = divide(caN, sumDNALength);
  ctP = divide(ctN, sumDNALength);
  ccP = divide(ccN, sumDNALength);
  cgP = divide(cgN, sumDNALength);

  gaP = divide(gaN, sumDNALength);
  gtP = divide(gtN, sumDNALength);
  gcP = divide(gcN, sumDNALength);
  ggP = divide(ggN, sumDNALength);
}

float DNAStats::calculateVariance(string filename)
{
  /*
  Takes the calculated mean and rereads the file.
  Each line is parsed for length and that length is subtracted from the mean.
  Those numbers are squared, and divided by the total number of lines
  this is the varience in the line length in the file
  */
  ifstream file(filename);
  string line;

  while(getline(file, line))
    {

      float thisLength = meanDNALength - line.length();
      sumSquaredDNALengths += thisLength * thisLength;
    }
   varianceDNALength = (float)sumSquaredDNALengths / numberDNALines;

   return varianceDNALength;
}

double DNAStats::calculateStandardGaussian()
{
  //Calculates a standard gaussian ratio.
  //does this version of G++ have the gaussian distribution random function?


  double uniformRandomA = rand() * (1.0 / RAND_MAX);
  double uniformRandomB = rand() * (1.0 / RAND_MAX);

  double pi = 3.141592653589793238463; //wasn't sure if MATH_PI exists in this c++ version

  double standardGaussian = sqrt(-2.0 * log(uniformRandomA)) * cos(2 * pi * uniformRandomB);
  return standardGaussian;
}

int DNAStats::newDNALength()
{
  //calculates a length based on a standard gaussian, the varience, and the mean DNA Length
  double doubleLength;
  doubleLength = meanDNALength + (sqrt(varianceDNALength) * calculateStandardGaussian());
  int intLength = round(doubleLength);
  if(intLength <= 0)
  {
    //makes it so the function can't return a 0 length line
    intLength = newDNALength();
  }
  return intLength;
}

string DNAStats::chooseNucleotide()
{
/*
Chooses a neucliotide based on the relative percents.
Plots the percents on a numberline from 0-99, proportionally split
A random number is chosen, and wherever it falls on that line, that number is selected.
Slight margin of error exists in accuracy in a percent accuracy greater than 2 numbers (.XX/1, or XX%)
*/


  float aPerc = aP * 100;
  float tPerc = tP * 100;
  float cPerc = cP * 100;
  float gPerc = gP * 100;

  tPerc += aPerc;
  cPerc += tPerc;
  gPerc += cPerc;

  float choice = rand() % 100;
  if(choice < aPerc)
  {
    return "A";
  }

  if(aPerc <= choice && choice < tPerc)
  {
    return "T";
  }

  if(tPerc <= choice && choice < cPerc)
  {
    return "C";
  }

  if(cPerc <= choice && choice < gPerc)
  {
    return "G";
  }

  return "G";
}

string DNAStats::createDNA()
{
  /*
  Chooses a Length, and fills a string with that many Neucleotides
  */
  int length = newDNALength();
  string newDNA = "";

  for(int i = 0; i < length; i++)
  {
    newDNA.append(chooseNucleotide());
  }
  return newDNA;
}
