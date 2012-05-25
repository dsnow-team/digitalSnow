#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

void comparison(std::string s, std::string t, unsigned int& n, unsigned int& n1, unsigned int& n2)
{
  if ( ( (s.compare("0")==0) && (t.compare("0")==0) )
       || ( (s.compare("1")==0) && (t.compare("1")==0) ) )
    {
      ++n; 
    }
  else
    {
      if ( (s.compare("0")==0) && (t.compare("1")==0) )
	{
	  ++n2;
	}
      else if ( (s.compare("1")==0) && (t.compare("0")==0) )
	{
	  ++n1; 
	}
      else
	{
	  std::cerr << "# Error. not equal to 0 or 1" << std::endl; 
	}
    }
}

int main(int argc, char** argv)
{
  
  std::cerr << "# Args:";
  for ( int i = 0; i < argc; ++i )
    std::cerr << " " << argv[ i ];
  std::cerr << std::endl;

  if (argc < 3) 
    {
      std::cerr << "# Usage: " << argv[0] << " resultsFile1 resultsFile2 configsFile" << std::endl;
      std::cerr << "# Compare line per line resultsFile1 and resultsFile2 with respect to configsFile" << std::endl;
      return 0;  
    }

  //////////////////////////////////////////////////////parameters
  //files
  std::ifstream f1( argv[1] ); 
  std::ifstream f2( argv[2] ); 
  std::ifstream f3( argv[3] ); 

  ///////////////////////////////////////////////////read and compare
  unsigned int nEq = 0; //number of lines that match
  unsigned int nDif1 = 0; //number of lines for which f1 contains 1 but f2 0
  unsigned int nDif2 = 0; //number of lines for which f2 contains 1 but f1 0

  if ( f1 && f2 && f3) 
  {
    std::string l1, l2, l3; // current lines of file1 and file2 and file3

    while ( std::getline( f1, l1 ) && std::getline( f2, l2 ) )
    {
      unsigned int nDifs = nDif1; 
      comparison(l1, l2, nEq, nDif1, nDif2);
      std::getline( f3, l3); 
      if (nDif1 != nDifs) std::cout << l3 << std::endl;
    }
  } 
  else 
  {
    std::cerr << "# Error. files not open" << std::endl; 
    return 0; 
  } 

  unsigned int n = nEq + nDif1 + nDif2; 
  std::cerr << "# equa: " << nEq  << " / " << n << std::endl; 
  std::cerr << "# diff: " << (nDif1+nDif2) << " / " << n << std::endl;
  std::cerr << "# true in " << argv[1] << " but false in " << argv[2] << ": " << nDif1 << std::endl; 
  std::cerr << "# true in " << argv[2] << " but false in " << argv[1] << ": " << nDif2 << std::endl; 

  return 1;
}

