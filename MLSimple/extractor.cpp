#include <iostream>
#include <sstream>
#include <fstream>
#include <string>


int main(int argc, char** argv)
{
  
  std::cerr << "# Args:";
  for ( int i = 0; i < argc; ++i )
    std::cerr << " " << argv[ i ];
  std::cerr << std::endl;

  if (argc < 4) 
    {
      std::cerr << "# Usage: " << argv[0] << " configsFile resultsFile value" << std::endl;
      std::cerr << "# write to the output all lines of configsFile whose corresponding lines in resultsFile match to value" << std::endl;
      return 0;  
    }

  //////////////////////////////////////////////////////parameters
  //files
  std::ifstream f1( argv[1] ); 
  std::ifstream f2( argv[2] ); 

  std::string value( argv[3] ); 
  ///////////////////////////////////////////////////read and compare

  if ( f1 && f2 ) 
  {
    std::string l1, l2; // current lines of file1 and file2

    int c = 0; 
    while ( std::getline( f1, l1 ) && std::getline( f2, l2 ) )
    {
      if (l2.compare(value) == 0) //equal 
      {
        std::cout << l1 << std::endl;
        ++c; 
      }
    }
    std::cerr << "# " << c << " lines " << std::endl; 
  } 
  else 
  {
    std::cerr << "# Error. files not open" << std::endl; 
    return 0; 
  } 

  return 1;
}

