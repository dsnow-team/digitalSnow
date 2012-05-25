#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>

#include "DGtal/topology/helpers/SimplePointHelper.h"

///////////////// helpers


template<DGtal::Dimension dimension, typename TCoordinate>
typename DGtal::PointVector<dimension, TCoordinate> 
makePoint(const TCoordinate& aValue) 
{
  typedef typename DGtal::PointVector<dimension, TCoordinate> Point; 
  Point p; 
  for (typename Point::Iterator i = p.begin(); i != p.end(); ++i)
    {
      *i = aValue; 
    }
  return p;
} 


template<typename TImage>
void labelsImageToString(const TImage& aImg, std::string& config)
{

  std::ostringstream oss;

  typename TImage::Domain d = aImg.domain(); 
  typename TImage::Domain::ConstIterator it = d.begin(); 
  typename TImage::Domain::ConstIterator itEnd = d.end();
  for ( ; it != itEnd; ++it)
    { 
      oss << aImg(*it);
    }
  config = oss.str();
}
///////////////// main functions



bool writeValidRandomConfig(std::ostream& out, const unsigned int& s, const unsigned int& n, const double& proba)
{

  //generate distinct labels as char
  std::set<char> labels; 
  for (unsigned int i = 0; i < n; ++i)
    {
      labels.insert( 97 + i ); 
    }

  //image creation 
  typedef char Label; 
  typedef SpaceND<3,int > Space; 
  typedef DGtal::ImageContainerBySTLVector<HyperRectDomain<Space >, Label >  LabelsImage; 
  int bound = s/2; 
  LabelsImage img( makePoint<3,int>(-bound), makePoint<3,int>(bound) );
  //generate a valid config
  if (SimplePointHelper<LabelsImage>::generateRandomConfiguration(img, labels, proba) )
    { 
      //write image as a string
      std::string config; 
      labelsImageToString(img, config); 

      //center goes into region of label l :
      std::string l = "b"; 

      out << config << " " << l << std::endl;
      return true; 
    }
  else 
    {
      return false; 
    }
}


/////////////////////////////////// main

int main(int argc, char** argv)
{

  std::srand ( std::time(NULL) );

  std::cerr << "# Args:";
  for ( int i = 0; i < argc; ++i )
    std::cerr << " " << argv[ i ];
  std::cerr << std::endl;

  //////////////////////////////////////////////////////parameters
  if (argc <= 3)
    {
      std::cerr << "# Usage: " << argv[0]; 
      std::cerr << " cubeSize distinctLabelsNumber configurationsNumber " << std::endl;  
    }
  else 
    {
      //nb of expected configurations 
      unsigned int nbC = 10; 
      std::istringstream sC( argv[3] ); 
      sC >> nbC;  
      //nb of distinct labels (should be less than 26, but greater than 2)
      unsigned int nbL = 2; 
      std::istringstream sL( argv[2] ); 
      sL >> nbL;
      if ( (nbL > 26) || (nbL < 2) ) 
	{
	  std::cerr << "Error. bad number of distinct labels" << std::endl; 
	  return 0; 
	}
      //size of the domain (odd number)
      unsigned int size = 1; 
      std::istringstream sS( argv[1] ); 
      sS >> size;
      size = size*2+1; 
  
      std::cerr << "# " << nbC << " cubes of length " << size; 
      std::cerr << " with " << nbL << " distinct labels" << std::endl; 

      ///////////////////////////////////////////////////writing configurations
      unsigned int i = 0; 
      unsigned int c = 0; 
      while(i < nbC)
	{
	  if (writeValidRandomConfig(std::cout, size, nbL, 0.3))
	    c++; 
	  i++;
	}
      std::cerr << "# " << c << " generated configurations" << std::endl; 

    }

  return 1;

}

