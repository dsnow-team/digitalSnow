#include "DGtal/topology/helpers/SimplePointHelper.h"
#include <iostream>
#include <exception>


///////////////// helpers
template<DGtal::Dimension dimension, typename TCoordinate>
typename DGtal::PointVector<dimension, TCoordinate> makePoint(const TCoordinate& aValue) 
{
  typedef typename DGtal::PointVector<dimension, TCoordinate> Point; 
  Point p; 
  for (typename Point::Iterator i = p.begin(); i != p.end(); ++i)
  {
    *i = aValue; 
  }
  return p;
} 


///////////////// main functions

bool basicTest(std::string& config, std::string& c)
{
  std::string::iterator itb = config.begin(); 
  std::string::iterator ite = config.end(); 
  int size = std::ceil(std::pow( (ite-itb), ((double) 1 / (double) 3)));
  if ( (ite-itb) != (size*size*size) )
  {
    std::cerr << "Error. Bad input (" << (ite-itb) << " not equal to "; 
    std::cerr << size << " ^3 = " << (size*size*size) << ")" << std::endl;
    throw std::exception(); 
  }  
  int size2 = size/2; 

  typedef std::iterator_traits<std::string::iterator>::value_type Label; 
  typedef DGtal::ImageContainerBySTLVector<HyperRectDomain<SpaceND<3, int> >, Label> Image; 
  Image img( makePoint<3,int>(-size2), makePoint<3,int>(size2) ); 
  DGtal::SimplePointHelper<Image >::readConfiguration(img, itb, ite); 
  DGtal::SimplePointHelper<Image > h( img ); 

  return ( h.isMLSimple(makePoint<3,int>(0), *(c.begin())) ); 
}

bool read(std::istream& in, std::string& config, std::string& c)
{
  std::string str;
  bool flag = getline( in, str );
  if ( ! in.good() ) return false;
  if ( ( str.size() > 0 ) && ( str.at(0) != '#' ) )
  {
    std::istringstream str_in( str );
    str_in >> config >> c;
  }

  return flag; 
}

int main(int argc, char** argv)
{
  
  trace.beginBlock ( "Testing class SimplePointHelper" );
  trace.info() << "# Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << std::endl;

  unsigned int counter = 0; 

  //read config in the standard input
  std::string config = ""; 
  std::string c = ""; 
  
  while ( read(std::cin, config, c) )
  { 

    //test
    unsigned int res = (basicTest(config,c))?1:0; 
    std::cout << res << std::endl; 
    counter += res; 

  }

  std::cerr << "# " << counter << std::endl; 
  
  trace.endBlock();
  return 1;

}

