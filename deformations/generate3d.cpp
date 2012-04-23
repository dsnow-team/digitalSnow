#include <sstream>
#include <iomanip>

/////////////////////
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

namespace po = boost::program_options;

/////////////////////
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>

using namespace Z3i; 


/////////////////////////// useful functions
#include "deformationFunctions.h"
#include "deformationDisplay3d.h"
#include "DGtal/io/readers/VolReader.h"

///////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  
  DGtal::trace.info() << "Generating 3d partitions using DGtal ";
  DGtal::trace.emphase() << "(version "<< DGTAL_VERSION << ")"<< std::endl;
  


  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("domainSize,d",  po::value<int>()->default_value(64), "Domain size (if default starting interface)" )
    ("shape,s", po::value<string>()->default_value("ball"), 
     "Generated shape: either <ball> (default) or <2balls> " )
    ("outputFiles,o",   po::value<string>()->default_value("interface"), "Output files basename" );

  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  po::notify(vm);    
  if(vm.count("help")||argc<=1)
    {
      trace.info()<< "Generating 3d partitions" << std::endl
		  << "Basic usage: "<<std::endl
		  << argv[0] << " [other options] -d <size> -s <shape> " 
		  << std::endl
		  << general_opt << "\n";
      return 0;
    }
  
  //Parse options
  //domain size
  int dsize; 
  if (!(vm.count("domainSize"))) trace.info() << "Domain size default value: 32" << std::endl; 
  dsize = vm["domainSize"].as<int>(); 

  //files
  std::string outputFiles; 
  if (!(vm.count("outputFiles"))) 
    trace.info() << "output files begin with : interface" << std::endl;
  outputFiles = vm["outputFiles"].as<std::string>();

  //image and implicit function
  typedef ImageContainerBySTLVector<Domain,short int> LabelImage;
  Point p = Point::diagonal(0);
  Point q = Point::diagonal(dsize); 
  Domain d = Domain(p,q); 
  LabelImage labelImage( Domain(p,q) ); 

  DGtal::trace.beginBlock("image reading..."); 
  if (vm.count("shape"))
    {
      if ( (vm["shape"].as<std::string>()) == "ball" )
        {
          Point c = Point::diagonal(dsize/2); 
          double r = (((double)dsize*3.0/5.0)/2.0); 
  Domain::ConstIterator cIt = d.begin(); 
  Domain::ConstIterator cItEnd = d.end(); 
  for ( ; cIt != cItEnd; ++cIt)
  { //for each domain point
    double dist = (*cIt-c).norm(Point::L_2) - r; 
    LabelImage::Value v = ( dist <= 0 )?255:0; 
    labelImage.setValue(*cIt, v);  
  }

        }

      if ( (vm["shape"].as<std::string>()) == "2balls" )
        {
          Point c1 = Point::diagonal(dsize/3 + 1); 
          Point c2 = Point::diagonal(dsize/3*2 - 1); 
          double r = (dsize/3); 
  Domain::ConstIterator cIt = d.begin(); 
  Domain::ConstIterator cItEnd = d.end(); 
  for ( ; cIt != cItEnd; ++cIt)
  { //for each domain point
    double dist1 = (*cIt-c1).norm(Point::L_2) - r; 
    double dist2 = (*cIt-c2).norm(Point::L_2) - r; 
    LabelImage::Value v = ( dist1 <= 0 )?255:0; 
    if (dist2 <= 0) v = 127; 
    labelImage.setValue(*cIt, v);  
  }
        }

/*      if ( (vm["shape"].as<std::string>()) == "pyramid" )
        {
          Point c = Point::diagonal(dsize/2); 
          initWithBallPredicate( labelImage, c, (dsize*3/5)/2 );
          inv( labelImage ); 
        }
*/

    } 
  else 
    trace.info() << "shape not known, use option -h" << std::endl; 

  trace.info() << "starting interface initialized with a " 
   << (vm["shape"].as<std::string>()) << std::endl;
  DGtal::trace.endBlock(); 

  //write into a vol file
  std::stringstream s; 
  s << outputFiles << ".vol";
  typedef GradientColorMap<LabelImage::Value, DGtal::CMAP_GRAYSCALE> ColorMap; 
  VolWriter<LabelImage,ColorMap>::exportVol( s.str(), labelImage, 0, 255 );

  return 1;
}

