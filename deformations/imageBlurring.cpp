/////////////////////
#include <iostream>
#include <fstream>
#include <sstream>

#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>

/////////////////////
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

namespace po = boost::program_options;

/////////////////////
using namespace std; 
using namespace DGtal; 
using namespace Z2i;

#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/PGMWriter.h"


#include "ExactDiffusionEvolver.h"
#include "WeickertKuhneEvolver.h"
#include "GrayscaleMapCast.h"

template< typename TImage >
void write ( const TImage& img, std::string basename, int min = 0, int max = 255 ) 
{

    typedef GrayscaleMapCast<typename TImage::Value> ColorMap;
    std::stringstream s; 
    s << basename << ".pgm"; 
    PGMWriter<TImage,ColorMap>::exportPGM( s.str(), img, ColorMap(min, max), true );
} 



////////////////////////////////////////
int main(int argc, char** argv)
{

  DGtal::trace.info() << "image blurring using DGtal ";
  DGtal::trace.emphase() << "(version "<< DGTAL_VERSION << ")"<< std::endl; 

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("inputImage,i",  po::value<string>(), "gray level image (.pgm)" )
    ("algo,a",  po::value<string>()->default_value("exact"), 
"can be: \n 'exact'  \n or 'weickert' " )
    ("sigma,s",  po::value<double>()->default_value(2), 
"controls the amount of blurring \n (between 0 and 50) " )
    ("outputFile,o",   po::value<string>()->default_value("output"), "Output file basename" );

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  po::notify(vm);    
  if(vm.count("help")||argc<=1)
    {
      trace.info()<< "Image-blurring" << std::endl
      << "Basic usage: "<<std::endl
      << argv[0] << " -a exact -s 2 " << std::endl
      << general_opt << "\n";
      return 0;
    }
  
  //Parse options
  //input
  std::string inputFilename; 
  if (!(vm.count("inputImage"))) 
    {
      trace.error() << "an input image is required" << std::endl;
      return 0; 
    }
  inputFilename = vm["inputImage"].as<std::string>();

  //algo
  std::string algo; 
  if (!(vm.count("algo"))) trace.info() << "default algorithm: exact" << std::endl; 
  algo = vm["algo"].as<string>(); 

  //sigma
  double sigma; 
  if (!(vm.count("sigma"))) trace.info() << "sigma default value: 2" << std::endl; 
  sigma = vm["sigma"].as<double>(); 

  //output file
  std::string outputFilename; 
  if (!(vm.count("outputFile"))) trace.info() << "output file nammed 'output'" << std::endl;
  outputFilename = vm["outputFile"].as<std::string>();

  /////////////////////////////////////////////////////////////////////////////////////////
  //reading image
  typedef ImageContainerBySTLVector<Domain,unsigned char> GrayImage; 
  GrayImage img = PGMReader<GrayImage>::importPGM( inputFilename ); 

  //Diffusion
  trace.info() << "sigma: " << sigma << std::endl; 

  if (algo.compare("exact")==0)
  {
    ExactDiffusionEvolver<GrayImage> e; 
    e.update(img,sigma);

  } else if (algo.compare("weickert")==0)
  {

    double tstep = (sigma*sigma)/2.0; 
    trace.info() << "time step: " << tstep << std::endl;

    //data functions
    ImageContainerBySTLVector<Domain,double> a(img.domain()); 
    std::fill(a.begin(),a.end(), 1.0 );  
    ImageContainerBySTLVector<Domain,double> b(img.domain()); 
    std::fill(b.begin(),b.end(), 1.0 );  
    ImageContainerBySTLVector<Domain,double> g(img.domain()); 
    std::fill(g.begin(),g.end(), 1.0 );  

    // pb of types unsigned char / double
    ImageContainerBySTLVector<Domain,double> img2(img.domain());
    std::copy(img.begin(), img.end(), img2.begin());
 
    //evolver
    WeickertKuhneEvolver<ImageContainerBySTLVector<Domain,double> > e(a,b,g,0,1);
    double tstepMax = 5.0;
    double q = std::floor( tstep / tstepMax );
    double r = tstep - q*tstepMax;  
    int k = (int) q; 
    for (unsigned int i = 0; i < k ; ++i ) 
      {
        trace.info() << "iteration #" << (i+1) << std::endl;
        e.update(img2, tstepMax);
      }
    if (r != 0) 
      {
        trace.info() << "iteration #" << (k+1) << std::endl;
        e.update(img2, r);
      }
    // pb of types unsigned char / double
    std::copy(img2.begin(), img2.end(), img.begin());

  } else trace.error() << "unknown algo. Try option -h to see the available algorithms " << std::endl;

  write(img, outputFilename); 


	return 1;
}
