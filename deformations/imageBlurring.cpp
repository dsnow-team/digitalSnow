/////////////////////
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>

/////////////////////
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

namespace po = boost::program_options;

/////////////////////

#include "DGtal/io/readers/PNMReader.h"
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/PNMWriter.h"

using namespace Z2i;

#include "ExactDiffusionEvolver.h"

template< typename TImage >
void write ( const TImage& img, std::string basename ) 
{

    typedef GradientColorMap<typename TImage::Value, DGtal::CMAP_GRAYSCALE> ColorMap; 

    std::stringstream s; 
    s << basename << ".pgm"; 
    PNMWriter<TImage,ColorMap>::exportPGM( s.str(), img, 0, 255, true );
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
    ("sigma,s",  po::value<double>()->default_value(0.01), "controls the amount of blurring (between 0.001 and 0.1)" )
    ("outputFile,o",   po::value<string>()->default_value("output"), "Output file basename" );

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  po::notify(vm);    
  if(vm.count("help")||argc<=1)
    {
      trace.info()<< "Image-blurring" << std::endl
      << "Basic usage: "<<std::endl
      << argv[0] << " -s 0.02" << std::endl
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

  //sigma
  double sigma; 
  if (!(vm.count("sigma"))) trace.info() << "sigma default value: 0.01" << std::endl; 
  sigma = vm["sigma"].as<double>(); 

  //output file
  std::string outputFilename; 
  if (!(vm.count("outputFile"))) trace.info() << "output file nammed 'output'" << std::endl;
  outputFilename = vm["outputFile"].as<std::string>();

  /////////////////////////////////////////////////////////////////////////////////////////
  //reading image
  typedef ImageContainerBySTLVector<Domain,unsigned char> GrayImage; 
  GrayImage img = PNMReader<GrayImage>::importPGMImage( inputFilename ); 

  //Diffusion
  trace.info() << "sigma: " << sigma << std::endl; 
  double tstep = (sigma*sigma)/2.0; 
  trace.info() << "time step: " << tstep << std::endl; 
  ExactDiffusionEvolver<GrayImage> e; 
  e.update(img,tstep);

  //writting image
  write(img, outputFilename); 

	return 1;
}
