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
#include "DGtal/io/readers/VolReader.h"

using namespace Z3i; 

//evolvers
//level set
#include "WeickertKuhneEvolver.h"

//phase field
#include "ExactDiffusionEvolver.h"
#include "ExactReactionEvolver.h"
#include "ExplicitReactionEvolver.h"
#include "LieSplittingEvolver.h"

/////////////////////////// useful functions
#include "deformationFunctions.h"
#include "deformationDisplay3d.h"


///////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  
  DGtal::trace.info() << "3d interface evolution using DGtal ";
  DGtal::trace.emphase() << "(version "<< DGTAL_VERSION << ")"<< std::endl;
  


  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("inputImage,i",  po::value<string>(), "Binary image to initialize the starting interface (vol format)" )
    ("domainSize,s",  po::value<int>()->default_value(32), "Domain size (if default starting interface)" )
    ("timeStep,t",  po::value<double>()->default_value(1.0), "Time step for the evolution" )
    ("displayStep,d",  po::value<int>()->default_value(1), "Number of time steps between 2 drawings" )
    ("stepsNumber,n",  po::value<int>()->default_value(1), "Maximal number of steps" )
    ("algo,a",  po::value<string>()->default_value("levelSet"), 
"can be: \n <levelSet>  \n or <phaseField> " )
    ("balloonForce,k",  po::value<double>()->default_value(0.0), "Balloon force" )
    ("epsilon,e",  po::value<double>()->default_value(3.0), "Interface width (only for phase fields)" )
    ("outputFiles,o",   po::value<string>()->default_value("interface"), "Output files basename" )
    ("outputFormat,f",   po::value<string>()->default_value("png"), 
"Output files format: either <png> (3d to 2d, default) or <vol> (3d)" )
    ("withVisu", "Enables interactive 3d visualization before and after evolution" );

  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  po::notify(vm);    
  if(vm.count("help")||argc<=1)
    {
      trace.info()<< "Evolution of a 3d interface" << std::endl
      << "Basic usage: "<<std::endl
      << argv[0] << " [other options] -t <time step> --withVisu" << std::endl
      << general_opt << "\n";
      return 0;
    }
  
  //Parse options
  //domain size
  int dsize; 
  if (!(vm.count("domainSize"))) trace.info() << "Domain size default value: 32" << std::endl; 
  dsize = vm["domainSize"].as<int>(); 

  //time step
  double tstep; 
  if (!(vm.count("timeStep"))) trace.info() << "time step default value: 1.0" << std::endl; 
  tstep = vm["timeStep"].as<double>(); 
    
  //iterations
  int step; 
  if (!(vm.count("displayStep"))) trace.info() << "number of steps between two drawings: 1 by default" << std::endl; 
  step = vm["displayStep"].as<int>(); 
  int max; 
  if (!(vm.count("stepsNumber"))) trace.info() << "maximal number of steps: 1 by default" << std::endl; 
  max = vm["stepsNumber"].as<int>(); 


  //files
  std::string outputFiles; 
  if (!(vm.count("outputFiles"))) 
    trace.info() << "output files beginning with : interface" << std::endl;
  outputFiles = vm["outputFiles"].as<std::string>();

  //files format
  std::string format; 
  if (!(vm.count("outputFormat"))) 
    trace.info() << "output files format is png (3d to 2d) " << std::endl;
  format = vm["outputFormat"].as<std::string>();
  if ((format != "png")&&(format != "vol")) 
    {
    trace.info() << "format is expected to be either png or vol " << std::endl;
    return 0; 
    }


  //image and implicit function
  Point p(0,0,0);
  Point q(dsize,dsize,dsize); 
  Point c(dsize/2,dsize/2,dsize/2); 
  ImageContainerBySTLVector<Domain,double> implicitFunction( Domain(p,q) ); 
  //initWithBall( implicitFunction, c, (dsize*3/5)/2);
  initWithFlower( implicitFunction, c, (dsize*3/5)/2, (dsize*1/5)/2, 5 ); 

  if (!(vm.count("inputImage"))) 
    trace.info() << "starting interface initialized with a flower shape" << std::endl;
  else
    { 
      string imageFileName = vm["inputImage"].as<std::string>();
      trace.emphase() << imageFileName <<std::endl; 
      typedef ImageContainerBySTLVector<Domain,unsigned char> BinaryImage; 
      BinaryImage img = VolReader<BinaryImage>::importVol( imageFileName);
      Domain d = img.domain(); 
      p = d.lowerBound(); q = d.upperBound(); 
      implicitFunction = ImageContainerBySTLVector<Domain,double>( Domain(p,q) ); 
      initWithDT( img, implicitFunction );

    }

  //3d to 2d display
  std::stringstream ss; 
  ss << outputFiles << "0001"; 
  writeImage( implicitFunction, ss.str(), format );

  //interactive display before the evolution
  if (vm.count("withVisu")) displayImage( argc, argv, implicitFunction ); 

  //balloon force
  double k; 
  if (!(vm.count("balloonForce"))) trace.info() << "balloon force default value: 0" << std::endl; 
  k = vm["balloonForce"].as<double>(); 

  //algo
  std::string algo; 
  if (!(vm.count("algo"))) trace.info() << "default algorithm: levelSet" << std::endl; 
  algo = vm["algo"].as<string>(); 

  if (algo.compare("levelSet")==0)
  {

    //data functions
    ImageContainerBySTLVector<Domain,double> a( Domain(p,q) ); 
    std::fill(a.begin(),a.end(), 1.0 );  
    ImageContainerBySTLVector<Domain,double> b( Domain(p,q) ); 
    std::fill(b.begin(),b.end(), 1.0 );  
    ImageContainerBySTLVector<Domain,double> g( Domain(p,q) ); 
    std::fill(g.begin(),g.end(), 1.0 );  

    //interface evolver
    WeickertKuhneEvolver<ImageContainerBySTLVector<Domain,double> > e(a,b,g,k,1); 

    for (unsigned int i = 1; i <= max; ++i) 
    {
      std::stringstream s0; 
      s0 << "iteration # " << i; 
      DGtal::trace.beginBlock( s0.str() );

      //update
      e.update(implicitFunction,tstep); 

      if ((i%step)==0) 
      {

         //3d to 2d display
         std::stringstream s; 
         s << outputFiles << setfill('0') << std::setw(4) << (i/step)+1; 
         writeImage( implicitFunction, s.str(), format );

      }
      DGtal::trace.endBlock();   
    }

    //interactive display after the evolution
    if (vm.count("withVisu")) displayImage( argc, argv, implicitFunction ); 

  } else if (algo.compare("phaseField")==0)
  {

    double epsilon = 3.0; 
    if (!(vm.count("epsilon"))) trace.info() << "epsilon default value: 3.0" << std::endl; 
    epsilon = vm["epsilon"].as<double>(); 
    if (epsilon <= 0) 
      {
        trace.error() << "epsilon should be greater than 0" << std::endl;
        return 0; 
      } 

    //computing the profile from the signed distance
    Profile p(epsilon); 
    std::transform(implicitFunction.begin(), implicitFunction.end(), implicitFunction.begin(), p); 

    ImageContainerBySTLVector<Domain,double> a( Domain( implicitFunction.domain() ) ); 
    std::fill(a.begin(),a.end(), 1.0 );  

    typedef ExactDiffusionEvolver<ImageContainerBySTLVector<Domain,double> > Diffusion; 
    typedef ExplicitReactionEvolver<ImageContainerBySTLVector<Domain,double>, 
      ImageContainerBySTLVector<Domain,double> > Reaction; 
    Diffusion diffusion; 
    Reaction reaction( epsilon, a, k );
    LieSplittingEvolver<Diffusion,Reaction> e(diffusion, reaction); 

    for (unsigned int i = step; i <= max; i += step) 
    {
      std::stringstream s0; 
      s0 << "iteration # " << i; 
      DGtal::trace.beginBlock( s0.str() );

      e.update( implicitFunction, (tstep*step) ); 

      //3d to 2d display
      std::stringstream s; 
      s << outputFiles << setfill('0') << std::setw(4) << (i/step)+1; 
      writeImage( implicitFunction, s.str(), format, 0.5 );

      DGtal::trace.endBlock();   
    }

    //interactive display after the evolution
    if (vm.count("withVisu")) displayImage( argc, argv, implicitFunction, 0.5 ); 

  } else trace.error() << "unknown algo. Try option -h to see the available algorithms " << std::endl;

  
  return 0;
}

