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

using namespace DGtal; 
using namespace Z3i; 
using namespace std; 

//evolvers
//level set
#include "WeickertKuhneEvolver.h"

//phase field
#include "ExactDiffusionEvolver.h"
#include "ExactReactionEvolver.h"
#include "ExplicitReactionEvolver.h"
#include "LieSplittingEvolver.h"

//local level set
#include "SimplePointHelper.h"
#include "PartitionEvolver.h"

/////////////////////////// useful functions
#include "deformationFunctions.h"
#include "deformationDisplay3d.h"
#include "DGtal/io/readers/VolReader.h"

// Qt
#include <QApplication>

///////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  
  DGtal::trace.info() << "3d interface evolution using DGtal ";
  DGtal::trace.emphase() << "(version "<< DGTAL_VERSION << ")"<< std::endl;
  
  // QApplication initialization with command-line parameters
  QApplication application(argc, argv);

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("inputImage,i",  po::value<string>(), "Binary image to initialize the starting interface (vol format)" )
    ("domainSize,d",  po::value<int>()->default_value(64), "Domain size (if default starting interface)" )
    ("shape,s", po::value<string>()->default_value("ball"), 
     "Generated shape: either <ball> (default) or <flower> " )
    ("timeStep,t",  po::value<double>()->default_value(0.25), "Time step for the evolution" )
    ("displayStep",  po::value<int>()->default_value(1), "Number of time steps between 2 drawings" )
    ("stepsNumber,n",  po::value<int>()->default_value(1), "Maximal number of steps" )
    ("algo,a",  po::value<string>()->default_value("levelSet"), 
     "can be: \n <levelSet>  \n or <phaseField> \n or <localLevelSet>" )
    ("balloonForce,k",  po::value<double>()->default_value(0.0), "Balloon force" )
    ("epsilon,e",  po::value<double>()->default_value(3.0), "Interface width (only for phase fields)" )
    ("outputFiles,o",   po::value<string>()->default_value("interface"), "Output files basename" )
    ("outputFormat,f",   po::value<string>()->default_value("vol"), 
     "Output files format: either <png> (3d to 2d with QGLViewer), <pngc> (3d to 2d with Cairo) or <vol> (3d, default)" )
    ("withVisu", "Enables interactive 3d visualization after evolution" );

  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  po::notify(vm);    
  if(vm.count("help")||argc<=1)
    {
      trace.info()<< "Evolution of a 3d interface" << std::endl
		  << "Basic usage: "<<std::endl
		  << argv[0] << " [other options] -t <time step> -n <number of steps> --withVisu" 
		  << std::endl
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
  if ((format != "png")&&(format != "vol")&&(format != "pngc")) 
    {
      trace.info() << "format is expected to be either <png>, <pngc> or <vol> " << std::endl;
      return 0; 
    }

  //image and implicit function
  typedef ImageContainerBySTLVector<Domain,short int> LabelImage;
  LabelImage* labelImage = NULL; 
 
  if (vm.count("inputImage")) 
    { 
      string imageFileName = vm["inputImage"].as<std::string>();
      trace.emphase() << imageFileName <<std::endl; 
      DGtal::trace.beginBlock("image reading..."); 
      LabelImage tmp = VolReader<LabelImage>::importVol( imageFileName );
      labelImage = new LabelImage( tmp ); 
      DGtal::trace.endBlock(); 
    }
  else
    {
      DGtal::trace.beginBlock("image reading..."); 
      Point p(0,0,0);
      Point q(dsize,dsize,dsize); 
      Point c(dsize/2,dsize/2,dsize/2); 
      labelImage = new LabelImage( Domain(p,q) ); 
      if (vm.count("shape"))
        {
          if ( (vm["shape"].as<std::string>()) == "flower" )
	    initWithFlowerPredicate( *labelImage, c, (dsize*3/5)/2, (dsize*1/5)/2, 5 );
          else 
	    initWithBallPredicate( *labelImage, c, (dsize*3/5)/2 );
        } 
      else 
        initWithBallPredicate( *labelImage, c, (dsize*3/5)/2 );
      trace.info() << "starting interface initialized with a " 
		   << (vm["shape"].as<std::string>()) << std::endl;
      DGtal::trace.endBlock(); 
    }

  //domain
  Domain d = Domain( labelImage->domain().lowerBound(), labelImage->domain().upperBound() );

  if (vm["outputFormat"].as<std::string>() == "png") //visu for choosing camera position
    displayPartition( *labelImage );     

  //algo
  std::string algo; 
  if (!(vm.count("algo"))) trace.info() << "default algorithm: levelSet" << std::endl; 
  algo = vm["algo"].as<string>(); 

  if (algo.compare("levelSet")==0)
    {

      //balloon force
      double k; 
      if (!(vm.count("balloonForce"))) trace.info() << "balloon force default value: 0" << std::endl; 
      k = vm["balloonForce"].as<double>(); 

      ImageContainerBySTLVector<Domain,double> implicitFunction( d ); 
      initWithDT( *labelImage, implicitFunction );

      //data functions
      ImageContainerBySTLVector<Domain,double> a( d ); 
      std::fill(a.begin(),a.end(), 1.0 );  
      ImageContainerBySTLVector<Domain,double> b( d ); 
      std::fill(b.begin(),b.end(), 1.0 );  
      ImageContainerBySTLVector<Domain,double> g( d ); 
      std::fill(g.begin(),g.end(), 1.0 );  

      //interface evolver
      WeickertKuhneEvolver<ImageContainerBySTLVector<Domain,double> > e(a,b,g,k,1); 

      DGtal::trace.beginBlock( "Deformation (Weickert's level set method)" );

      double sumt = 0; 
      for (unsigned int i = 1; i <= max; ++i) 
	{
	  DGtal::trace.info() << "iteration # " << i << std::endl; 

	  //update
	  e.update(implicitFunction,tstep); 

	  if ((i%step)==0) 
	    {

	      //display
	      std::stringstream s; 
	      s << outputFiles << setfill('0') << std::setw(4) << (i/step);
	      updateLabelImage( *labelImage, implicitFunction ); 
	      writePartition( *labelImage, s.str(), format );

	    }
	  sumt += tstep; 
	  DGtal::trace.info() << "Time spent: " << sumt << std::endl;    
	}

      DGtal::trace.endBlock();

      //interactive display after the evolution
      updateLabelImage( *labelImage, implicitFunction, 0 ); 
      if (vm.count("withVisu")) 
	displayImageWithInfo( *labelImage, implicitFunction, a, b ); 

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

      ImageContainerBySTLVector<Domain,double> implicitFunction( d ); 
      initWithDT( *labelImage, implicitFunction );

      //computing the profile from the signed distance
      Profile p(epsilon); 
      std::transform(implicitFunction.begin(), implicitFunction.end(), implicitFunction.begin(), p); 

      typedef ExactDiffusionEvolver<ImageContainerBySTLVector<Domain,double> > Diffusion; 
      typedef ExactReactionEvolver<ImageContainerBySTLVector<Domain,double> > Reaction; 
      Diffusion diffusion; 
      Reaction reaction( epsilon );
      // typedef ExplicitReactionEvolver<ImageContainerBySTLVector<Domain,double>, 
      //   ImageContainerBySTLVector<Domain,double> > Reaction; 
      // Diffusion diffusion; 
      // ImageContainerBySTLVector<Domain,double> a( Domain( implicitFunction.domain() ) ); 
      // std::fill(a.begin(),a.end(), 1.0 );  
      // Reaction reaction( epsilon, a, k );
      LieSplittingEvolver<Diffusion,Reaction> e(diffusion, reaction); 

      DGtal::trace.beginBlock( "Deformation (phase field)" );

      double sumt = 0; 
      for (unsigned int i = 1; i <= max; ++i) 
	{
	  DGtal::trace.info() << "iteration # " << i << std::endl; 

	  //update
	  e.update(implicitFunction,tstep); 

	  if ((i%step)==0) 
	    {

	      //display
	      std::stringstream s; 
	      s << outputFiles << setfill('0') << std::setw(4) << (i/step); 
	      updateLabelImage( *labelImage, implicitFunction, 0.5 ); 
	      writePartition( *labelImage, s.str(), format );

	    }
	  sumt += tstep; 
	  DGtal::trace.info() << "Time spent: " << sumt << std::endl;    
	}

      DGtal::trace.endBlock();

      //interactive display after the evolution
      updateLabelImage( *labelImage, implicitFunction, 0.5 ); 
      if (vm.count("withVisu")) displayPartition( *labelImage ); 

    } else if (algo.compare("localLevelSet")==0)
    {
      //space
      KSpace ks;
      ks.init( d.lowerBound(), d.upperBound(), true ); 
   
      //distance image...
      typedef ImageContainerBySTLVector<Domain,double> DistanceImage; 
      //and extern data...
      DistanceImage g( d );
      std::fill(g.begin(), g.end(), 1.0); 

      // topological predicate
      typedef SimplePointHelper<LabelImage> TopologicalPredicate; 
      TopologicalPredicate topologicalPredicate(*labelImage); 
      
      //frontier evolver
      PartitionEvolver<KSpace, LabelImage, DistanceImage, DistanceImage, 
	TopologicalPredicate> 
	e(ks, *labelImage, g, topologicalPredicate); 

      DGtal::trace.info() << e << std::endl; 
      
      DGtal::trace.beginBlock( "Deformation (narrow band with topological control)" );

      double sumt = 0; 
      for (unsigned int i = 1; i <= max; ++i) 
	{
	  DGtal::trace.info() << "iteration # " << i << std::endl; 

	  //update
	  e.update(tstep); 

	  if ((i%step)==0) 
	    {
	      //display
	      std::stringstream s; 
	      s << outputFiles << setfill('0') << std::setw(4) << (i/step); 
	      writePartition( *labelImage, s.str(), format );
	    }
	  sumt += tstep; 
	  DGtal::trace.info() << "Time spent: " << sumt << std::endl;    
      }

      DGtal::trace.endBlock();

      //interactive display after the evolution
      if (vm.count("withVisu")) displayPartition( *labelImage ); 
      

    } else trace.error() << "unknown algo. Try option -h to see the available algorithms " << std::endl;

  //free
  delete( labelImage ); 

  return 1;
}

