#include <sstream>
#include <iomanip>
#include <cstddef>

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
#include "MultiPhaseField.h"

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

  // Default options
  size_t dsize = 64;          // Domain size
  double tstep = 0.25;        // Time step
  size_t disp_step  = 1;      // Display step
  size_t max_step = 1;        // Maximum number of steps
  string shape = "ball";      // Generated shape
  string algo  = "levelSet";  // Algorithm used for evolution
  double balloon = 0.;        // Balloon Force
  double epsilon = 3.;        // Interface width

  string outputFiles  = "interface";  // Output files basename
  string outputFormat = "vol";        // Output files format
  bool flagWithCstVol = false;  // Volume Conservation
  bool flagWithVisu   = false;  // Interactive 3d visualization after evolution

  // Parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h",          "display this message")
    ("inputImage,i",    po::value<string>(), "Binary image to initialize the starting interface (vol format)" )
    ("domainSize,d",    po::value<size_t>(&dsize)->default_value(dsize), "Domain size (if default starting interface)" )
    ("shape,s",         po::value<string>(&shape)->default_value(shape), 
        "Generated shape: either <ball> or <flower> " )
    ("timeStep,t",      po::value<double>(&tstep)->default_value(tstep), "Time step for the evolution" )
    ("displayStep",     po::value<size_t>(&disp_step)->default_value(disp_step), "Number of time steps between 2 drawings" )
    ("stepsNumber,n",   po::value<size_t>(&max_step)->default_value(max_step), "Maximal number of steps" )
    ("algo,a",          po::value<string>(&algo)->default_value(algo), 
        "can be: \n <levelSet>  \n or <phaseField> \n or <multiPhaseField> \n or <localLevelSet>" )
    ("balloonForce,k",  po::value<double>(&balloon)->default_value(balloon), "Balloon force" )
    ("epsilon,e",       po::value<double>(&epsilon)->default_value(epsilon), "Interface width (only for phase fields)" )
    ("withCstVol",      po::bool_switch(&flagWithCstVol), "with volume conservation (only for phase fields)" )
    ("outputFiles,o",   po::value<string>(&outputFiles)->default_value(outputFiles), "Output files basename" )
    ("outputFormat,f",  po::value<string>(&outputFormat)->default_value(outputFormat), 
        "Output files format: either <png> (3d to 2d with QGLViewer), <pngc> (3d to 2d with Cairo) or <vol> (3d)" )
    ("withVisu",        po::bool_switch(&flagWithVisu), "Enables interactive 3d visualization after evolution" );


  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  po::notify(vm);    
  if(vm.count("help") || argc<=1)
    {
      trace.info()<< "Evolution of a 3d interface" << std::endl
        << "Basic usage: "<<std::endl
        << argv[0] << " [other options] -t <time step> -n <number of steps> --withVisu" 
        << std::endl
        << general_opt << "\n";
      return 1;
    }

  // Options validity
  // Files format
  if ( (outputFormat != "png") && (outputFormat != "vol") && (outputFormat != "pngc") ) 
    {
      trace.info() << "format is expected to be either <png>, <pngc> or <vol> " << std::endl;
      return 1; 
    }

  // Generated shape
  if ( vm.count("inputImage") == 0 && shape != "ball" && shape != "flower" )
    {
      trace.info() << "if no input file is specified, shape is expected to be either <ball> or <flower> " << std::endl;
      return 1;
    }

  // Evolution algorithm
  if ( algo != "levelSet" && algo != "phaseField" && algo != "multiPhaseField" && algo != "localLevelSet" )
    {
      trace.info() << "algo is expected to be either <levelSet>, <phaseField>, <multiPhaseField> or <localLevelSet> " << std::endl;
      return 1;
    }


  // Image and implicit function
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
      Point p(0, 0, 0);
      Point q(dsize, dsize, dsize); 
      Point c(dsize/2, dsize/2, dsize/2); 
      labelImage = new LabelImage( Domain(p,q) ); 

      if ( (vm["shape"].as<std::string>()) == "flower" )
        initWithFlowerPredicate( *labelImage, c, (dsize*3/5)/2, (dsize*1/5)/2, 5 );
      else 
        initWithBallPredicate( *labelImage, c, (dsize*3/5)/2 ); 
      
      trace.info() << "starting interface initialized with a " << shape << std::endl;

      DGtal::trace.endBlock(); 
    }

  // Domain
  Domain d = Domain( labelImage->domain().lowerBound(), labelImage->domain().upperBound() );

  // Pre-visualization for choosing camera position
  if (outputFormat == "png")
    displayPartition( *labelImage );     

  // Algorithm dispatch
  if ( algo == "levelSet" )
    {

      // Distance function
      ImageContainerBySTLVector<Domain,double> implicitFunction( d ); 
      initWithDT( *labelImage, implicitFunction );

      // Data functions
      ImageContainerBySTLVector<Domain,double> a( d ); 
      std::fill(a.begin(),a.end(), 1.0 );  
      ImageContainerBySTLVector<Domain,double> b( d ); 
      std::fill(b.begin(),b.end(), 1.0 );  
      ImageContainerBySTLVector<Domain,double> g( d ); 
      std::fill(g.begin(),g.end(), 1.0 );  

      // Interface evolver
      WeickertKuhneEvolver<ImageContainerBySTLVector<Domain,double> > e(a, b, g, balloon, 1); 

      DGtal::trace.beginBlock( "Deformation (Weickert's level set method)" );

      // Time integration
      double sumt = 0; 
      for (unsigned int i = 1; i <= max_step; ++i) 
        {
          DGtal::trace.info() << "iteration # " << i << std::endl; 

          // Update
          e.update(implicitFunction,tstep); 

          // Display
          if ((i % disp_step) == 0) 
            {
              std::stringstream s; 
              s << outputFiles << setfill('0') << std::setw(4) << (i/disp_step);
              updateLabelImage( *labelImage, implicitFunction ); 
              writePartition( *labelImage, s.str(), outputFormat );
            }

          sumt += tstep; 
          DGtal::trace.info() << "Time spent: " << sumt << std::endl;    
        }

      updateLabelImage( *labelImage, implicitFunction, 0 );
      DGtal::trace.info() << "Volume: " << getSize(*labelImage, 0) << std::endl;
      DGtal::trace.endBlock();

      // Interactive display after the evolution
      if ( flagWithVisu ) 
        displayImageWithInfo( *labelImage, implicitFunction, a, b ); 

    } 
  else if ( algo == "phaseField" )
    {
      // Option's validity
      if (epsilon <= 0) 
        {
          trace.error() << "epsilon should be greater than 0" << std::endl;
          return 1; 
        } 

      // Distance function
      ImageContainerBySTLVector<Domain,double> implicitFunction( d ); 
      initWithDT( *labelImage, implicitFunction );

      // Computing the profile from the signed distance
      Profile p(epsilon); 
      std::transform(implicitFunction.begin(), implicitFunction.end(), implicitFunction.begin(), p); 

      // Diffusion evolver
      typedef ExactDiffusionEvolver<ImageContainerBySTLVector<Domain,double> > Diffusion; 
      Diffusion diffusion;

      // Exact reaction evolver (no possible volume conservation)
      typedef ExactReactionEvolver<ImageContainerBySTLVector<Domain,double> > Reaction; 
      Reaction reaction( epsilon );

      // Explicit reaction evolver (with possible volume conservation)
      /*
      typedef ExplicitReactionEvolver<
          ImageContainerBySTLVector<Domain,double>, 
          ImageContainerBySTLVector<Domain,double> 
        > Reaction; 
      ImageContainerBySTLVector<Domain,double> a( Domain( implicitFunction.domain() ) ); 
      std::fill(a.begin(), a.end(), 1.0 );  
      Reaction reaction( epsilon, a, balloon, flagWithCstVol );
      */

      // Lie splitting
      LieSplittingEvolver<Diffusion,Reaction> e(diffusion, reaction); 

      DGtal::trace.beginBlock( "Deformation (phase field)" );

      // Initial state export
      std::stringstream s; 
      s << outputFiles << setfill('0') << std::setw(4) << 0; 
      writePartition( *labelImage, s.str(), outputFormat );
      
      // Time integration
      double sumt = 0; 
      for (unsigned int i = 1; i <= max_step; ++i) 
        {
          DGtal::trace.info() << "iteration # " << i << std::endl; 

          // Update
          e.update(implicitFunction, tstep); 

          // Display
          if ((i % disp_step)==0) 
            {
              std::stringstream s; 
              s << outputFiles << setfill('0') << std::setw(4) << (i/disp_step); 
              updateLabelImage( *labelImage, implicitFunction, 0.5 ); 
              writePartition( *labelImage, s.str(), outputFormat );
            }

          sumt += tstep; 
          DGtal::trace.info() << "Time spent: " << sumt << std::endl;    
        }

      // Post-processing
      updateLabelImage( *labelImage, implicitFunction, 0.5 );
      DGtal::trace.info() << "Volume: " << getSize(*labelImage, 0.5) << std::endl;
      DGtal::trace.endBlock();

      // Interactive display after the evolution
      if ( flagWithVisu ) 
        displayPartition( *labelImage ); 

    } 
  else if ( algo == "multiPhaseField" )
    {
      // Options's validity
      if (epsilon <= 0) 
        {
          trace.error() << "epsilon should be greater than 0" << std::endl;
          return 1; 
        }

      // Field image
      typedef ImageContainerBySTLVector<Domain, double> FieldImage;

      // Diffusion evolver
      typedef ExactDiffusionEvolver< FieldImage > Diffusion; 
      Diffusion diffusion;
     
      // Exact Reaction Evolver (no possible volume conservation)
      typedef ExactReactionEvolver < FieldImage > Reaction; 
      Reaction reaction( epsilon );
      
      // Explicit Reaction Evolver (with possible volume conservation)
      /*
      typedef ExplicitReactionEvolver<
          ImageContainerBySTLVector<Domain,double>, 
          ImageContainerBySTLVector<Domain,double> 
        > Reaction; 
      ImageContainerBySTLVector<Domain,double> a( d ); 
      std::fill(a.begin(), a.end(), 1.0 );  
      Reaction reaction( epsilon, a, balloon, flagWithCstVol );
      */

      // Lie splitting
      LieSplittingEvolver< Diffusion, Reaction > phaseEvolver(diffusion, reaction); 

      // Multi phase-field
      MultiPhaseField< LabelImage, FieldImage, LieSplittingEvolver<Diffusion,Reaction> > evolver(*labelImage, phaseEvolver);
      
      DGtal::trace.beginBlock( "Deformation (multi phase field)" );

      // Initial state export
      std::stringstream s; 
      s << outputFiles << setfill('0') << std::setw(4) << 0; 
      writePartition( *labelImage, s.str(), outputFormat );
      
      // Time integration
      double sumt = 0; 
      for (unsigned int i = 1; i <= max_step; ++i) 
        {
          DGtal::trace.info() << "iteration # " << i << std::endl; 

          // Update
          evolver.update( tstep ); 

          // Display
          if ( (i % disp_step) == 0 ) 
            {
              std::stringstream s; 
              s << outputFiles << setfill('0') << std::setw(4) << (i/disp_step); 
              writePartition( *labelImage, s.str(), outputFormat );

            }

          sumt += tstep; 
          DGtal::trace.info() << "Time spent: " << sumt << std::endl;    
        }

      //TODO: Volume of each label ?
      //DGtal::trace.info() << "Volume: " << getSize(*labelImage, 1.5) << std::endl;
      
      DGtal::trace.endBlock();

      // Interactive display after the evolution
      if ( flagWithVisu ) 
        displayPartition( *labelImage ); 

    }
  else if ( algo == "localLevelSet" )
    {
      // Space
      KSpace ks;
      ks.init( d.lowerBound(), d.upperBound(), true ); 

      // Distance image...
      typedef ImageContainerBySTLVector<Domain,double> DistanceImage; 
      // ... and extern data.
      DistanceImage g( d );
      std::fill(g.begin(), g.end(), 1.0); 

      // Topological predicate
      typedef SimplePointHelper<LabelImage> TopologicalPredicate; 
      TopologicalPredicate topologicalPredicate(*labelImage); 

      // Frontier evolver
      PartitionEvolver<KSpace, LabelImage, DistanceImage, DistanceImage, 
        TopologicalPredicate> 
          e(ks, *labelImage, g, topologicalPredicate); 

      DGtal::trace.info() << e << std::endl; 

      DGtal::trace.beginBlock( "Deformation (narrow band with topological control)" );

      // Initial state export
      std::stringstream s; 
      s << outputFiles << setfill('0') << std::setw(4) << 0; 
      writePartition( *labelImage, s.str(), outputFormat );
      
      // Time integration
      double sumt = 0; 
      for (unsigned int i = 1; i <= max_step; ++i) 
        {
          DGtal::trace.info() << "iteration # " << i << std::endl; 

          // Update
          e.update(tstep); 

          // Display
          if ((i % disp_step) == 0) 
            {
              std::stringstream s; 
              s << outputFiles << setfill('0') << std::setw(4) << (i/disp_step); 
              writePartition( *labelImage, s.str(), outputFormat );
            }

          sumt += tstep; 
          DGtal::trace.info() << "Time spent: " << sumt << std::endl;    
        }

      DGtal::trace.endBlock();

      // Interactive display after the evolution
      if ( flagWithVisu ) displayPartition( *labelImage ); 

    } 
  else 
    {
      trace.error() << "unknown algo. Try option -h to see the available algorithms " << std::endl;
      return 1;
    }

  // Free
  delete( labelImage ); 

  return 0;
}

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

