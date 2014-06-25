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
#include "DGtal/io/readers/PNMReader.h"

using namespace DGtal; 
using namespace Z2i;
using namespace std; 
//evolvers
//level set
#include "WeickertKuhneEvolver.h"

//phase field
#include "ExactDiffusionEvolver.h"
#include "ExactReactionEvolver.h"
#include "LieSplittingEvolver.h"


/////////////////////////// useful functions
#include "deformationFunctions.h"
#include "deformationDisplay2d.h"

template< typename TEvolver, typename TImage >
void evolution (TEvolver& e, TImage& img, 
		const int& n, const double& tstep, 
		const double& R0, const double& threshold = 0.0) 
{

  int i = 0; 
  int area = setSize( img, threshold ); 
  while ( area > 0 && (i < n) ) 
    {
      ++i; 
      std::stringstream s0; 
      s0 << "iteration # " << i; 
      DGtal::trace.beginBlock( s0.str() );

      e.update( img, tstep); 

      DGtal::trace.endBlock(); 

      area = setSize( img, threshold ); 
      std::cout << "# time   expected area   computed area" << std::endl; 
      double t = i*tstep; 
      double Rt = std::sqrt(R0*R0 - (2*t)); 
      std::cout << t << " " << (M_PI*Rt*Rt) << " " << area << std::endl; 
    }

}

///////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  
  DGtal::trace.info() << "2d interface evolution testing using DGtal ";
  DGtal::trace.emphase() << "(version "<< DGTAL_VERSION << ")"<< std::endl; 

  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("domainSize,s",  po::value<int>()->default_value(64), "Domain size" )
    ("timeStep,t",  po::value<double>()->default_value(1.0), "Time step for the evolution" )
    ("iterationsNumber,n",  po::value<int>()->default_value(10), "Maximal number of iterations" )
    ("algo,a",  po::value<string>()->default_value("levelSet"), 
     "can be: \n <levelSet>  \n or <phaseField> " )
    ("epsilon,e",  po::value<double>()->default_value(3.0), "Interface width (only for phase fields)" )
    ("withFunction",   po::value<string>(), "Output pgm file basename, where the starting implicit function is stored" );

  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  po::notify(vm);    
  if(vm.count("help")||argc<=1)
    {
      trace.info()<< "MCM test on disks" << std::endl
		  << "Basic usage: "<<std::endl
		  << argv[0] << " [other options] -s 128" << std::endl
		  << general_opt << "\n";
      return 0;
    }
  
  //Parse options
  //domain size
  int dsize; 
  if (!(vm.count("domainSize"))) trace.info() << "Domain size default value: 64" << std::endl; 
  dsize = vm["domainSize"].as<int>(); 

  //time step
  double tstep; 
  if (!(vm.count("timeStep"))) trace.info() << "time step default value: 1.0" << std::endl; 
  tstep = vm["timeStep"].as<double>(); 
    
  int n = 10; 
  if (!(vm.count("iterationsNumber"))) trace.info() << "maximal number of iteration: 10" << std::endl; 
  n = vm["iterationsNumber"].as<int>(); 
    
  //image and implicit function
  Point p(0,0);
  Point q(dsize,dsize); 
  Point c(dsize/2,dsize/2); 
  ImageContainerBySTLVector<Domain,double> implicitFunction(Domain(p,q));
  double R0 = (dsize*3/5)/2; 
  initWithBall( implicitFunction, c, R0 ); 
  trace.info() << "# starting interface initialized with a disk of radius "; 
  trace.info() << R0 << std::endl;
 

  //algo
  std::string algo; 
  if (!(vm.count("algo"))) trace.info() << "default algorithm: levelSet" << std::endl; 
  algo = vm["algo"].as<string>(); 

  if (algo.compare("levelSet")==0)
    {


      if (vm.count("withFunction")) 
	drawFunction( implicitFunction, vm["withFunction"].as<string>() ); 


      //data functions
      ImageContainerBySTLVector<Domain,double> f(Domain(p,q)); 
      std::fill(f.begin(),f.end(), 1.0 );  

      //evolution
      WeickertKuhneEvolver<ImageContainerBySTLVector<Domain,double> > e(f,f,f,0,1); 
      evolution(e, implicitFunction, n, tstep, R0); 

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

      if (vm.count("withFunction")) 
        drawFunction( implicitFunction, vm["withFunction"].as<string>() ); 


      typedef ExactDiffusionEvolver<ImageContainerBySTLVector<Domain,double> > Diffusion; 
      typedef ExactReactionEvolver<ImageContainerBySTLVector<Domain,double> > Reaction; 
      Diffusion diffusion; 
      Reaction reaction( epsilon );
      LieSplittingEvolver<Diffusion,Reaction> e(diffusion, reaction); 
      evolution(e, implicitFunction, n, tstep, R0, 0.5); 

    } else trace.error() << "unknown algo. Try option -h to see the available algorithms " << std::endl;


  return 1;
}

