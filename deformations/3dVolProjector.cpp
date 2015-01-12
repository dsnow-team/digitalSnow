/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 * @file visuDistanceTransform.cpp
 * @ingroup Examples
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2011/01/04
 *
 * An example file named qglViewer.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <sstream>
#include <iomanip>

// Qt
#include <QApplication>
#include <QCoreApplication>

#include "DGtal/base/Common.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"


using namespace std;
using namespace DGtal;
using namespace Z3i;

#include "deformationDisplay3d.h"

/////////////////////
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

namespace po = boost::program_options;

int displayOneFile( 
		   const string& inputFilename, 
		   const string& outputBasename,
		   const int& offset = 0, const double& step = 0)
{
  //image reading
  typedef ImageSelector<Domain, unsigned char>::Type Image;
  Image image = VolReader<Image>::importVol( inputFilename );
 
  QCoreApplication* application = QCoreApplication::instance();
  Viewer3D<> viewer;
  viewer.show();
 
  //display
  //  displayPartition(viewer, image); 
  Domain domain(image.domain());

  GradientColorMap<long> colorMap( 0, 510 );
  colorMap.addColor(Color::Yellow);
  colorMap.addColor(Color::Blue);
  colorMap.addColor(Color::Red);
  colorMap.addColor(Color::Green);

  for(Domain::ConstIterator it = domain.begin(), itend=domain.end(); it!=itend; ++it){
    unsigned char  val= image( *it );     
    Color c = colorMap( val );
    if(val > 0){
      viewer << CustomColors3D(c, c);     
      viewer << *it;     
    }     
  }


  viewer << Viewer3D<>::updateDisplay;

  if (!viewer.restoreStateFromFile())
    {
      string s = viewer.stateFileName().toStdString(); 
      trace.emphase() << " file " << s 
		      << " not found " 
		      << std::endl;
    }
  viewer.updateGL(); 

  if ( (offset != 0)&&(step != 0) )
    {
      viewer.camera()->setOrientation( -(offset*step), 0.0);
      viewer.showEntireScene(); 
    }

  viewer.setSnapshotFileName(outputBasename.c_str());  
  viewer.setSnapshotFormat("PNG");  
  viewer.saveSnapshot(true, true); 

  {//rename snapshot
    std::stringstream olds;
    olds << viewer.snapshotFileName().toStdString()
	 << "-" << setfill('0') << std::setw(4) 
	 << (viewer.snapshotCounter()-1) << ".png"; 
    string oldf = olds.str(); 
    std::stringstream news; 
    news << outputBasename << ".png";
    string newf = news.str();  
    if (rename (oldf.c_str(), newf.c_str()) == -1) 
      trace.info() << "renaming " << oldf << " into " 
		   << newf << " failed " << std::endl; 
  }

  {//rename state file
    string oldf = viewer.stateFileName().toStdString();
    std::stringstream news; 
    news << ".qglviewer" << (QGLViewer::QGLViewerIndex(&viewer)+1) << ".xml";
    string newf = news.str();  
    if (rename (oldf.c_str(), newf.c_str()) == -1) 
      trace.info() << "renaming " << oldf << " into " 
		   << newf << " failed " << std::endl; 
  }

  application->exit();

  return 0; 
}

int main( int argc, char** argv )
{
  // Init Qt with command-line parameters
  QApplication application(argc, argv);
    
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("output-file,o",   po::value<string>()->default_value("interface"), "output file(s) basename" )
    ("input-file,i", po::value<std::string>(), "volume file" )
    ("multi-input,mi", po::value<std::string>(), "volume files basename " )
    ("start,s",  po::value<int>()->default_value(1), "starting number (for option -mi)" )
    ("end,e",  po::value<int>()->default_value(2), "ending number, not included (for option -mi)" )
    ("angle-step,a",  po::value<int>()->default_value(360), "angle step as a fraction of 2PI (0 disables this feature)" )
    ("number,n",  po::value<int>()->default_value(90), "number of angle steps when moving camera" );
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  po::notify(vm);    
  if(vm.count("help")||argc<=1)
    {
      std::cout << "Usage: " << argv[0] << " [input-file]\n"
		<< "Display volume file as a set of digital frontiers"
		<< general_opt << "\n";
      return 0;
    }
  
  //files
  if (!(vm.count("output-file"))) 
    trace.info() << "output file begin with : interface" << std::endl;
  string outputBasename = vm["output-file"].as<std::string>();


  if(vm.count("input-file"))
    {
      if (vm.count("multi-input"))
	{
	  trace.error() << " Cannot use both input options in the same time " 
			<< endl;      
	  return 1;
	}
      string inputFilename = vm["input-file"].as<std::string>();

      if(vm.count("angle-step"))
	{
	  int den = vm["angle-step"].as<int>(); 
	  double angleStep = (den == 0)?0.0:( (2.0*M_PI)/((double) den) );
	  int n = vm["number"].as<int>(); 
	  for (int i = 0; i < n; ++i)
	    {
	      trace.info() << i << std::endl; 
	      std::stringstream so;
	      so << outputBasename << setfill('0') << std::setw(4) 
		 << i; 

	      displayOneFile(inputFilename, so.str(), i+1, angleStep); 
	    }

	  {//rename state file
	    std::stringstream olds;
	    olds << ".qglviewer" << n << ".xml"; 
	    string oldf = olds.str();
	    string newf = ".qglviewer.xml";  
	    if (rename (oldf.c_str(), newf.c_str()) == -1) 
	      trace.info() << "renaming " << oldf << " into " 
			   << newf << " failed " << std::endl; 
	  }

	}
      else 
	{
	  displayOneFile(inputFilename, outputBasename); 


	  {//rename state file
	    string oldf = ".qglviewer1.xml";
	    string newf = ".qglviewer.xml";  
	    if (rename (oldf.c_str(), newf.c_str()) == -1) 
	      trace.info() << "renaming " << oldf << " into " 
			   << newf << " failed " << std::endl; 
	  }
	}
    }
  else
    {
      if (vm.count("multi-input"))
	{
	  string inputBasename = vm["multi-input"].as<std::string>();

	  int den = vm["angle-step"].as<int>(); 
	  double angleStep = (den == 0)?0.0:( (2.0*M_PI)/((double) den) );

	  int start = vm["start"].as<int>();
	  int end = vm["end"].as<int>();
	  for (int i = start; i != end; ++i)
	    {
	      std::stringstream si;
	      si << inputBasename << setfill('0') << std::setw(4) 
		 << i << ".vol"; 
	      trace.info() << si.str() << std::endl; 
	      std::stringstream so;
	      so << outputBasename << setfill('0') << std::setw(4) 
		 << i; 

	      displayOneFile(si.str(), so.str(), i+1, angleStep); 

	    }


	  {//rename state file
	    std::stringstream olds;
	    olds << ".qglviewer" << (end-start) << ".xml"; 
	    string oldf = olds.str();
	    string newf = ".qglviewer.xml";  
	    if (rename (oldf.c_str(), newf.c_str()) == -1) 
	      trace.info() << "renaming " << oldf << " into " 
			   << newf << " failed " << std::endl; 
	  }

	}
      else
	{
	  trace.error() << " No file name defined" << endl;      
	  return 1;
	}
    }


  return 0; 
}
