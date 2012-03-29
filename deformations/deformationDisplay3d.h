#include <fstream>

//display 3D
// static 
  #ifdef WITH_CAIRO
#include "DGtal/io/boards/Board3DTo2D.h"
  #endif

#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/VolWriter.h"

template< typename TImage >
bool writeImage(const TImage& img, string filename, string format, const double& threshold = 0)
{

  if (format.compare("png")==0)
  {
  #ifdef WITH_CAIRO
    Board3DTo2D viewer;
    
    Domain d = img.domain(); 
    Domain::ConstIterator cIt = d.begin(); 
    Domain::ConstIterator cItEnd = d.end(); 
    for ( ; cIt != cItEnd; ++cIt)
    { 
      if (img(*cIt) <= threshold) 
	      viewer << *cIt; 
    }

  //reading camera configuration for the 3d to 2d projection
  std::ifstream file(".camera", ios::in); 
  if(file) 
  {       

    std::vector<std::vector<double> > p; //position/direciton
    double znear, zfar; 

    std::string line;
    getline(file, line); //skip the first line
    //for the three following lines
    unsigned int nbLines = 1; 
    while( std::getline(file, line) && (nbLines <= 3) ) 
    {
      p.push_back( std::vector<double>(3) ); 

      std::istringstream isline( line );
      std::string word;
      std::getline( isline, word, ' ' ); //skip the first word
      unsigned int k = 0; 
      while ( std::getline( isline, word, ' ' ) && (k <= 3) )
      {
          std::istringstream isword( word.substr (0,word.size()-1) );
          isword >> p[nbLines-1][k];
          ++k; 
      }
      ++nbLines; 
    }
    //for the last line
    if ( std::getline(file, line) ) 
    { 
      std::istringstream isline( line );
      std::string word;
      std::getline( isline, word, ' ' ); //skip the first word
      std::getline( isline, word, ' ' );
      std::istringstream is1( word );
      is1 >> znear;
      std::getline( isline, word, ' ' ); //skip the third word
      std::getline( isline, word, ' ' ); //skip the fourth word
      std::getline( isline, word, ' ' );
      std::istringstream is2( word );
      is1 >> zfar;
    }
    file.close();

    //setting camera configuration
    viewer << CameraPosition(p[0][0], p[0][1], p[0][2])
	   << CameraDirection(p[1][0], p[1][1], p[1][2]) 
     << CameraUpVector(p[2][0], p[2][1], p[2][2])
     << CameraZNearFar(znear, zfar);

  }
  else  
  {
    trace.emphase() << "Failed to read '.camera'. Default camera configuration" << std::endl;

    //default config
    typename TImage::Vector v = img.extent(); 
    viewer << CameraPosition(v.at(0)/2, v.at(1)/2, 2*v.at(2))
	   << CameraDirection(0, 0, -1) 
     << CameraUpVector(0, 1, 0)
     << CameraZNearFar(v.at(2)/2, 3*v.at(2));
  }

    int size = img.extent().at(0); 
    std::stringstream s; 
    s << filename << ".png";
    viewer.saveCairo(s.str().c_str(),Board3DTo2D::CairoPNG,3*size/2,3*size/2 ); 
    return true; 
  #else
    trace.emphase() << "Failed to use Cairo 3d to 2d (not installed)" << std::endl;
    return false; 
  #endif

  } else if (format.compare("vol")==0)
  {

    //create a label image from the implicit function
    typedef ImageContainerBySTLVector<Domain,int> LabelImage; 
    LabelImage labelImage( img.domain() ); 
    Domain d(labelImage.domain()); 
    Domain::ConstIterator cIt = d.begin(); 
    Domain::ConstIterator cItEnd = d.end(); 
    for ( ; cIt != cItEnd; ++cIt)
    { 
      if (img(*cIt) <= threshold) 
	       labelImage.setValue(*cIt,255);
      else  
	       labelImage.setValue(*cIt,0);
    }
    //write it into a vol file
    std::stringstream s; 
    s << filename << ".vol";
    typedef GradientColorMap<typename LabelImage::Value, DGtal::CMAP_GRAYSCALE> ColorMap; 
    VolWriter<LabelImage,ColorMap>::exportVol( s.str(), labelImage, 0, 255 );

    return true; 

 } else return false; 
  
}


// interactive
#include <QtGui/qapplication.h>
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/CDrawableWithDisplay3D.h"

#include "DGtal/images/imagesSetsUtils/SimpleThresholdForegroundPredicate.h"

#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/helpers/Surfaces.h"

template< typename TImage >
bool displayImage(int argc, char** argv, const TImage& img, const double& threshold = 0)
{

  //KhalimskySpace
  Domain d = img.domain(); 
  KSpace K;
  K.init(d.lowerBound(), d.upperBound(), true);
  //adjacency  
  SurfelAdjacency<3> SAdj( true );
  std::vector<std::vector<SCell> > vectConnectedSCell;
  //predicate
  SimpleThresholdForegroundPredicate<TImage> predicate(img,threshold);
  //tracking 
  Surfaces<KSpace>::extractAllConnectedSCell(vectConnectedSCell,K, SAdj, predicate, true);
  

  #ifdef WITH_VISU3D_QGLVIEWER
  QApplication application(argc,argv);
  Viewer3D viewer;
  viewer.show();

  for(unsigned int i=0; i< vectConnectedSCell.size();i++){
    for(unsigned int j=0; j< vectConnectedSCell.at(i).size();j++){
      viewer << vectConnectedSCell.at(i).at(j);
    }    
  }

  viewer << Viewer3D::updateDisplay;

  return application.exec();
#else
  return false; 
#endif
}



