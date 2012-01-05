//display 3D
// static 
#include "DGtal/io/boards/Board3DTo2D.h"

#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/VolWriter.h"

template< typename TImage >
bool writeImage(const TImage& img, string filename, string format)
{

  if (format.compare("png")==0)
  {
    Board3DTo2D viewer;
    
    Domain d(img.lowerBound(), img.upperBound()); 
    Domain::ConstIterator cIt = d.begin(); 
    Domain::ConstIterator cItEnd = d.end(); 
    for ( ; cIt != cItEnd; ++cIt)
    { 
      if (img(*cIt) <= 0) 
	      viewer << *cIt; 
    }

    viewer << CameraPosition(-10.0, -10.0, -10.0)
	   << CameraDirection(1.000000, 1.000000, 1.000000); 
    // << CameraUpVector(0.000000, 1.000000, 0.000000)
    // << CameraZNearFar(4.578200, 22.578199);

    int size = img.extent().at(0); 
    std::stringstream s; 
  #ifdef WITH_CAIRO
    s << filename << ".png";
    viewer.saveCairo(s.str().c_str(),Board3DTo2D::CairoPNG,3*size/2,3*size/2 ); 
    return true; 
  #else
    return false; 
  #endif

  } else if (format.compare("vol")==0)
  {

    //create a label image from the implicit function
    typedef ImageContainerBySTLVector<Domain,int> LabelImage; 
    LabelImage labelImage( img.lowerBound(), img.upperBound() ); 
    Domain d(labelImage.lowerBound(), labelImage.upperBound()); 
    Domain::ConstIterator cIt = d.begin(); 
    Domain::ConstIterator cItEnd = d.end(); 
    for ( ; cIt != cItEnd; ++cIt)
    { 
      if (img(*cIt) <= 0) 
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
bool displayImage(int argc, char** argv, const TImage& img)
{

  //KhalimskySpace
  Domain d = img.domain(); 
  KSpace K;
  K.init(d.lowerBound(), d.upperBound(), true);
  //adjacency  
  SurfelAdjacency<3> SAdj( true );
  std::vector<std::vector<SCell> > vectConnectedSCell;
  //predicate
  SimpleThresholdForegroundPredicate<TImage> predicate(img,0);
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



