//display 3D
// static 
#include "DGtal/io/boards/Board3DTo2D.h"

template< typename TImage >
void displayImage(const TImage& img, string filename)
{

  Board3DTo2D viewer;
  
  Z3i::Domain d(img.lowerBound(), img.upperBound()); 
  Z3i::Domain::ConstIterator cIt = d.begin(); 
  Z3i::Domain::ConstIterator cItEnd = d.end(); 
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
#endif

}


// interactive
#include <QtGui/qapplication.h>
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/CDrawableWithDisplay3D.h"

template< typename TImage >
bool displayImage(int argc, char** argv, const TImage& img)
{

  bool flag = true;    
  #ifdef WITH_VISU3D_QGLVIEWER

  QApplication application(argc,argv);
  Viewer3D viewer;
  viewer.show();

  
  Z3i::Domain d(img.lowerBound(), img.upperBound()); 
  Z3i::Domain::ConstIterator cIt = d.begin(); 
  Z3i::Domain::ConstIterator cItEnd = d.end(); 
  for ( ; cIt != cItEnd; ++cIt)
  { 
    if (img(*cIt) <= 0) 
	viewer << *cIt; 
  }
  viewer << Viewer3D::updateDisplay;

  return application.exec();
#else
  return false; 
#endif
}



