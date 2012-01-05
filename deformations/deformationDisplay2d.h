//////////////////////////////////// display
//display
#include <DGtal/io/boards/Board2D.h>
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/PNMWriter.h"

template< typename TImage >
bool drawContour(const TImage& img, std::string filename, std::string format)
{

  if (format.compare("vector")==0)
  {

  Color green( 0, 192, 0 );
  Color blue( 0, 0, 192 );
  
  Board2D b; 
  Z2i::Domain d(img.lowerBound(), img.upperBound()); 
  Z2i::Domain::ConstIterator cIt = d.begin(); 
  Z2i::Domain::ConstIterator cItEnd = d.end(); 
  b << img.lowerBound() << img.upperBound();
  for ( ; cIt != cItEnd; ++cIt)
  { 
    if (img(*cIt) <= 0) 
      b << CustomStyle( (*cIt).className(), new CustomFillColor(blue) );
    else
      b << CustomStyle( (*cIt).className(), new CustomFillColor(green) );
    b << *cIt; 
  }

  #ifdef WITH_CAIRO
  std::stringstream s; 
  s << filename << ".png"; 
  b.saveCairo(s.str().c_str(),Board2D::CairoPNG);
  #else
  std::stringstream s; 
  s << filename << ".eps"; 
  b.saveEPS(s.str().c_str());
  #endif
  return true; 

  } else if (format.compare("raster")==0)
  {

    //create a label image from the implicite function
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
    //write it into a pgm file
    std::stringstream s; 
    s << filename << ".pgm";
    typedef GradientColorMap<typename LabelImage::Value, DGtal::CMAP_GRAYSCALE> ColorMap; 
    PNMWriter<LabelImage,ColorMap>::exportPGM( s.str(), labelImage, 0, 255, true );

    return true; 
 } else return false; 
}



