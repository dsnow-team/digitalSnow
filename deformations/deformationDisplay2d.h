//////////////////////////////////// display
//display
#include <DGtal/io/boards/Board2D.h>
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/PNMWriter.h"


template< typename TImage >
bool drawContour(const TImage& img, std::string filename, std::string format, 
		 const double& threshold = 0)
{

  if (format.compare("vector")==0)
  {

  Color green( 0, 192, 0 );
  Color blue( 0, 0, 192 );
  
  Board2D b; 
  Z2i::Domain d(img.domain()); 
  Z2i::Domain::ConstIterator cIt = d.begin(); 
  Z2i::Domain::ConstIterator cItEnd = d.end(); 
  b << d.lowerBound() << d.upperBound();
  for ( ; cIt != cItEnd; ++cIt)
  { 
    if (img(*cIt) <= threshold) 
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

    //create a label image from the implicit function
    typedef ImageContainerBySTLVector<Domain,int> LabelImage; 
    LabelImage labelImage( Domain( img.domain() ) ); 
    Domain d = labelImage.domain(); 
    Domain::ConstIterator cIt = d.begin(); 
    Domain::ConstIterator cItEnd = d.end(); 
    for ( ; cIt != cItEnd; ++cIt)
    { 
      if (img(*cIt) <= threshold) 
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

template< typename TImage >
bool drawContours(const TImage& img, std::string filename, std::string format, 
		 const double& threshold = 0)
{

  if (format.compare("vector")==0)
  {

    typedef GradientColorMap<typename TImage::Value, DGtal::CMAP_GRAYSCALE> ColorMap; 
    ColorMap colormap(0,255); 
  
    Board2D b; 
    Z2i::Domain d(img.domain()); 
    Z2i::Domain::ConstIterator cIt = d.begin(); 
    Z2i::Domain::ConstIterator cItEnd = d.end(); 
    for ( ; cIt != cItEnd; ++cIt)
      { 
	b << CustomStyle( (*cIt).className(), new CustomFillColor( colormap( img(*cIt) ) ) );
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

    //write it into a pgm file
    std::stringstream s; 
    s << filename << ".pgm";
    typedef GradientColorMap<typename TImage::Value, DGtal::CMAP_GRAYSCALE> ColorMap; 
    PNMWriter<TImage,ColorMap>::exportPGM( s.str(), img, 0, 255, true );

    return true; 
 } else return false; 
}


template< typename TImage >
void drawFunction( const TImage& img, std::string basename) 
{
    typedef GradientColorMap<typename TImage::Value, DGtal::CMAP_GRAYSCALE> ColorMap; 
    typename TImage::Value min = *min_element( img.begin(),img.end() ); 
    typename TImage::Value max = *max_element( img.begin(),img.end() );
    //std::cerr << min << " " << max << std::endl;  
    std::stringstream s; 
    s << basename << ".pgm"; 
    PNMWriter<TImage,ColorMap>::exportPGM( s.str(), img, min, max, true );
} 

