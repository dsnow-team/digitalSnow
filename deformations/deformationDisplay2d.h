//////////////////////////////////// display
//display
#include <DGtal/io/boards/Board2D.h>

template< typename TImage >
void drawContour(const TImage& img, std::string filename)
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

  std::stringstream s; 
#ifdef WITH_CAIRO
  s << filename << ".png";
  b.saveCairo(s.str().c_str(),Board2D::CairoPNG ); 
#else
  s << filename << ".eps"; 
  b.saveEPS(s.str().c_str());
#endif
}



