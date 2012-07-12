#if defined(deformationFunctions_RECURSES)
#error Recursive header files inclusion detected in deformationFunctions.h
#else // defined(deformationFunctions_RECURSES)
/** Prevents recursive inclusion of headers. */
#define deformationFunctions_RECURSES

#if !defined deformationFunctions_h
/** Prevents repeated inclusion of headers. */
#define deformationFunctions_h


//images
#include <DGtal/images/ImageContainerBySTLVector.h>

/////////////////////////// useful functions
template< typename TImage >
int getSize(TImage& img, const double& threshold = 0)
{
 
  int c = 0; //counter

  typename TImage::Domain d = img.domain(); 
  typename TImage::Domain::ConstIterator cIt = d.begin(); 
  typename TImage::Domain::ConstIterator cItEnd = d.end(); 
  for ( ; cIt != cItEnd; ++cIt)
  { //for each domain point

    typedef typename TImage::Point Point; 
    Point p( *cIt ); //point p

    if (img(p) <= threshold) ++c; 
  }

  return c; 
}

template< typename TLabelImage, typename TDistanceImage >
void updateLabelImage(TLabelImage& limg, const TDistanceImage& dimg, const double& threshold = 0)
{
 
  typename TLabelImage::Domain d = limg.domain(); 
  typename TLabelImage::Domain::ConstIterator cIt = d.begin(); 
  typename TLabelImage::Domain::ConstIterator cItEnd = d.end(); 
  for ( ; cIt != cItEnd; ++cIt)
  { //for each domain point
    if (dimg(*cIt) <= threshold) 
      limg.setValue(*cIt, 255);   
    else 
      limg.setValue(*cIt, 0);  
  }

}

template< typename TImage >
void inv(TImage& img, const double& threshold = 0)
{
 
  typename TImage::Domain d = img.domain(); 
  typename TImage::Domain::ConstIterator cIt = d.begin(); 
  typename TImage::Domain::ConstIterator cItEnd = d.end(); 
  for ( ; cIt != cItEnd; ++cIt)
  { //for each domain point
    if (img(*cIt) <= threshold) 
      img.setValue( *cIt, 1 );  
    else 
      img.setValue( *cIt, 0 );   
  }

}


template< typename TImage >
void initWithBall(TImage& img, const typename TImage::Point& c, const double& r)
{
 
  typename TImage::Domain d = img.domain(); 
  typename TImage::Domain::ConstIterator cIt = d.begin(); 
  typename TImage::Domain::ConstIterator cItEnd = d.end(); 
  for ( ; cIt != cItEnd; ++cIt)
  { //for each domain point

    typedef typename TImage::Point Point; 
    Point p( *cIt ); //point p

    double dist = (p-c).norm(Point::L_2); 
    dist = dist-r;

    img.setValue(p, (typename TImage::Value) dist);  
  }

}

template< typename TImage >
void initWithBallPredicate(TImage& img, const typename TImage::Point& c, const double& r)
{
 
  typename TImage::Domain d = img.domain(); 
  typename TImage::Domain::ConstIterator cIt = d.begin(); 
  typename TImage::Domain::ConstIterator cItEnd = d.end(); 
  for ( ; cIt != cItEnd; ++cIt)
  { //for each domain point

    typedef typename TImage::Point Point; 
    Point p( *cIt ); //point p

    double dist = (p-c).norm(Point::L_2); 
    dist = dist-r;

    img.setValue(p, (typename TImage::Value)( dist <= 0 ) );  
  }

}


template< typename TImage >
void initWithFlower(TImage& img, const typename TImage::Point& c, double r, double v, double k)
{
 
  typename TImage::Domain d = img.domain(); 
  
  typename TImage::Domain::ConstIterator cIt = d.begin(); 
  typename TImage::Domain::ConstIterator cItEnd = d.end(); 
  for ( ; cIt != cItEnd; ++cIt)
  { //for each domain point


    typedef typename TImage::Point Point; 
    Point p( *cIt ); //point p

    //distance au centre calcule
    double rho = r; 
    double deviation = 0; 
    typedef typename TImage::Dimension Dimension; 
    for (Dimension i = 1; i < TImage::dimension; ++i)
      {
	double t = std::abs(std::atan2((p[i]-c[i]),(p[0]-c[0])));
	deviation += std::cos(k*t);
      }
    rho += v*deviation; 

    //distance au centre
    double dist = (p-c).norm(Point::L_2); 
    dist = dist - rho;

    img.setValue(p, (typename TImage::Value) dist);  
  }
}

template< typename TImage >
void initWithFlowerPredicate(TImage& img, const typename TImage::Point& c, double r, double v, double k)
{
 
  typename TImage::Domain d = img.domain(); 
  
  typename TImage::Domain::ConstIterator cIt = d.begin(); 
  typename TImage::Domain::ConstIterator cItEnd = d.end(); 
  for ( ; cIt != cItEnd; ++cIt)
  { //for each domain point

    typedef typename TImage::Point Point; 
    Point p( *cIt ); //point p

    //distance au centre calcule
    double rho = r; 
    double deviation = 0; 
    typedef typename TImage::Dimension Dimension; 
    for (Dimension i = 1; i < TImage::dimension; ++i)
      {
	double t = std::abs(std::atan2((p[i]-c[i]),(p[0]-c[0])));
	deviation += std::cos(k*t);
      }
    rho += v*deviation; 

    //distance au centre
    double dist = (p-c).norm(Point::L_2); 
    dist = dist - rho;

    img.setValue(p, (typename TImage::Value)(dist <= 0) );  
  }
}


#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"

template< typename TValue >
TValue aFunction(const TValue& v)
{
  if (v == 0) return 1; 
  else return 0;
}

template< typename TImage >
void initWithDT(const TImage& inputImage, ImageContainerBySTLVector<typename TImage::Domain,double>& outputImage)
{

  //inv
  typename TImage::Domain d = inputImage.domain(); 
  TImage rInputImage(d); 
  std::transform(inputImage.begin(), inputImage.end(),
                 rInputImage.begin(), aFunction<typename TImage::Value> ); 
  //DT 
  typedef  DistanceTransformation<TImage, 2> DT;
  DT dt;
  typename DT::OutputImage dtImage = dt.compute ( inputImage );
  typename DT::OutputImage rDtImage = dt.compute ( rInputImage );

  //deduce the signed distance function 
  typename TImage::Domain::ConstIterator cIt = d.begin(); 
  typename TImage::Domain::ConstIterator cItEnd = d.end(); 
  for ( ; cIt != cItEnd; ++cIt)
  { //for each domain point

    typedef typename TImage::Point Point; 
    Point p( *cIt ); //point p

    double dist = 0; //signed distance
    if ( dtImage( p ) == 0 )
      {
        if ( rDtImage( p ) == 0 )
          {
            std::cerr << "Error init with DT" << std::endl;  
          }
        else
          {
            dist = std::sqrt( (double) rDtImage( p ) ) - 0.5;  
          }
      }
    else
      {
        dist = - ( std::sqrt( (double) dtImage( p ) ) - 0.5 );  
      }

    outputImage.setValue(p, (double) dist);  
  }
}

class Profile {
private:  
  double myEpsilon; 
public: 
  Profile (const double& anEpsilon) : myEpsilon( anEpsilon ) {}
  double operator()( const double& v )
  {
   return 0.5 - 0.5*std::tanh(-v/(2*myEpsilon)); 
  }
}; 

#endif // !defined deformationFunctions_h

#undef deformationFunctions_RECURSES
#endif // else defined(deformationFunctions_RECURSES)
