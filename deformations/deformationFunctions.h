//images
#include <DGtal/images/ImageContainerBySTLVector.h>

/////////////////////////// useful functions
template< typename TImage >
int setSize(TImage& img, const double& threshold = 0)
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

#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"

#include <functional>

template< typename TValue >
TValue aFunction(const TValue& v)
{
  if (v == 0) return 1; 
  else return 0;
}

struct Predicate 
{
  template< typename TValue >
  bool operator()(const TValue& v)
  {
    if (v == 0) return true; 
    else return false;
  }
}; 

template< typename TImage >
void initWithDT(const TImage& inputImage, ImageContainerBySTLVector<typename TImage::Domain,double>& outputImage)
{
  typedef typename TImage::Domain::Space Space; 

  //domain
  typename TImage::Domain d = inputImage.domain(); 

  //metric
  typedef ExactPredicateLpSeparableMetric<Space, 2> L2Metric;
  L2Metric l2; 

  //Foreground 
  //helpers
  typedef Thresholder<typename TImage::Value, true, true> Binarizer; 
  Binarizer binarizer(0);
  typedef PointFunctorPredicate<TImage,Binarizer> PredicateOnPoints; 
  PredicateOnPoints predicate(inputImage, binarizer); 
  //DT
  typedef DistanceTransformation<Space, PredicateOnPoints, L2Metric > DT;
  DT dt( d, predicate, l2);

  //Background
  //helpers
  typedef Thresholder<typename TImage::Value, false, false> Binarizer2; 
  Binarizer2 binarizer2(0);
  typedef PointFunctorPredicate<TImage,Binarizer2> PredicateOnPoints2; 
  PredicateOnPoints2 predicate2(inputImage, binarizer2); 
  //DT
  typedef DistanceTransformation<Space, PredicateOnPoints2, L2Metric > DT2;
  DT2 dt2( d, predicate2, l2);

  //Signed distance 
  typename ImageContainerBySTLVector<typename TImage::Domain,double>::OutputIterator out; 
  for ( typename DT::ConstRange::ConstIterator 
	  it = dt.constRange().begin(), itend = dt.constRange().end();
	it != itend; ++it)
    {
      double dist = (*it); 
      if (dist != 0) 
	dist = std::sqrt( dist ) - 0.5;  
      *out++ = dist; 
      // outputImage.setValue(p, dist);  
    }

  for ( typename DT2::ConstRange::ConstIterator 
	  it = dt2.constRange().begin(), itend = dt2.constRange().end();
	it != itend; ++it)
    {
      double dist = (*it); 
      if (dist != 0) 
	{ //only write if dist != 0
	  dist = std::sqrt( dist ) - 0.5;  
	  *out++ = dist; 
	  // outputImage.setValue(p, dist);  
	}
      else 
	out++; 
    }

  // typedef  DistanceTransformation<TImage, 2> DT;
  // DT dt;
  // typename DT::OutputImage dtImage = dt.compute ( inputImage );

  // //inv
  // TImage rInputImage(d); 
  // std::transform(inputImage.begin(), inputImage.end(),
  //                rInputImage.begin(), aFunction<typename TImage::Value> ); 
  // typename DT::OutputImage rDtImage = dt.compute ( rInputImage );

  // //deduce the signed distance function 
  // typename TImage::Domain::ConstIterator cIt = d.begin(); 
  // typename TImage::Domain::ConstIterator cItEnd = d.end(); 
  // for ( ; cIt != cItEnd; ++cIt)
  // { //for each domain point

  //   typedef typename TImage::Point Point; 
  //   Point p( *cIt ); //point p

  //   double dist = 0; //signed distance
  //   if ( dtImage( p ) == 0 )
  //     {
  //       if ( rDtImage( p ) == 0 )
  //         {
  //           std::cerr << "Error init with DT" << std::endl;  
  //         }
  //       else
  //         {
  //           dist = std::sqrt( (double) rDtImage( p ) ) - 0.5;  
  //         }
  //     }
  //   else
  //     {
  //       dist = - ( std::sqrt( (double) dtImage( p ) ) - 0.5 );  
  //     }

  //   outputImage.setValue(p, (double) dist);  
  // }
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
