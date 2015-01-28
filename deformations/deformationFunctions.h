#if defined(deformationFunctions_RECURSES)
#error Recursive header files inclusion detected in deformationFunctions.h
#else // defined(deformationFunctions_RECURSES)
/** Prevents recursive inclusion of headers. */
#define deformationFunctions_RECURSES

#if !defined deformationFunctions_h
/** Prevents repeated inclusion of headers. */
#define deformationFunctions_h

#include <numeric> // Algorithms

//images
#include <DGtal/images/ImageContainerBySTLVector.h>

using namespace DGtal::functors;

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

/** Calculate histogram of an image (ie the size of each label) 
 *
 * @tparam  TImage  type of the image to analyse
 * @tparam  TMap    type of the map (std) in which to store the result
 *
 * @param   img     the image
 * @param   map     the map 
*/
template <
  typename TImage, 
  typename TMap
>
void calcHistogram( TImage const& img, TMap & map )
{
  typedef typename TImage::Domain TDomain;
  typedef typename TMap::key_type T;
  
  TDomain d = img.domain();
  typename TDomain::ConstIterator cIt    = d.begin();
  typename TDomain::ConstIterator cItEnd = d.end();

  for ( ; cIt != cItEnd; ++cIt )
    {
      const T label = static_cast<T>( img(*cIt) );
      if ( map.count(label) == 0 )
        map[label] = 1;
      else
        ++(map[label]);
    }
}

/** Calculate area/volume of a phase from his implicit function
 *
 * It simply integrates the function over the full domain.
 *
 * @tparam T      type of the result.
 * @tparam TField type of the implicit function (must be iterable).
 * @param field the implicit function.
 */
template < typename T, typename TField >
inline
T getVolume( TField const& field )
{
  return std::accumulate( std::begin(field), std::end(field), T(0) );
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
  typename ImageContainerBySTLVector<typename TImage::Domain,double>::Range::Iterator 
    out = outputImage.range().begin(); 
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

  out = outputImage.range().begin(); 
  for ( typename DT2::ConstRange::ConstIterator 
	  it = dt2.constRange().begin(), itend = dt2.constRange().end();
	it != itend; ++it)
    {
      double dist = (*it); 
      if (dist != 0) 
	{ //only write if dist != 0
	  dist = -std::sqrt( dist ) + 0.5;  
	  *out++ = dist; 
	  // outputImage.setValue(p, dist);  
	}
      else 
	out++; 
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

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

