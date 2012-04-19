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

#pragma once

/**
 * @file FrontierEvolver.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/03/01
 *
 * Header file for module FrontierEvolver.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(FrontierEvolver_RECURSES)
#error Recursive header files inclusion detected in FrontierEvolver.h
#else // defined(FrontierEvolver_RECURSES)
/** Prevents recursive inclusion of headers. */
#define FrontierEvolver_RECURSES

#if !defined FrontierEvolver_h
/** Prevents repeated inclusion of headers. */
#define FrontierEvolver_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/images/CImage.h"
#include "DGtal/kernel/CPointFunctor.h"
#include "DGtal/kernel/CPointPredicate.h"

//predicates
#include "DGtal/base/BasicBoolFunctions.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "PointPredicates.h"



// set
#include "DGtal/kernel/sets/DigitalSetFromMap.h"
#include "DGtal/kernel/sets/DigitalSetBySTLSet.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/images/ImageHelper.h"

// FMM
#include "FMM.h"

// frontier
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/helpers/FrontierPredicate.h"
#include "DGtal/topology/LightExplicitDigitalSurface.h"

//partition
#include "PartitionEvolver.h"

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{


  //------------------------------------------------------------------------------
  template<typename I, typename D, typename V>
struct SetFromImageDomainValueTraits
  { 
  public: 
    typedef DigitalSetBySTLSet<D> Set; 
  public: 
  public: 
    static Set get(I& aImage)
    {
      return Set(aImage.domain()); 
    }
  }; 
  //Partial specialization
  template<typename D, typename V>
struct SetFromImageDomainValueTraits<
    ImageContainerBySTLMap<D,V>, 
    D, V >
  {
  public: 
    typedef DigitalSetFromMap<ImageContainerBySTLMap<D,V> > Set; 
  public: 
    static Set get(ImageContainerBySTLMap<D,V>& aImage)
    {
      return Set(aImage); 
    }
  }; 
  //------------------------------------------------------------------------------
  template<typename I>
struct SetFromImageSelector
  { 
  public: 
    BOOST_CONCEPT_ASSERT(( CImage<I> )); 
    typedef typename I::Domain Domain; 
    typedef typename I::Value Value; 

    typedef typename SetFromImageDomainValueTraits<I, Domain, Value>::Set Set;

  public: 
    static Set get(I& i) 
    {
      return SetFromImageDomainValueTraits<I, Domain, Value>::get(i); 
    }
  }; 

  //------------------------------------------------------------------------------
  namespace details
  {
    class CompareSecondElement {
    public: 
      template <typename T>
      bool operator()(const T& a, const T& b) 
      {
	return ( a.second < b.second ); 
      }
    };

  }


  /////////////////////////////////////////////////////////////////////////////
  // template class FrontierEvolver
  /**
   * Description of template class 'FrontierEvolver' <p>
   * \brief Aim: This class is a way of deforming an image of labels
   * around a connected contact surface between two regions, 
   * according to a speed field, whose computation is 
   * delegated to a point functor. 
   *   
   * At each step, a implicit function is extended from the known
   * values at the boundary points adjacent to the contact surface. 
   * The points lying around the interface are sorted according to
   * their time of zero-crossing (ie. their distance to the interface 
   * divided by their speed) so that they are flipped from one region
   * to another, one by one and in order, until a time greater than 
   * a given threshold is reached, but only if a given predicate 
   * based on topological properties returns true. 
   *
   * @tparam TKSpace a model of CCellularGridSpaceND
   * @tparam TLabelImage a model of CImage (storing labels)
   * @tparam TDistanceImage a model of CImage (storing distance values)
   * @tparam TFunctor a model of CPointFunctor
   * @tparam TTopoPredicate a model of topological predicate
   */
  template <typename TKSpace, 
	    typename TLabelImage, typename TDistanceImage, 
	    typename TFunctor, typename TTopoPredicate>
  class FrontierEvolver
  {


    BOOST_CONCEPT_ASSERT(( CImage<TLabelImage> )); 
    BOOST_CONCEPT_ASSERT(( CImage<TDistanceImage> )); 
    BOOST_STATIC_ASSERT
    (( ConceptUtils::SameType< typename TKSpace::Point,
       typename TLabelImage::Point>::value ));
    BOOST_STATIC_ASSERT
    (( ConceptUtils::SameType< typename TKSpace::Point,
       typename TDistanceImage::Point>::value ));

    BOOST_CONCEPT_ASSERT(( CPointFunctor<TFunctor> )); 
    BOOST_STATIC_ASSERT
    (( ConceptUtils::SameType< typename TKSpace::Point,
       typename TFunctor::Point>::value ));

    // BOOST_CONCEPT_ASSERT(( CPointPredicate<TTopoPredicate> )); 
    // BOOST_STATIC_ASSERT
    // (( ConceptUtils::SameType< typename TKSpace::Point,
    //    typename TTopoPredicate::Point>::value ));
    //TODO testing TTopoPredicate as a binary predicate on points and labels

    // ----------------------- Types ------------------------------
  public:

    /// Image of labels
    typedef TLabelImage LImage;
    typedef typename LImage::Value Label;
    typedef typename LImage::Domain Domain;
    typedef typename Domain::Point Point;

    /// Image of distance values
    typedef TDistanceImage DImage;
    typedef typename DImage::Value Distance;

    /// Point functor for the mapping points-speed
    typedef TFunctor Functor; 
    typedef typename Functor::Value Speed; 
    typedef std::pair<Distance, Speed> DistanceSpeed; 
    typedef std::pair<Point,double> PointTime; 

    /// Topological predicate 
    typedef TTopoPredicate TopoPredicate;

    /// Set of points where the distance values are known
    typedef typename SetFromImageSelector<DImage>::Set PointSet; 

    /// Khalimsky space
    typedef TKSpace KSpace;

    /// Frontier
    typedef FrontierPredicate<KSpace, LImage> SurfelPredicate;
    typedef LightExplicitDigitalSurface<KSpace, SurfelPredicate> Frontier;
    /// Surfel 
    typedef typename Frontier::Surfel Surfel;
    typedef typename Frontier::SurfelConstIterator SurfelIterator;


    /// Partition
    typedef PartitionEvolver<KSpace, LImage, DImage, 
			     typename Functor::ExternField, 
			     TopoPredicate> Partition; 

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param aK khalimsky space where the digital frontier is defined 
     * @param aI an image of labels
     * @param aD an image of distance values
     * @param aS a surfel lying between two regions of @a aI
     * @param aF a point functor mapping a speed to points 
     * @param aP any topological predicate
     */
    FrontierEvolver(const KSpace& aK, LImage& aI, DImage& aD, const Surfel& aS, 
		    const Functor& aF, const TopoPredicate& aP, 
		    Partition* aPartitionPtr = NULL);

    /**
     * Destructor. Does nothing.
     */
    ~FrontierEvolver();

     /**
     * Deform the image of labels during @a aT
     *
     * @param aT time step
     * @return time spent during the deformation 
     * (equal to aT). 
     */
    double update(const double& aT);


    /**
     * @return starting surfel of the digital frontier.
     */
    Surfel surfel() const;

    /**
     * @param aSurfel new starting surfel of the digital frontier.
     */
    void setSurfel(const Surfel& aSurfel);

    /**
     * @return begin iterator on the surfels of the digital frontier.
     */
    SurfelIterator begin() const;

    /**
     * @return end iterator on the surfels of the digital frontier.
     */
    SurfelIterator end() const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;


    // ------------------------- Protected Datas ------------------------------
  protected:
    // ------------------------- Private Datas --------------------------------
  private:

    /**
     * Constant reference on the khalimsky space
     */
    const KSpace& myKSpace; 
    /**
     * Reference on the image of labels
     */
    LImage& myLImage; 
    /**
     * Reference on the image of distance values
     */
    DImage& myDImage; 
    /**
     * Set of points where the distance values are known
     */
    PointSet myPointSet; 
    /**
     * Starting surfel of the digital frontier 
     */
    Surfel mySurfel; 
    /**
     * Constant reference on the functor
     */
    const Functor& myFunctor; 
    /**
     * Constant reference on the topological predicate
     */
    const TopoPredicate& myTopoPred; 
    /**
     * Label of the inner region 
     */
    Label myInnerLabel; 
    /**
     * Label of the outer region 
     */
    Label myOuterLabel; 
    /**
     * Surfel predicate telling whether a given surfel
     * belongs to the frontier or not
     */
    SurfelPredicate mySurfelPred; 
    /**
     * (implicit) digital frontier
     */
    const Frontier* myFrontier; 
    /**
     * Aliasing pointer on the partition the frontier belongs to
     */
    Partition* myPartitionPtr; 

    // ------------------------- Hidden services ------------------------------
  protected:


    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    FrontierEvolver ( const FrontierEvolver & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    FrontierEvolver & operator= ( const FrontierEvolver & other );

  private:



    // ------------------------- Internals ------------------------------------
  private:

     /**
     * Return through @a res the points
     * of the narrow band
     * for which the distance value has been
     * computed and stored in @a myDImage
     *
     * @tparam TOutputIterator a model of output iterator
     *
     * @param res an output iterator
     */
    template <typename TOutputIterator>
    void initNarrowBand(TOutputIterator res);

     /**
     * Return through @a res the points
     * that can be flipped into the 
     * adjacent region
     *
     * @param itb begin iterator on points
     * @param ite end iterator on points
     * @param aDistanceSpeedIto an output iterator on DistanceSpeed pairs
     * @param aCandidateIto an output iterator on candidates
     *
     * @tparam TInputIterator a model of input iterator on points
     * @tparam TOutputIterator1 a model of output iterator on DistanceSpeed pairs
     * @tparam TOutputIterator2 a model of output iterator on candidates
     *
     */
    template <typename TInputIterator, 
	      typename TOutputIterator1, typename TOutputIterator2>
    void initCandidates( const TInputIterator& itb, const TInputIterator& ite, 
			 TOutputIterator1 aDistanceSpeedIto, TOutputIterator2 aCandidateIto );


    /**
     * Checks whether the other frontiers (if any)
     * can be modified by the flip of @a aPoint
      * @param aPoint any (digital) point
     */
    void checkPartition ( const Point& aPoint );

     /**
     * Update @a myLImage for each candidates of the 
     * range [@a itb , @a ite ) 
     *
     * The points that are not allowed to flip
     * (because of the topological predicate) 
     * are stored through @a ito into a container
     *
     * @param itb begin iterator on candidates
     * @param ite end iterator on candidates
     * @param ito output iterator on points
     * @param aP (returned) last point
     * @param aT (returned) zero-crossing time of the last point
     * @param aTMax maximal accepted time (equal to the evolution time step)
     *
     * @tparam TInputIterator a model of input iterator on candidates
     * @tparam TOutputIterator a model of input iterator on points
     *
     * @return number of flipped points
     */
    template <typename TInputIterator, typename TOutputIterator>
    int updateLabelImage( const TInputIterator& itb, const TInputIterator& ite,
			  TOutputIterator ito, 
			  Point& aP, double &aT, const double& aTMax );

     /**
     * Update @a myDImage at each point of the 
     * range [@a itb , @a ite ) knowing their
     * distance and speed returning by @a aDistanceSpeedIt
     *
     * The points that belongs to @a aSet are not allowed to flip. 
     *
     * @param itb begin iterator on points
     * @param ite end iterator on points
     * @param aDistanceSpeedIt an input iterator on DistanceSpeed pairs
     * @param aSet any set of points
     * @param t evolution time step
     *
     * @tparam TInputIterator1 a model of input iterator on points
     * @tparam TInputIterator2 a model of input iterator on DistanceSpeed pairs
     * @tparam TSet STL set like type 
     *
     */
    template <typename TInputIterator1, typename TInputIterator2, typename TSet>
    void updateDistanceImage( const TInputIterator1& itb, const TInputIterator1& ite,
			      TInputIterator2 aDistanceSpeedIt, 
			      const TSet& aSet, const double& t );


    /**
     * Update the starting surfel @a mySurfel
     * of the digital frontier from point @a p
     * @param p any (digital) point
     */
    void updateFrontier ( const Point& p );

    /**
     * Get inner point.
     * @param s a surfel
     * @return the inner point
     */
    Point getInnerPoint ( const Surfel& s ) const ;

    /**
     * Get outer point.
     * @param s a surfel
     * @return the outer point
     */
    Point getOuterPoint ( const Surfel& s ) const ;


  }; // end of class FrontierEvolver


  /**
   * Overloads 'operator<<' for displaying objects of class 'FrontierEvolver'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'FrontierEvolver' to write.
   * @return the output stream after the writing.
   */
  template <typename TKSpace, typename TLabelImage, typename TDistanceImage, 
	    typename TFunctor, typename TTopoPredicate>
  std::ostream&
  operator<< ( std::ostream & out, 
	       const FrontierEvolver<TKSpace, TLabelImage, TDistanceImage, 
	       TFunctor, TTopoPredicate> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "FrontierEvolver.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined FrontierEvolver_h

#undef FrontierEvolver_RECURSES
#endif // else defined(FrontierEvolver_RECURSES)
