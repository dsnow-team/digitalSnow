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
   * according to a velocity field, whose computation is 
   * delegated to a point functor. 
   *   
   * At each step, a signed distance function is built. 
   * The points are sorted according to their time of zero-crossing
   * (ie. their distance to the interface divided by their velocity)
   * so that they are flipped from a region to another, one by one and 
   * in order, until a time greater than a threshold is reached or 
   * until a point predicate (possibly based on topological properties)
   * returns false. 
   *
   * @tparam TKSpace a model of CCellularGridSpaceND
   * @tparam TLabelImage a model of CImage (storing labels)
   * @tparam TDistanceImage a model of CImage (storing distance values)
   * @tparam TFunctor a model of CPointFunctor
   * @tparam TPredicate a model of CPointPredicate
   */
  template <typename TKSpace, 
	    typename TLabelImage, typename TDistanceImage, 
	    typename TFunctor, typename TPredicate>
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

    // BOOST_CONCEPT_ASSERT(( CPointPredicate<TPredicate> )); 
    // BOOST_STATIC_ASSERT
    // (( ConceptUtils::SameType< typename TKSpace::Point,
    //    typename TPredicate::Point>::value ));
    //TODO testing TPredicate as a binary predicate on points and labels

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


    /// Point functor for the mapping points-velocity
    typedef TFunctor Functor; 
    typedef typename Functor::Value Velocity; 
    /// Topological predicate 
    typedef TPredicate Predicate;
 

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param aK khalimsky space where the digital frontier is defined 
     * @param aI an image of labels
     * @param aD an image of distance values
     * @param aS a surfel lying between two regions of @a aI
     * @param aF a point functor mapping a velocity to points 
     * @param aP any point predicate
     * @param aW maximal width of the deformation band (1.0 by default)
     */
    FrontierEvolver(const KSpace& aK, LImage& aI, DImage& aD, Surfel& aS, 
		    const Functor& aF, const Predicate& aP, const double& aW = 1.0);

    /**
     * Destructor. Does nothing.
     */
    ~FrontierEvolver();

    //  /**
    //  * Deform the image of labels around the digital frontier
    //  *
    //  * @return time spent during the deformation. 
    //  */
    // double update();

     /**
     * Deform the image of labels during @a aT
     *
     * @param aT time step
     * @return time spent during the deformation 
     * (equal to aT). 
     */
    double update(const double& aT);

     /**
     * Return through @a out the points
     * for which the distance value has been
     * computed and stored in @a myDImage ,
     * which are candidate to the flip
     *
     * @tparam TOutputIterator a model of output iterator
     *
     * @param out an output iterator
     */
    template <typename TOutputIterator>
    void init(const TOutputIterator& out);

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
     * Constant reference on the predicate
     * TODO: rename it into myTopoPred
     */
    const Predicate& myPointPred; 
    /**
     * Maximal width of the band whithin which
     * the frontier is moving
     */
    double myW; 
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

    /**
     * Update the starting surfel @a mySurfel
     * of the digital frontier from point @a p
     * @param p any (digital) point
     */
    void updateFrontier ( const Point& p );

  }; // end of class FrontierEvolver


  /**
   * Overloads 'operator<<' for displaying objects of class 'FrontierEvolver'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'FrontierEvolver' to write.
   * @return the output stream after the writing.
   */
  template <typename TKSpace, typename TLabelImage, typename TDistanceImage, 
	    typename TFunctor, typename TPredicate>
  std::ostream&
  operator<< ( std::ostream & out, 
	       const FrontierEvolver<TKSpace, TLabelImage, TDistanceImage, 
	       TFunctor, TPredicate> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "FrontierEvolver.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined FrontierEvolver_h

#undef FrontierEvolver_RECURSES
#endif // else defined(FrontierEvolver_RECURSES)
