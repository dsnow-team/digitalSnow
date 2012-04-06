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
 * @file PartitionEvolver.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/04/06
 *
 * Header file for module PartitionEvolver.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(PartitionEvolver_RECURSES)
#error Recursive header files inclusion detected in PartitionEvolver.h
#else // defined(PartitionEvolver_RECURSES)
/** Prevents recursive inclusion of headers. */
#define PartitionEvolver_RECURSES

#if !defined PartitionEvolver_h
/** Prevents repeated inclusion of headers. */
#define PartitionEvolver_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/images/CImage.h"
#include "DGtal/kernel/CPointFunctor.h"
#include "DGtal/kernel/CPointPredicate.h"

#include "DGtal/base/CowPtr.h"
#include "DGtal/base/CountedPtr.h"

#include "LocalMCM.h"
#include "FrontierEvolver.h"

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class PartitionEvolver
  /**
   * Description of template class 'PartitionEvolver' <p>
   * \brief Aim: This class is a way of deforming an image of labels
   * around all the connected contact surfaces between two adjacent regions, 
   * according to a speed field, whose computation is delegated to a point functor. 
   *   
   *
   * @tparam TKSpace a model of CCellularGridSpaceND
   * @tparam TLabelImage a model of CImage (storing labels)
   * @tparam TDistanceImage a model of CImage (storing distance values)
   * @tparam TExternImage a model of CImage (storing an extern scalar field)
   * @tparam TTopoPredicate
   */
  template <typename TKSpace, 
	    typename TLabelImage, typename TDistanceImage, 
	    typename TExternImage, typename TTopoPredicate>
  class PartitionEvolver
  {


    BOOST_CONCEPT_ASSERT(( CImage<TLabelImage> )); 
    BOOST_CONCEPT_ASSERT(( CImage<TDistanceImage> )); 
    BOOST_STATIC_ASSERT
    (( ConceptUtils::SameType< typename TKSpace::Point,
       typename TLabelImage::Point>::value ));
    BOOST_STATIC_ASSERT
    (( ConceptUtils::SameType< typename TKSpace::Point,
       typename TDistanceImage::Point>::value ));

    BOOST_CONCEPT_ASSERT(( CImage<TExternImage> )); 
    BOOST_STATIC_ASSERT
    (( ConceptUtils::SameType< typename TKSpace::Point,
       typename TExternImage::Point>::value ));

    // BOOST_CONCEPT_ASSERT(( CPoinTTopoPredicate<TTopoPredicate> )); 
    // BOOST_STATIC_ASSERT
    // (( ConceptUtils::SameType< typename TKSpace::Point,
    //    typename TTopoPredicate::Point>::value ));
    //TODO testing TTopoPredicate as a binary TopoPredicate on points and labels

    // ----------------------- Types ------------------------------
  public:

    /// Khalimsky space
    typedef TKSpace KSpace;
    typedef typename KSpace::Cell Cell;
    typedef typename KSpace::SCell SCell;

    /// Image of labels
    typedef TLabelImage LabelImage;
    typedef typename LabelImage::Value Label;
    typedef typename LabelImage::Domain Domain;
    typedef typename Domain::Point Point;

    /// Images of distance values
    typedef TDistanceImage DistanceImage;
    typedef CountedPtr<DistanceImage> DistanceImagePtr;

    /// Extern image
    typedef TExternImage ExternImage;

    /// Point functor for the mapping points-speed
    typedef LocalMCM<DistanceImage, DistanceImage > Functor;
    typedef CountedPtr<Functor> FunctorPtr; 

    /// Topological predicate 
    typedef TTopoPredicate TopoPredicate;

    /// Frontier evolver
    typedef FrontierEvolver<KSpace, LabelImage, DistanceImage, Functor, TopoPredicate > Evolver; 
    typedef CountedPtr<Evolver> FrontierEvolverPtr; 
 

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param aK khalimsky space where the digital frontiers are defined 
     * @param aI an image of labels
     * @param aF any extern data field 
     * @param aP any point TopoPredicate
     */
    PartitionEvolver(const KSpace& aK, LabelImage& aI, const ExternImage& aF, const TopoPredicate& aP);

    /**
     * Destructor. Does nothing.
     */
    ~PartitionEvolver();


    /**
     * Deform the image of labels during @a aT
     *
     * @param aT time step
     * @return time spent during the deformation 
     * (equal to aT). 
     */
    double update(const double& aT);

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


    // ------------------------- References --------------------------------
    /**
     * Constant reference on the khalimsky space
     */
    const KSpace& myKSpace; 
    /**
     * Reference on the image of labels
     */
    LabelImage& myLabelImage; 
    /**
     * Constant reference on the extern image
     */
    const ExternImage& myExternImage; 
    /**
     * Constant reference on the topological predicate
     */
    const TopoPredicate& myTopoPred; 


    // ------------------------- Data --------------------------------
    /**
     * Set of smart owning pointers on images (of distance values)
     */
    std::vector<DistanceImagePtr> myImages; 
    /**
     * Set of smart owning pointers on functors, 
     * which locally computes the differential estimations
     * and displacement speed
     */
    std::vector<FunctorPtr> myFunctors; 
    /**
     * Set of smart owning pointers on frontier evolvers
     */
    std::vector<FrontierEvolverPtr> myEvolvers; 

    // ------------------------- Hidden services ------------------------------
  protected:


    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    PartitionEvolver ( const PartitionEvolver & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    PartitionEvolver & operator= ( const PartitionEvolver & other );

  private:



    // ------------------------- Internals ------------------------------------
  private:

    /**
     * Insert bels, ie cells lying between two points
     * of different labels in @a aImg, 
     * in @a aSet
     * @param aSet (returned) set of bels.
     * @param aKSpace khalimsky space.
     * @param aImg label image.
     * @param aLowerBound.
     * @param aUpperBound.
     */
    void getBels ( std::set<Cell>& aSet, 
		   const KSpace & aKSpace,
		   const LabelImage & aImg,
		   const Point & aLowerBound, 
		   const Point & aUpperBound  );

  }; // end of class PartitionEvolver


  /**
   * Overloads 'operator<<' for displaying objects of class 'PartitionEvolver'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'PartitionEvolver' to write.
   * @return the output stream after the writing.
   */
  template <typename TKSpace, typename TLabelImage, typename TDistanceImage, 
	    typename TExternImage, typename TTopoPredicate>
  std::ostream&
  operator<< ( std::ostream & out, 
	       const PartitionEvolver<TKSpace, TLabelImage, TDistanceImage, 
	       TExternImage, TTopoPredicate> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "PartitionEvolver.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined PartitionEvolver_h

#undef PartitionEvolver_RECURSES
#endif // else defined(PartitionEvolver_RECURSES)
