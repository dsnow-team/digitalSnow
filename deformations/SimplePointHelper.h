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
 * @file SimplePointHelper.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en LabelsImage et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/04/05
 *
 * @brief Header file for module SimplePointHelper.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(SimplePointHelper_RECURSES)
#error Recursive header files inclusion detected in SimplePointHelper.h
#else // defined(SimplePointHelper_RECURSES)
/** Prevents recursive inclusion of headers. */
#define SimplePointHelper_RECURSES

#if !defined SimplePointHelper_h
/** Prevents repeated inclusion of headers. */
#define SimplePointHelper_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions

#include <iostream>
#include <vector>
#include <queue>
#include <set>
#include <functional>
#include <cstdlib>

#include "DGtal/base/Common.h"
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/base/CLabel.h"

//simplicity
#include "DGtal/topology/MetricAdjacency.h"
#include "DGtal/topology/DomainMetricAdjacency.h"
#include "DGtal/topology/DomainAdjacency.h"
#include "DGtal/topology/DigitalTopology.h"
#include "DGtal/topology/Object.h"

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

#include "SimplePointHelperDetails.ih"

  /////////////////////////////////////////////////////////////////////////////
  // template class SimplePointHelper
  /**
   * \brief Aim: Tests whether a given voxel v,
   * which belongs to a region V,
   * is ML-simple with respect to a given region R.
   *
   * In short, v can belong to R if it is a simple point 
   * for all groups of up to three regions including either V or R. 
   *
   * [Bertrand, 1994] A voxel v is simple for a region X 
   * iff #C6 [G6(v, X)] = #C18[G18(v, X^c)] = 1, 
   * where #Ck [Y] denotes the number of k-connected components of a set Y.
   *
   * [Bazin et. al., 2007] A voxel v is ML-simple for R iff: 
   * -1 v is simple for V 
   * -2 v is simple for R U v
   * -3 for each region O (distinct from V and R) 
   * that is adjacent to v: 
   *  -a v is simple for O U V 
   *  -b v is simple for O U R U v
   * -4 for each region Q (distinct from V and R) 
   * that is adjacent to v: 
   *  -a for all O, v is simple for Q U O U V
   *  -b for all O, v is simple for Q U O U V
   *
   *
   * In order to test if a point is ML-simple for a given label, 
   * it is enough to call the `operator()` or the `isMLSimple` method: 
   *
   * @code 

   DGtal::SimplePointHelper<TImage> h(anImage); 
   trace.info() << h.isMLsimple(aPoint, aLabel) << endl;
   trace.info() << h(aPoint, aLabel) << endl;

   * @endcode
   *
   * @tparam TImage  type of image.
   */
  template <typename TImage>
  class SimplePointHelper
  {

    BOOST_STATIC_ASSERT( (TImage::dimension<=3)&&(TImage::dimension>=1) );

    // ----------------------- Types ------------------------------
  public:
    

    typedef TImage LabelsImage; 
    typedef typename TImage::Value Label; 
    BOOST_CONCEPT_ASSERT ((CLabel<Label>));

    typedef typename TImage::Domain::Space Space; 
    typedef typename TImage::Domain Domain; 
    typedef typename Domain::Point Point; 
    typedef typename Domain::Point Vector;
    typedef typename Point::Coordinate Coordinate;  

    /// domain should be rectangular
    BOOST_STATIC_ASSERT ((boost::is_same< Domain, 
			  HyperRectDomain<SpaceND<Domain::dimension, Coordinate> > >::value));


    // ----------------------- Standard services ------------------------------
  public:


    /**
     * Main constructor.
     *
     * @param aImg any image 
     * @param aInftyRegionLabel label of the infinite region 
     * located outside the image domain 
     * (default value if not provided)
     */
    SimplePointHelper(LabelsImage& aImg, 
		      const Label& aInftyRegionLabel);

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    SimplePointHelper ( const SimplePointHelper & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    SimplePointHelper & operator= ( const SimplePointHelper & other );

    /**
     * Destructor. Does nothing.
     */
    ~SimplePointHelper();


    /**
     * Checks if @a aPoint is ML-simple with respect 
     * to the region having label @a aLabel.
     *
     * @param aPoint any point
     * @param aLabel  any label
     * @return 'true' if ML-simple, 'false' otherwise.
     */
    bool isMLSimple(const Point& aPoint, const Label& aLabel) const;

    /**
     * Checks if @a aPoint is ML-simple with respect 
     * to the region having label @a aLabel.
     *
     * @param aPoint any point
     * @param aLabel  any label
     * @return 'true' if ML-simple, 'false' otherwise.
     *
     * @see isMLSimple
     */
    bool operator()(const Point& aPoint, const Label& aLabel) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 
     */
    bool isValid() const;


    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;


    // ------------------------- Static methods ------------------------------

    // ------------------------- Protected Datas ------------------------------
  protected:
    // ------------------------- Private Datas --------------------------------
  private:


    /**
     * Aliasing pointer on an image of labels.
     */
    LabelsImage* myImgPtr; 

    /**
     * Label of the infinite region located outside the image domain.
     */
    const Label myInftyRegionLabel; 

    // ------------------------- Hidden services ------------------------------
  protected:



  private:


    // ------------------------- Internals ------------------------------------
  private:
    
    // ------------------------- Main algorithms ------------------------------------


    /**
     * Checks if @a aPoint is ML-simple 
     * for the region of label @a aLabel
     *
     * @param aImg any image
     * @param aPoint any point of the image domain
     * @param aLabel label to assign to @a aPoint
     * @return 'true' if ML-simple, 'false' otherwise.
     */
    bool isMLSimple(LabelsImage& aImg, const Point& aPoint, const Label& aLabel) const;

    /**
     * Checks if @a aPoint, assumed to belong to region X, 
     * is simple for X, where X is defined as the 
     * set of voxels that make @a aPredicate be true.
     *
     * NB: implies the simplicity with the (dim-1,dim-2) adjacency
     *
     * @param aPoint any point to check
     * @tparam TPredicate  type of predicate
     * @param aPredicate  any predicate
     * @return 'true' if simple, 'false' otherwise.
     */
    template <typename TPredicate>
    bool isSimple(const Point& aPoint, const TPredicate& aPredicate) const;



    /**
     * Return the hyper-rectangular domain around @aPoint
     * corresponding to its <r-neighborhood
     * @param aPoint any point
     * @return  the domain
     */
    static Domain getLocalDomain(const Point& aPoint);

  
  
  }; // end of class SimplePointHelper


  /**
   * Overloads 'operator<<' for displaying objects of class 'SimplePointHelper'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'SimplePointHelper' to write.
   * @return the output stream after the writing.
   */
  template <typename TImage>
  std::ostream&
  operator<< ( std::ostream & out, const SimplePointHelper<TImage> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "SimplePointHelper.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined SimplePointHelper_h

#undef SimplePointHelper_RECURSES
#endif // else defined(SimplePointHelper_RECURSES)
