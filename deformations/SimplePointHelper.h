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
 * @date 2011/11/03
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
#include "DGtal/images/ImageContainerBySTLVector.h"

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  #include "SimplePointHelperDetails.ih"

  /////////////////////////////////////////////////////////////////////////////
  // template class SimplePointHelper
  /**
   * \brief Aim: Tests whether a given voxel v
   * is ML-simple with respect to a given region.
   *
   * Two voxels u and w are r-adjacent (0 <= r <= dim) iff 
   * exactly (dim-r) coordinates are equal and r differ by one. 
   * The r-neighborhood of a voxel u is the set of voxels that are r-adjacent to u. 
   * The <r-neighborhood of a voxel u is the set of voxels that are k-adjacent to u
   * for all k between 0 and r. 
   *
   *
   * A voxel v of the region V is ML-simple with respect to the region R if:
   * -1) v is simple for V 
   * -2) v is simple for R U v
   * -3) for each region O (distinct from V and R), 
   * that is (dim-1)-adjacent to v: 
   *  -a) v is simple for O U V 
   *  -b) v is simple for O U R U v
   *
   * The notion of simplicity used here implies the simplicity of Bertrand based on 
   * the (dim-1,dim-2)- of (dim-1,dim-3)- adjacency pair. The reverse is not true, 
   * because some configurations implying the disappearance of folds on the surface 
   * are rejected. 
   *
   * In order to test if a point is ML-simple for a given label, 
   * it is enough to call method isMLSimple(): 
   *
   * @code 

   DGtal::SimplePointHelper<TImage> h(anImage); 
   trace.info() << h.isMLsimple(aPoint, aLabel) << endl;

   * @endcode
   *
   * It is also possible to test if a point is simple: 
   *
   * @code 

   DGtal::SimplePointHelper<TImage> h(anImage); 
   trace.info() << h.isMLsimple(aPoint) << endl;

   * @endcode
   *
   * NB: The method isSimple() returns 'true' if the center voxel v 
   * is simple for the region V it belongs to (and the voxels that
   * does not belong to V are assumed to belong to the complementary set). 
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

    typedef typename TImage::Domain Domain; 
    typedef typename Domain::Point Point; 
    typedef typename Domain::Point Vector;
    typedef typename Point::Coordinate Coordinate;  

  private: 
    
    typedef ImageContainerBySTLVector<Domain, bool> BinaryImage; 

    // ----------------------- Standard services ------------------------------
  public:


    /**
     * Main constructor.
     *
     * @param aImg any image 
     * @param aLabel label of the infinite region 
     * located outside the image domain 
     * (default value if not provided)
     */
    SimplePointHelper(LabelsImage& aImg, const Label& aLabel);

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
     * Checks if @a aPoint is simple with respect 
     * to the region @a aPoint belongs to.
     *
     * @param aPoint any point
     * @return 'true' if simple, 'false' otherwise.
     */
    bool isSimple(const Point& aPoint) const;

    /**
     * Checks if the underlying image is labeled, 
     * ie. sets of voxels of same label are (dim-1)-connected.
     *
     * @return 'true' if the image is labeled, 'false' otherwise.
     *
     * NB: trivial algorithm in O(|D|.|L|), 
     * where |D| is the size of the domain and 
     * |L| is the number of distinct labels in the image.
     */
    bool isLabeled() const;


    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 
     * ie. the underlying image is labeled, 
     * 'false' otherwise.
     *
     * @see isLabeled
     */
    bool isValid() const;


    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;


    // ------------------------- Static methods ------------------------------

    /**
     * Read labels from the range [ @a itb , @a ite ) for filling @a aImg.
     *
     * @param aImg labeled image to fill
     * @param itb begin iterator
     * @param ite end iterator
     * @tparam TIterator  type of iterator on labels.
     */
    template <typename TIterator>
    static void readConfiguration ( LabelsImage& aImg, 
				    const TIterator& itb, const TIterator& ite );

    /**
     * Randomly fill the image @a aImg with the labels contained in @a labels
     * following a region growing strategy from the center voxel.  
     *
     * NB: The number of distinct adjacent regions to the center voxel
     * depends of the number of distinct labels, but also depends on 
     * the probability @a prob of filling the neighbors of a given voxel
     * (the less the probability is, the more there are distincts regions). 
     *
     * @param aImg labeled image to fill
     * @param labels set of distinct labels
     * @param prob probability of growing from a given voxel to its neighbor
     * @return 'true' is the computed configuration is valid, 
     * 'false' otherwise
     */
    static bool generateRandomConfiguration ( LabelsImage& aImg, const std::set<Label>& labels, 
					      const double& prob );

    // ------------------------- Protected Datas ------------------------------
  protected:
    // ------------------------- Private Datas --------------------------------
  private:


    /**
     * Reference on an image of labels.
     */
    LabelsImage& myImg; 
    
    /**
     * Label of the infinite region located outside the image domain.
     */
    Label myLabel; 

    // ------------------------- Hidden services ------------------------------
  protected:



  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    SimplePointHelper ( const SimplePointHelper & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    SimplePointHelper & operator= ( const SimplePointHelper & other );






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
    bool isSimple12(const Point& aPoint, const TPredicate& aPredicate) const;

    /**
     * Checks if @a aPoint, assumed to belong to region X, 
     * is simple for X, where X is defined as the 
     * set of voxels that make @a aPredicate be true.
     *
     * NB: implies the simplicity with the (dim-1,dim-3) adjacency
     *
     * @param aPoint any point to check
     * @tparam TPredicate  type of predicate
     * @param aPredicate  any predicate
     * @return 'true' if simple, 'false' otherwise.
     */
    template <typename TPredicate>
    bool isSimple13(const Point& aPoint, const TPredicate& aPredicate) const;

    /**
     * Checks if @a aPoint, assumed to belong to region X, 
     * is simple for X, where X is defined as the 
     * set of voxels that make @a aPredicate be true.
     *
     * NB: implies the simplicity with the (dim-1,dim-2) adjacency,
     * plus reject disappearance of folds on the surface of X. 
     *
     * @param aPoint any point to check
     * @tparam TPredicate  type of predicate
     * @param aPredicate  any predicate
     * @return 'true' if simple, 'false' otherwise.
     */
    template <typename TPredicate>
    bool isSimple12WithoutFolds(const Point& aPoint, const TPredicate& aPredicate) const;

    /**
     * Checks if @a aPoint, assumed to belong to region X, 
     * is simple for X, where X is defined as the 
     * set of voxels that make @a aPredicate be true.
     *
     * NB: implies the simplicity with the (dim-1,dim-3) adjacency
     * plus reject disappearance of folds on the surface of X. 
     *
     * @param aPoint any point to check
     * @tparam TPredicate  type of predicate
     * @param aPredicate  any predicate
     * @return 'true' if simple, 'false' otherwise.
     */
    template <typename TPredicate>
    bool isSimple13WithoutFolds(const Point& aPoint, const TPredicate& aPredicate) const;

    // ------------------------- small internal helpers ------------------------------------    


    /**
     * Tests wether the set of points that are marked true in @a aImg
     * is (dim-1)-connected or not. 
     *
     * @param aImg  any binary image.
     * @param aN  the total number of points marked as true.
     * @param aStartingPoint  any point marked as true.
     * @param aNeighborhoodFunctor  a functor returning the neighborhood of a point. 
     * @tparam TFunctor type of a neighborhood functor
     *
     * @return true if connected, false otherwise
     */
    template <typename TFunctor>
    static bool isConnected( BinaryImage& aImg, unsigned int& aN, const Point& aStartingPoint, 
		      const TFunctor& aNeighborhoodFunctor);


    /**
     * Return the hyper-rectangular domain around @aPoint
     * corresponding to its <r-neighborhood
     * @param aPoint any point
     * @return  the domain
     */
    static Domain getLocalDomain(const Point& aPoint);


    /**
     * Return a random number between 0 and @a aN
     * @param aN any number
     * @return generated random number
     *
     * @see generateRandomConfiguration
     */
    static unsigned int randomInt(const unsigned int& aN);
  
  
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
