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
 * @file PointPredicates.h
 *
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 * @date 2012/02/02
 *
 * This files contains several basic classes representing binary predicates
 *
 * This file is part of the DGtal library.
 */

#if defined(PointPredicates_RECURSES)
#error Recursive header files inclusion detected in PointPredicates.h
#else // defined(PointPredicates_RECURSES)
/** Prevents recursive inclusion of headers. */
#define PointPredicates_RECURSES

#if !defined PointPredicates_h
/** Prevents repeated inclusion of headers. */
#define PointPredicates_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class OneLabelPredicate
  /**
   * \brief Aim: Predicate returning true or false
   * according to the comparison of a label with a reference label. 
   *
   * @tparam TLabel  type of labels 
   (must be default-constructible and equally-comparable)
   * @tparam TBinaryPredicate  type of binary predicate
   (must have a () operator taking two input parameters
   and returning a boolean). 
   *
   */
  template <typename TImage, typename TBinaryPredicate>
  class OneLabelPredicate
  {
  public: 
    typedef typename TImage::Value Value; 
    typedef typename TImage::Point Point; 

  public: 
    /**
     * Constructor 
     * @param aImg  any image
     * @param aLabel any label
     * @param aFunctor  any binary functor
     */
    OneLabelPredicate(const TImage& aImg, const Value& aLabel, 
		      const TBinaryPredicate& aFunctor = TBinaryPredicate() )
      : myImg(aImg), myLabel(aLabel), myF(aFunctor)
    { 
    }
    
  private: 
    /**
     * Reference on an image 
     */
    const TImage& myImg; 
    /**
     * Value to compare with
     */
    Value myLabel; 
    /**
     * Comparison method
     */
    TBinaryPredicate myF; 

  public: 
    /**
     * Compare @a aLabel to @a myLabel with @a myF 
     * @param aPoint any point of the image domain
     * whose label has to be compared to @a myLabel
     * @return true or false
     */
    bool operator()(const Point& aPoint) const
    {
      return myF(myLabel, myImg(aPoint) ); 
    }
    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const 
    { 
      return true;
    }
    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const 
    {
      out << "(" << myLabel << ")"; 
    }  

  }; 


  /////////////////////////////////////////////////////////////////////////////
  // template class TwoLabelsPredicate
  /**
   * \brief Aim: Predicate returning true or false
   * according to the comparison of a label with 
   * two reference labels. 
   *
   * @tparam TLabel  type of labels 
   (must be default-constructible and equally-comparable)
   * @tparam TBinaryPredicate  type of binary predicate
   (must have a () operator taking two input parameters
   and returning a boolean). 
   *
   */
  template <typename TImage, typename TBinaryPredicate>
  class TwoLabelsPredicate
  {

  public: 
    typedef typename TImage::Value Value; 
    typedef typename TImage::Point Point; 

  public: 
    /**
     * Constructor 
     * @param aImg any image
     * @param aLabel1  any value
     * @param aLabel2  any value
     * @param aFunctor  any binary functor
     */
    TwoLabelsPredicate(const TImage& aImg, 
		       const Value& aLabel1, const Value& aLabel2, 
		       const TBinaryPredicate& aFunctor = TBinaryPredicate() )
      : myImg(aImg), myLabel1(aLabel1), myLabel2(aLabel2), myF(aFunctor)
    { 
    }
    
  private: 
    /**
     * Reference on an image 
     */
    const TImage& myImg; 
    /**
     * Value to compare with
     */
    Value myLabel1;
    /**
     * Value to compare with
     */
    Value myLabel2;   
    /**
     * Comparison method
     */
    TBinaryPredicate myF; 

  public: 
    /**
     * Compare the label of @a aPoint to @a myLabel1 and @a myLabel2 
     * with @a myF 
     * @param aPoint any point whose label has to be compared
     * @return true or false
     */
    bool operator()(const Point& aPoint) const
    {
      return ( myF(myLabel1, myImg(aPoint)) || myF(myLabel2, myImg(aPoint)) ); 
    }
    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const 
    { 
      return true;
    }
    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const 
    {
      out << "(" << myLabel1 << " U " << myLabel2 << ")"; 
    }
  }; 


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
//#include "PointPredicates.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined PointPredicates_h

#undef PointPredicates_RECURSES
#endif // else defined(PointPredicates_RECURSES)
