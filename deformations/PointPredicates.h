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
  // template class CascadingPointPredicate
  /**
   * Description of template class 'CascadingPointPredicate' <p> \brief
   * Aim: The predicate returns true when the given binary functor
   * returns true for the two PointPredicates given at construction.
   *
   * @tparam PointPredicate1 the left predicate type.
   * @tparam PointPredicate2 the right predicate type.
   */
  template <typename TPointPredicate1, typename TPointPredicate2>
  struct CascadingPointPredicate
  {
    typedef TPointPredicate1 PointPredicate1;
    typedef TPointPredicate2 PointPredicate2;
    typedef typename PointPredicate1::Point Point;
    // should be the same.
    BOOST_STATIC_ASSERT ((boost::is_same< Point, typename PointPredicate2::Point >::value)); 
    typedef typename PointPredicate2::Point Point2;

    /**
       Constructor from predicates
       @param pred1 the left predicate.
       @param pred2 the right predicate.
     */
    CascadingPointPredicate( const PointPredicate1 & pred1,
        const PointPredicate2 & pred2 );

    /**
       Copy constructor.
       @param other the object to copy
      */
    CascadingPointPredicate(  const CascadingPointPredicate& other );

    /**
       Assignement
       @param other the object to copy
       @return reference to the current object
     */
    CascadingPointPredicate& operator=( const CascadingPointPredicate& other );

    /**
       Destructor
     */
    ~CascadingPointPredicate();

    /**
     * @param p any point.
     * @return the value of the predicate at this point.
     */
    bool operator()( const Point & p ) const;

    /// aliasing pointer to the left predicate.
    const PointPredicate1* myPred1;
    /// aliasing pointer to the right predicate.
    const PointPredicate2* myPred2;
  };

//------------------------------------------------------------------------------
template <typename TPointPredicate1, typename TPointPredicate2>
inline
DGtal::CascadingPointPredicate<TPointPredicate1,TPointPredicate2>
::CascadingPointPredicate( const PointPredicate1 & pred1,
      const PointPredicate2 & pred2 )
  : myPred1( &pred1 ), myPred2( &pred2 )
{
}
//------------------------------------------------------------------------------
template <typename TPointPredicate1, typename TPointPredicate2>
inline
DGtal::CascadingPointPredicate<TPointPredicate1,TPointPredicate2>
::CascadingPointPredicate( const CascadingPointPredicate& other )
  : myPred1( other.pred1 ), myPred2( other.pred2 )
{
}
//------------------------------------------------------------------------------
template <typename TPointPredicate1, typename TPointPredicate2>
inline
DGtal::CascadingPointPredicate<TPointPredicate1,TPointPredicate2>&
DGtal::CascadingPointPredicate<TPointPredicate1,TPointPredicate2>
::operator=( const CascadingPointPredicate& other )
{
  if (this != &other)
    {
      myPred1 = other.myPred1; 
      myPred2 = other.myPred2; 
    }
}
//------------------------------------------------------------------------------
template <typename TPointPredicate1, typename TPointPredicate2>
inline
DGtal::CascadingPointPredicate<TPointPredicate1,TPointPredicate2>
::~CascadingPointPredicate()
{
}
//------------------------------------------------------------------------------
template <typename TPointPredicate1, typename TPointPredicate2>
inline
bool 
DGtal::CascadingPointPredicate<TPointPredicate1,TPointPredicate2>
::operator()( const Point & p ) const
{
  if ( myPred1->operator()( p ) ) 
		return myPred2->operator()( p );
  else return false; 
}


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
