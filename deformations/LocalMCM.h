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
 * @file LocalMCM.h
 *
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 * @date 2012/03/01
 *
 * This files contains several basic classes representing Functors
 * on points.
 *
 * This file is part of the DGtal library.
 */

#if defined(LocalMCM_RECURSES)
#error Recursive header files inclusion detected in LocalMCM.h
#else // defined(LocalMCM_RECURSES)
/** Prevents recursive inclusion of headers. */
#define LocalMCM_RECURSES

#if !defined LocalMCM_h
/** Prevents repeated inclusion of headers. */
#define LocalMCM_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include "DGtal/kernel/CPointFunctor.h"

#include "DGtal/images/DifferentialOperators.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class LocalMCM
  /**
   * Description of template class 'LocalMCM' <p>
   * \brief Aim: Functor that maps a point 
   * to a velocity computed from an implicit function
   * and some extern data (two extern scalar fields 
   * for weighting). 
   *
   * @tparam TFunction type of implicit function, 
   * which a model of point functor
   * @tparam TExternField type of the extern field, 
   * which a model of point functor too
   */
  template <typename TFunction, typename TExternField >
  struct LocalMCM
  {
    typedef TFunction PointFunctor;
    typedef TExternField ExternField; 

    BOOST_CONCEPT_ASSERT(( CPointFunctor<PointFunctor> )); 
    BOOST_CONCEPT_ASSERT(( CPointFunctor<ExternField> )); 
    BOOST_STATIC_ASSERT
    (( ConceptUtils::SameType< typename PointFunctor::Point,
       typename ExternField::Point>::value ));

    typedef typename PointFunctor::Point Point; 
    typedef double Value; 
    
    typedef typename Gradient<CentralDifference<PointFunctor> >::OutputValue Normal; 
    typedef typename Divergence<WeightedDifference2<PointFunctor> >::OutputValue Curvature; 

    /**
     * Constructor
     */
    LocalMCM(PointFunctor& aF1, const ExternField& aA, const ExternField& aB);

    /**
     * Main operator
     * @param aPoint any point.
     * @tparam TInputPoint type of point
     * @return the velocity at @a aPoint.
     */
    template<typename TInputPoint>
    double operator()( const TInputPoint& aPoint ) const;

    /**
     * Gradient.
     * @param aPoint any point.
     * @tparam TInputPoint type of point
     * @return the gradient at @a aPoint.
     */
    template<typename TInputPoint>
    Normal getNormal( const TInputPoint& aPoint ) const;

    /**
     * Curvature.
     * @param aPoint any point.
     * @tparam TInputPoint type of point
     * @return curvature at @a aPoint.
     */
    template<typename TInputPoint>
    Curvature getCurvature( const TInputPoint& aPoint ) const;

   private: 
    /**
     * Aliasing pointer to the implicit function
     */
    const PointFunctor* myFuncPtr;  
    /**
     * Constant aliasing pointer to the extern field A
     */
    const ExternField* myAPtr;  
    /**
     * Constant aliasing pointer to the extern field B
     */
    const ExternField* myBPtr;  

  }; // end of class LocalMCM


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "LocalMCM.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined LocalMCM_h

#undef LocalMCM_RECURSES
#endif // else defined(LocalMCM_RECURSES)
