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
 * @file LocalBalloonForce.h
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

#if defined(LocalBalloonForce_RECURSES)
#error Recursive header files inclusion detected in LocalBalloonForce.h
#else // defined(LocalBalloonForce_RECURSES)
/** Prevents recursive inclusion of headers. */
#define LocalBalloonForce_RECURSES

#if !defined LocalBalloonForce_h
/** Prevents repeated inclusion of headers. */
#define LocalBalloonForce_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include "DGtal/kernel/CPointFunctor.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class LocalBalloonForce
  /**
   * Description of template class 'LocalBalloonForce' <p>
   * \brief Aim: Functor that maps a point 
   * to a velocity computed from an implicit function
   * and some extern data (a balloon force and 
   * an extern scalar field for weighting). 
   *
   * @tparam TFunction type of implicit function, 
   * which a model of point functor
   * @tparam TExternField type of the extern field, 
   * which a model of point functor too
   */
  template <typename TFunction, typename TExternField >
  struct LocalBalloonForce
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

    /**
     * Constructor
     */
    LocalBalloonForce(PointFunctor& aF1, const ExternField& aF2, const double& aK = 0);

    /**
     * Main operator
     * @param aPoint any point.
     * @tparam TInputPoint type of point
     * @return the velocity at @a aPoint.
     */
    template<typename TInputPoint>
    double operator()( const TInputPoint& aPoint ) const;

   private: 
    /**
     * Aliasing pointer to the implicit function
     */
    PointFunctor* myFuncPtr;  
    /**
     * Constant aliasing pointer to the extern field
     */
    const ExternField* myFieldPtr;  
    /**
     * Balloon force
     */
    double myK; 

  }; // end of class LocalBalloonForce


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "LocalBalloonForce.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined LocalBalloonForce_h

#undef LocalBalloonForce_RECURSES
#endif // else defined(LocalBalloonForce_RECURSES)
