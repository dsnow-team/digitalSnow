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
 * @file BinaryPredicates.h
 *
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 * @date 2012/02/02
 *
 * This files contains several basic classes representing binary predicates
 *
 * This file is part of the DGtal library.
 */

#if defined(BinaryPredicates_RECURSES)
#error Recursive header files inclusion detected in BinaryPredicates.h
#else // defined(BinaryPredicates_RECURSES)
/** Prevents recursive inclusion of headers. */
#define BinaryPredicates_RECURSES

#if !defined BinaryPredicates_h
/** Prevents repeated inclusion of headers. */
#define BinaryPredicates_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class ConstantBinaryPredicate
  /**
   * Description of template class 'ConstantPointPredicate' <p>
   * \brief Aim: The predicate that returns always the same value boolCst
   *
   * @tparam boolCst any boolean value
   */
  template <bool boolCst>
  struct ConstantBinaryPredicate
  {
    /**
     * @param t1 first argument
     * @param t2 second argument
     * @tparam T1 type of the first argument
     * @tparam T2 type of the second argument
     * @return the value of the predicate 
     */
    template <typename T1, typename T2>
    bool operator()( const T1& t1, const T2& t2 ) const
      {
        return boolCst; 
      }

  }; // end of class ConstantBinaryPredicate

  typedef ConstantBinaryPredicate<true> TrueBinaryPredicate; 
  typedef ConstantBinaryPredicate<false> FalseBinaryPredicate; 

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
//#include "BinaryPredicates.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined BinaryPredicates_h

#undef BinaryPredicates_RECURSES
#endif // else defined(BinaryPredicates_RECURSES)
