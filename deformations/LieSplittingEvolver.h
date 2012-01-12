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
 * @file LieSplittingEvolver.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/01/09
 *
 * Header file for module LieSplittingEvolver.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(LieSplittingEvolver_RECURSES)
#error Recursive header files inclusion detected in LieSplittingEvolver.h
#else // defined(LieSplittingEvolver_RECURSES)
/** Prevents recursive inclusion of headers. */
#define LieSplittingEvolver_RECURSES

#if !defined LieSplittingEvolver_h
/** Prevents repeated inclusion of headers. */
#define LieSplittingEvolver_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/images/CImageContainer.h"


//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class LieSplittingEvolver
  /**
   * Description of template class 'LieSplittingEvolver' <p>
   * \brief Aim: Computes the combination of two evolution terms 
   * following a Lie splitting scheme.  
   * Actually, to transform a function f at time t0 into a next one
   * at time t0 + delta, f is first evolved during a time step equal to 
   * delta according to the first term and then evolved during the same
   * time step according to the second term. 
   *
   * @tparam TEvolver1, first type of evolver 
   * @tparam TEvolver2, second type of evolver 
   */
  template <typename TEvolver1, typename TEvolver2>
  class LieSplittingEvolver
  {

    //concepts to write
    //ASSERT
    //BOOST_CONCEPT_ASSERT(( CImplicitEvolver<TEvolver1> )); 
    //BOOST_CONCEPT_ASSERT(( CImplicitEvolver<TEvolver2> )); 

    // ----------------------- Types ------------------------------
  public:

    typedef TEvolver1 Evolver1;
    typedef TEvolver2 Evolver2;

    typedef typename Evolver1::Image Image; 
    BOOST_STATIC_ASSERT ((boost::is_same< Image, typename Evolver2::Image >::value)); 

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param anEvolver1 first evolver
     * @param anEvolver2 second evolver
     */
    LieSplittingEvolver(Evolver1& anEvolver1, Evolver2& anEvolver2);

    /**
     * Destructor. Does nothing.
     */
    ~LieSplittingEvolver();

    /**
     * Transform function @a aF at time t0 into a next one
     * at time t0 + @a aT
     *
     * @param aF any implicite function
     * @param aT any time step
     */
    void update ( Image& aF, const double& aT );


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
     * reference on the first evolver
     */
    Evolver1& myEvolver1; 
    /**
     * reference on the second evolver
     */
    Evolver2& myEvolver2; 
    // ------------------------- Hidden services ------------------------------
  protected:


    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    LieSplittingEvolver ( const LieSplittingEvolver & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    LieSplittingEvolver & operator= ( const LieSplittingEvolver & other );

  private:



    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class LieSplittingEvolver


  /**
   * Overloads 'operator<<' for displaying objects of class 'LieSplittingEvolver'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'LieSplittingEvolver' to write.
   * @return the output stream after the writing.
   */
   template <typename TEvolver1, typename TEvolver2>
  std::ostream&
  operator<< ( std::ostream & out, const LieSplittingEvolver<TEvolver1, TEvolver2> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "LieSplittingEvolver.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined LieSplittingEvolver_h

#undef LieSplittingEvolver_RECURSES
#endif // else defined(LieSplittingEvolver_RECURSES)
