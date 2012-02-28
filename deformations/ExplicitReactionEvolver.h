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
 * @file ExplicitReactionEvolver.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/01/09
 *
 * Header file for module ExplicitReactionEvolver.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ExplicitReactionEvolver_RECURSES)
#error Recursive header files inclusion detected in ExplicitReactionEvolver.h
#else // defined(ExplicitReactionEvolver_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ExplicitReactionEvolver_RECURSES

#if !defined ExplicitReactionEvolver_h
/** Prevents repeated inclusion of headers. */
#define ExplicitReactionEvolver_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/images/CImage.h"


//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class ExplicitReactionEvolver
  /**
   * Description of template class 'ExplicitReactionEvolver' <p>
   * \brief Aim: Computes the reaction term following an 
   * explicit time discretisation scheme 
   * for a double-well function W(s) = 0.5*s^2*(1-s)^2. 
   *
   * @tparam TImage, the type of image 
   */
  template <typename TImage>
  class ExplicitReactionEvolver
  {

    //ASSERT
    BOOST_CONCEPT_ASSERT(( CImage<TImage> )); 

    // ----------------------- Types ------------------------------
  public:

    typedef TImage Image;
    typedef typename Image::Value Value; 
    typedef typename Image::Point Point;
    typedef typename Image::Vector Vector;  
    typedef typename Image::Domain Domain;

    typedef typename Image::Dimension Dimension;
    static const typename Image::Dimension dimension = Image::dimension;
    

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param anEpsilon width of the interface in the phase field equation
     */
    ExplicitReactionEvolver(const double& anEpsilon);

    /**
     * Destructor. Does nothing.
     */
    ~ExplicitReactionEvolver();

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

    double myEpsilon; 
    // ------------------------- Hidden services ------------------------------
  protected:


    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    ExplicitReactionEvolver ( const ExplicitReactionEvolver & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    ExplicitReactionEvolver & operator= ( const ExplicitReactionEvolver & other );

  private:

    /**
     * Derivative function W' of the double well function W.
     * @param aV any value
     * @return the returned value of W' from @a aValue
     */
    double derivative ( const double & aV ) const;

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class ExplicitReactionEvolver


  /**
   * Overloads 'operator<<' for displaying objects of class 'ExplicitReactionEvolver'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'ExplicitReactionEvolver' to write.
   * @return the output stream after the writing.
   */
   template <typename TImage>
  std::ostream&
  operator<< ( std::ostream & out, const ExplicitReactionEvolver<TImage> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "ExplicitReactionEvolver.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ExplicitReactionEvolver_h

#undef ExplicitReactionEvolver_RECURSES
#endif // else defined(ExplicitReactionEvolver_RECURSES)
