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
 * @file ExactReactionEvolver.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/01/09
 *
 * Header file for module ExactReactionEvolver.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ExactReactionEvolver_RECURSES)
#error Recursive header files inclusion detected in ExactReactionEvolver.h
#else // defined(ExactReactionEvolver_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ExactReactionEvolver_RECURSES

#if !defined ExactReactionEvolver_h
/** Prevents repeated inclusion of headers. */
#define ExactReactionEvolver_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/images/CImage.h"


//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class ExactReactionEvolver
  /**
   * Description of template class 'ExactReactionEvolver' <p>
   * \brief Aim: Computes the exact reaction term at any time t 
   * for a double-well function W(s) = 0.5*s^2*(1-s)^2. 
   *
   * @tparam TImage, the type of image 
   */
  template <typename TImage>
  class ExactReactionEvolver
  {

    //ASSERT
    BOOST_CONCEPT_ASSERT(( concepts::CImage<TImage> )); 

    // ----------------------- Types ------------------------------
  public:

    typedef TImage Image;
    typedef typename Image::Value Value; 
    typedef typename Image::Point Point;
    typedef typename Image::Vector Vector;  
    typedef typename Image::Domain Domain;

    typedef typename Image::iterator Iterator;
    typedef typename Image::const_iterator ConstIterator;

    typedef typename Image::Dimension Dimension;
    static const typename Image::Dimension dimension = Image::dimension;
    

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param anEpsilon width of the interface in the phase field equation
     */
    ExactReactionEvolver(const double& anEpsilon);

    /**
     * Destructor. Does nothing.
     */
    ~ExactReactionEvolver();

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
    ExactReactionEvolver ( const ExactReactionEvolver & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    ExactReactionEvolver & operator= ( const ExactReactionEvolver & other );

  private:



    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class ExactReactionEvolver


  /**
   * Overloads 'operator<<' for displaying objects of class 'ExactReactionEvolver'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'ExactReactionEvolver' to write.
   * @return the output stream after the writing.
   */
   template <typename TImage>
  std::ostream&
  operator<< ( std::ostream & out, const ExactReactionEvolver<TImage> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "ExactReactionEvolver.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ExactReactionEvolver_h

#undef ExactReactionEvolver_RECURSES
#endif // else defined(ExactReactionEvolver_RECURSES)
