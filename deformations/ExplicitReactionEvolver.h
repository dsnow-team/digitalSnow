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

#include "DGtal/images/DifferentialOperators.h"

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
   * @tparam TImage the type of image used as the implicit function
   * @tparam TFunctor the type of functor used as an extern field
   */
  template <typename TImage, typename TFunctor>
  class ExplicitReactionEvolver
  {

    //ASSERT
    BOOST_CONCEPT_ASSERT(( CImage<TImage> )); 

    // ----------------------- Types ------------------------------
  public:

    typedef TImage Image;
    typedef TFunctor ExternImage; 
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
     * @param anEps width of the interface in the phase field equation
     * @param aF extern scalar field
     * @param aG any signed value standing for the balloon force (default 0)
     * @param aFlag boolean equal to 'true' if the volume need to remain constant
     * and equal to 'false' otherwise (default)
     */
    ExplicitReactionEvolver(const double& anEps, const ExternImage& aF, 
			    const double& aG = 0.0, 
			    bool aFlag = false );

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

    /**
     * Width of the interface
     */
    double myEpsilon;

    /**
     * Const aliasing pointer on 
     * the extern scalar field
     */
    const ExternImage* myExternField; 

    /**
     * Balloon force
     */
    double myG; 

    /**
     * Flag indicating whether the volume 
     * of the characteristic function need 
     * to be constant or not
     */
    bool myWithVolumeConservation; 

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
     * Double well function W.
     * @param aV any value
     * @return the value returned by W from @a aValue
     */
    double function ( const double & aV ) const;

    /**
     * Derivative function W' of the double well function W.
     * @param aV any value
     * @return the value returned by W' from @a aValue
     */
    double derivative ( const double & aV ) const;

    /**
     * Balloon force needed for the volume conservation
     * during the evolution of @a aF
     * @param aF implicit function to evolve
     * @return the force
     */
    double force ( Image& aF ) const;

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class ExplicitReactionEvolver


  /**
   * Overloads 'operator<<' for displaying objects of class 'ExplicitReactionEvolver'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'ExplicitReactionEvolver' to write.
   * @return the output stream after the writing.
   */
  template <typename TImage, typename TFunctor>
  std::ostream&
  operator<< ( std::ostream & out, const ExplicitReactionEvolver<TImage,TFunctor> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "ExplicitReactionEvolver.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ExplicitReactionEvolver_h

#undef ExplicitReactionEvolver_RECURSES
#endif // else defined(ExplicitReactionEvolver_RECURSES)
