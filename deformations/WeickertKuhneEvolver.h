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
 * @file WeickertKuhneEvolver.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2011/09/05
 *
 * Header file for module WeickertKuhneEvolver.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(WeickertKuhneEvolver_RECURSES)
#error Recursive header files inclusion detected in WeickertKuhneEvolver.h
#else // defined(WeickertKuhneEvolver_RECURSES)
/** Prevents recursive inclusion of headers. */
#define WeickertKuhneEvolver_RECURSES

#if !defined WeickertKuhneEvolver_h
/** Prevents repeated inclusion of headers. */
#define WeickertKuhneEvolver_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <list>
#include <cstdlib>
#include <cmath>
#include "DGtal/base/Common.h"
#include "DGtal/images/CImageContainer.h"


//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class WeickertKuhneEvolver
  /**
   * Description of template class 'WeickertKuhneEvolver' <p>
   * \brief Aim: Evolves an implicit contour. 
   *
   * ﻿@note J. Weickert and G. Kühne,
   Fast Methods for Implicit Active Contour Models, 
    in Geometric Level Set Methods in Imaging, Vision, and Graphics, 
  Springer New York,
  43-57, 2003.
   *
   * @tparam TImage, the type of image 
   */
  template <typename TImage>
  class WeickertKuhneEvolver
  {

    //ASSERT
    BOOST_CONCEPT_ASSERT(( CImageContainer<TImage> )); 

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
    
    typedef std::vector<Value> Values; 


    


    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param aAImage Reference on a function  
     * @param aBImage Reference on a function 
     * @param aGImage Reference on a stopping function (values between 0 and 1) 
     * @param aK constant force terme (dilation when positive, erosion when negative, default 0)
     * @param aGridStep grid step size (default 1)
     *
     * NB: 
     * - a = g, b = 1 yields the geometric model 
     * - a = 1, b = g yields the geodesic model
     * - a = b = 1 and k = 0 yields the mean curvature motion
     */
    WeickertKuhneEvolver( const Image& aAImage, const Image& aBImage, const Image& aGImage, 
                          const double& aK = 0, 
			  const double& aGridStep = 1.0);

    /**
     * Destructor. Does nothing.
     */
    ~WeickertKuhneEvolver();

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
     * Reference on the  a function
     */
    const Image& myA; 
    /**
     * Reference on the  b function
     */
    const Image& myB; 
    /**
     * Reference on the   g stopping function 
     */
    const Image& myG; 
     /**
     * Constant force term
     */
    double myK;   
     /**
     * grid step
     */
    double myH; 


    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    WeickertKuhneEvolver();

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    WeickertKuhneEvolver ( const WeickertKuhneEvolver & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    WeickertKuhneEvolver & operator= ( const WeickertKuhneEvolver & other );

  private:


    /**
     * Transform function @a aU at time t0 into a next one
     * at time t0 + @a aTimeStep, but according to axis @a aDim.  
     *
     * @param aF any implicite function
     * @param aT any time step
     * @param aDim any axis
     */
    void diffusion(Image& aF, const double& aT, const Dimension& aDim); 

    
    /**
     * Returns the velocity term from four values
     * @param bi  value b function at i
     * @param gi  value of the gradient modulus of the implicit function at i
     * @param bj  value b function at j
     * @param gj  value of the gradient modulus of the implicit function at j
     * return the velocity term
     */
    Value velocity(const Value& bi, const Value& gi, const Value& bj, const Value& gj);

    /**
     * Thomas algorithm for trilinear diagonally dominant system:
     * we wand to find u, knowing B and d such that B.u = d. 
     * B is given by its 3 diagonals: 
     * - alpha  the main diagonal 
     * - beta  the upper diagonal
     * - gamma  the lower diagonal  
     * @param aU the returned vector of values 
     * @param aAlpha values vector standing for the main diagonal of the matrix
     * @param aBeta values vector  standing for the upper diagonal of the matrix
     * @param aGamma values vector  standing for the lower diagonal of the matrix
     * @param aD input values vector 
     */
    void thomasAlgorithm(Values& aU, 
      Values& aAlpha, Values& aBeta, Values& aGamma, const Values& aD); 


    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class WeickertKuhneEvolver


  /**
   * Overloads 'operator<<' for displaying objects of class 'WeickertKuhneEvolver'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'WeickertKuhneEvolver' to write.
   * @return the output stream after the writing.
   */
   template <typename TImage>
  std::ostream&
  operator<< ( std::ostream & out, const WeickertKuhneEvolver<TImage> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "WeickertKuhneEvolver.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined WeickertKuhneEvolver_h

#undef WeickertKuhneEvolver_RECURSES
#endif // else defined(WeickertKuhneEvolver_RECURSES)
