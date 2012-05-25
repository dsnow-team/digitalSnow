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
 * @file ExactDiffusionEvolver.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/01/09
 *
 * Header file for module ExactDiffusionEvolver.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ExactDiffusionEvolver_RECURSES)
#error Recursive header files inclusion detected in ExactDiffusionEvolver.h
#else // defined(ExactDiffusionEvolver_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ExactDiffusionEvolver_RECURSES

#if !defined ExactDiffusionEvolver_h
/** Prevents repeated inclusion of headers. */
#define ExactDiffusionEvolver_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/images/CImage.h"

#include "FFT.h"
#include "IFFT.h"

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class ExactDiffusionEvolver
  /**
   * Description of template class 'ExactDiffusionEvolver' <p>
   * \brief Aim: Computes the exact diffusion at any time t of
   * an implicit function (i.e. computes its convolution 
   * with the heat kernel in the fourier space). 
   *
   * NB: the heat kernel is a gaussian kernel s.t sigma = sqrt(2t)
   *
   * @tparam TImage, the type of image 
   */
  template <typename TImage>
  class ExactDiffusionEvolver
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

    typedef typename Image::iterator Iterator;
    typedef typename Image::const_iterator ConstIterator;

    typedef typename Image::Dimension Dimension;
    static const typename Image::Dimension dimension = Image::dimension;
    

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     *
     */
    ExactDiffusionEvolver();

    /**
     * Destructor. Does nothing.
     */
    ~ExactDiffusionEvolver();

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

    // ------------------------- Hidden services ------------------------------
  protected:


    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    ExactDiffusionEvolver ( const ExactDiffusionEvolver & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    ExactDiffusionEvolver & operator= ( const ExactDiffusionEvolver & other );

  private:



    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class ExactDiffusionEvolver


  /**
   * Overloads 'operator<<' for displaying objects of class 'ExactDiffusionEvolver'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'ExactDiffusionEvolver' to write.
   * @return the output stream after the writing.
   */
   template <typename TImage>
  std::ostream&
  operator<< ( std::ostream & out, const ExactDiffusionEvolver<TImage> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "ExactDiffusionEvolver.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ExactDiffusionEvolver_h

#undef ExactDiffusionEvolver_RECURSES
#endif // else defined(ExactDiffusionEvolver_RECURSES)
