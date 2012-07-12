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
 * @file MultiPhaseField.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/07/12
 *
 * Header file for module MultiPhaseField.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(MultiPhaseField_RECURSES)
#error Recursive header files inclusion detected in MultiPhaseField.h
#else // defined(MultiPhaseField_RECURSES)
/** Prevents recursive inclusion of headers. */
#define MultiPhaseField_RECURSES

#if !defined MultiPhaseField_h
/** Prevents repeated inclusion of headers. */
#define MultiPhaseField_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/images/CImage.h"

#include "DGtal/base/CowPtr.h"
#include "DGtal/base/CountedPtr.h"

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{


  /////////////////////////////////////////////////////////////////////////////
  // template class MultiPhaseField
  /**
   * Description of template class 'MultiPhaseField' <p>
   * \brief Aim: This class is a way of deforming an image of labels. 
   * Each region (ie. set of points having a same label) is viewed as 
   * the set of points having a value greater than 0.5 for a given phase field. 
   * Each region is evolved through its phase field. 
   *
   * @tparam TLabelImage a model of CImage (storing labels)
   * @tparam TFieldImage a model of CImage (storing phase field values)
   * @tparam TEvolver a model of phase field evolver
   */
  template <typename TLabelImage, typename TFieldImage, typename TEvolver>
  class MultiPhaseField
  {

    // ----------------------- Types check -----------------------

    BOOST_CONCEPT_ASSERT(( CImage<TLabelImage> )); 
    BOOST_CONCEPT_ASSERT(( CImage<TFieldImage> )); 
    BOOST_STATIC_ASSERT
    (( ConceptUtils::SameType< typename TLabelImage::Point,
       typename TFieldImage::Point>::value ));


    // ----------------------- Types ------------------------------
  public:

    /// Image of labels
    typedef TLabelImage LabelImage;
    typedef typename LabelImage::Value Label;
    typedef typename LabelImage::Domain Domain;
    typedef typename Domain::Point Point;

    /// Images of phase field values
    typedef TFieldImage FieldImage;
    typedef CowPtr<FieldImage> FieldImagePtr; 


    /// Phase field evolver
    typedef TEvolver Evolver; 
  
    // ------------------------- Protected Datas ------------------------------
  protected:
    // ------------------------- Private Datas --------------------------------
  private:


    // ------------------------- References --------------------------------
    /**
     * Reference on the image of labels
     */
    LabelImage& myLabelImage; 
    /**
     * Constant reference on the evolver
     */
    const Evolver& myEvolver; 


    // ------------------------- Data --------------------------------
    /**
     * List of smart owning pointers on phase fields
     */
    std::vector<FieldImagePtr> myFields; 
    /**
     * List of labels
     */
    std::vector<Label> myLabels; 

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param aI an image of labels
     * @param aE any phase field evolver
     */
    MultiPhaseField(LabelImage& aI, const Evolver& aE);

    /**
     * Destructor. Does nothing.
     */
    ~MultiPhaseField();


    /**
     * Deform the image of labels during @a aT
     *
     * @param aT time step
     * @return time spent during the deformation 
     * (equal to aT). 
     */
    double update(const double& aT);


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




    // ------------------------- Hidden services ------------------------------
  protected:


    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    MultiPhaseField ( const MultiPhaseField & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    MultiPhaseField & operator= ( const MultiPhaseField & other );

  private:



    // ------------------------- Internals ------------------------------------
  private:

    /**
     * Init @a aImage by a signed distance to the frontier
     * of the region having @a aLabel as label
     * @param aLabel region id
     * @param aImage image to initialize
     */
    void getSignedDistance(const Label& aLabel, FieldImage& aImage); 

  }; // end of class MultiPhaseField


  /**
   * Overloads 'operator<<' for displaying objects of class 'MultiPhaseField'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'MultiPhaseField' to write.
   * @return the output stream after the writing.
   */
  template <typename TLabelImage, typename TFieldImage, typename TEvolver>
  std::ostream&
  operator<< ( std::ostream & out, 
	       const DGtal::MultiPhaseField<TLabelImage, TFieldImage, TExternImage, TEvolver>& object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "MultiPhaseField.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined MultiPhaseField_h

#undef MultiPhaseField_RECURSES
#endif // else defined(MultiPhaseField_RECURSES)
