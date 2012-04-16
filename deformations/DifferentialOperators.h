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
 * @file DifferentialOperators.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 201/12/19
 *
 * Header file for module DifferentialOperators.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(DifferentialOperators_RECURSES)
#error Recursive header files inclusion detected in DifferentialOperators.h
#else // defined(DifferentialOperators_RECURSES)
/** Prevents recursive inclusion of headers. */
#define DifferentialOperators_RECURSES

#if !defined DifferentialOperators_h
/** Prevents repeated inclusion of headers. */
#define DifferentialOperators_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/images/CConstImage.h"

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{


  /////////////////////////////////////////////////////////////////////////////
  // template class ForwardDifference
  /**
   * Description of template class 'ForwardDifference' <p>
   * \brief Aim: Computes the forward difference at a point. 
   *
   *
   * @tparam TFonctor model of CPointFunctor 
   * @tparam TPointPredicate model of CPointPredicate
   * @tparam TOutputValue type of returned value (default TFunctor::Value)
   */
  template <typename TImage, typename TOutputValue = typename TImage::Value >
  class ForwardDifference
  {

    BOOST_CONCEPT_ASSERT(( CConstImage<TImage> )); 

    // ----------------------- Types ------------------------------
  public:

    typedef TImage Image;
    typedef TOutputValue OutputValue; 
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
     *
     * @param aStartingImage  any image of signed values
     * @param aGridStep any length (=1 by default)
     */
    ForwardDifference( Image& aStartingImage, const OutputValue& aGridStep = 1);

    /**
     * Destructor. Does nothing.
     */
    ~ForwardDifference() {}

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const {return true;}

    /**
     * Difference.
     *
     * @param aPoint the point where the derivative is computed
     * @param aDim the axis along which the derivative is computed
     * @return first derivative along axis @a aDim at @ aPoint
     */
    OutputValue operator() ( const Point& aPoint, const Dimension& aDim ) const; 

    // ------------------------- Private Datas --------------------------------
  private:

    /**
     * Reference on an image
     */
    Image& myU; 

    /**
     * Grid step
     */
    OutputValue myH; 

  }; 

  /////////////////////////////////////////////////////////////////////////////
  // template class BackwardDifference
  /**
   * Description of template class 'BackwardDifference' <p>
   * \brief Aim: Computes the backward difference at a point
   * of an image. 
   *
   * @code 
   * @endcode
   *
   * @tparam TImage type of image 
   * @tparam TOutputValue type of returned value (default TImage::Value)
   */
  template <typename TImage, typename TOutputValue = typename TImage::Value >
  class BackwardDifference
  {

    BOOST_CONCEPT_ASSERT(( CConstImage<TImage> )); 

    // ----------------------- Types ------------------------------
  public:

    typedef TImage Image;
    typedef TOutputValue OutputValue; 
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
     *
     * @param aStartingImage  any image of signed values
     * @param aGridStep any length (=1 by default)
     */
    BackwardDifference( Image& aStartingImage, const OutputValue& aGridStep = 1);

    /**
     * Destructor. Does nothing.
     */
    ~BackwardDifference() {}

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const {return true;}

    /**
     * Difference.
     *
     * @param aPoint the point where the derivative is computed
     * @param aDim the axis along which the derivative is computed
     * @return first derivative along axis @a aDim at @ aPoint
     */
    OutputValue operator() ( const Point& aPoint, const Dimension& aDim ) const; 

    // ------------------------- Private Datas --------------------------------
  private:

    /**
     * Reference on an image
     */
    Image& myU; 

    /**
     * Grid step
     */
    OutputValue myH; 

  }; 

  /////////////////////////////////////////////////////////////////////////////
  // template class CentralDifference
  /**
   * Description of template class 'CentralDifference' <p>
   * \brief Aim: Computes the backward difference at a point
   * of an image. 
   *
   * @code 
   * @endcode
   *
   * @tparam TImage type of image 
   * @tparam TOutputValue type of returned value (default TImage::Value)
   */
  template <typename TImage, typename TOutputValue = typename TImage::Value >
  class CentralDifference
  {

    BOOST_CONCEPT_ASSERT(( CConstImage<TImage> )); 

    // ----------------------- Types ------------------------------
  public:

    typedef TImage Image;
    typedef TOutputValue OutputValue; 
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
     *
     * @param aStartingImage  any image of signed values
     * @param aGridStep any length (=1 by default)
     */
    CentralDifference( Image& aStartingImage, const OutputValue& aGridStep = 1);

    /**
     * Destructor. Does nothing.
     */
    ~CentralDifference() {}

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const {return true;}

    /**
     * Difference.
     *
     * @param aPoint the point where the derivative is computed
     * @param aDim the axis along which the derivative is computed
     * @return first derivative along axis @a aDim at @ aPoint
     */
    OutputValue operator() ( const Point& aPoint, const Dimension& aDim ) const; 

    // ------------------------- Private Datas --------------------------------
  private:

    /**
     * Reference on an image
     */
    Image& myU; 

    /**
     * Grid step
     */
    OutputValue myH; 

  }; 


  /////////////////////////////////////////////////////////////////////////////
  // template class Difference2
  /**
   * Description of template class 'Difference2' <p>
   * \brief Aim: Computes the second difference at a point
   * of an image. 
   *
   * @code 
   * @endcode
   *
   * @tparam TImage type of image 
   * @tparam TOutputValue type of returned value (default TImage::Value)
   */
  template <typename TImage, typename TOutputValue = typename TImage::Value >
  class Difference2
  {

    BOOST_CONCEPT_ASSERT(( CConstImage<TImage> )); 

    // ----------------------- Types ------------------------------
  public:

    typedef TImage Image;
    typedef TOutputValue OutputValue; 
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
     *
     * @param aStartingImage  any image of signed values
     * @param aGridStep any length (=1 by default)
     */
    Difference2( Image& aStartingImage, const OutputValue& aGridStep = 1);

    /**
     * Destructor. Does nothing.
     */
    ~Difference2() {}

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const {return true;}

    /**
     * Second forward/backward difference.
     *
     * @param aPoint the point where the derivative is computed
     * @param aDim the axis along which the derivative is computed
     * @return first derivative along axis @a aDim at @ aPoint
     */
    OutputValue operator() ( const Point& aPoint, const Dimension& aDim ) const; 

    // ------------------------- Private Datas --------------------------------
  private:

    /**
     * Reference on an image
     */
    Image& myU; 

    /**
     * Grid step
     */
    OutputValue myH; 

  }; 

  /////////////////////////////////////////////////////////////////////////////
  // template class NormalizedDifference2
  /**
   * Description of template class 'NormalizedDifference2' <p>
   * \brief Aim: Computes the second difference at a point
   * of an image. 
   *
   * @code 
   * @endcode
   *
   * @tparam TImage type of image 
   * @tparam TOutputValue type of returned value (default TImage::Value)
   */
  template <typename TImage, typename TOutputValue = typename TImage::Value >
  class NormalizedDifference2
  {

    BOOST_CONCEPT_ASSERT(( CConstImage<TImage> )); 

    // ----------------------- Types ------------------------------
  public:

    typedef TImage Image;
    typedef TOutputValue OutputValue; 
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
     *
     * @param aStartingImage  any image of signed values
     * @param aGridStep any length (=1 by default)
     */
    NormalizedDifference2( Image& aStartingImage, const OutputValue& aGridStep = 1);

    /**
     * Destructor. Does nothing.
     */
    ~NormalizedDifference2() {}

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const {return true;}

    /**
     * Second forward/backward difference
     * normalized by the inverse of the gradient modulus
     *
     * @param aPoint the point where the derivative is computed
     * @param aDim the axis along which the derivative is computed
     * @return first derivative along axis @a aDim at @ aPoint
     */
    OutputValue operator() ( const Point& aPoint, const Dimension& aDim ) const; 

    // ------------------------- Private Datas --------------------------------
  private:

    /**
     * Reference on an image
     */
    Image& myU; 

    /**
     * Grid step
     */
    OutputValue myH; 

    // ------------------------- Internals --------------------------------
  private:

    /**
     * Return the harmonic average of the inverse of @a aV1 and @a aV2 
     *
     * @param aV1 a first value
     * @param aV2 a second value
     * @return the average 
     */
    double average ( const Value& aV1, const Value& aN2 ) const; 

  }; 

  /////////////////////////////////////////////////////////////////////////////
  // template class WeightedDifference2
  /**
   * Description of template class 'WeightedDifference2' <p>
   * \brief Aim: Computes the second difference at a point
   * of an image. 
   *
   * @code 
   * @endcode
   *
   * @tparam TImage type of image 
   * @tparam TOutputValue type of returned value (default TImage::Value)
   */
  template <typename TImage, typename TOutputValue = typename TImage::Value >
  class WeightedDifference2
  {

    BOOST_CONCEPT_ASSERT(( CConstImage<TImage> )); 

    // ----------------------- Types ------------------------------
  public:

    typedef TImage Image;
    typedef TOutputValue OutputValue; 
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
     *
     * @param aImage  any image
     * @param aWImage  any image of weights
     * @param aGridStep any length (=1 by default)
     */
    WeightedDifference2( Image& aImage, Image& aWImage, const OutputValue& aGridStep = 1);

    /**
     * Destructor. Does nothing.
     */
    ~WeightedDifference2() {}

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const {return true;}

    /**
     * Second forward/backward differences on @a aImage
     * weighted by @a aWImage and normalized 
     * by the inverse of the gradient modulus
     *
     * @param aPoint the point where the derivative is computed
     * @param aDim the axis along which the derivative is computed
     * @return second derivative of @a aImage along axis @a aDim at @ aPoint
     */
    OutputValue operator() ( const Point& aPoint, const Dimension& aDim ) const; 

    // ------------------------- Private Datas --------------------------------
  private:

    /**
     * Reference on an image
     */
    Image& myImage; 
    /**
     * Reference on the weights image
     */
    Image& myWImage; 

    /**
     * Grid step
     */
    OutputValue myH; 

    // ------------------------- Internals --------------------------------
  private:

    /**
     * Return the harmonic average of @a aV1 and @a aV2, 
     * @a aV1 and @a aV2 being normalized by @a aN1 and @a aN2.
     *
     * @param aV1 a first value
     * @param aN1 any value dividing @a aV1
     * @param aV2 a second value
     * @param aN2 any value dividing @a aV2
     * @return the normalized harmonic average of @a aV1 and @a aV2
     */
    double average ( const Value& aV1, const double& aN1, 
		     const Value& aV2, const double& aN2 ) const; 

  }; 

  ///////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////
  // template class Gradient
  /**
   * Description of template class 'Gradient' <p>
   * \brief Aim: Computes the gradient at a point
   * of an image. 
   *
   * @code 
   * @endcode
   *
   * @tparam TDifference  type of directionnal differential operator 
   */
  template <typename TDifference>
  class Gradient
  {


    // ----------------------- Types ------------------------------
  public:

    typedef TDifference FiniteDifference;
    typedef typename FiniteDifference::Dimension Dimension; 
    static const typename FiniteDifference::Dimension dimension = FiniteDifference::dimension;
    typedef typename FiniteDifference::Image Image; 
    typedef typename FiniteDifference::Point Point; 
   
    typedef DGtal::PointVector<dimension,typename FiniteDifference::OutputValue> OutputValue;
 
    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     *
     * @param aD a finite difference operator
     */
    Gradient( const FiniteDifference& aD): myD(aD) {}


    /**
     * Constructor.
     *
     * @param aStartingImage  any image of signed values
     * @param aGridStep any length (=1 by default)
     */
    Gradient( typename FiniteDifference::Image& aStartingImage, 
	      const typename FiniteDifference::OutputValue& aGridStep = 1);

    /**
     * Destructor. Does nothing.
     */
    ~Gradient() {}

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const {return true;}

    /**
     * Gradient
     *
     * @param aPoint the point where the gradient is computed
     * @return gradient of @a myU at @ aPoint
     */
    OutputValue operator() ( const Point& aPoint ) const; 

    // ------------------------- Private Datas --------------------------------
  private:

    /**
     * Finite difference
     */
    FiniteDifference myD; 

  }; 

  /////////////////////////////////////////////////////////////////////////////
  // template class UpwindGradient
  /**
   * Description of template class 'UpwindGradient' <p>
   * \brief Aim: Computes the gradient at a point
   * of an image. 
   *
   * @code 
   * @endcode
   *
   * @tparam TDifference  type of directionnal differential operator 
   */
  template <typename TImage, typename TOutputValue = typename TImage::Value>
  class UpwindGradient
  {


    // ----------------------- Types ------------------------------
  public:

    BOOST_CONCEPT_ASSERT(( CConstImage<TImage> )); 

    // ----------------------- Types ------------------------------
  public:

    typedef TImage Image;
    typedef typename Image::Value Value;
 
    typedef typename Image::Point Point;
    typedef typename Image::Vector Vector;  
    typedef typename Image::Domain Domain;

    typedef typename Image::Dimension Dimension;
    static const typename Image::Dimension dimension = Image::dimension;

    typedef DGtal::PointVector<dimension,TOutputValue> OutputValue;
 
    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     *
     * @param aStartingImage  any image of signed values
     * @param aGridStep any length (=1 by default)
     */
    UpwindGradient( Image& aStartingImage, 
	      const TOutputValue& aGridStep = 1);

    /**
     * Destructor. Does nothing.
     */
    ~UpwindGradient() {}

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const {return true;}

    /**
     * Gradient
     *
     * @param aPoint the point where the gradient is computed
     * @param aVector displacement vector of the interface
     * @return gradient of @a myU at @ aPoint
     */
    OutputValue operator() ( const Point& aPoint, const Vector& aVector ) const; 

    // ------------------------- Private Datas --------------------------------
  private:

    /**
     * Forward difference
     */
    ForwardDifference<Image,TOutputValue> myF; 
    /**
     * Backward difference
     */
    BackwardDifference<Image,TOutputValue> myB; 

  }; 

  /////////////////////////////////////////////////////////////////////////////
  // template class GodunovGradient
  /**
   * Description of template class 'GodunovGradient' <p>
   * \brief Aim: Computes the gradient at a point
   * of an image. 
   *
   * @code 
   * @endcode
   *
   * @tparam TDifference  type of directionnal differential operator 
   */
  template <typename TImage, typename TOutputValue = typename TImage::Value>
  class GodunovGradient
  {


    // ----------------------- Types ------------------------------
  public:

    BOOST_CONCEPT_ASSERT(( CConstImage<TImage> )); 

    // ----------------------- Types ------------------------------
  public:

    typedef TImage Image;
    typedef typename Image::Value Value;
 
    typedef typename Image::Point Point;
    typedef typename Image::Vector Vector;  
    typedef typename Image::Domain Domain;

    typedef typename Image::Dimension Dimension;
    static const typename Image::Dimension dimension = Image::dimension;

    typedef DGtal::PointVector<dimension,TOutputValue> OutputValue;
 
    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     *
     * @param aStartingImage  any image
     * @param isPositive flag equal to 'true' if the 
     * interface moves in the direction of the normal
     * and 'false' if it moves in the opposite direction 
     * (='true' by default)
     * @param aGridStep any length (=1 by default)
     */
    GodunovGradient( Image& aStartingImage, bool isPositive = true, 
	      const TOutputValue& aGridStep = 1);

    /**
     * Destructor. Does nothing.
     */
    ~GodunovGradient() {}

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const {return true;}

    /**
     * Gradient
     *
     * @param aPoint the point where the gradient is computed
     * @return gradient of @a myU at @ aPoint
     */
    OutputValue operator() ( const Point& aPoint ) const; 

    /**
     * Gradient
     *
     * @param aPoint the point where the gradient is computed
     * @param isPositive flag equal to 'true' if the 
     * interface moves in the direction of the normal
     * and 'false' if it moves in the opposite direction 
     * (='true' by default)
     * @return gradient of @a myU at @ aPoint
     */
    OutputValue operator() ( const Point& aPoint, bool isPositive ) const; 

    // ------------------------- Private Datas --------------------------------
  private:

    /**
     * Forward difference
     */
    ForwardDifference<Image,TOutputValue> myF; 
    /**
     * Backward difference
     */
    BackwardDifference<Image,TOutputValue> myB; 
    /**
     * Flag
     */
    bool myIsPositive; 
  }; 

  /////////////////////////////////////////////////////////////////////////////
  // template class Divergence
  /**
   * Description of template class 'Divergence' <p>
   * \brief Aim: Computes the divergence at a point
   * of an image. 
   *
   * @code 
   * @endcode
   *
   * @tparam TDifference  type of directionnal differential operator 
   */
  template <typename TDifference>
  class Divergence
  {


    // ----------------------- Types ------------------------------
  public:

    typedef TDifference FiniteDifference;
    typedef typename FiniteDifference::Dimension Dimension; 
    static const typename FiniteDifference::Dimension dimension = FiniteDifference::dimension;
    typedef typename FiniteDifference::Image Image; 
    typedef typename FiniteDifference::Point Point; 
    typedef typename FiniteDifference::OutputValue OutputValue;
 
    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     *
     * @param aD a finite difference operator
     */
    Divergence( const FiniteDifference& aD ): myD(aD) {}

    /**
     * Constructor.
     *
     * @param aStartingImage  any image of signed values
     * @param aGridStep any length (=1 by default)
     */
    Divergence( typename FiniteDifference::Image& aStartingImage, 
	      const typename FiniteDifference::OutputValue& aGridStep = 1);

    /**
     * Destructor. Does nothing.
     */
    ~Divergence() {}

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const {return true;}

    /**
     * Divergence
     *
     * @param aPoint the point where the divergence is computed
     * @return gradient of @a myU at @ aPoint
     */
    OutputValue operator() ( const Point& aPoint ) const; 

    // ------------------------- Private Datas --------------------------------
  private:

    /**
     * Finite difference
     */
    FiniteDifference myD; 

  }; 

  /////////////////////////////////////////////////////////////////////////////
  // template class GradientModulus
  /**
   * Description of template class 'GradientModulus' <p>
   * \brief Aim: Computes the gradient modulus at a point
   * of an image. 
   *
   * @code 
   * @endcode
   *
   * @tparam TGradient type of gradient 
   */
  template <typename TGradient>
  class GradientModulus
  {


    // ----------------------- Types ------------------------------
  public:

    typedef TGradient Gradient;
    typedef typename Gradient::Dimension Dimension; 
    static const typename Gradient::Dimension dimension = Gradient::dimension;
    typedef typename Gradient::Point Point; 
   
    typedef double OutputValue;
 
    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     *
     * @param aG any gradient operator
     */
    GradientModulus( const Gradient& aG) : myG(aG) {}

    /**
     * Constructor.
     *
     * @param aStartingImage  any image of signed values
     * @param aGridStep any length (=1 by default)
     */
    GradientModulus( typename Gradient::Image& aStartingImage, 
		     const typename Gradient::OutputValue::Component& aGridStep = 1);


    /**
     * Destructor. Does nothing.
     */
    ~GradientModulus() {}

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const {return true;}

    /**
     * Gradient modulus
     *
     * @param aPoint the point where the gradient modulus is computed
     * @return gradient modulus of @a myU at @ aPoint
     */
    OutputValue operator() ( const Point& aPoint ) const; 

    // ------------------------- Private Datas --------------------------------
  private:

    /**
     * Gradient
     */
    Gradient myG; 

  }; 


  ///////////////////////////////////////////////////////////////////



} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DifferentialOperators.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined DifferentialOperators_h

#undef DifferentialOperators_RECURSES
#endif // else defined(DifferentialOperators_RECURSES)
