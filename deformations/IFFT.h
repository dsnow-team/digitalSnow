#ifndef IFFT_H
#define IFFT_H

#include <fftw3.h>

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <DGtal/images/CImage.h>

namespace DGtal{

  /////////////////////////////////////////////////////////////////////////////
  // template class IFFT
  /**

   * @tparam TImage, the type of image 
   */
  template <typename TImage>
  class IFFT
  {

    //ASSERT
    BOOST_CONCEPT_ASSERT(( concepts::CImage<TImage> )); 
    BOOST_STATIC_ASSERT((boost::is_same< typename TImage::Value, 
			 std::complex<double> >::value));

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
     * @param anImage any image
     */
    IFFT( const Image& anImage);

    /**
     * Destructor. Does nothing.
     */
    ~IFFT();

    /**
     * Computes the transform
     *
     * @param anImage the returned output image
     * @tparam TOutputImage type of returned image
     */
    template < typename TOutputImage >
    void compute (TOutputImage& anImage);


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
     * Reference on the image
     */
    const Image& myImage; 


    // ------------------------- Hidden services ------------------------------
  protected:


  private:



    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class IFFT

} //end namespace DGtal 




///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "IFFT.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif

