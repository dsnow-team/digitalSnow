#ifndef FFT_H
#define FFT_H

#include <fftw3.h>

#include <complex>

#include <DGtal/images/CImage.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

namespace DGtal{

  /////////////////////////////////////////////////////////////////////////////
  // template class FFT
  /**

   * @tparam TImage, the type of image 
   */
  template <typename TImage>
  class FFT
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

    typedef typename Image::Dimension Dimension;
    static const typename Image::Dimension dimension = Image::dimension;
    
    typedef std::complex<double> Complex; 
    typedef ImageContainerBySTLVector<  Domain, Complex > ComplexImage;

    


    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     * @param anImage any image
     */
    FFT( const Image& anImage);

    /**
     * Destructor. Does nothing.
     */
    ~FFT();

    /**
     * Computes the transform
     *
     * @param anImage the returned output image
     */
    void compute (ComplexImage& anImage);


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

  }; // end of class FFT

} //end namespace DGtal 




///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "FFT.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif

