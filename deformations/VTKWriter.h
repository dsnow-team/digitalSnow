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
 * @file VTKWriter.h
 * @author Roland Denis (\c roland.denis@univ-savoie.fr )
 * LAboratory of MAthematics - LAMA (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2015/01/22
 *
 * This file is part of the DGtal library.
 */

#if defined(VTKWriter_RECURSES)
#error Recursive header files inclusion detected in VTKWriter.h
#else // defined(VTKWriter_RECURSES)
/** Prevents recursive inclusion of headers. */
#define VTKWriter_RECURSES

#if !defined VTKWriter_h
/** Prevents repeated inclusion of headers. */
#define VTKWriter_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <string>
#include <fstream>

#include <boost/concept/assert.hpp>

#include <DGtal/base/Exceptions.h>
#include <DGtal/kernel/domains/HyperRectDomain.h>
#include <DGtal/images/CConstImage.h>
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class VTKWriter
  /**
   * VTK writer for DGtal Images
   *
   * @tparam TDomain Type of the domain used for export
   * @tparam Binary 'true' for BINARY format, 'false' for ASCII format
   *
   * @see http://www.vtk.org/VTK/img/file-formats.pdf
   */
  template <
    typename TDomain,
    bool Binary = true
  >
  class VTKWriter;

  template <typename TSpace, bool Binary>
  class VTKWriter< HyperRectDomain<TSpace>, Binary >
  {
    // ----------------------- Standard services ------------------------------
  public:

    typedef HyperRectDomain<TSpace> Domain;

    /**
     * Constructor
     *
     * @param filename  name of the file ...
     * @param domain    domain of the datas to be exported.
     */
    VTKWriter(std::string const& filename, Domain domain);

    /**
     * Destructor
     */
    ~VTKWriter();

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    /**
     * Write the VTK header.
     *
     * It is automatically done at the first data export
     * @return a reference to the writer instance
     */
    VTKWriter<Domain,Binary> & init();


    /**
     * Set the name of the next field
     *
     * @param fieldname the name of the next field.
     * @return a reference to the writer instance.
     */
    VTKWriter<Domain,Binary> & operator<< ( std::string const& fieldname );
    VTKWriter<Domain,Binary> & operator<< ( const char* fieldname );

    /**
     * Write a field
     *
     * @tparam TImage type of the image.
     * @param field   the image.
     * @return a reference to the writer instance.
     *
     * @todo enable only if TImage is a model a concepts::CImage !!
     */
    template <typename TImage>
    BOOST_CONCEPT_REQUIRES( (( concepts::CConstImage<TImage> )),
    ( VTKWriter<HyperRectDomain<TSpace>,Binary> & ))
    operator<< ( TImage const& field );

    /**
     * Write a field, given his name and optionally specifying export type.
     *
     * By specifying T template, it is possible to cast values to a specified
     * type before export. For example, it can export double values as float,
     * to save disk space.
     *
     * @tparam TImage type of the image.
     * @tparam T      type of values.
     * @return a reference to the write instance.
     */
    template <
      typename TImage,
      typename T = typename TImage::Value
    >
    VTKWriter<Domain,Binary> & write( std::string const& fieldname, TImage const& field );
    

    /**
     * Close the file
     */
    void close();

    // ------------------------- Protected Datas ------------------------------
  protected:
    // ------------------------- Private Datas --------------------------------
  private:
    Domain m_domain;
    std::string m_fieldname;
    std::ofstream m_fstream;
    bool m_header;

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    VTKWriter();

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    VTKWriter ( const VTKWriter & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    VTKWriter & operator= ( const VTKWriter & other );

    // ------------------------- Internals ------------------------------------
  private:

    /**
     * Internal data stream
     */
    struct DataStream
      {
        DataStream( std::ofstream & fstream );

        std::string data_format();

        template <typename TValue>
        DataStream& operator<< (TValue const& value);

        void separator();
        
        std::ofstream & m_fstream;
      };

    DataStream m_dataStream;

  }; // end of class VTKWriter

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "VTKWriter.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined VTKWriter_h

#undef VTKWriter_RECURSES
#endif // else defined(VTKWriter_RECURSES)

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

