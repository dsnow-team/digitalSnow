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
 * @file FrontierEvolverHelpers.h
 *
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 * @date 2012/04/19
 *
 * This file is part of the DGtal library.
 */

#if defined(FrontierEvolverHelpers_RECURSES)
#error Recursive header files inclusion detected in FrontierEvolverHelpers.h
#else // defined(FrontierEvolverHelpers_RECURSES)
/** Prevents recursive inclusion of headers. */
#define FrontierEvolverHelpers_RECURSES

#if !defined FrontierEvolverHelpers_h
/** Prevents repeated inclusion of headers. */
#define FrontierEvolverHelpers_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////


//#define WITHINFO

namespace DGtal
{


  //------------------------------------------------------------------------------
  template<typename I, typename D, typename V>
struct SetFromImageDomainValueTraits
  { 
  public: 
    typedef DigitalSetBySTLSet<D> Set; 
  public: 
  public: 
    static Set get(I& aImage)
    {
      return Set(aImage.domain()); 
    }
  }; 
  //Partial specialization
  template<typename D, typename V>
struct SetFromImageDomainValueTraits<
    ImageContainerBySTLMap<D,V>, 
    D, V >
  {
  public: 
    typedef DigitalSetFromMap<ImageContainerBySTLMap<D,V> > Set; 
  public: 
    static Set get(ImageContainerBySTLMap<D,V>& aImage)
    {
      return Set(aImage); 
    }
  }; 
  //------------------------------------------------------------------------------
  template<typename I>
struct SetFromImageSelector
  { 
  public: 
    BOOST_CONCEPT_ASSERT(( concepts::CImage<I> )); 
    typedef typename I::Domain Domain; 
    typedef typename I::Value Value; 

    typedef typename SetFromImageDomainValueTraits<I, Domain, Value>::Set Set;

  public: 
    static Set get(I& i) 
    {
      return SetFromImageDomainValueTraits<I, Domain, Value>::get(i); 
    }
  }; 

  //------------------------------------------------------------------------------
  namespace details
  {
    class CompareSecondElement {
    public: 
      template <typename T>
      bool operator()(const T& a, const T& b) 
      {
	return ( a.second < b.second ); 
      }
    };

  }

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
//#include "FrontierEvolverHelpers.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined FrontierEvolverHelpers_h

#undef FrontierEvolverHelpers_RECURSES
#endif // else defined(FrontierEvolverHelpers_RECURSES)
