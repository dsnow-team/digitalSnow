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
 * @file CPointFunctorAdpter.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/03/01
 *
 * Header file for concept CPointFunctorAdpter.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(CPointFunctorAdpter_RECURSES)
#error Recursive header files inclusion detected in CPointFunctorAdpter.h
#else // defined(CPointFunctorAdpter_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CPointFunctorAdpter_RECURSES

#if !defined CPointFunctorAdpter_h
/** Prevents repeated inclusion of headers. */
#define CPointFunctorAdpter_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"

#include "DGtal/kernel/CPointFunctor.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class CPointFunctorAdpter
  /**
     Description of \b concept '\b CPointFunctorAdpter' <p>
     @ingroup Concepts
     @brief Aim: Defines the concept describing a point functor, 
    adapted from another one, which can be returned by the base() method
    and can be set the attach() method. 

     <p> Refinement : 
     - \t PointFunctor
   
     <p> Nested types : 
     - \t PointFunctor : a model of CPointFunctor 
  
     <p> Notation
     - \t X : a model of CPointFunctorAdpter
     - \t x : object of type X
  
     <p> Definitions
    
     <p> Valid expressions and semantics <br>
     <table> 
      <tr> 
        <td class=CName> \b Name </td> 
        <td class=CExpression> \b Expression </td>
        <td class=CRequirements> \b Type requirements </td> 
        <td class=CReturnType> \b Return type </td>
        <td class=CPrecondition> \b Precondition </td> 
        <td class=CSemantics> \b Semantics </td> 
        <td class=CPostCondition> \b Postcondition </td> 
        <td class=CComplexity> \b Complexity </td>
      </tr>
      <tr> 
        <td class=CName> underlying point functor  </td> 
        <td class=CExpression> x.base()     </td>
        <td class=CRequirements>    </td> 
        <td class=CReturnType> const PointFunctor&     </td>
        <td class=CPrecondition>    </td> 
        <td class=CSemantics> returns the underlying point functor </td> 
        <td class=CPostCondition>       </td> 
        <td class=CComplexity>   </td>
      </tr>
      <tr> 
        <td class=CName> underlying point functor  </td> 
        <td class=CExpression> x.attach(const PointFunctor&)     </td>
        <td class=CRequirements>    </td> 
        <td class=CReturnType>     </td>
        <td class=CPrecondition>    </td> 
        <td class=CSemantics> set the underlying point functor </td> 
        <td class=CPostCondition>       </td> 
        <td class=CComplexity>   </td>
      </tr>
     </table>
    
     <p> Invariants <br>
    
     <p> Models <br>
     
    
     <p> Notes <br>

     @tparam X the type that should be a model of CPointFunctorAdpter.
   */
  template <typename X> 
  struct CPointFunctorAdpter : CPointFunctor<X> 
  {
    // ----------------------- Concept checks ------------------------------
  public:
    // Inner types
    typedef typename X::PointFunctor PointFunctor;
    BOOST_CONCEPT_ASSERT(( CPointFunctor <PointFunctor> ));
 
    // Methods
    BOOST_CONCEPT_USAGE( CPointFunctorAdpter )
    {
      ConceptUtils::sameType( myPF, myX.base() );
      myX.attach( myPF ); 
    }
    // ------------------------- Private Datas --------------------------------
  private:
    X myX;
    const PointFunctor& myPF;  
  
    // ------------------------- Internals ------------------------------------
  private:
    
  }; // end of concept CPointFunctorAdpter
  
} // namespace DGtal

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CPointFunctorAdpter_h

#undef CPointFunctorAdpter_RECURSES
#endif // else defined(CPointFunctorAdpter_RECURSES)
