/***
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that neither the name of Emin
 * Martinian nor the names of any contributors are be used to endorse or
 * promote products derived from this software without specific prior
 * written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "misc.h"

/***********************************************************************/
/*  FUNCTION:  void Assert(int assertion, char* error)  */
/**/
/*  INPUTS: assertion should be a predicated that the programmer */
/*  assumes to be true.  If this assumption is not true the message */
/*  error is printed and the program exits. */
/**/
/*  OUTPUT: None. */
/**/
/*  Modifies input:  none */
/**/
/*  Note:  If DEBUG_ASSERT is not defined then assertions should not */
/*         be in use as they will slow down the code.  Therefore the */
/*         compiler will complain if an assertion is used when */
/*         DEBUG_ASSERT is undefined. */
/***********************************************************************/
void Assert(int assertion, char const *error) {
  if(!assertion) {
    std::string buf("Assertion Failed: ");
    throw Rcpp::exception(buf.append(error).append("\n").c_str());
  }
}

/***********************************************************************/
/*  FUNCTION:  SafeMalloc */
/**/
/*    INPUTS:  size is the size to malloc */
/**/
/*    OUTPUT:  returns pointer to allocated memory if succesful */
/**/
/*    EFFECT:  mallocs new memory.  If malloc fails, prints error message */
/*             and terminates program. */
/**/
/*    Modifies Input: none */
/**/
/***********************************************************************/
void * SafeMalloc(size_t size) {
  void * result;

  if ( (result = malloc(size)) ) { /* assignment intentional */
    return(result);
  } else {
    throw Rcpp::exception("Memory overflow: malloc failed in SafeMalloc.");
    return(0);
  }
}
/*  NullFunction does nothing it is included so that it can be passed */
/*  as a function to RBTreeCreate when no other suitable function has */
/*  been defined */
void NullFunction(void * junk) { }
