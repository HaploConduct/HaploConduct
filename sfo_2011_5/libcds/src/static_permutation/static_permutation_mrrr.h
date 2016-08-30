/* static_permutation_mrrr.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * Permutation
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef _STATIC_PERMUTATION_MRRR_H
#define _STATIC_PERMUTATION_MRRR_H

#include <basics.h>
#include <static_permutation.h>
#include <perm.h>

/** Wrapper for Diego Arroyuelo's implementation of Munro et al.'s permutations.
 *  @author Francisco Claude
 */
class static_permutation_mrrr : public static_permutation {
  public:
    static_permutation_mrrr(uint * elems, uint nelems, uint t, static_bitsequence_builder * bmb);
    virtual ~static_permutation_mrrr();
    /** Computes the i-th element of the permutation */
    virtual uint pi(uint i);
    /** Computes the inverse of i */
    virtual uint rev_pi(uint i);
    /** Saves the permutation to fp, returns 0 in case of success */
    virtual uint save(FILE *fp);
    /** Returns the size of the permutation */
    virtual uint size();
    /** Loads a static_permutation from fp */
    static static_permutation_mrrr * load(FILE *fp);
  protected:
    perm permutation;
    static_permutation_mrrr();
};

#endif
