/* static_permutation.h
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

#ifndef _STATIC_PERMUTATION_H
#define _STATIC_PERMUTATION_H

#include <basics.h>

#define STATIC_PERMUTATION_MRRR_HDR 2

/** Base class for static permutations
 *  @author Francisco Claude
 */
class static_permutation {
  public:
    virtual ~static_permutation() {}
    /** Computes the i-th element of the permutation */
    virtual uint pi(uint i)=0;
    /** Computes the inverse of i */
    virtual uint rev_pi(uint i)=0;
    /** Saves the permutation to fp, returns 0 in case of success */
    virtual uint save(FILE *fp)=0;
    /** Returns the size of the permutation */
    virtual uint size()=0;
    /** Loads a static_permutation from fp */
    static static_permutation * load(FILE *fp);
};

#include <static_permutation_mrrr.h>

#endif
