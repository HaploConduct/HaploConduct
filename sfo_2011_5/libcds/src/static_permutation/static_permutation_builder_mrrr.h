/* static_permutation_builder_mrrr.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * Permutation builder
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

#ifndef _STATIC_PERMUTATION_BUILDER_MRRR_H
#define _STATIC_PERMUTATION_BUILDER_MRRR_H

#include <static_permutation_builder.h>
#include <static_bitsequence.h>
#include <static_bitsequence_builder.h>

/** Base class for static permutation builders
 *  @author Francisco Claude
 */
class static_permutation_builder_mrrr : public static_permutation_builder {
  public:
    static_permutation_builder_mrrr(uint t, static_bitsequence_builder * bmb);
    virtual ~static_permutation_builder_mrrr();
    /** Returns a new permutation build for perm */
    virtual static_permutation * build(uint * perm, uint len);
    
   protected:
    uint t;
    static_bitsequence_builder * bmb;
};

#include <static_permutation_builder_mrrr.h>

#endif
