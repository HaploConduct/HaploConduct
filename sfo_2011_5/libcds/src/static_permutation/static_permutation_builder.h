/* static_permutation_builder.h
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

#ifndef _STATIC_PERMUTATION_BUILDER_H
#define _STATIC_PERMUTATION_BUILDER_H

#include <basics.h>
#include <static_permutation.h>

/** Base class for static permutation builders
 *  @author Francisco Claude
 */
class static_permutation_builder {
  public:
    virtual ~static_permutation_builder() {}
    /** Returns a new permutation build for perm */
    virtual static_permutation * build(uint * perm, uint len)=0;
};

#include <static_permutation_builder_mrrr.h>

#endif
