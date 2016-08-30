/* static_sequence_builder.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * Sequence builder
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

#ifndef _STATIC_SEQUENCE_BUILDER_H
#define _STATIC_SEQUENCE_BUILDER_H

#include <basics.h>
#include <static_sequence.h>

/** Base class for static sequence builders
 *  @author Francisco Claude
 */
class static_sequence_builder {
  public:
    virtual ~static_sequence_builder() {}
    /** Returns a new sequence build for seq */
    virtual static_sequence * build(uint * seq, uint len)=0;
};

#include <static_sequence_builder_wvtree.h>
#include <static_sequence_builder_wvtree_noptrs.h>
#include <static_sequence_builder_gmr.h>
#include <static_sequence_builder_gmr_chunk.h>

#endif
