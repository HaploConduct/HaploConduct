/* static_bitsequence_builder.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_bitsequence_builder definition
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

#ifndef _STATIC_BITSEQUENCE_BUILDER_SDARRAY_H
#define _STATIC_BITSEQUENCE_BUILDER_SDARRAY_H

#include <static_bitsequence.h>
#include <static_bitsequence_builder.h>

class static_bitsequence_builder_sdarray : public static_bitsequence_builder {
  public:
		static_bitsequence_builder_sdarray() {}
    virtual ~static_bitsequence_builder_sdarray() {}
    /** Builds a static_bitsequence for the bitmap bitsequence of length len */
    virtual static_bitsequence * build(uint * bitsequence, uint len);
};

#endif /* _STATIC_BITSEQUENCE_BUILDER_H */
