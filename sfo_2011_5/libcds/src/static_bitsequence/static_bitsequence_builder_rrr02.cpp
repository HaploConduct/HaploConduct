/* static_bitsequence_builder_rrr02.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_bitsequence_builder_rrr02 definition
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

#include <static_bitsequence_builder_rrr02.h>

static_bitsequence_builder_rrr02::static_bitsequence_builder_rrr02(uint sampling) {
  sample_rate=sampling;
}

static_bitsequence * static_bitsequence_builder_rrr02::build(uint * bitsequence, uint len) {
  return new static_bitsequence_rrr02(bitsequence,len,sample_rate);
}
