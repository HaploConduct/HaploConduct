/* static_sequence_builder_gmr.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * Sequence builder gmr 
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

#include <static_sequence_builder_gmr.h>

static_sequence_builder_gmr::static_sequence_builder_gmr(uint chunk_length, static_bitsequence_builder *bmb, static_sequence_builder *ssb) {
  this->chunk_length = chunk_length;
  this->bmb = bmb;
  this->ssb = ssb;
}

static_sequence * static_sequence_builder_gmr::build(uint * seq, uint len) {
  return new static_sequence_gmr(seq,len,chunk_length,bmb,ssb);
}
