/* static_sequence_builder_wvtree.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * Sequence builder wavelet tree
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

#ifndef _STATIC_SEQUENCE_BUILDER_WVTREE_H
#define _STATIC_SEQUENCE_BUILDER_WVTREE_H

#include <basics.h>
#include <wt_coder.h>
#include <alphabet_mapper.h>
#include <static_sequence_wvtree.h>
#include <static_sequence_builder.h>

/** Builder for wavelet trees
 *  @author Francisco Claude
 */
class static_sequence_builder_wvtree : public static_sequence_builder {
  public:
    static_sequence_builder_wvtree(wt_coder * wc, static_bitsequence_builder *bmb, alphabet_mapper * am);
    virtual ~static_sequence_builder_wvtree() {}
    virtual static_sequence * build(uint * seq, uint len);
    
  protected:
    alphabet_mapper * am;
    wt_coder * wc;
    static_bitsequence_builder *bmb;
};

#endif
