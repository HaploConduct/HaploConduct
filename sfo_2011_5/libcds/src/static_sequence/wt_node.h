/* wt_node.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * wt_node
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
 
#ifndef wt_node_h
#define wt_node_h

#include <basics.h>
#include <wt_coder.h>
#include <vector>

#define WT_NODE_NULL_HDR 0
#define WT_NODE_INTERNAL_HDR 2
#define WT_NODE_LEAF_HDR 3

/** Base clase for nodes in the wavelet tree  
 * 
 *  @author Francisco Claude
 */
class wt_node {
	public:
    virtual ~wt_node() {} 
		virtual uint rank(uint symbol, uint pos, uint l, wt_coder * c)=0;
		virtual uint rankLessThan(uint &symbol, uint pos) = 0;
		virtual uint select(uint symbol, uint pos, uint l, wt_coder * c)=0;
		virtual uint access(uint pos)=0;
		virtual uint access(uint pos, uint &rank)
                {
                    assert(0); // Implemented only in wt_node_internal
                    return -1;
                }
                virtual void access(std::vector<int> &result, uint i, uint j, uint min, uint max, uint l, uint pivot)=0;
                virtual void access(std::vector<int> &result, uint i, uint j)=0;
                virtual uint access(uint i, uint j, uint min, uint max, uint l, uint pivot)=0;
		virtual uint size()=0;
    virtual uint save(FILE *fp)=0;
    static wt_node * load(FILE *fp);
};

#include <wt_node_internal.h>
#include <wt_node_leaf.h>

#endif
