/* wt_node_internal.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * wt_node_internal
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
 
#ifndef wt_node_internal_h
#define wt_node_internal_h

#include <wt_node.h>
#include <wt_node_leaf.h>
#include <wt_coder.h>
#include <basics.h>
#include <static_bitsequence.h>
#include <static_bitsequence_builder.h>
#include <cassert>

/** Clase for representing internal nodes  
 * 
 *  @author Francisco Claude
 */
class wt_node_internal: public wt_node {
	public:
		wt_node_internal(uint * seq, uint n, uint l, wt_coder * c, static_bitsequence_builder * bmb);
		wt_node_internal(uchar * seq, uint n, uint l, wt_coder * c, static_bitsequence_builder * bmb, uint, uint *);
		virtual ~wt_node_internal();
		virtual uint rank(uint symbol, uint pos, uint level, wt_coder * c);
		virtual uint rankLessThan(uint &symbol, uint pos);
		virtual uint select(uint symbol, uint pos, uint level, wt_coder * c);
		virtual uint access(uint pos);
                virtual uint access(uint pos, uint &rank);
                virtual void access(vector<int> &result, uint i, uint j, uint min, uint max, uint l, uint pivot);
                virtual void access(vector<int> &result, uint i, uint j);
                virtual uint access(uint i, uint j, uint min, uint max, uint l, uint pivot);
		virtual uint size();
    virtual uint save(FILE *fp);
    static wt_node_internal * load(FILE *fp);

		
	protected:
    wt_node_internal();
		wt_node *left_child, *right_child;
		static_bitsequence * bitmap;
		//uint length;
};

#endif
