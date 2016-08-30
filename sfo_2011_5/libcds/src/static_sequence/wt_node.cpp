/* wt_node.cpp
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

#include <wt_node.h>

wt_node * wt_node::load(FILE *fp) {
  uint rd;
  if(fread(&rd,sizeof(uint),1,fp)!=1) return NULL;
  if(rd==WT_NODE_NULL_HDR) return NULL;
  fseek(fp,-sizeof(uint),SEEK_CUR);
  switch(rd) {
    case WT_NODE_INTERNAL_HDR: return wt_node_internal::load(fp);
    case WT_NODE_LEAF_HDR: return wt_node_leaf::load(fp);
  }
  return NULL;
}
