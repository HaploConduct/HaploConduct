/* wt_node_internal.cpp
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
 
#include <wt_node_internal.h>

wt_node_internal::wt_node_internal(uint * symbols, uint n, uint l, wt_coder * c, static_bitsequence_builder * bmb) {
	uint * ibitmap = new uint[n/W+1];
	for(uint i=0;i<n/W+1;i++)
		ibitmap[i]=0;
	for(uint i=0;i<n;i++) 
		if(c->is_set(symbols[i],l))
			bitset(ibitmap,i);
	bitmap = bmb->build(ibitmap, n);
  delete [] ibitmap;
	uint count_right = bitmap->rank1(n-1);
	uint count_left = n-count_right+1;
	uint * left = new uint[count_left+1];
	uint * right = new uint[count_right+1];
	count_right = count_left = 0;
	bool match_left = true, match_right = true;
	for(uint i=0;i<n;i++) {
		if(bitmap->access(i)) {
			right[count_right++]=symbols[i];
			if(count_right>1)
				if(right[count_right-1]!=right[count_right-2])
					match_right = false;
		}
		else {
			left[count_left++]=symbols[i];
			if(count_left>1)
				if(left[count_left-1]!=left[count_left-2])
					match_left = false;
		}
	}
	if(count_left>0) {
		if(match_left/* && c->done(left[0],l+1)*/)
			left_child = new wt_node_leaf(left[0], count_left);
		else
			left_child = new wt_node_internal(left, count_left, l+1, c, bmb);
	} else {
		left_child = NULL;
	}
	if(count_right>0) {
		if(match_right/* && c->done(right[0],l+1)*/)
			right_child = new wt_node_leaf(right[0], count_right);
		else
			right_child = new wt_node_internal(right, count_right, l+1, c, bmb);
	} else {
		right_child = NULL;
	}
	delete [] left;
	delete [] right;
}

wt_node_internal::wt_node_internal(uchar * symbols, uint n, uint l, wt_coder * c, static_bitsequence_builder * bmb, uint left, uint *done) {
	uint * ibitmap = new uint[n/W+1];
	for(uint i=0;i<n/W+1;i++)
		ibitmap[i]=0;
	for(uint i=0;i<n;i++) 
		if(c->is_set((uint)symbols[i + left],l))
			bitset(ibitmap,i);
	bitmap = bmb->build(ibitmap, n);
        delete [] ibitmap;

	uint count_right = bitmap->rank1(n-1);
	uint count_left = n-count_right;

        for (uint i=0;i<n;i++)
            set_field(done, 1, i+left, 0);

        for (uint i = 0; i < n; ) 
        {
            uint j = i;
            uchar swap = symbols[j+left];
            while (!get_field(done, 1, j+left)) { // swapping
                ulong k = j; 
                if (!c->is_set(swap,l)) 
                    j = bitmap->rank0(k)-1;
                else 
                    j = count_left + bitmap->rank1(k)-1;
                uchar temp = symbols[j+left];
                symbols[j+left] = swap;
                swap = temp;
                set_field(done,1,k+left,1);
            }

            while (get_field(done,1,i+left))
                   ++i;
        }

	bool match_left = true, match_right = true;
        for (uint i=1; i < count_left; i++)
            if (symbols[i+left] != symbols[i+left-1])
                match_left = false;
        for (uint i=count_left + 1; i < n; i++)
            if (symbols[i+left] != symbols[i+left-1])
                match_right = false;


	if(count_left>0) {
		if(match_left/* && c->done(left[0],l+1)*/)
                    left_child = new wt_node_leaf((uint)symbols[left], count_left);
		else
                    left_child = new wt_node_internal(symbols, count_left, l+1, c, bmb, left, done);
	} else {
		left_child = NULL;
	}
	if(count_right>0) {
		if(match_right/* && c->done(right[0],l+1)*/)
                    right_child = new wt_node_leaf((uint)symbols[left+count_left], count_right);
		else 
                    right_child = new wt_node_internal(symbols, count_right, l+1, c, bmb, left+count_left, done);
	} else {
		right_child = NULL;
	}
}


wt_node_internal::wt_node_internal() { }

wt_node_internal::~wt_node_internal() {
	delete bitmap;
	if(right_child!=NULL) delete right_child;
	if(left_child!=NULL) delete left_child;
}

uint wt_node_internal::rank(uint symbol, uint pos, uint l, wt_coder * c) {
	bool is_set = c->is_set(symbol,l);
	if(!is_set) {
		if(left_child==NULL) return 0;
		return left_child->rank(symbol, bitmap->rank0(pos)-1,l+1,c);
	}
	else {
		if(right_child==NULL) return 0;
		return right_child->rank(symbol, bitmap->rank1(pos)-1,l+1,c);
	}
}

// return value is rank of symbol (less or equal to the given symbol) that has rank > 0, 
// the parameter symbol is updated accordinly
uint wt_node_internal::rankLessThan(uint &symbol, uint pos) 
{
    uint result = -1;
    using std::cout;
    using std::endl;
//    cout << "pos = " << pos << ", symbol = " << symbol << endl;
    
    if (pos == (uint)-1)
        return (uint)-1;
    if(right_child!=NULL)
        result = right_child->rankLessThan(symbol, bitmap->rank1(pos)-1);
    if(result == (uint)-1 && left_child!=NULL)
        return left_child->rankLessThan(symbol, bitmap->rank0(pos)-1);
    return result;
}


uint wt_node_internal::select(uint symbol, uint pos, uint l, wt_coder * c) {
	bool is_set = c->is_set(symbol, l);
	uint ret = 0;
	if(!is_set) {
		if(left_child==NULL)
			return (uint)(-1);
		uint new_pos = left_child->select(symbol, pos, l+1,c);
		if(new_pos+1==0) return (uint)(-1);
		ret = bitmap->select0(new_pos)+1;
	} else {
		if(right_child==NULL)
			return (uint)(-1);
		uint new_pos = right_child->select(symbol, pos, l+1,c);
		if(new_pos+1==0) return (uint)(-1);
		ret = bitmap->select1(new_pos)+1;
	}
	if(ret==0) return (uint)-1;
	return ret;
}

uint wt_node_internal::access(uint pos) {
	bool is_set = bitmap->access(pos);
	if(!is_set) {
		assert(left_child!=NULL);
		return left_child->access(bitmap->rank0(pos)-1);
	} else {
		assert(right_child!=NULL);
		return right_child->access(bitmap->rank1(pos)-1);
	}
}

// Returns the value at given position and its rank
uint wt_node_internal::access(uint pos, uint &rank) 
{
    bool is_set = bitmap->access(pos);
    if(!is_set)
    {
        // recurse left
        pos = bitmap->rank0(pos)-1;
        return left_child->access(pos, rank);
    } 
    else 
    {
        // recurse right
        pos = bitmap->rank1(pos)-1;
        return right_child->access(pos, rank);
    }
}


void wt_node_internal::access(vector<int> &result, uint i, uint j, uint min, uint max, uint l, uint pivot)
{
    uint symbol = pivot | (1 << l);
//    std::cout << "At l = " << l << ", [" << i << ", " << j  << "], [" << min << ", " << max << "], symbol = " << symbol << std::endl;

    if (j < i || max < min)
        return;

    if (min < symbol)
    {
        // Recurse left
        uint newi = 0;
        if (i > 0)
            newi = bitmap->rank0(i - 1);
        uint newj = bitmap->rank0(j);

        uint newmax = max < symbol - 1 ? max : symbol - 1;
        if (left_child != NULL && newj > 0)
            left_child->access(result, newi, newj-1, min, newmax, l-1, pivot);
    }
    
    if (max >= symbol)
    {
        // Recurse right
        uint newi = 0;
        if (i > 0)
            newi = bitmap->rank1(i - 1);
        uint newj = bitmap->rank1(j);

        uint newmin = min > symbol ? min : symbol;
        if (right_child != NULL && newj > 0)
            right_child->access(result, newi, newj-1, newmin, max, l-1, symbol);
    }
}

void wt_node_internal::access(vector<int> &result, uint i, uint j)
{
//    std::cout << "At l = " << l << ", [" << i << ", " << j  << "], [" << min << ", " << max << "], symbol = " << symbol << std::endl;

    if (j < i)
        return;

    {
        // Recurse left
        uint newi = 0;
        if (i > 0)
            newi = bitmap->rank0(i - 1);
        uint newj = bitmap->rank0(j);

        if (left_child != NULL && newj > 0)
            left_child->access(result, newi, newj-1);
    }
    
    {
        // Recurse right
        uint newi = 0;
        if (i > 0)
            newi = bitmap->rank1(i - 1);
        uint newj = bitmap->rank1(j);

        if (right_child != NULL && newj > 0)
            right_child->access(result, newi, newj-1);
    }
}

// Count
uint wt_node_internal::access(uint i, uint j, uint min, uint max, uint l, uint pivot)
{
    uint count = 0;
    uint symbol = pivot | (1 << l);
//    std::cout << "At l = " << l << ", [" << i << ", " << j  << "], [" << min << ", " << max << "], symbol = " << symbol << std::endl;

    if (j < i || max < min)
        return 0;

    if (min < symbol)
    {
        // Recurse left
        uint newi = 0;
        if (i > 0)
            newi = bitmap->rank0(i - 1);
        uint newj = bitmap->rank0(j);

        uint newmax = max < symbol - 1 ? max : symbol - 1;
        if (left_child != NULL && newj > 0)
            count += left_child->access(newi, newj-1, min, newmax, l-1, pivot);
    }
    
    if (max >= symbol)
    {
        // Recurse right
        uint newi = 0;
        if (i > 0)
            newi = bitmap->rank1(i - 1);
        uint newj = bitmap->rank1(j);

        uint newmin = min > symbol ? min : symbol;
        if (right_child != NULL && newj > 0)
            count += right_child->access(newi, newj-1, newmin, max, l-1, symbol);
    }
    return count;
}


uint wt_node_internal::size() {
	uint s = bitmap->size()+sizeof(wt_node_internal);
	if(left_child!=NULL)
		s += left_child->size();
	if(right_child!=NULL)
		s += right_child->size();
	return s;
}

uint wt_node_internal::save(FILE *fp) {
  uint wr = WT_NODE_INTERNAL_HDR;
  wr = fwrite(&wr,sizeof(uint),1,fp);
  if(wr!=1) return 1;
  if(bitmap->save(fp)) return 1;
  if(left_child!=NULL) {
    if(left_child->save(fp)) return 1;
  } else {
    wr = WT_NODE_NULL_HDR;
    wr = fwrite(&wr,sizeof(uint),1,fp);
    if(wr!=1) return 1;
  }
  if(right_child!=NULL) {
    if(right_child->save(fp)) return 1;
  } else {
    wr = WT_NODE_NULL_HDR;
    wr = fwrite(&wr,sizeof(uint),1,fp);
    if(wr!=1) return 1;
  }
  return 0;
}

wt_node_internal * wt_node_internal::load(FILE *fp) {
  uint rd;
  if(fread(&rd,sizeof(uint),1,fp)!=1) return NULL;
  if(rd!=WT_NODE_INTERNAL_HDR) return NULL;
  wt_node_internal * ret = new wt_node_internal();
  ret->bitmap = static_bitsequence::load(fp);
  ret->left_child = wt_node::load(fp);
  ret->right_child = wt_node::load(fp);
  return ret;
}
