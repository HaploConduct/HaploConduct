/* huff.cpp
   Copyright (C) 2008, Gonzalo Navarro, all rights reserved.

   Canonical Huffman

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
// implements canonical Huffman 

#include <huff.h>

typedef struct
   { uint freq;
     uint symb;
     union
       { int prev;
	 uint depth;
       } h;
     int ch1,ch2;
   } Ttree;
  
static void sort (Ttree *tree, int lo, int up)

   { uint i, j;
     Ttree temp;
     while (up>lo) 
        { i = lo;
	  j = up;
	  temp = tree[lo]; 
	  while (i<j) 
	     { while (tree[j].freq > temp.freq) j--;
	       tree[i] = tree[j]; 
	       while (i<j && tree[i].freq <= temp.freq) i++;
	       tree[j] = tree[i];
	     }
	  tree[i] = temp;
	  if (i-lo < up-i) { sort(tree,lo,i-1); lo = i+1; }
	  else { sort(tree,i+1,up); up = i-1; }
	}
   }

static void setdepths (Ttree *tree, uint node, int depth)

   { if (tree[node].ch1 == -1) // leaf
	{ tree[node].h.depth = depth; 
	  return;
	}
     setdepths (tree,tree[node].ch1,depth+1);
     setdepths (tree,tree[node].ch2,depth+1);
   }

THuff createHuff (uint *freq, uint lim)

   { THuff H;
     int i,j,d;
     Ttree *tree;
     uint ptr,last,fre;
         // remove zero frequencies
     H.max = lim;
     tree = (Ttree*)malloc((2*(lim+1)-1)*sizeof(Ttree));
     j = 0;
     for (i=0;i<=(int)lim;i++) 
	{ if (freq[i]>0) 
	     { tree[j].freq = freq[i];   
               tree[j].symb = i;   
	       j++;
	     }
	}
     H.lim = lim = j-1;
     	// now run Huffman algorithm
     sort (tree,0,lim);
     for (i=0;i<=(int)lim;i++) 
         { tree[i].h.prev = i+1;
	   tree[i].ch1 = tree[i].ch2 = -1;
	 }
     tree[lim].h.prev = -1;
        // last = next node to process, ptr = search point, fre = next free cell
	// leaves are in 0..lim in decreasing freq order
	// internal nodes are in lim+1.. 2*lim, created in incr. fre order
     last=0; ptr = 0; fre = lim+1;
     for (i=0;i<(int)lim;i++)
	{ tree[fre].ch1 = last;
	  last = tree[last].h.prev;
	  tree[fre].ch2 = last;
	  tree[fre].freq = tree[tree[fre].ch1].freq+tree[tree[fre].ch2].freq;
	  while ((tree[ptr].h.prev != -1) &&
		 (tree[tree[ptr].h.prev].freq <= tree[fre].freq))
	      ptr = tree[ptr].h.prev;
	  tree[fre].h.prev = tree[ptr].h.prev;
	  tree[ptr].h.prev = fre;
	  last = tree[last].h.prev;
	  fre++;
	}
        // now assign depths recursively
     setdepths (tree,2*lim,0);
     H.s.spos = (uint*)malloc ((H.max+1)*sizeof(uint));
     for (i=0;i<=(int)H.max;i++) H.s.spos[i] = ~0;
     H.num = (uint*)malloc ((lim+1)*sizeof(uint)); // max possible depth
     d=0;
     for (i=lim;i>=0;i--)
	{ H.s.spos[tree[i].symb] = i;
	  while ((int)tree[i].h.depth > d) 
	      { H.num[d] = i+1; d++; }
	}
     H.num[d] = 0;
     H.depth = d;
     for (d=H.depth;d>0;d--) H.num[d] = H.num[d-1] - H.num[d];
     H.num[0] = (lim == 0);
     H.num = (uint*)realloc(H.num,(H.depth+1)*sizeof(uint));
     H.total = 0;
     for (i=0;i<=(int)lim;i++) 
       H.total += freq[tree[i].symb] * tree[i].h.depth;
     free (tree);
     return H;
   }

void bitzero (register uint *e, register uint p, 
         register uint len)

   { e += p/W; p %= W;
     if (p+len >= W)
  { *e &= ~((1<<p)-1);
    len -= p;
    e++; p = 0;
  }
     while (len >= W)
  { *e++ = 0;
    len -= W;
  }
     if (len > 0)
  *e &= ~(((1<<len)-1)<<p);
   }
   
ulong encodeHuff (THuff H, uint symb, uint *stream, ulong ptr)

   { uint pos;
     uint code;
     uint d;
     pos = H.s.spos[symb];
     code = 0;
     d = H.depth; 
     while (pos >= H.num[d])
       { code = (code + H.num[d]) >> 1;
         pos -= H.num[d--];
       }
     code += pos;
     if (d > W) { bitzero(stream,ptr,d-W); ptr += d-W; d = W; }
     while (d--)
        { if ((code >> d) & 1) bitset(stream,ptr);
	  else bitclean(stream,ptr);
	  ptr++;
	}
     return ptr;
   } 

ulong decodeHuff (THuff H, uint *symb, uint *stream, ulong ptr)

   { uint pos;
     uint d;
     pos = 0;
     d = 0;
     while (pos < H.fst[d])
        { pos = (pos << 1) | bitget(stream,ptr);
          ptr++; d++;
        }
     *symb = H.s.symb[H.num[d]+pos-H.fst[d]];
     return ptr;
   }

/*   { uint pos;  // This "improved" code is actually slower!
     int d;
     uint wrd,off;
     stream += ptr/W;
     off = ptr & (W-1);
     wrd = *stream >> off;
     pos = 0;
     d = 0;
     while (pos < H.fst[d])
       { pos = (pos << 1) | (wrd & 1);
	 d++; wrd >>= 1; off++;
	 if (off == W) { wrd = *++stream; off = 0; }
       }
     *symb = H.s.symb[H.num[d]+pos-H.fst[d]];
     return ptr+d;
   } 
*/
void saveHuff (THuff H, FILE *f)

   { uint *symb = new uint[H.lim+1];
     uint i;
		 for(i=0;i<(H.lim+1);i++) symb[i] = 0;
     for (i=0;i<=H.max;i++) 
	 		 if (H.s.spos[i] != (uint)~0) symb[H.s.spos[i]] = i;
     uint l=fwrite (&H.max,sizeof(uint),1,f); 
     l += fwrite (&H.lim,sizeof(uint),1,f); 
     l += fwrite (&H.depth,sizeof(uint),1,f); 
     l += fwrite (symb,sizeof(uint),H.lim+1,f); 
     l += fwrite (H.num,sizeof(uint),H.depth+1,f); 
     delete [] (symb);
   }

uint sizeHuff (THuff H)

   { return (4+(H.lim+1)+(H.depth+1))*sizeof(uint);
   }

void freeHuff (THuff H)

  { free (H.s.spos); free (H.num);
  }

THuff loadHuff (FILE *f, int enc)

   { THuff H;
     uint *symb;
     //uint *num;
     uint i,d,dold,dact;
     uint l = fread (&H.max,sizeof(uint),1,f); 
     l += fread (&H.lim,sizeof(uint),1,f); 
     l += fread (&H.depth,sizeof(uint),1,f); 
     symb = (uint*)malloc ((H.lim+1)*sizeof(uint));
     l += fread (symb,sizeof(uint),H.lim+1,f); 
     if (enc) 
          { H.s.spos = (uint*)malloc ((H.max+1)*sizeof(uint));
            for (i=0;i<=H.max;i++) H.s.spos[i] = (uint)~0;
            for (i=0;i<=H.lim;i++) H.s.spos[symb[i]] = i;
            free (symb);
	  }
     else H.s.symb = symb;
     H.num = (uint*)malloc ((H.depth+1)*sizeof(uint));
     l += fread (H.num,sizeof(uint),H.depth+1,f); 
     if (!enc) 
          { H.fst = (uint*)malloc ((H.depth+1)*sizeof(uint));
            H.fst[H.depth] = 0; dold = 0;
            for (d=H.depth-1;d>=0;d--)
	        { dact = H.num[d+1];
		  H.fst[d] = (H.fst[d+1]+dact) >> 1;
		  H.num[d+1] = dold;
		  dold += dact;
		}
	    H.num[0] = dold;
	  }
     return H;
   }
