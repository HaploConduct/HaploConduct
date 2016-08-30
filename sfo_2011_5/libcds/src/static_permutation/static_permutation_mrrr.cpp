/* static_permutation_mrrr.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * Permutation
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


#include <static_permutation_mrrr.h>

static_permutation_mrrr::static_permutation_mrrr(uint * elems, uint nelems, uint t, static_bitsequence_builder * bmb) {
  permutation = createPerm(elems, nelems, t, bmb);
}

static_permutation_mrrr::static_permutation_mrrr() {
}

static_permutation_mrrr::~static_permutation_mrrr() {
  destroyPerm(permutation);
}

uint static_permutation_mrrr::size() {
  return sizeof(static_permutation)+sizeofPerm(permutation);
}

uint static_permutation_mrrr::pi(uint i) {
  return getelemPerm(permutation,i);
}

uint static_permutation_mrrr::rev_pi(uint i) {
  return inversePerm(permutation,i);
}

uint static_permutation_mrrr::save(FILE *fp) {
	uint wr = STATIC_PERMUTATION_MRRR_HDR;
	wr = fwrite(&wr,sizeof(uint),1,fp);
	if(wr!=1) return 1;
  return savePerm(permutation,fp);
}

static_permutation_mrrr * static_permutation_mrrr::load(FILE *fp) {
  uint rd;
  if(fread(&rd,sizeof(uint),1,fp)!=1) return NULL;
  if(rd!=STATIC_PERMUTATION_MRRR_HDR) return NULL;
  static_permutation_mrrr * ret = new static_permutation_mrrr();
  ret->permutation = loadPerm(fp);
  return ret;
}
