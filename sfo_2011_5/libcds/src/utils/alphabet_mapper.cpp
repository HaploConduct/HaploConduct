/* alphabet_mapper.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_bitsequence definition
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

#include <alphabet_mapper.h>

alphabet_mapper::alphabet_mapper() {
	user_count=0;
}

void alphabet_mapper::use() {
	user_count++;
}

void alphabet_mapper::unuse() {
	user_count--;
	if(user_count==0)
		delete this;
}

alphabet_mapper * alphabet_mapper::load(FILE *fp) {
  uint rd;
  if(fread(&rd,sizeof(uint),1,fp)!=1) return NULL;
  fseek(fp,-1*sizeof(uint),SEEK_CUR);
  switch(rd) {
    case ALPHABET_MAPPER_NONE_HDR: return alphabet_mapper_none::load(fp);
    case ALPHABET_MAPPER_CONT_HDR: return alphabet_mapper_cont::load(fp);
  }
  return NULL;
}
