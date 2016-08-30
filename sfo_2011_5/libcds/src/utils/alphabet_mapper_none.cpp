/* alphabet_mapper_none.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * alphabet_mapper definition
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
 
#include <alphabet_mapper_none.h>

alphabet_mapper_none::alphabet_mapper_none() { }

uint alphabet_mapper_none::map(uint s) {return s;}

uint alphabet_mapper_none::unmap(uint s) {return s;}

uint alphabet_mapper_none::size() { return sizeof(alphabet_mapper_none); }

uint alphabet_mapper_none::save(FILE *fp) {
  uint wr = ALPHABET_MAPPER_NONE_HDR;
  wr = fwrite(&wr,sizeof(uint),1,fp);
  if(wr!=1) return 1;
  return 0;
}

alphabet_mapper_none * alphabet_mapper_none::load(FILE * fp) {
  uint rd;
  if(fread(&rd,sizeof(uint),1,fp)!=1) return NULL;
  if(rd!=ALPHABET_MAPPER_NONE_HDR) return NULL;
  return new alphabet_mapper_none();
}
