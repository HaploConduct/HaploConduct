/* test_to_int.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * text_to_int
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

#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <basics.h>
#include <static_bitsequence.h>
#include <static_bitsequence_builder.h>
#include <alphabet_mapper.h>
#include <static_sequence.h>
#include <static_sequence_builder.h>
#include "static_sequence_tester.h"

int main(int argc, char ** argv) {
  if(argc!=3) {
    cout << "Usage: " << argv[0] << " <file> <output>" << endl;
    return 0;
  }
  char * fname = argv[1];
  char * oname = argv[2];

  FILE * fp = fopen(fname,"r");
  if(fp==NULL) {
    cout << "could not open " << fname << endl;
    return 1;
  }
  struct stat text_info;
  if(stat(fname,&text_info)<0) {
    cout << "could not stat: " << fname << endl;
    return 1;
  }
  
  uint n= (uint)text_info.st_size;
  uint * text = new uint[n];

  for(uint i=0;i<n;i++) {
    uchar c;
    text[i] = fread(&c,sizeof(uchar),1,fp);
    text[i] = (uint)c;
  }
	fclose(fp);
  
  FILE * out = fopen(oname,"w");
  if(out==NULL) {
    cout << "could not open " << oname << endl;
    return 1;
  }
  
  n = fwrite(text,sizeof(uint),n,out);
  fclose(out);
  
  delete [] text;
	return 0;
}

