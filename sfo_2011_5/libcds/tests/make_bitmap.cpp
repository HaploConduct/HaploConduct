/* make_bitmap.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * make_bitmap
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

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <basics.h>
#include <cmath>

using namespace std;

int main(int argc, char ** argv) {
	if(argc!=4) {
		cout << "usage: " << argv[0] << " <bitmap> <length> <ones>" << endl;
		return 0;
	}
	char * fname = argv[1];
	uint len = atoi(argv[2]);
	uint ones = atoi(argv[3]);
	uint * bm = new uint[uint_len(len,1)];
	for(uint i=0;i<uint_len(len,1);i++)
		bm[i] = 0;
	for(uint i=0;i<ones;i++) {
		while(true) {
			uint p = rand()%len;
			if(!bitget(bm,p)) {
				bitset(bm,p);
				break;
			}
		}
	}
	FILE * fp = fopen(fname,"w");
	uint l = fwrite(&len,sizeof(uint),1,fp);
	l += fwrite(bm,sizeof(uint),uint_len(len,1),fp);
	fclose(fp);
	delete [] bm;
	return 0;
}

