
#include <iostream>
#include <basics.h>
#include <static_bitsequence.h>

using namespace std;

int main(int argc, char ** argv) {
	if(argc!=2) {
		cout << "usage: " << argv[0] << " <bitmap_file>" << endl;
		return 0;
	}
	FILE * fp = fopen(argv[1],"r");
	if(fp==NULL) {
		cout << "Error opening " << argv[1] << endl;
		return -1;
	}
	uint *bitseq, len;
	uint l = fread(&len, sizeof(uint), 1, fp);
	bitseq = new uint[uint_len(len,1)];
	l += fread(bitseq, sizeof(uint), uint_len(len,1), fp);
	static_bitsequence * bs = new static_bitsequence_naive(bitseq,len);
	cout << "Bitmap length: " << len << " =? " << bs->length() << endl;
	uint ones = 0;
	for(uint i=0;i<len;i++) {
		if(bitget(bitseq,i)) ones++;
		if(bs->rank1(i) != ones) {
			cout << "Rank1 error for position " << i << endl;
			cout << " got: " << bs->rank1(i) << " expected: " << ones << endl;
		}
		if(bitget(bitseq,i) && bs->select1(ones) != i) {
			cout << "Select1 error for position " << i << endl;
			cout << " got: " << bs->select1(ones) << " expected: " << i << endl;
		}
		if(bs->rank0(i) != i+1-ones) {
			cout << "Rank0 error for position " << i << endl;
			cout << " got: " << bs->rank0(i) << " expected: " << ones << endl;
		}
		if(!bitget(bitseq,i) && bs->select0(i+1-ones) != i) {
			cout << "Select0 error for position " << i << endl;
			cout << " got: " << bs->select0(ones) << " expected: " << i << endl;
		}
	}
	delete bs;
	fclose(fp);
	cout << "Test completed." << endl;
	return 0;
}

