
#include <iostream>
#include <sstream>
#include <basics.h>
#include <static_bitsequence.h>

using namespace std;

int main(int argc, char ** argv) {
	if(argc!=3) {
		cout << "usage: " << argv[0] << " <bitmap_file> <sample_rate>" << endl;
		return 0;
	}
	FILE * fp = fopen(argv[1],"r");
	if(fp==NULL) {
		cout << "Error opening " << argv[1] << endl;
		return -1;
	}
	uint *bitseq, len;
	uint l=fread(&len, sizeof(uint), 1, fp);
	bitseq = new uint[uint_len(len,1)];
	l+=fread(bitseq, sizeof(uint), uint_len(len,1), fp);
	fclose(fp);

	static_bitsequence * bs = new static_bitsequence_rrr02(bitseq,len);
	
	char * fname = new char[string(argv[1]).length()+10];
	sprintf(fname,"%s.rrr",argv[1]);
	
	fp = fopen(fname,"w");
	cout << "Save: " << bs->save(fp) << endl;
	fclose(fp);
	delete bs;

	fp = fopen(fname,"r");
	bs = static_bitsequence::load(fp);
	uint sample_rate;
	stringstream ss(argv[2]);
	ss >> sample_rate;
	((static_bitsequence_rrr02*)bs)->create_sampling(sample_rate);
	fclose(fp);
	delete [] fname;

	cout << "Bitmap length: " << len << " =? " << bs->length() << endl;
  cout << "Ones: " << bs->count_one() << endl;
	cout << "Bitmap size: " << bs->size() << endl;
	/*for(uint i=0;i<64;i++) {
		if(i%15==0) cout << " ";
		cout << (bs->access(i)?"1":"0");
	}
	cout << endl;*/
	uint ones = 0;
	for(uint i=0;i<len;i++) {
		//cout << "i="<<i<< endl;
		if(i%max(1,(bs->length()/100))==0) cout << i/max(1,(bs->length()/100)) << "%" << endl;
		if(bitget(bitseq,i)) ones++;
		if(bs->access(i) != (bitget(bitseq,i)!=0)) {
			cout << "Access error for position " << i << endl;
			cout << " got: " << bs->access(i) << " expected: " << (bitget(bitseq,i)!=0) << endl;
		}
		if(bs->rank1(i) != ones) {
			cout << "Rank1 error for position " << i << endl;
			cout << " got: " << bs->rank1(i) << " expected: " << ones << endl;
		} 
		if(bitget(bitseq,i) && bs->select1(ones) != i) {
			cout << "Select1 error for position " << i << " ones:" << ones << endl;
			cout << " got: " << bs->select1(ones) << " expected: " << i << endl;
		}
		if(bs->rank0(i) != i+1-ones) {
			cout << "Rank0 error for position " << i << endl;
			cout << " got: " << bs->rank0(i) << " expected: " << ones << endl;
		}
		if(!bitget(bitseq,i) && bs->select0(i+1-ones) != i) {
			cout << "Select0 error for position " << i << endl;
			cout << " got: " << bs->select0(i+1-ones) << " expected: " << i << endl;
		}
	}
	delete bs;
  //static_bitsequence_rrr02::delete_E();
  delete [] bitseq;
	cout << "Test completed." << endl;
	return 0;
}

