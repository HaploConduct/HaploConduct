
#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <basics.h>
#include <static_bitsequence.h>
#include <alphabet_mapper.h>
#include <static_sequence.h>

using namespace std;

void test_static_sequence(uint * symbols, uint n, static_sequence * ss) {
  cout << "Size: " << ss->size() << endl;
  uint max_v=0;
  for(uint i=0;i<n;i++) 
    max_v = max(max_v,symbols[i]);
  uint * occ = new uint[max_v+1];
  for(uint i=0;i<=max_v;i++)
    occ[i] = 0;
  bool error = false;
  for(uint i=0;i<n && !error;i++) {
    if(i!=0 && i%max(1,(n-1)/100)==0) { cout << "."; cout.flush(); }
    if(i!=0 && i%max(1,(n-1)/10)==0) cout << endl;
    occ[symbols[i]]++;
    uint a = ss->access(i);
    uint r = ss->rank(symbols[i],i);
    uint s = ss->select(symbols[i],occ[symbols[i]]);
    uint rM1 = (i==0)?0:ss->rank(symbols[i],i-1);
    if(r!=occ[symbols[i]]) {
      cout << "Error in rank for symbol " << symbols[i] << " at position " << i << endl;
      cout << "value: " << r << endl;
      cout << "Expected: " << occ[symbols[i]] << endl;
      error = true;
    }
    if(s!=i) {
      cout << "Error in select for symbol " << symbols[i] << " at position " << occ[symbols[i]] << endl;
      cout << "value: " << s << endl;
      cout << "Expected: " << i << endl;
      error = true;
    }
    if(a!=symbols[i]) {
      cout << "Error in access at position " << i << endl;
      cout << "value: " << a << endl;
      cout << "Expected: " << symbols[i] << endl;
      error = true;
    }
    if(rM1!=occ[symbols[i]]-1) {
      cout << "Error in rankM1 for symbol " << symbols[i] << " at position " << i-1 << endl;
      cout << "value: " << rM1 << endl;
      cout << "Expected: " << occ[symbols[i]]-1 << endl;
      error = true;
    }
  }
  if(!error)
    cout << "Test OK! It works :)" << endl;
  delete [] occ;  
}

int main(int argc, char ** argv) {
  if(argc!=3) {
    cout << "usage: " << argv[0] << " <file> <samp>" << endl;
    return 0;
  }
  struct stat text_info;
  if(stat(argv[1],&text_info)<0) {
    cout << "could not stat: " << argv[1] << endl;
    return -1;
  }
  
  stringstream ss;
  ss << argv[2];
  uint samp;
  ss >> samp;

  uint n= (uint)text_info.st_size/4;
  uint * text = new uint[n+1];
  FILE * fp = fopen(argv[1],"r");
  if(fp==NULL) {
    cout << "could not open " << argv[1] << endl;
    return -1;
  }

  cout << "File: " << argv[1] << endl;
  cout << "Length: " << n << endl;

  uint max_symbol = 0;
  for(uint i=0;i<n;i++) {
    uint c;
    uint read=fread(&c,sizeof(uint),1,fp);
    //assert(read==1);
    text[i] = (uint)c;
    c += read;
    max_symbol = max(max_symbol,text[i]);
  }
  max_symbol++;

  fclose(fp);

  /*uint *occ = new uint[max_symbol];
  for(uint i=0;i<max_symbol;i++)
    occ[i] = 0;
  for(uint i=0;i<n;i++)
    occ[text[i]]++;*/

  alphabet_mapper * am = new alphabet_mapper_none();
  static_bitsequence_builder * bmb = new static_bitsequence_builder_rrr02(samp);
	cout << "Building Huffman table..."; cout.flush();
	wt_coder * wtc = new wt_coder_huff(text,n,am);
	cout << "done" << endl; cout.flush();
  cout << "Building static_sequence..."; cout.flush();
  static_sequence * wt = new static_sequence_wvtree(text,n,wtc,bmb,am);
  cout << "done" << endl; cout.flush();
  delete bmb;
  
  char * fname = new char[10+string(argv[1]).length()];
  sprintf(fname,"%s.wt",argv[1]);
  fp = fopen(fname,"w");
  wt->save(fp);
  fclose(fp);
  delete wt;
  fp = fopen(fname,"r");
  wt = static_sequence::load(fp);
  fclose(fp);
  delete [] fname;
  
  test_static_sequence(text,n,wt);
  
  cout << "WT Size: " << wt->size() << endl;
  cout << "ft = " << 1.*wt->size()/(bits(max_symbol-1)*n/8) << endl;

  delete [] text;
  delete wt;

}
