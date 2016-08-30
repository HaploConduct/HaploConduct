
#include <static_bitsequence_sdarray.h>

static_bitsequence_sdarray::static_bitsequence_sdarray(uint * buff, uint len) {
  uint * tmp_seq = new uint[uint_len(len,1)+1];
  ones = 0;
  for(uint i=0;i<uint_len(len,1)+1;i++)
    tmp_seq[i] = 0;
  for(uint i=0;i<len;i++)
  if(bitget(buff,i)) {
    __setbit(tmp_seq,i,1);
    ones++;
  }
  if(ones)
    selects3_construct(&sd,len,tmp_seq);
  this->len = len;
	//sd.lasti=(uint)-3;
  //this->ones = sd.m;
  delete [] tmp_seq;
}


static_bitsequence_sdarray::static_bitsequence_sdarray() {make___selecttbl();}

static_bitsequence_sdarray::~static_bitsequence_sdarray() {
  if(ones)
    selects3_free(&sd);
}


uint static_bitsequence_sdarray::rank1(uint i) {
  if(i>=len) return -1;
  if(ones)
    return selects3_rank(&sd,i);
  else
    return 0;
}


uint static_bitsequence_sdarray::select1(uint i) {
  if(i>ones || i==0) return -1;
  if(ones)
    return selects3_select(&sd,i);
  else
    return (uint)-1;
}


uint static_bitsequence_sdarray::select_next1(uint i) {
  return selects3_selectnext(&sd,i);
}


uint static_bitsequence_sdarray::size() {
  return sizeof(static_bitsequence_sdarray)+(ones?(sd.size + sd.sd0->size + sd.sd1->size):0);
}


int static_bitsequence_sdarray::save(FILE * fp) {
  uint wr = SDARRAY_HDR;
  wr = fwrite(&wr,sizeof(uint),1,fp);
  wr += fwrite(&len,sizeof(uint),1,fp);
  wr += fwrite(&ones,sizeof(uint),1,fp);
  if(wr!=3 || (ones?(selects3_save(&sd,fp)):false))
    return 1;
  return 0;
}


static_bitsequence_sdarray * static_bitsequence_sdarray::load(FILE * fp) {
  uint id;
  if(fread(&id,sizeof(uint),1,fp)!=1) return NULL;
  if(id!=SDARRAY_HDR) return NULL;
  static_bitsequence_sdarray * ret = new static_bitsequence_sdarray();
  id = fread(&ret->len,sizeof(uint),1,fp);
  id += fread(&ret->ones,sizeof(uint),1,fp);
  if(ret->ones && selects3_load(&ret->sd,fp)) {
    delete ret;
    return NULL;
  }
  return ret;
}
