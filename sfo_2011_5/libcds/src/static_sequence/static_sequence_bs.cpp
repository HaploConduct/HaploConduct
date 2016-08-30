
#include <static_sequence_bs.h>

static_sequence_bs::static_sequence_bs(uint * seq, uint n, alphabet_mapper * am, static_bitsequence_builder * bmb) {
	sigma = 0;
	len = n;
	this->am = am;
	am->use();
	for(uint i=0;i<n;i++) sigma=max(sigma,am->map(seq[i]));
	sigma++;
	uint * occ = new uint[sigma+1];
	for(uint i=0;i<=sigma;i++) occ[i] = 0;
	for(uint i=0;i<n;i++) occ[am->map(seq[i])+1]++;
	for(uint i=1;i<sigma;i++) occ[i] += occ[i-1];
	uint * pos = new uint[n];
	for(uint i=0;i<n;i++) pos[i] = 0;
	for(uint i=0;i<n;i++) pos[occ[am->map(seq[i])]++]=i;
	bitmaps = new static_bitsequence*[sigma];
	uint * bm = new uint[uint_len(n,1)];
	uint pp=0;
	for(uint i=0;i<sigma;i++) {
		for(uint j=0;j<uint_len(n,1);j++) bm[j]=0;
		while(pp<occ[i]) {
			bitset(bm,pos[pp]);
			pp++;
		}
		bitmaps[i] = bmb->build(bm,len);
	}
	delete [] bm;
	delete [] occ;
	delete [] pos;
}

static_sequence_bs::static_sequence_bs() {
	len = 0;
	sigma = 0;
	bitmaps = NULL;
	am = NULL;
}

static_sequence_bs::~static_sequence_bs() {
	if(bitmaps!=NULL) {
		for(uint i=0;i<sigma;i++) {
			if(bitmaps[i]!=NULL) delete bitmaps[i];
		}
		delete [] bitmaps;
	}
	if(am!=NULL) am->unuse();
}

uint static_sequence_bs::rank(uint c, uint i) {
	if(am->map(c)>=sigma) return (uint)-1;
	return bitmaps[am->map(c)]->rank1(i);
}

uint static_sequence_bs::select(uint c, uint i) {
	if(am->map(c)>=sigma) return (uint)-1;
	return bitmaps[am->map(c)]->select1(i);
}

uint static_sequence_bs::select_next(uint c, uint i) {
	if(am->map(c)>=sigma) return (uint)-1;
	return bitmaps[am->map(c)]->select_next1(i);
}

uint static_sequence_bs::access(uint i) {
	for(uint j=0;j<sigma;j++) {
		if(bitmaps[j]->access(i)) return am->unmap(j);
	}
	return (uint)-1;
}

uint static_sequence_bs::size() {
	uint size = sizeof(static_sequence_bs)+am->size();
	for(uint i=0;i<sigma;i++) 
		size += bitmaps[i]->size();
	return size;
}

uint static_sequence_bs::save(FILE * fp) {
	uint wr = BS_HDR;
	wr = fwrite(&wr,sizeof(uint),1,fp);
	wr += fwrite(&len,sizeof(uint),1,fp);
	wr += fwrite(&sigma,sizeof(uint),1,fp);
	if(wr!=3) return 1;
	for(uint i=0;i<sigma;i++)
		if(bitmaps[i]->save(fp)) return 2;
	if(am->save(fp)) return 3;
	return 0;
}

static_sequence_bs * static_sequence_bs::load(FILE * fp) {
	uint rd = 0;
	uint type = 0;
	rd += fread(&type,sizeof(uint),1,fp);
	static_sequence_bs * ret = new static_sequence_bs();
	rd += fread(&ret->len,sizeof(uint),1,fp);
	rd += fread(&ret->sigma,sizeof(uint),1,fp);
	if(rd!=3 || type != BS_HDR) {
		delete ret;
		return NULL;
	}
	ret->bitmaps = new static_bitsequence*[ret->sigma];
	for(uint i=0;i<ret->sigma;i++)
		ret->bitmaps[i] = NULL;
	for(uint i=0;i<ret->sigma;i++)
		if((ret->bitmaps[i]=static_bitsequence::load(fp))==NULL) {
			delete ret;
			return NULL;
		}
	if((ret->am = alphabet_mapper::load(fp))==NULL) {
		delete ret;
		return NULL;
	}
	ret->am->use();
	return ret;
}

