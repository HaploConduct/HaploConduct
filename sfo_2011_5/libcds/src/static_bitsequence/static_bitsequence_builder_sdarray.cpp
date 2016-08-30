
#include <static_bitsequence_builder_sdarray.h>

static_bitsequence * static_bitsequence_builder_sdarray::build(uint * buff, uint len) {
	return new static_bitsequence_sdarray(buff,len);
}

