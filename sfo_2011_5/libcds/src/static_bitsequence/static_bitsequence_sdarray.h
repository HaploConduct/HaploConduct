
#ifndef _STATIC_BITSEQUENCE_SDARRAY_H
#define _STATIC_BITSEQUENCE_SDARRAY_H

#include <basics.h>
#include <static_bitsequence.h>
#include <sdarray.h>

class static_bitsequence_sdarray: public static_bitsequence {
	public:
		static_bitsequence_sdarray(uint * buff, uint len);
		virtual ~static_bitsequence_sdarray();
		virtual uint select1(uint i);
		virtual uint rank1(uint i);
		virtual uint select_next1(uint i);
		virtual uint size();
		virtual int save(FILE * fp);
		static static_bitsequence_sdarray * load(FILE * fp);

	protected:
		selects3 sd;
		static_bitsequence_sdarray();

};

#endif

