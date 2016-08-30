#ifndef _BLOCK_ARRAY_H_
#define _BLOCK_ARRAY_H_
#include "Tools.h"
#include <iostream>
#include <stdexcept>

class BlockArray
{
private:
    ulong* data;
    ulong n;
    ulong index;
    ulong blockLength;
public:
    BlockArray(ulong len, ulong blockLen) {
       n = len;
       blockLength = blockLen;
       data = new ulong[n*blockLength/W +1];
       for (ulong i = 0; i < n*blockLength/W +1; ++i)
           data[i] = 0;
    }
    ~BlockArray() {
       delete [] data;
    }

    BlockArray& operator[](ulong i)  {
       index = i;
       return *this;
    }  

    void operator=(const ulong x) {
       Tools::SetField(data,blockLength,index,x);
    }  

    BlockArray& operator=(const BlockArray& ba) {
       if (this == &ba) return *this;
       ulong value = Tools::GetField(ba.data, ba.blockLength, ba.index);
       Tools::SetField(data,blockLength,index,value);
       return *this;
    }  
        
    operator ulong() {
       return Tools::GetField(data,blockLength,index);
    }
    
    ulong spaceInBits() {
        return n*blockLength+W; // plus 4 ulong's
    }

    /**
     * Saving data fields:
     *     ulong n;
     *     ulong blockLength;
     *     ulong* data;
     */
    void Save(FILE *file) const
    {
        if (std::fwrite(&(this->n), sizeof(ulong), 1, file) != 1)
            throw std::runtime_error("BlockArray::Save(): file write error (n).");
        if (std::fwrite(&(this->blockLength), sizeof(ulong), 1, file) != 1)
            throw std::runtime_error("BlockArray::Save(): file write error (blockLength).");
    
        if (std::fwrite(this->data, sizeof(ulong), n*blockLength/W+1, file) != n*blockLength/W+1)
            throw std::runtime_error("BlockArray::Save(): file write error (data).");
    }

    /**
     * Load from file
     */
    BlockArray(FILE *file)
    {
        if (std::fread(&(this->n), sizeof(ulong), 1, file) != 1)
            throw std::runtime_error("BlockArray::Load(): file read error (n).");
        if (std::fread(&(this->blockLength), sizeof(ulong), 1, file) != 1)
            throw std::runtime_error("BlockArray::Load(): file read error (blockLength).");

        data = new ulong[n*blockLength/W+1];
        if (std::fread(this->data, sizeof(ulong), n*blockLength/W+1, file) != n*blockLength/W+1)
            throw std::runtime_error("BlockArray::Load(): file read error (data).");
    }
};

#endif

