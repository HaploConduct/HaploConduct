/*
 * Collection of basic tools and defines
 */

#include "Tools.h"
#include <cstdio>


time_t Tools::startTime;

void Tools::StartTimer()
{
    startTime = time(NULL);
}

double Tools::GetTime()
{
    time_t stopTime = time(NULL);
    return difftime( stopTime, startTime );
}

uchar * Tools::GetRandomString(unsigned min, unsigned max, unsigned &alphabetSize)
{
    unsigned len = std::rand() % (max - min) + min;
    alphabetSize = std::rand() % 26 + 1;
    uchar* temp = new uchar[len + 2];
    for (unsigned i = 0; i < len; i++)
        temp[i] = 97 + std::rand() % alphabetSize;
    temp[len] = 0u ;temp[len+1] = '\0';
    return temp;
}


void Tools::PrintBitSequence(ulong *A, ulong len)
{
    for(ulong i = 0; i < len; i++)
        if (GetField(A, 1, i))
            std::cout << "1";
        else
            std::cout << "0";
    std::cout << "\n";
}

unsigned Tools::FloorLog2(ulong i)
{
    unsigned b = 0;
    if (i == 0)
        return 0;
    while (i)
    { 
        b++; 
        i >>= 1; 
    }
    return b - 1;
}

unsigned Tools::CeilLog2(ulong i)
{
    unsigned j = FloorLog2(i);
    if ((ulong)(1lu << j) != i)
        return j + 1;
        
    return j;
}

uchar * Tools::GetFileContents(char *filename, ulong maxSize)
{
    std::ifstream::pos_type posSize;
    std::ifstream file ((char *)filename, std::ios::in|std::ios::binary|std::ios::ate);
    if (file.is_open())
    {
        posSize = file.tellg();
        ulong size = posSize;
        if (maxSize != 0 && size > maxSize)
            size = maxSize;
        char *memblock = new char [size + 1];
        file.seekg (0, std::ios::beg);
        file.read (memblock, size);
        memblock[size] = '\0';
        file.close();
	return (uchar *)memblock;
    }
    else
        return 0;
}

unsigned Tools::bits (ulong n)

   { unsigned b = 0;
     while (n)
    { b++; n >>= 1; }
     return b;
   }


