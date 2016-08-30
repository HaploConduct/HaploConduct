/******************************************************************************
 *   Copyright (C) 2006-2008 by Veli Mäkinen and Niko Välimäki                *
 *                                                                            *
 *                                                                            *
 *   This program is free software; you can redistribute it and/or modify     *
 *   it under the terms of the GNU Lesser General Public License as published *
 *   by the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                      *
 *                                                                            *
 *   This program is distributed in the hope that it will be useful,          *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *   GNU Lesser General Public License for more details.                      *
 *                                                                            *
 *   You should have received a copy of the GNU Lesser General Public License *
 *   along with this program; if not, write to the                            *
 *   Free Software Foundation, Inc.,                                          *
 *   51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.            *
 *****************************************************************************/
#include "TCImplementation.h"

//#define DEBUG_MEMUSAGE
#ifdef DEBUG_MEMUSAGE
#include "HeapProfiler.h" // FIXME remove
#endif

#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cassert>
#include <cstring> // For strlen()
using std::vector;
using std::pair;
using std::make_pair;
using std::map;

namespace SXSI
{

// Save file version info
const uchar TCImplementation::versionFlag = 7;

/**
 * Constructor inits an empty dynamic FM-index.
 * Samplerate defaults to TEXTCOLLECTION_DEFAULT_SAMPLERATE.
 */
TCImplementation::TCImplementation(uchar * bwt, ulong length, unsigned samplerate_, 
                                   unsigned numberOfTexts_, ulong maxTextLength_, ulong numberOfSamples_, char tsType, bool colorcodes)
    : n(length), samplerate(samplerate_), alphabetrank(0), color(colorcodes), sampled(0), suffixes(0), 
      suffixDocId(0), numberOfTexts(numberOfTexts_), maxTextLength(maxTextLength_), Doc(0)
{
    this->setColorCoded(color);

    makewavelet(bwt); // Deletes bwt!
    bwt = 0;
 
    // Make sampling tables
    maketables(numberOfSamples_, tsType);
}

bool TCImplementation::EmptyText(DocId k) const
{
    assert(k < (DocId)numberOfTexts); 
    return false; // Empty texts are not indexed
}

uchar * TCImplementation::GetText(DocId k) const
{
    assert(k < (DocId)numberOfTexts);

    return textStorage->GetText(k);
}

/**
 * Save index to a file handle
 *
 * Throws a std::runtime_error exception on i/o error.
 * First byte that is saved represents the version number of the save file.
 * In version 2 files, the data fields are:
 *  uchar versionFlag;
    TextPosition n;
    unsigned samplerate;
    unsigned C[256];
    TextPosition bwtEndPos;
    static_sequence * alphabetrank;
    BSGAP *sampled; 
    BlockArray *suffixes;
    BlockArray *suffixDocId;
    unsigned numberOfTexts;
    ulong maxTextLength;
    static_sequence *docId;
 */
void TCImplementation::Save(FILE *file) const
{
    // Saving version info:
    if (std::fwrite(&versionFlag, 1, 1, file) != 1)
        throw std::runtime_error("TCImplementation::Save(): file write error (version flag).");

    if (std::fwrite(&(this->n), sizeof(TextPosition), 1, file) != 1)
        throw std::runtime_error("TCImplementation::Save(): file write error (n).");
    if (std::fwrite(&(this->samplerate), sizeof(unsigned), 1, file) != 1)
        throw std::runtime_error("TCImplementation::Save(): file write error (samplerate).");

    for(ulong i = 0; i < 256; ++i)
        if (std::fwrite(this->C + i, sizeof(unsigned), 1, file) != 1)
            throw std::runtime_error("TCImplementation::Save(): file write error (C table).");

    if (std::fwrite(&(this->bwtEndPos), sizeof(TextPosition), 1, file) != 1)
        throw std::runtime_error("TCImplementation::Save(): file write error (bwt end position).");
    
    alphabetrank->save(file);
    sampled->save(file);
    suffixes->Save(file);
    suffixDocId->Save(file);
    
     if (std::fwrite(&(this->color), sizeof(bool), 1, file) != 1)
        throw std::runtime_error("TCImplementation::Save(): file write error (color).");
    if (std::fwrite(&(this->numberOfTexts), sizeof(unsigned), 1, file) != 1)
        throw std::runtime_error("TCImplementation::Save(): file write error (numberOfTexts).");
    if (std::fwrite(&(this->maxTextLength), sizeof(ulong), 1, file) != 1)
        throw std::runtime_error("TCImplementation::Save(): file write error (maxTextLength).");

    Doc->save(file);
    textStorage->Save(file);
    fflush(file);
}


/**
 * Load index from a file handle
 *
 * Throws a std::runtime_error exception on i/o error.
 * For more info, see TCImplementation::Save().
 */
TCImplementation::TCImplementation(FILE *file, unsigned samplerate_)
    : n(0), samplerate(samplerate_), alphabetrank(0), color(0), sampled(0), suffixes(0), 
      suffixDocId(0), numberOfTexts(0), maxTextLength(0), Doc(0)
{
    uchar verFlag = 0;
    if (std::fread(&verFlag, 1, 1, file) != 1)
        throw std::runtime_error("TCImplementation::Load(): file read error (version flag).");
    if (verFlag != TCImplementation::versionFlag)
        throw std::runtime_error("TCImplementation::Load(): invalid save file version.");

    if (std::fread(&(this->n), sizeof(TextPosition), 1, file) != 1)
        throw std::runtime_error("TCImplementation::Load(): file read error (n).");
    if (std::fread(&samplerate, sizeof(unsigned), 1, file) != 1)
        throw std::runtime_error("TCImplementation::Load(): file read error (samplerate).");
// FIXME samplerate can not be changed during load.
//    if (this->samplerate == 0)
//        this->samplerate = samplerate;

    for(ulong i = 0; i < 256; ++i)
        if (std::fread(this->C + i, sizeof(unsigned), 1, file) != 1)
            throw std::runtime_error("TCImplementation::Load(): file read error (C table).");

    if (std::fread(&(this->bwtEndPos), sizeof(TextPosition), 1, file) != 1)
        throw std::runtime_error("TCImplementation::Load(): file read error (bwt end position).");

    //std::cout << "Loading alphabet rank (" << Tools::GetTime() << " s)." << std::endl;
    alphabetrank = static_sequence::load(file);
    //std::cout << "Loading samples (" << Tools::GetTime() << " s)." << std::endl;
    sampled = static_bitsequence::load(file);
    suffixes = new BlockArray(file);
    suffixDocId = new BlockArray(file);

    if (std::fread(&(this->color), sizeof(bool), 1, file) != 1)
        throw std::runtime_error("TCImplementation::Load(): file read error (color).");
    this->setColorCoded(color);

    if (std::fread(&(this->numberOfTexts), sizeof(unsigned), 1, file) != 1)
        throw std::runtime_error("TCImplementation::Load(): file read error (numberOfTexts).");
    if (std::fread(&(this->maxTextLength), sizeof(ulong), 1, file) != 1)
        throw std::runtime_error("TCImplementation::Load(): file read error (maxTextLength).");

    Doc = new ArrayDoc(file);

    textStorage = TextStorage::Load(file);


    // FIXME Construct data structures with new samplerate
    //maketables(); 
}


TCImplementation::~TCImplementation() {
    delete alphabetrank;       
    delete sampled;
    delete suffixes;
    delete suffixDocId;
    delete Doc;
    delete textStorage;
}

void TCImplementation::makewavelet(uchar *bwt)
{
    ulong i, min = 0,
             max;
    for (i=0;i<256;i++)
        C[i]=0;
    for (i=0;i<n;++i)
        C[(int)bwt[i]]++;
    for (i=0;i<256;i++)
        if (C[i]>0) {min = i; break;}          
    for (i=255;i>=min;--i)
        if (C[i]>0) {max = i; break;}
    
    ulong prev=C[0], temp;
    C[0]=0;
    for (i=1;i<256;i++) {          
        temp = C[i];
        C[i]=C[i-1]+prev;
        prev = temp;
    }

#ifdef DEBUG_MEMUSAGE
    std::cerr << "max heap usage before WT: " << HeapProfiler::GetMaxHeapConsumption()/(1024*1024) << " Mbytes" << std::endl;
    HeapProfiler::ResetMaxHeapConsumption(); 
#endif

    static_bitsequence_builder * bmb = new static_bitsequence_builder_brw32(8);
    alphabet_mapper * am = new alphabet_mapper_cont(bwt, n, bmb); 
    delete bmb;
    bmb = new static_bitsequence_builder_brw32(8);
    wt_coder * wtc = new wt_coder_binary(bwt,n,am); // Note: wt_coder_huff is not thread-safe
    alphabetrank = new static_sequence_wvtree(bwt,n,wtc,bmb,am);
    delete bmb;
    bwt = 0; // already deleted
   
#ifdef DEBUG_MEMUSAGE
    std::cerr << "heap usage after WT: " << HeapProfiler::GetHeapConsumption()/(1024*1024) << " Mbytes" << std::endl;
    std::cerr << "max heap usage after WT: " << HeapProfiler::GetMaxHeapConsumption()/(1024*1024) << " Mbytes" << std::endl;
#endif
}

void TCImplementation::maketables(ulong sampleLength, char tsType)
{
    // Calculate BWT end-marker position (of last inserted text)
    {
        ulong i = 0;
        uint alphabetrank_i_tmp = 0;
        uchar c  = alphabetrank->access(i, alphabetrank_i_tmp);
        while (c != '\0')
        {
            i = C[c]+alphabetrank_i_tmp-1;
            c = alphabetrank->access(i, alphabetrank_i_tmp);
        }

        this->bwtEndPos = i;
    }

#ifdef DEBUG_MEMUSAGE
    std::cerr << "heap usage before BWT traverse: " << HeapProfiler::GetHeapConsumption()/(1024*1024) << " / " << HeapProfiler::GetMaxHeapConsumption()/(1024*1024) << " Mbytes, " <<  HeapProfiler::GetHeapConsumption() << " / " <<  HeapProfiler::GetMaxHeapConsumption() << std::endl;
    HeapProfiler::ResetMaxHeapConsumption();
#endif

    // Build up array for text starting positions
//    BlockArray* textStartPos = new BlockArray(numberOfTexts, Tools::CeilLog2(this->n));
//    (*textStartPos)[0] = 0; 

    // Mapping from end-markers to doc ID's:
    unsigned logNumberOfTexts = Tools::CeilLog2(numberOfTexts);
//    uint *endmarkerDocId = new uint[(numberOfTexts * logNumberOfTexts)/(8*sizeof(uint)) + 1];
    BlockArray *endmarkerDocId = new BlockArray(numberOfTexts, logNumberOfTexts);

    BlockArray* positions = new BlockArray(sampleLength, Tools::CeilLog2(this->n));
    uint *sampledpositions = new uint[n/(sizeof(uint)*8)+1];
    for (ulong i = 0; i < n / (sizeof(uint)*8) + 1; i++)
        sampledpositions[i] = 0;
    
    ulong x,p=bwtEndPos;
    ulong sampleCount = 0;
    // Keeping track of text position of prev. end-marker seen
    ulong posOfSuccEndmarker = n-1;
    DocId textId = numberOfTexts;
    ulong ulongmax = 0;
    ulongmax--;
    uint alphabetrank_i_tmp =0;

    TextStorageBuilder tsbuilder(n);
    Tools::StartTimer();

    for (ulong i=n-1;i<ulongmax;i--) {
        // i substitutes SA->GetPos(i)
        x=(i==n-1)?0:i+1;

        uchar c = alphabetrank->access(p, alphabetrank_i_tmp);
        tsbuilder[i] = c;

        if ((posOfSuccEndmarker - i) % samplerate == 0 && c != '\0')
        {
            set_field(sampledpositions,1,p,1);
            (*positions)[sampleCount] = p;
            sampleCount ++;
        }

        if (c == '\0')
        {
            --textId;
            
            // Record the order of end-markers in BWT:
            ulong endmarkerRank = alphabetrank_i_tmp - 1;
            (*endmarkerDocId)[endmarkerRank] = (textId + 1) % numberOfTexts;
            

            // Store text length and text start position:
            if (textId < (DocId)numberOfTexts - 1)
            {
//                (*textStartPos)[textId + 1] = x;  // x-1 is text position of end-marker.
                posOfSuccEndmarker = i;
            }

            // LF-mapping from '\0' does not work with this (pseudo) BWT (see details from Wolfgang's thesis).
            p = textId; // Correct LF-mapping to the last char of the previous text.
        }
        else // Now c != '\0', do LF-mapping:
            p = C[c]+alphabetrank_i_tmp-1;
    }
    assert(textId == 0);

#ifdef DEBUG_MEMUSAGE
    std::cerr << "heap usage before tsbuilder init: " << HeapProfiler::GetHeapConsumption()/(1024*1024) << " / " << HeapProfiler::GetMaxHeapConsumption()/(1024*1024) << " Mbytes, " <<  HeapProfiler::GetHeapConsumption() << " / " <<  HeapProfiler::GetMaxHeapConsumption() << std::endl;
    HeapProfiler::ResetMaxHeapConsumption();
#endif

    textStorage = tsbuilder.InitTextStorage(tsType);

#ifdef DEBUG_MEMUSAGE
    std::cerr << "heap usage after tsbuilder init: " << HeapProfiler::GetHeapConsumption()/(1024*1024) << " / " << HeapProfiler::GetMaxHeapConsumption()/(1024*1024) << " Mbytes, " <<  HeapProfiler::GetHeapConsumption() << " / " <<  HeapProfiler::GetMaxHeapConsumption() << std::endl;
    HeapProfiler::ResetMaxHeapConsumption();
#endif

    sampled = new static_bitsequence_brw32(sampledpositions, n, 16);
    delete [] sampledpositions;
    assert(sampleCount == sampleLength);
    assert(sampled->rank1(n-1) == sampleLength);

#ifdef DEBUG_MEMUSAGE
    std::cerr << "heap usage after sampled bit vector: " << HeapProfiler::GetHeapConsumption()/(1024*1024) << " / " << HeapProfiler::GetMaxHeapConsumption()/(1024*1024) << " Mbytes, " <<  HeapProfiler::GetHeapConsumption() << " / " <<  HeapProfiler::GetMaxHeapConsumption() << std::endl;
    HeapProfiler::ResetMaxHeapConsumption();
#endif

    // Suffixes store an offset from the text start position
    suffixes = new BlockArray(sampleLength, Tools::CeilLog2(maxTextLength));
    suffixDocId = new BlockArray(sampleLength, Tools::CeilLog2(numberOfTexts));

    x = n - 2;
    posOfSuccEndmarker = n-1;
    for(ulong i=0; i<sampleLength; i++) {
        // Find next sampled text position
        while ((posOfSuccEndmarker - x) % samplerate != 0)
        {
            --x;
            assert(x != ~0lu);
            if (textStorage->IsEndmarker(x))
                posOfSuccEndmarker = x--;
        }
        assert((*positions)[i] < n);
        ulong j = sampled->rank1((*positions)[i]);

        assert(j != 0); // if (j==0) j=sampleLength;
        
        TextPosition textPos = (x==n-1)?0:x+1;
        (*suffixDocId)[j-1] = textStorage->DocIdAtTextPos(textPos);

        assert((*suffixDocId)[j-1] < numberOfTexts);
        // calculate offset from text start:
        (*suffixes)[j-1] = textPos - textStorage->TextStartPos((*suffixDocId)[j-1]);
        --x;
        if (x != ~0lu && textStorage->IsEndmarker(x))
            posOfSuccEndmarker = x--;
    }

    delete positions;

#ifdef DEBUG_MEMUSAGE
    std::cerr << "heap usage after sampled arrays: " << HeapProfiler::GetHeapConsumption()/(1024*1024) << " / " << HeapProfiler::GetMaxHeapConsumption()/(1024*1024) << " Mbytes, " <<  HeapProfiler::GetHeapConsumption() << " / " <<  HeapProfiler::GetMaxHeapConsumption() << std::endl;
    HeapProfiler::ResetMaxHeapConsumption();
#endif

#ifdef DEBUG_MEMUSAGE
    std::cerr << "max heap usage before Doc: " << HeapProfiler::GetMaxHeapConsumption()/(1024*1024) << " Mbytes" << std::endl;
    HeapProfiler::ResetMaxHeapConsumption();
#endif

    Doc = new ArrayDoc(endmarkerDocId);

#ifdef DEBUG_MEMUSAGE
    std::cerr << "max heap usage after Doc: " << HeapProfiler::GetMaxHeapConsumption()/(1024*1024) << " Mbytes" << std::endl;
#endif
}


/**
 * Finds document identifier for given text position
 *
 * Starting text position of the document is stored into second parameter.
 * Binary searching on text starting positions. 
 */
TextCollection::DocId TCImplementation::DocIdAtTextPos(BlockArray* textStartPos, TextPosition i) const
{
    assert(i < n);

    DocId a = 0;
    DocId b = numberOfTexts - 1;
    while (a < b)
    {
        DocId c = a + (b - a)/2;
        if ((*textStartPos)[c] > i)
            b = c - 1;
        else if ((*textStartPos)[c+1] > i)
            return c;
        else
            a = c + 1;
    }

    assert(a < (DocId)numberOfTexts);
    assert(i >= (*textStartPos)[a]);
    assert(i < (a == (DocId)numberOfTexts - 1 ? n : (*textStartPos)[a+1]));
    return a;
}


} // namespace SXSI

