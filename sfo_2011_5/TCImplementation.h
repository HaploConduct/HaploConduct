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

#ifndef _TCImplementation_H_
#define _TCImplementation_H_
#include "BitRank.h"
#include "TextCollection.h"
#include "BlockArray.h"
#include "ArrayDoc.h"

// Include  from XMLTree/libcds
#include <basics.h> // Defines W == 32
#include <static_bitsequence.h>
#include <alphabet_mapper.h>
#include <static_sequence.h>
#include <static_sequence_wvtree_noptrs.h>

// Re-define word size to ulong:
#undef W
#if __WORDSIZE == 64
#   define W 64
#else
#   define W 32
#endif
#undef bitset
#undef bitget

#include "TextStorage.h"
#include <set>

namespace SXSI 
{

class MismatchOverlap;
class EditDistanceOverlap;
class SuffixFilterOverlap;
class SuffixFilterOverlapEd;

/**
 * Implementation of the TextCollection interface
 *
 */
class TCImplementation : public SXSI::TextCollection {
    // Evil friend hack: Queries access C[] and alphabetrank
    friend class MismatchOverlap;
    friend class EditDistanceOverlap;
    friend class SuffixFilterOverlap;
    friend class SuffixFilterOverlapEd;
public:
    TCImplementation(uchar *, ulong, unsigned, unsigned, ulong, ulong, char, bool);
    ~TCImplementation();

    bool EmptyText(DocId) const;

    /**
     * Extracting one text.
     *
     * Call DeleteText() for each pointer returned by GetText()
     * to avoid possible memory leaks.
     */
    uchar * GetText(DocId) const;
    void DeleteText(uchar *text) const
    { textStorage->DeleteText(text); }

    /**
     * Returns a pointer to the beginning of texts i, i+1, ..., j.
     * Texts are separated by a '\0' byte.
     *
     * Call DeleteText() for each pointer returned by GetText()
     * to avoid possible memory leaks.
     */
    uchar * GetText(DocId i, DocId j) const
    { return textStorage->GetText(i, j); }

    bool isColorCoded() const
    { return color; }
    unsigned GetNumberOfTexts() const
    { return numberOfTexts; }
      
    // Index from/to disk
    TCImplementation(FILE *, unsigned);
    void Save(FILE *) const;


    /**
     * Enumerate document+position pairs (full_result) of
     * each suffix in given interval.
     */
    typedef std::vector<std::pair<DocId, TextPosition> > full_result;
    inline void EnumeratePositions(full_result &result, TextPosition sp, TextPosition ep) const
    {
        result.reserve(result.size() + ep-sp+1);

        uint tmp_rank_c = 0; // Cache rank value of c.
        for (; sp <= ep; ++sp)
        {
            TextPosition i = sp;
            TextPosition dist = 0;
            uchar c = alphabetrank->access(i, tmp_rank_c);
            while (c != '\0' && !sampled->access(i))
            {
                i = C[c]+tmp_rank_c-1; //alphabetrank->rank(c,i)-1;
                c = alphabetrank->access(i, tmp_rank_c);
                ++ dist;
            }
            if (c == '\0')
            {
                // Rank among the end-markers in BWT
                unsigned endmarkerRank = tmp_rank_c-1; //alphabetrank->rank(0, i) - 1;
                DocId docId = Doc->access(endmarkerRank);
                result.push_back(make_pair(docId, dist)); 
            }
            else
            {
                TextPosition textPos = (*suffixes)[sampled->rank1(i)-1] + dist;
                DocId docId = (*suffixDocId)[sampled->rank1(i)-1];

                result.push_back(make_pair(docId, textPos));
            }
        }
    }


private:
    typedef std::vector<std::pair<ulong, ulong> > suffix_range_vector;

    static const uchar versionFlag;
    TextPosition n;
    unsigned samplerate;
    unsigned C[256];
    TextPosition bwtEndPos;
    static_sequence * alphabetrank;

    bool color;

    // Sample structures for texts longer than samplerate
    static_bitsequence * sampled;
    BlockArray *suffixes;
    BlockArray *suffixDocId;

    // Total number of texts in the collection
    unsigned numberOfTexts;
    // Length of the longest text
    ulong maxTextLength;

    // Array of document id's in the order of end-markers in BWT
    ArrayDoc *Doc;

    // Text storage for fast extraction
    TextStorage * textStorage;

    // Following methods are not part of the public API
    uchar * BWT(uchar *);
    void makewavelet(uchar *);
    void maketables(ulong, char);
    DocId DocIdAtTextPos(BlockArray*, TextPosition) const;


}; // class TCImplementation

} // namespace SXSI

#endif
