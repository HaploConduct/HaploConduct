/******************************************************************************
 *   Copyright (C) 2009 by Niko Valimaki <nvalimak@cs.helsinki.fi>            *
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
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                *
 ******************************************************************************/ 

#ifndef _SXSI_EditDistanceOverlap_h_
#define _SXSI_EditDistanceOverlap_h_

#include "EditDistance.h"
#include "OverlapQuery.h"
#include "TCImplementation.h"

namespace SXSI
{
/**
 * Overlaps with k-mismatches
 * Note: N's are not wildchars.
 */
class EditDistanceOverlap : public OverlapQuery
{
public:
    EditDistanceOverlap(unsigned t, unsigned k)
        : result(full_dist_result()), tc(0), pattern(0), patlen(0), threshold(t), maxk(k), ed(0)
    {
        if (threshold == 0)
            std::cerr << "Threshold must be > 0!" << endl;
    }

    void KErrorSuffixPrefix(ulong sp, ulong ep, unsigned j)
    {
        if(j >= threshold){        
            unsigned l = 0;
            if (sp != 0)
                l = tc->alphabetrank->rank(0, sp - 1);
            unsigned r = tc->alphabetrank->rank(0, ep);
	    const unsigned sl = ed->getSuffixLength();
	    const unsigned min = ed->getMinimum();
            if (r > 0)
                for (unsigned i = l; i <= r-1; ++i)
		  result.push_back(dist_result(tc->Doc->access(i), sl, j, 0, min));

	    if (ed->dist() <= maxk)
	    {
	        TCImplementation::full_result result;
		tc->EnumeratePositions(result, sp, ep);
		MyersEditDistance * e = new MyersEditDistance(patlen, maxk);
		for (TCImplementation::full_result::iterator it = result.begin(); it != result.end(); ++it)
		  {
		    e->init(pattern);
		    unsigned docId = (*it).first;
		    ulong offset = (*it).second;
		    if (offset == 0)
		      continue; // Do not output suffix/prefix matches here, just inclusions (offset > 0)

		    uchar *text = tc->GetText(docId) + offset;
		    int candlen = strlen((char*)text) < patlen+maxk ? strlen((char *)text) : patlen+maxk;
		    pair<unsigned, unsigned> dist = e->prefixDist(text, candlen);
		    if (dist.first <= maxk)
		      this->result.push_back(dist_result(docId, patlen, dist.second, offset, dist.first));
		  }
	    }
        }

        if (ed->end())
            return;

        for (const char *cp = TextCollection::ALPHABET; 
             cp < TextCollection::ALPHABET + TextCollection::ALPHABET_SIZE; ++cp)
        {
            ulong spnew = tc->C[(int)*cp]+tc->alphabetrank->rank(*cp,sp-1);
            ulong epnew = tc->C[(int)*cp]+tc->alphabetrank->rank(*cp,ep)-1;
            if (spnew>epnew) 
                continue; // No match

//	    if (j <= 3) // Debug print
//		std::cerr << "At depth j = " << j << ", push = " << *cp << std::endl;

            unsigned min = ed->pushChar(*cp);
            if (min <= maxk)
                KErrorSuffixPrefix(spnew, epnew, j+1);
            ed->popChar();
        }
    }


    full_dist_result SuffixPrefixMatch(TCImplementation const *tc, uchar const *pattern)
    {
        result.clear();
        this->tc = tc;
        this->pattern = pattern;

        this->patlen = strlen((char *)pattern);
        if (patlen < threshold)
            return full_dist_result();

        uchar *tmp = new uchar[patlen+1]; // Reverse
        for (unsigned i = 0; i < patlen; i++) {
            tmp[i] = pattern[patlen - i - 1];
        }
        tmp[patlen] = 0;

        ed = new MyersEditDistanceIncremental(tmp, patlen, maxk);
        KErrorSuffixPrefix(0, tc->n - 1, 0);
        delete ed;
        ed = 0;
        delete [] tmp;
        return result;
    }

    private:
        full_dist_result result;
        TCImplementation const *tc;
        uchar const * pattern;
        unsigned patlen;
        unsigned threshold;
        unsigned maxk;
        
        MyersEditDistanceIncremental *ed;
    };
} // Namespace
#endif


