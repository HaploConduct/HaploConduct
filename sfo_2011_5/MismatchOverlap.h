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

#ifndef _SXSI_MismatchOverlap_h_
#define _SXSI_MismatchOverlap_h_

#include "OverlapQuery.h"
#include "TCImplementation.h"

namespace SXSI
{
/**
 * Overlaps with k-mismatches
 * Note: N's are not wildchars.
 */
class MismatchOverlap : public OverlapQuery
{
public:
    MismatchOverlap(unsigned t, unsigned k)
        : result(full_dist_result()), tc(0), pattern(0), patlen(0), threshold(t), maxk(k)
    {
        if (threshold == 0)
            std::cerr << "Threshold must be > 0!" << endl;
    }

    
    void ExactSuffixPrefix(ulong sp, ulong ep, ulong i)
    {
        int c;
        while (sp<=ep && i >= 1 && (patlen - i) < threshold) 
        {
            c = (int)pattern[--i];
            sp = tc->C[c]+tc->alphabetrank->rank(c,sp-1);
            ep = tc->C[c]+tc->alphabetrank->rank(c,ep)-1;
        }    
    
        //After the threadhold, count the ones that finish inside the suffix of 'pattern'
        while (sp<=ep && i >= 1) 
        {
            //document_result dr = tc->EnumerateEndmarkers(sp, ep);
            // Map to endmarkers in Doc
            unsigned l = 0;
            if (sp != 0)
                l = tc->alphabetrank->rank(0, sp - 1);
            unsigned r = tc->alphabetrank->rank(0, ep);
            if (r > 0)
                for (unsigned j = l; j <= r-1; ++j)
                    result.push_back(dist_result(tc->Doc->access(j), (patlen - i), (patlen - i), 0));
        
            c = (int)pattern[--i];
            sp = tc->C[c]+tc->alphabetrank->rank(c,sp-1);
            ep = tc->C[c]+tc->alphabetrank->rank(c,ep)-1;
        }

        if (sp<=ep)
        {
/*            document_result dr = tc->EnumerateEndmarkers(sp, ep);
            for (unsigned ii=0; ii < dr.size(); ii++){
                result.push_back(make_pair(dr[ii], make_pair(0, (patlen - i))));
                }*/
            // Map to endmarkers in Doc
            if (sp != 0)
                sp = tc->alphabetrank->rank(0, sp - 1);
            ep = tc->alphabetrank->rank(0, ep);
            if (ep == 0)
                return;
            --ep;
            
            for (unsigned j = sp; j <= ep; ++j)
                result.push_back(dist_result(tc->Doc->access(j), (patlen - i), (patlen - i), 0));
        }
    }

    void KMismatchSuffixPrefix(ulong sp, ulong ep, ulong j, unsigned k)
    {
        if((patlen - j) >= threshold){        
/*            document_result dr = tc->EnumerateEndmarkers(sp, ep);
            for (unsigned ii=0; ii < dr.size(); ii++){
                result.push_back(make_pair(dr[ii], make_pair(0, (patlen - j))));
                }*/
            unsigned l = 0;
            if (sp != 0)
                l = tc->alphabetrank->rank(0, sp - 1);
            unsigned r = tc->alphabetrank->rank(0, ep);
            if (r > 0)
                for (unsigned i = l; i <= r-1; ++i)
                    result.push_back(dist_result(tc->Doc->access(i), (patlen - j), (patlen - j), 0));
        }

        if (j == 0)
            return;

        if (k == 0)
        {
            ExactSuffixPrefix(sp, ep, j);
            return;
        }


        for (const char *cp = TextCollection::ALPHABET; 
             cp < TextCollection::ALPHABET + TextCollection::ALPHABET_SIZE; ++cp)
        {
            ulong spnew = tc->C[(int)*cp]+tc->alphabetrank->rank(*cp,sp-1);
            ulong epnew = tc->C[(int)*cp]+tc->alphabetrank->rank(*cp,ep)-1;
            if (spnew>epnew) 
                continue; // No match

            if (*cp == pattern[j-1]) 
                KMismatchSuffixPrefix(spnew, epnew, j-1, k);
            else
                KMismatchSuffixPrefix(spnew, epnew, j-1, k-1);
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

        KMismatchSuffixPrefix(0, tc->n - 1, patlen, maxk);
        return result;
    }

    private:
        full_dist_result result;
        TCImplementation const *tc;
        uchar const * pattern;
        unsigned patlen;
        unsigned threshold;
        unsigned maxk;
    };
} // Namespace
#endif


