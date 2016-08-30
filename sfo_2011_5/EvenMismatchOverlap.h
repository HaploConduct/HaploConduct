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

#ifndef _SXSI_EvenMismatchOverlap_h_
#define _SXSI_EvenMismatchOverlap_h_

#include "OverlapQuery.h"
#include "TCImplementation.h"

namespace SXSI
{
/**
 * Overlaps with evenly distributed k-mismatches
 *
 * Note: N's are not wildchars.
 */
class EvenMismatchOverlap : public OverlapQuery
{
public:
    EvenMismatchOverlap(unsigned t, unsigned r)
        : result(full_dist_result()), tc(0), pattern(0), patlen(0), threshold(t), errorrate(r)
    {
        if (threshold == 0)
            std::cerr << "Threshold must be > 0!" << endl;

        if (errorrate == 0)
            std::cerr << "Error-rate must be > 0!" << endl;
        if (errorrate < 10)
            std::cerr << "Warning: Error-rate was < 10!" << endl;
    }

    
    void ExactSuffixPrefix(ulong sp, ulong ep, ulong i)
    {
        int c;
        while (sp<=ep && i >= 1 && (patlen - i) < threshold 
               && (patlen - i) % errorrate) 
        {
            c = (int)pattern[--i];
            sp = tc->C[c]+tc->alphabetrank->rank(c,sp-1);
            ep = tc->C[c]+tc->alphabetrank->rank(c,ep)-1;
        }    
    
        //After the threadhold, count the ones that finish inside the suffix of 'pattern'
        while (sp<=ep && i >= 1 && (patlen - i) % errorrate) 
        {
	  TextCollection::document_result dr = tc->EnumerateEndmarkers(sp, ep);
            for (unsigned ii=0; ii < dr.size(); ii++){
	      result.push_back(dist_result(dr[ii], patlen-i, patlen-i, 0));
            }
        
            c = (int)pattern[--i];
            sp = tc->C[c]+tc->alphabetrank->rank(c,sp-1);
            ep = tc->C[c]+tc->alphabetrank->rank(c,ep)-1;

        }

        if (i > 0 && (patlen - i) % errorrate == 0)
            KMismatchSuffixPrefix(sp, ep, i, 0);
        else if (sp<=ep)
        {
            TextCollection::document_result dr = tc->EnumerateEndmarkers(sp, ep);
            for (unsigned ii=0; ii < dr.size(); ii++){
	      result.push_back(dist_result(dr[ii], patlen-i, patlen-i, 0));
            }
        }
    }

    void KMismatchSuffixPrefix(ulong sp, ulong ep, ulong j, unsigned k)
    {
        if((patlen - j) >= threshold){        
            TextCollection::document_result dr = tc->EnumerateEndmarkers(sp, ep);
            for (unsigned ii=0; ii < dr.size(); ii++){
	      result.push_back(dist_result(dr[ii], patlen-j, patlen-j, 0));
            }
        }

        if (j == 0)
            return;

        if ((patlen - j) % errorrate == 0)
            ++k;

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

        KMismatchSuffixPrefix(0, tc->n - 1, patlen, 0);
        return result;
    }

    private:
        full_dist_result result;
        TCImplementation const *tc;
        uchar const * pattern;
        unsigned patlen;
        unsigned threshold;
        unsigned errorrate;
    };
} // Namespace
#endif


