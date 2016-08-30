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

#ifndef _SXSI_SuffixFilterOverlap_h_
#define _SXSI_SuffixFilterOverlap_h_

#include "OverlapQuery.h"
#include "TCImplementation.h"
#include "EditDistance.h"

namespace SXSI
{
/**
 * Overlaps with suffix filters (k-mismatch)
 * 
 * FIXME Does not output alignments inside reads
 * 
 * Note: N's are not wildchars.
 */
class SuffixFilterOverlap : public OverlapQuery
{
public:
    SuffixFilterOverlap(unsigned t, unsigned r)
      : result(full_dist_result()), tc(0), pattern(0), patlen(0), threshold(t), errorrate(r)
    {
        if (threshold == 0)
            std::cerr << "Threshold must be > 0!" << endl;

        if (errorrate == 0)
            std::cerr << "Error-rate must be > 0!" << endl;
        if (errorrate < 10)
            std::cerr << "Warning: Error-rate was < 10!" << endl;
    }

    void checkMatch(ulong sp, ulong ep, ulong i)
    {
        assert((patlen - i) >= threshold);

        unsigned l = patlen - i;
        unsigned maxk = l % errorrate == 0 ? l / errorrate : l / errorrate + 1;
	//        cout << "checkMatch at " << i << ", overlap length = " << l << ", maxk = " << maxk << ", errorrate = " << errorrate << endl;

	// Map to endmarkers in Doc
	if (sp != 0)
            sp = tc->alphabetrank->rank(0, sp - 1);
        ep = tc->alphabetrank->rank(0, ep);
        if (ep == 0)
	  return;
        --ep;

        for (unsigned j = sp; j <= ep; ++j)
        {
	    uchar *text = tc->GetText(tc->Doc->access(j));
//            cout << "comparing: " << endl << pattern << endl << text << endl;
            uchar const *tmp = pattern + i;
            
            unsigned k = maxk + 1;
            while (*tmp != 0 && *text != 0 && k > 0)
                if (*tmp++ != *text++)
                    --k;

            if (*tmp == 0 && k > 0)
	      result.push_back(dist_result(tc->Doc->access(j), l, l, 0, maxk - k + 1));
        }
    }

    void ExactSuffixPrefix(ulong sp, ulong ep, unsigned i, unsigned split)
    {
        int c;
        while (sp<=ep && i != splitPos[split])
        {
            c = (int)pattern[--i];
            sp = tc->C[c]+tc->alphabetrank->rank(c,sp-1);
            ep = tc->C[c]+tc->alphabetrank->rank(c,ep)-1;
        }    
    
        if (sp > ep)
            return;

        // After the threshold, check endmarkers
        if ((patlen - i) >= threshold)
            checkMatch(sp, ep, i);
        
        if (split > 0)
            suffixFilterNextStep(sp, ep, i, 1, split-1);
    }

    void suffixFilterNextStep(ulong sp, ulong ep, ulong j, unsigned k, unsigned split)
    {
      // Have to check after each step here
        if ((patlen - j) >= threshold)
	    checkMatch(sp, ep, j);

        if (j == splitPos[split]) 
        {
            if (j == 0)
                return;

            ++k; // Next piece
            --split;
        }
        else if (k == 0)
        {
            ExactSuffixPrefix(sp, ep, j, split);
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
                suffixFilterNextStep(spnew, epnew, j-1, k, split);
            else
                suffixFilterNextStep(spnew, epnew, j-1, k-1, split);
        }
    }

    void ExactSecondPhase(ulong sp, ulong ep, unsigned i, unsigned split)
    {
        int c;
        while (sp<=ep && i != splitPos[split])
        {
            c = (int)pattern[--i];
            sp = tc->C[c]+tc->alphabetrank->rank(c,sp-1);
            ep = tc->C[c]+tc->alphabetrank->rank(c,ep)-1;

	    // After the threshold, check endmarkers
	    if ((patlen - i) >= threshold)
	      checkMatch(sp, ep, i);
        }    
    
        if (sp > ep)
            return;
        
        if (split > 0)
            suffixFilterSecondPhase(sp, ep, i, 1, split-1);
    }

    void suffixFilterSecondPhase(ulong sp, ulong ep, ulong j, unsigned k, unsigned split)
    {
        if ((patlen - j) >= threshold)
	    checkMatch(sp, ep, j);

        if (j == splitPos[split]) 
        {
            if (j == 0)
                return;

            ++k; // Next piece
            --split;
        }

	if (k == 0)
        {
            ExactSecondPhase(sp, ep, j, split);
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
                suffixFilterSecondPhase(spnew, epnew, j-1, k, split);
            else
                suffixFilterSecondPhase(spnew, epnew, j-1, k-1, split);
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

        unsigned klimit = threshold%errorrate == 0 ? threshold/errorrate : threshold/errorrate + 1;
        psize = threshold%(klimit+1) == 0 ? threshold / (klimit + 1) : threshold / (klimit + 1) + 1;
//        cout << "orig p = " << psize << ", ";
        for (unsigned i = threshold + 1; i < patlen; ++i)
        {
            unsigned tmp = i%errorrate == 0 ? i/errorrate : i/errorrate + 1;
            if (psize > (i%(tmp+1) == 0 ? i / (tmp + 1) : i / (tmp + 1) + 1))
                psize = i%(tmp+1) == 0 ? i / (tmp + 1) : i / (tmp + 1) + 1;
        }
//        cout << "final psize = " << psize << endl;


        // uniform partition of klimit pieces  
        klimit = patlen%psize == 0 ? patlen/psize : patlen/psize + 1;
        assert(klimit < 1000);

        splitPos[klimit] = patlen;
	//cout << "splitPos: 0";
        for (unsigned i = klimit-1; i > 0; --i)
        {
            splitPos[i] = splitPos[i+1] - psize;
	    // cout << ", " << splitPos[i];
        }
        splitPos[0]=0;
        //cout << ", " << splitPos[klimit] << endl;


     
	// PHASE 1
	unsigned i = 2; // Skip the very first piece (check it in next phase)
	for (; i < klimit; i++) {// Last filter (splitPos[klimit]) is checked in second phase
	  //std::cout << "new try: " << splitPos[i] << "\n";
	  ExactSuffixPrefix(0, tc->n - 1, splitPos[i], i-1);
	}

	// PHASE 2
	// Start with one mismatch
	suffixFilterSecondPhase(0, tc->n - 1, patlen, 1, klimit-1);

        return result;
    }

private:
    full_dist_result result;
    TCImplementation const *tc;
    uchar const * pattern;
    unsigned patlen;
    unsigned threshold;
    unsigned errorrate;
    unsigned psize;
    unsigned splitPos[1000];
};
} // Namespace
#endif


