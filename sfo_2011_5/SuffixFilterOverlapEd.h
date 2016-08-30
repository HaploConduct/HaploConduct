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

#ifndef _SXSI_SuffixFilterOverlapEd_h_
#define _SXSI_SuffixFilterOverlapEd_h_

#include "OverlapQuery.h"
#include "TCImplementation.h"
#include "EditDistance.h"

namespace SXSI
{
/**
 * Overlaps with suffix filters (k-errors)
 * 
 * Note: N's are not wildchars.
 */
class SuffixFilterOverlapEd : public OverlapQuery
{
public:
    SuffixFilterOverlapEd(unsigned t, unsigned r)
        : result(full_dist_result()), tc(0), pattern(0), patlen(0), threshold(t), errorrate(r), ed(0)
    {
        if (threshold == 0)
            std::cerr << "Threshold must be > 0!" << endl;

        if (errorrate == 0)
            std::cerr << "Error-rate must be > 0!" << endl;
        if (errorrate < 10)
            std::cerr << "Warning: Error-rate was < 10!" << endl;
    }

    void checkMatch(ulong sp, ulong ep)
    {
        assert(patlen >= threshold);
        
        unsigned l = patlen;
        unsigned maxk = l % errorrate == 0 ? l / errorrate : l / errorrate + 1;
//        cout << "overlap length = " << l << ", maxk = " << maxk << ", errorrate = " << errorrate << endl;

        TCImplementation::full_result result;
        tc->EnumeratePositions(result, sp, ep);

//        cout << "suffixlen = " << l << endl;
//        cout << "suffix = " << pattern << endl;

        MyersEditDistance * e = new MyersEditDistance(l, maxk);
        for (TCImplementation::full_result::iterator it = result.begin(); it != result.end(); ++it)
        {
//            MyersEditDistance * e = new MyersEditDistance(pattern, l, maxk);
            e->init(pattern);
            unsigned docId = (*it).first;
            ulong offset = (*it).second;
	    if (offset == 0)
	      continue; // Do not output suffix/prefix matches here, just inclusions (offset > 0)

	    uchar *text = tc->GetText(docId) + offset;
//            cout << "comparing: " << endl << pattern << endl << text << endl;

            int candlen = strlen((char*)text) < l+maxk ? strlen((char *)text) : l+maxk;
//            cout << "cand = " << text << ", len = " << candlen << endl;
            pair<unsigned, unsigned> dist = e->prefixDist(text, candlen);
//            cout << "dist.second = " << dist.second << ", " << l << ", distance = " << dist.first << "(max " << maxk << ")" << endl;
            if (dist.first <= maxk)
	      this->result.push_back(dist_result(docId, l, dist.second, offset, dist.first));
        }
        delete e;
    }


    void checkMatch(ulong sp, ulong ep, ulong i)
    {
        if (outputInsideAlignments)
        { // Check "included" overlaps
	    // Max number of errors for the end block
	    unsigned maxk = ed->length%errorrate == 0 ? ed->length/errorrate-1 : ed->length/errorrate;
	    if (ed->length == patlen)
	      ++ maxk; // +1 in the phase 2
	    if (ed->dist() <= maxk)
	    {
	      checkMatch(sp, ep);
	      return;
	    }
	}

        assert((patlen - i) >= threshold);
        
        unsigned l = patlen - i;
//        cout << "overlap length = " << l << ", maxk = " << maxk << ", errorrate = " << errorrate << endl;

	// Map to endmarkers in Doc
	if (sp != 0)
            sp = tc->alphabetrank->rank(0, sp - 1);
        ep = tc->alphabetrank->rank(0, ep);
        if (ep == 0)
            return;
            --ep;

        unsigned suffixlen = ed->getSuffixLength(); // Pick the suffix having smallest distance
        suffixlen += patlen-ed->length;

	unsigned maxk = suffixlen % errorrate == 0 ? suffixlen / errorrate : suffixlen / errorrate + 1;
//        cout << "suffixlen = " << suffixlen << endl;

        uchar const *tmp = pattern + (patlen - suffixlen);
//        cout << "suffix = " << tmp << endl;
        MyersEditDistance * e = new MyersEditDistance(suffixlen, maxk);
        for (unsigned j = sp; j <= ep; ++j)
        {
            e->init(tmp);
	    uchar *text = tc->GetText(tc->Doc->access(j));
//            cout << "comparing: " << endl << pattern << endl << text << endl;

            int candlen = strlen((char*)text) < l+maxk ? strlen((char *)text) : l+maxk;
//            cout << "cand = " << text << ", len = " << candlen << endl;
            pair<unsigned, unsigned> dist = e->prefixDist(text, candlen);
//            cout << "dist.second = " << dist.second << ", " << l << ", distance = " << dist.first << "(max " << maxk << ")" << endl;
            if (dist.first <= maxk)
                this->result.push_back(dist_result(tc->Doc->access(j), suffixlen, dist.second, 0, dist.first));

        }
        delete e;
    }
    
    void ExactSuffixPrefix(ulong sp, ulong ep, unsigned i, unsigned split)
    {
        int c;
        while (sp<=ep && i != splitPos[split])
        {
            c = (int)pattern[--i];
            sp = tc->C[c]+tc->alphabetrank->rank(c,sp-1);
            ep = tc->C[c]+tc->alphabetrank->rank(c,ep)-1;
	    ed->pushChar(c);
        }    
    
        if (sp > ep)
            return;

        if (split > 0)
            suffixFilterNextStep(sp, ep, i, split);
	else
	  if ((patlen - i) >= threshold)
            checkMatch(sp, ep, i);
    }

    void suffixFilterNextStep(ulong sp, ulong ep, ulong j, const unsigned split)
    {
        // Have to check at each step
        if ((patlen - j) >= threshold)
	    checkMatch(sp, ep, j);
	if (j == 0)
	  ++j; // Avoid ~0u
        
        if (ed->end())
            return;

        for (const char *cp = TextCollection::ALPHABET; 
             cp < TextCollection::ALPHABET + TextCollection::ALPHABET_SIZE; ++cp)
        {
            ulong spnew = tc->C[(int)*cp]+tc->alphabetrank->rank(*cp,sp-1);
            ulong epnew = tc->C[(int)*cp]+tc->alphabetrank->rank(*cp,ep)-1;
            if (spnew>epnew) 
                continue; // No match

	    const unsigned colmin = ed->pushChar(*cp);
            // Edit distance is for reversed pattern, compute prefix lenght:
            const unsigned prefixl = ed->length - ed->getSuffixLength();
            // How many errors are we allowing on a prefix of given length?
            unsigned i = 1;
            while (splitPos[i] < prefixl)
                ++i;
            
            assert(i-1 < (patlen%psize==0 ? patlen/psize : patlen/psize+1));
            assert(split > i-1);
            const unsigned maxk = split - (i-1);

	    if (colmin <= maxk) 
                suffixFilterNextStep(spnew, epnew, j-1, split);
	    ed->popChar();
        }
    }

    void suffixFilterStartOne(ulong sp, ulong ep, ulong j, const unsigned split)
    {
	if ((patlen - j) >= threshold)
	  checkMatch(sp, ep, j);

        if (ed->end())
            return;

        for (const char *cp = TextCollection::ALPHABET; 
             cp < TextCollection::ALPHABET + TextCollection::ALPHABET_SIZE; ++cp)
        {
            ulong spnew = tc->C[(int)*cp]+tc->alphabetrank->rank(*cp,sp-1);
            ulong epnew = tc->C[(int)*cp]+tc->alphabetrank->rank(*cp,ep)-1;
            if (spnew>epnew) 
                continue; // No match

	    const unsigned colmin = ed->pushChar(*cp);
            // Edit distance is for reversed pattern, compute prefix lenght:
            const unsigned prefixl = ed->length - ed->getSuffixLength();

            // How many errors are we allowing on a prefix of given length?
            unsigned i = 1;
            while (splitPos[i] < prefixl)
                ++i;
            
            assert(i-1 < (patlen%psize==0 ? patlen/psize : patlen/psize+1));
            assert(split >= i-1);
            const unsigned maxk = split - (i-1) + 1; // +1 to start with one error

	    //	    std::cout << "push = " << *cp << ", colmin = " << colmin << " (" << maxk << "), prefixl = " << prefixl << ", split = " << s << " (" << split << "), j = " << j << std::endl;

	    if (colmin <= maxk) 
                suffixFilterStartOne(spnew, epnew, j-1, split);
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

	// Find the largest psize we are allowed to choose
        unsigned klimit = threshold%errorrate == 0 ? threshold/errorrate : threshold/errorrate + 1;
        psize = threshold%(klimit+1) == 0 ? threshold / (klimit + 1) : threshold / (klimit + 1) + 1;
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
	    // sanity check
	    if (splitPos[i+1] < psize)
	      {
		std::cerr << "error: splitPos[" << i << "] < 0, given: patlen = " << patlen << ", psize = " << psize << std::endl;
		std::exit(0);
	      }

            splitPos[i] = splitPos[i+1] - psize;
	    //  	    cout << ", " << splitPos[i];
        }
        splitPos[0]=0;
	//cout << ", " << splitPos[klimit] << endl;
     
	uchar *reversepat = new uchar[patlen+1];
	
	for (unsigned i = 2; i < klimit; i++) {  // Last filter (splitPos[klimit]) is checked in second phase
	  //                std::cout << "new try: " << splitPos[i] << "\n";
	  for (unsigned jj = 0; jj < splitPos[i]; jj++)
	    reversepat[jj] = pattern[splitPos[i]-jj-1];
	  reversepat[splitPos[i]] = 0;
	  ed = new MyersEditDistanceIncremental(reversepat, splitPos[i], i-1);
	  
	  ExactSuffixPrefix(0, tc->n - 1, splitPos[i], i-1);
	  
	  delete ed; ed = 0;
	}

	// Second phase: start the search with 1-error
	for (unsigned jj = 0; jj < splitPos[klimit]; jj++)
	    reversepat[jj] = pattern[splitPos[klimit]-jj-1];
	reversepat[patlen] = 0;
	ed = new MyersEditDistanceIncremental(reversepat, patlen, klimit);
	  
	suffixFilterStartOne(0, tc->n - 1, splitPos[klimit], klimit-1);
	  
	delete ed; ed = 0;

	delete [] reversepat;
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

    MyersEditDistanceIncremental * ed;
};
} // Namespace
#endif
