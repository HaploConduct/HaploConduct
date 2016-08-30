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

#ifndef _SXSI_OverlapQuery_h_
#define _SXSI_OverlapQuery_h_

#include "TextCollection.h"
#include "TCImplementation.h"
#include "Tools.h" // Defines ulong and uchar.
#include <vector>
#include <utility> // Defines std::pair.

namespace SXSI
{
/**
 * General interface for overlap queries
 *
 * Class is virtual
 */

class OverlapQuery
{
public:
    // Type of document identifier
    typedef TextCollection::DocId DocId;
    // Type for text position
    typedef TextCollection::TextPosition TextPosition;

    // Data type for results
    struct dist_result
    {
        DocId idB;
        TextPosition OLA;
        TextPosition OLB;
        TextPosition offsetB;
        unsigned k;
    dist_result(DocId id, TextPosition ola, TextPosition olb, TextPosition off, unsigned k_ = 255)
        : idB(id), OLA(ola), OLB(olb), offsetB(off), k(k_)
        { }
    };
    typedef std::vector<struct dist_result> full_dist_result;

    struct FullDistResultComp
    {
      bool operator()(const struct dist_result & s, const struct dist_result & e)
      {
	return s.idB < e.idB;
      }
    };

    virtual full_dist_result SuffixPrefixMatch(TCImplementation const *, uchar const *) = 0;

    /**
     * Virtual destructor
     */
    virtual ~OverlapQuery() { };
        
    void setInsideAlignments(bool t)
    {
      outputInsideAlignments = t;
    }

protected:
    // Output also alignments that go inside reads (not only suffix/prefix overlaps)
    bool outputInsideAlignments;

    // Protected constructor
    OverlapQuery() { };

    // No copy constructor or assignment
    OverlapQuery(OverlapQuery const&);
    OverlapQuery& operator = (OverlapQuery const&);
};
}
#endif


