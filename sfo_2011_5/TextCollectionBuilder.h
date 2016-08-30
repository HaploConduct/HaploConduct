/******************************************************************************
 *   Copyright (C) 2009 by Niko Valimaki <nvalimak@cs.helsinki.fi>            *
 *   Text collection interface for an in-memory XQuery/XPath engine           *
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

#ifndef _SXSI_TextCollectionBuilder_h_
#define _SXSI_TextCollectionBuilder_h_

#include "TextCollection.h"
#include "TextStorage.h"
#include "Tools.h" // Defines ulong and uchar.
#include <vector>
#include <utility> // Defines std::pair.
#include <cstring> // Defines std::strlen, added by Kim

// Un-comment to compare BWT against a BWT generated from class dynFMI:
//#define TCB_TEST_BWT

// Default samplerate for suffix array samples
#define TEXTCOLLECTION_DEFAULT_SAMPLERATE 16

// Default input length, used to calculate the buffer size.
#define TEXTCOLLECTION_DEFAULT_INPUT_LENGTH (950 * 1024 * 1024)


namespace SXSI
{
    struct TCBuilderRep; // Pimpl
    
    /**
     * Build an instance of the TextCollection class.
     */
    class TextCollectionBuilder
    {
    public:
        explicit TextCollectionBuilder(bool color, unsigned samplerate = TEXTCOLLECTION_DEFAULT_SAMPLERATE, 
                                       ulong estimatedInputLength =  TEXTCOLLECTION_DEFAULT_INPUT_LENGTH);
        ~TextCollectionBuilder();
        
        /** 
         * Insert text
         *
         * Must be a zero-terminated string from alphabet [1,255].
         * Can not be called after makeStatic().
         * The i'th text insertion gets an identifier value i-1.
         * In other words, document identifiers start from 0.
         */
        void InsertText(uchar const *);
        /**
         * Make static
         *
         * Convert to a static collection.
         * New texts can not be inserted after this operation.
         *
         * TextStorage type defaults to TYPE_PLAIN_TEXT, another
         * possible type is TYPE_LZ_INDEX.
         */
        TextCollection * InitTextCollection(char type = TextStorage::TYPE_PLAIN_TEXT);
        
    private:
        struct TCBuilderRep * p_;

        // No copy constructor or assignment
        TextCollectionBuilder(TextCollectionBuilder const&);
        TextCollectionBuilder& operator = (TextCollectionBuilder const&);
    };
}
#endif
