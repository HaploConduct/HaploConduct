/******************************************************************************
 *   Copyright (C) 2008 by Niko Valimaki <nvalimak@cs.helsinki.fi>            *
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

#ifndef _SXSI_TextCollection_h_
#define _SXSI_TextCollection_h_

#include "Tools.h" // Defines ulong and uchar.
#include <vector>
#include <utility> // Defines std::pair.

namespace SXSI
{
    /**
     * General interface for a text collection
     *
     * Class is virtual, make objects by calling 
     * the static method InitTextCollection().
     */
    class TextCollection
    {
    public:
        // Type of document identifier
        typedef unsigned DocId;
        // Type for text position (FIXME ulong or long?)
        typedef ulong TextPosition;

        // Static A,C,G,T,N alphabet for overlap searches
        static char const *ALPHABET;
        static const unsigned ALPHABET_SIZE = 5;
        static const char ALPHABET_DNA[];
        static const char ALPHABET_SOLID[];

        static const std::string REVERSE_EXTENSION;
        static const std::string FMINDEX_EXTENSION;

        /**
         * Load from a file
         *
         * New samplerate can be given, otherwise will use the one specified in the save file!
         * Note: This is not a static method; call InitTextCollection() first to get the object handle.
         *
         * Throws an exception if std::fread() fails.
         * 
         */
        static TextCollection* Load(FILE *, unsigned samplerate = 0);

        /**
         * Save data structure into a file
         * 
         * Throws an exception if std::fwrite() fails.
         */
        virtual void Save(FILE *) const = 0;

        /**
         * Virtual destructor
         */
        virtual ~TextCollection() { };

	/**
	 * Tests if the string pointed to by DocId is empty
         */
	virtual bool EmptyText(DocId) const = 0;

        /**
         * Displaying content
         *
         * Returns the i'th text in the collection.
         * The numbering starts from 0.
         *
         * Call DeleteText() for each pointer returned by GetText()
         * to avoid possible memory leaks.
         */
        virtual uchar* GetText(DocId) const = 0;
        virtual void DeleteText(uchar *text) const = 0;

        /**
         * Returns a pointer to the beginning of texts i, i+1, ..., j.
         * Texts are separated by a '\0' byte.
         *
         * Call DeleteText() for each pointer returned by GetText()
         * to avoid possible memory leaks.
         */
        virtual uchar * GetText(DocId i, DocId j) const = 0;

    protected:
        // Protected constructor; use TextCollectionBuilder
        TextCollection() { };

        void setColorCoded(bool color)
        {
            if (color)
                TextCollection::ALPHABET = TextCollection::ALPHABET_SOLID;
            else
                TextCollection::ALPHABET = TextCollection::ALPHABET_DNA;
        }

        // No copy constructor or assignment
        TextCollection(TextCollection const&);
        TextCollection& operator = (TextCollection const&);
    };
}
#endif
