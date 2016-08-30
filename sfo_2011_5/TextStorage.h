/******************************************************************************
 *   Copyright (C) 2009 Niko Välimäki                                         *
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

#ifndef _TextStorage_H_
#define _TextStorage_H_

#include "TextCollection.h"
#include "Tools.h"
//#include "incbwt/bits/deltavector.h"
#include <cassert>
#include <stdexcept>


// Include  from XMLTree/libcds
#include <static_bitsequence.h>

// Re-define word size to ulong:
#undef W
#if __WORDSIZE == 64
#   define W 64
#else
#   define W 32
#endif
#undef bitset
#undef bitget


namespace SXSI 
{

class TextStorageBuilder;
class TextStoragePlainText;
//class TextStorageLzIndex;

/**
 * Text collection that supports fast extraction.
 * Defines an abstact interface class.
 * See subclasses TextStorageLzIndex and TextStoragePlainText
 * below.
 */
class TextStorage
{
public:
    // Storage type
    const static char TYPE_PLAIN_TEXT = 0;
    const static char TYPE_LZ_INDEX = 1;

    // Call DeleteText() for each pointer returned by GetText()
    // to avoid possible memory leaks.
    virtual uchar * GetText(TextCollection::DocId docId) const = 0;
    virtual uchar * GetText(TextCollection::DocId i, TextCollection::DocId j) const = 0;
    virtual void DeleteText(uchar *) const = 0;

    static TextStorage * Load(FILE *file);
    virtual void Save(FILE *file) const = 0;

    virtual ~TextStorage()
    {
        delete offsets_;
        offsets_ = 0;
    }

    TextCollection::DocId DocIdAtTextPos(TextCollection::TextPosition i) const
    {
        assert(i < n_);
        return offsets_->rank1(i)-1;
    }

    TextCollection::TextPosition TextStartPos(TextCollection::DocId i) const
    {
        assert(i < (TextCollection::DocId)numberOfTexts_);
        return offsets_->select1(i+1);
    }

    bool IsEndmarker(TextCollection::TextPosition i) const
    {
        assert(i < n_);
        if (i >= n_ - 1)
            return true;
        return offsets_->access(i+1);
    }


protected:
    // Define a shortcut
    typedef TextCollection::TextPosition TextPosition;
    // Block size in DeltaVector
    //    const static CSA::usint DV_BLOCK_SIZE = 32;

    TextStorage(uchar const * text, TextPosition n)
        : n_(n), offsets_(0), numberOfTexts_(0)
    { 
        // Delta encoded bitvector of text offsets.
      /*        CSA::DeltaEncoder encoder(DV_BLOCK_SIZE);
        encoder.setBit(0); // Start of the first text.
      */
      uint *startpos = new uint[n/(sizeof(uint)*8)+1];
      for (unsigned long i = 0; i < n / (sizeof(uint)*8) + 1; i++)
        startpos[i] = 0;

        // Read offsets by finding text end positions:
      set_field(startpos,1,0,1);
        for (TextPosition i = 0; i < n_ - 1; ++i)
            if (text[i] == '\0')
	      set_field(startpos,1,i+1,1);
//encoder.setBit(i+1);
        

//        offsets_ = new CSA::DeltaVector(encoder, n_);

	offsets_ = new static_bitsequence_brw32(startpos, n, 16);
	delete [] startpos;
	

        for (ulong i = 0; i < n_-1; ++i)
            if ((text[i] == '\0') != IsEndmarker(i))
                std::cout << "misplaced endmarker at i = " << i << std::endl;

        numberOfTexts_ = offsets_->rank1(n_ - 1);
    }
    
    TextStorage(std::FILE *);
    void Save(FILE *file, char type) const;

    TextPosition n_;
    //CSA::DeltaVector *offsets_;
    static_bitsequence * offsets_ ;
    TextPosition numberOfTexts_;
};

/******************************************************************
 * Plain text collection.
 */
class TextStoragePlainText : public TextStorage
{
public:
    TextStoragePlainText(uchar *text, TextPosition n)
        : TextStorage(text, n), text_(text)
    { }

    TextStoragePlainText(FILE *file)
        : TextStorage(file), text_(0)
    {
        text_ = new uchar[n_];
        if (std::fread(this->text_, sizeof(uchar), n_, file) != n_)
            throw std::runtime_error("TextStorage::Load(): file read error (text_).");
    }

    void Save(FILE *file) const
    {
        TextStorage::Save(file, TYPE_PLAIN_TEXT);

        if (std::fwrite(this->text_, sizeof(uchar), n_, file) != n_)
            throw std::runtime_error("TextStorage::Save(): file write error (text_).");
    }

    ~TextStoragePlainText()
    {
        delete [] text_;
        text_ = 0;
        n_ = 0;
    }

    uchar * GetText(TextCollection::DocId docId) const
    {
        assert(docId < (TextCollection::DocId)numberOfTexts_);

        TextPosition offset = offsets_->select1(docId+1);
        return &text_[offset];
    }

    uchar * GetText(TextCollection::DocId i, TextCollection::DocId j) const
    {
        assert(i < (TextCollection::DocId)numberOfTexts_);
        assert(j < (TextCollection::DocId)numberOfTexts_);

        TextPosition offset = offsets_->select1(i+1);
        return &text_[offset];
    }

    // No operation, since text is a pointer to this->text_ 
    void DeleteText(uchar *text) const
    { }

private:
    uchar *text_;
}; // class TextStorage


/**
 * Builder for TextStorage class
 */
class TextStorageBuilder
{
public:
    // Define a shortcut
    typedef TextCollection::TextPosition TextPosition;

    // Build up simple uchar array
    explicit TextStorageBuilder(TextPosition n)
        : n_(n), text_(new uchar [n]), freeText(true)
    { }
    
    ~TextStorageBuilder()
    {
        if (freeText)
            delete [] text_;
        text_ = 0;
        n_ = 0;
    }
    
    // Write access to text[]
    uchar& operator[] (TextPosition i)
    {
        return text_[i];
    }

    // Init TextStorage
    // Type defaults to plain text.
    TextStorage * InitTextStorage(char type = TextStorage::TYPE_PLAIN_TEXT)
    {
        freeText = false; // Passing text to TextStorage.
        switch(type)
        {
        case (TextStorage::TYPE_PLAIN_TEXT):
            return new TextStoragePlainText(text_, n_);
//        case (TextStorage::TYPE_LZ_INDEX):
//            return new TextStorageLzIndex(text_, n_);
        default:
            std::cerr << "TextStorageBuilder: Unknown type given!" << std::endl;
            exit(1);
        }
    }

private:
    TextPosition n_;
    uchar *text_; // FIXME Replace with a succinct representation.
    bool freeText;
}; // class TextStorageBuilder

} // namespace SXSI

#endif // #ifndef _TextStorage_H_
