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

#include "TextStorage.h"

#undef W
#undef bitsW
#undef bitset
#undef bits

// Re-define word size to ulong:
#undef W
#if __WORDSIZE == 64
#   define W 64
#else
#   define W 32
#endif
#undef bitset
#undef bitget
#undef bits


namespace SXSI
{

/******************************************************************
 * Class TextStorage
 */

TextStorage * TextStorage::Load(std::FILE *file)
{
    char type = 0;
    if (std::fread(&type, sizeof(char), 1, file) != 1)
        throw std::runtime_error("TextStorage::Load(): file read error (type).");

    switch(type)
    {
    case (TYPE_PLAIN_TEXT):
        return new TextStoragePlainText(file);
//    case (TYPE_LZ_INDEX):
//        return new TextStorageLzIndex(file);
    default:
        std::cerr << "TextStorage::Load(): Unknown type in save file!" << std::endl;
        exit(1);
    }
}

TextStorage::TextStorage(std::FILE * file)
    : n_(0), offsets_(0), numberOfTexts_(0)
{
   if (std::fread(&(this->n_), sizeof(TextPosition), 1, file) != 1)
        throw std::runtime_error("TextStorage::Load(): file read error (n_).");

    if (std::fread(&(this->numberOfTexts_), sizeof(TextPosition), 1, file) != 1)
        throw std::runtime_error("TextStorage::Load(): file read error (numberOfTexts_).");

    offsets_ = static_bitsequence::load(file); //new CSA::DeltaVector(file);   
}

void TextStorage::Save(FILE *file, char type) const
{
    if (std::fwrite(&type, sizeof(char), 1, file) != 1)
        throw std::runtime_error("TextStorage::Save(): file write error (type).");
        
    if (std::fwrite(&(this->n_), sizeof(TextPosition), 1, file) != 1)
        throw std::runtime_error("TextStorage::Save(): file write error (n_).");
    
    if (std::fwrite(&(this->numberOfTexts_), sizeof(TextPosition), 1, file) != 1)
        throw std::runtime_error("TextStorage::Save(): file write error (n_).");

    offsets_->save(file);
}



} // namespace SXSI

