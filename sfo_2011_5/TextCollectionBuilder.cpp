#include "incbwt/rlcsa_builder.h"
#include "TextCollectionBuilder.h"
#include "TCImplementation.h"

namespace SXSI
{

struct TCBuilderRep
{
    unsigned samplerate;
    CSA::RLCSABuilder * sa;

    ulong n;
    // Total number of texts in the collection
    unsigned numberOfTexts;
    // Length of the longest text
    ulong maxTextLength;
    ulong numberOfSamples;
    bool color;
};

/**
 * Init text collection
 *
 */
    TextCollectionBuilder::TextCollectionBuilder(bool color, unsigned samplerate, ulong estimatedInputLength)
    : p_(new struct TCBuilderRep())
{
    p_->n = 0;
    p_->samplerate = samplerate;
    if (samplerate == 0)
        p_->samplerate = TEXTCOLLECTION_DEFAULT_SAMPLERATE;

    p_->numberOfTexts = 0;
    p_->numberOfSamples = 0;
    p_->color = color;
    
    // Current params: 8 bytes, no samples, buffer size n/10 bytes.
    // Buffer size is always at least 15MB:
    if (estimatedInputLength < TEXTCOLLECTION_DEFAULT_INPUT_LENGTH)
        estimatedInputLength = TEXTCOLLECTION_DEFAULT_INPUT_LENGTH;
    p_->sa = new CSA::RLCSABuilder(8, 0, estimatedInputLength/10);
    assert(p_->sa->isOk());
}

TextCollectionBuilder::~TextCollectionBuilder()
{
    delete p_->sa;
    delete p_;
}

void TextCollectionBuilder::InsertText(uchar const * text)
{
    TextCollection::TextPosition m = std::strlen((char *)text) + 1;
    if (m > p_->maxTextLength)
        p_->maxTextLength = m; // Store length of the longest text seen so far.

    if (m > 1)
    {
        p_->n += m;
        p_->numberOfTexts ++;
        p_->numberOfSamples += (m-1)/p_->samplerate;

        p_->sa->insertSequence((char*)text, m-1, 0);
        assert(p_->sa->isOk());
    }
    else
    {
        // FIXME indexing empty texts
        std::cerr << "TextCollectionBuilder::InsertText() error: can not index empty texts!" << std::endl;
        exit(1);
    }
}


TextCollection * TextCollectionBuilder::InitTextCollection(char type)
{
    uchar * bwt = 0;
    CSA::usint length = 0;
    if (p_->numberOfTexts == 0)
    {
        p_->numberOfTexts ++; // Add one empty text
        bwt = new uchar[2];
        bwt[0] = '\0';
        bwt[1] = '\0';
        length = 1;
        p_->maxTextLength = 1;
    }
    else
    {
        bwt = (uchar *)p_->sa->getBWT(length);
        delete p_->sa;
        p_->sa = 0;

        assert(length == p_->n);
    }

    TextCollection *result = new TCImplementation(bwt, (ulong)length, 
           p_->samplerate, p_->numberOfTexts, p_->maxTextLength, p_->numberOfSamples, type, p_->color);
    return result;
}


} // namespace SXSI
