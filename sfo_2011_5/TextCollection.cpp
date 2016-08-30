#include "TextCollection.h"
#include "TCImplementation.h"


namespace SXSI
{

const string TextCollection::REVERSE_EXTENSION = ".reverse";
const string TextCollection::FMINDEX_EXTENSION = ".fmi";

const char TextCollection::ALPHABET_DNA[] = {'A', 'C', 'G', 'T', 'N'};
const char TextCollection::ALPHABET_SOLID[] = {'0', '1', '2', '3', '.'};
char const *TextCollection::ALPHABET = ALPHABET_DNA;

    /**
     * Init text collection from a file
     *
     * See TCImplementation.h for more details.
     */
    TextCollection * TextCollection::Load(FILE *fp, unsigned samplerate)
    {
        TextCollection *result = new TCImplementation(fp, samplerate);
        return result;
    }
}
