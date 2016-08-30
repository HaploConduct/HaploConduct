/******************************************************************************
 *   Copyright (C) 2009-2011 Niko Välimäki (orig. Susana Ladra)               *
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

#include "TextCollectionBuilder.h"
#include "EditDistanceOverlap.h"
#include "MismatchOverlap.h"
#include "SuffixFilterOverlap.h"
#include "SuffixFilterOverlapEd.h"
#include "Tools.h"
#include "EditDistance.h"

#ifdef PARALLEL_SUPPORT
#include <omp.h>
#endif

#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <getopt.h>
using SXSI::TextCollection;
using SXSI::TCImplementation;
using SXSI::TextCollectionBuilder;
using SXSI::TextStorage;
using namespace std;

/**
 * Definitions for parsing command line options
 */
enum parameter_t { long_opt_indels = 256, long_opt_mismatch, long_opt_skip, long_opt_nreads };

void print_usage(char const *name)
{
    cerr << "usage: " << name << " [options] <input>" << endl
         << "Check README or `" << name << " --help' for more information." << endl;
}

void print_help(char const *name)
{
    cerr << "./sfoverlap [options] <index>" << endl
         << endl
         << "Input file:" << endl
         << "  <index>   Index filename, see `builder --help' for more information" << endl
         << "            about constructing indexes. By default, all reads in the index are" << endl
         << "            matched against all other reads in the index. Options --skip and" << endl
         << "            --nreads can be used to define the subset of reads to be searched." << endl
         << endl
         << "Basic alignment modes:" << endl
         << " -e <int>                    Uses error-rate 1/<int>, e.g. -e20 equals the" << endl
         << "                             error-rate of 1/20 = 0.05. The number of errors" << endl
         << "                             allowed in the alignment depens on the length of" << endl
         << "                             the overlap, that is, overlap length divided by " << endl
         << "                             <int> (round up). This is the recommended mode." << endl
         << endl
         << " -k <int>                    Uses a fixed number of errors for all overlap" << endl
         << "                             lengths. Usable only for small values of <int>" << endl
         << "                             (less than 4)." << endl
         << endl
         << "Alignment options (default is --indels):" << endl
         << " --indels                    Allow mismatches, insertions and deletions in" << endl
         << "                             the aligment." << endl
         << " --mismatch                  Allow only mismatches in the alignment." << endl
         << endl
         << "General options:" << endl
         << " -t <int>, --treshold <int>  Minimum overlap length threshold. Only overlaps" << endl
         << "                             longer than or equal to <int> are outputted" << endl
         << "                             (default is 40)." << endl
         << " --skip <int>                Skip first <int> reads the set (default is 0)." << endl
         << " --nreads <int>              Align <int> reads from the set (after skipping)" << endl
         << "                             (default is all reads)." << endl
         << " -c, --color                 Reads are in SOLiD color codes." << endl
         << " -v, --verbose               Verbose mode, uses the standard error output." << endl
         << " -h, --help                  Display command line options." << endl;
#ifdef PARALLEL_SUPPORT
    cerr << " -P <int>, --parallel <int>  Number of parallel threads to use (default is one," << endl
         << "                             give argument 0 to use all available cores)." << endl;
#else
    cerr << " -P <int>, --parallel <int>  Parallel computation requires OpenMP, see README." << endl;
#endif
}


int atoi_min(char const *value, int min, char const *parameter, char const *name)
{
    std::istringstream iss(value);
    int i;
    char c;
    if (!(iss >> i) || iss.get(c))
    {
        cerr << name << ": argument of " << parameter << " must be of type <int>, and greater than or equal to " << min << endl
             << "Check README or `" << name << " --help' for more information." << endl;
        std::exit(1);
    }

    if (i < min)
    {
        cerr << name << ": argument of " << parameter << " must be greater than or equal to " << min << endl
             << "Check README or `" << name << " --help' for more information." << endl;
        std::exit(1);
    }
    return i;
}

int main(int argc, char** argv)
{
    MyersEditDistanceIncremental::initMyersFourRussians();

    /**
     * Parse command line parameters
     */
    if (argc <= 1)
    {
        print_usage(argv[0]);
        return 0;
    }

#ifdef PARALLEL_SUPPORT
    unsigned parallel = 1;
#endif
    bool mismatch = false;
    bool color = false;
    bool suffixfilters = true;
    bool verbose = false;
    unsigned threshold = 40, maxerrors = ~0u;
    unsigned skipReads = 0, nReads = 0;

    static struct option long_options[] =
        {
            {"indels",    no_argument,       0, long_opt_indels},
            {"mismatch",  no_argument,       0, long_opt_mismatch},
            {"threshold", required_argument, 0, 't'},
            {"skip",      required_argument, 0, long_opt_skip},
            {"nreads",    required_argument, 0, long_opt_nreads},
            {"color",     no_argument,       0, 'c'},
            {"verbose",   no_argument,       0, 'v'},
            {"help",      no_argument,       0, 'h'},
            {"parallel",  required_argument, 0, 'P'},
            {0, 0, 0, 0}
        };
    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "e:k:t:cvhP:",
                            long_options, &option_index)) != -1) 
    {
        switch(c) 
        {
        case 'e':
            if (maxerrors != ~0u) { cerr << "error: options -e and -k are mutually exclusive." << endl; std::abort(); }
            suffixfilters = true;
            maxerrors = atoi_min(optarg, 10, "-e", argv[0]);
            break;
        case 'k':
            if (maxerrors != ~0u) { cerr << "error: options -e and -k are mutually exclusive." << endl; std::abort(); }
            suffixfilters = false;
            maxerrors = atoi_min(optarg, 0, "-k", argv[0]);
            if (maxerrors > 3) { cerr << "Warning: High values of -k <int> can be slooow!" << endl; }
            break;
        case long_opt_indels:
            if (mismatch) { cerr << "error: options --indels and --mismatch are mutually exclusive." << endl; std::abort(); }
            mismatch = false;
            break;
        case long_opt_mismatch:
            mismatch = true; break;
        case 't':
            threshold = atoi_min(optarg, 10, "-t, --threshold", argv[0]); break;
        case long_opt_skip:
            skipReads = atoi_min(optarg, 0, "--skip", argv[0]); break;
        case long_opt_nreads:
            nReads = atoi_min(optarg, 1, "--nreads", argv[0]); break;
        case 'c':
            color = true; break;
        case 'v':
            verbose = true; break;
        case 'h':
            print_help(argv[0]);
            return 0;
        case 'P':
#ifdef PARALLEL_SUPPORT
            parallel = atoi_min(optarg, 0, "-P, --parallel", argv[0]); break;
#else
            cerr << argv[0] << ": Parallel processing not currently available!" << endl 
                 << "Please recompile with parallel support; see README for more information." << endl;
            return 1;
#endif
        case '?':
            print_usage(argv[0]);
            return 1;
        default:
            print_usage(argv[0]);
            std::abort ();
        }
    }

    if (maxerrors == ~0u) 
    { 
        cerr << "error: Give either the option -e <int>  or  -k <int>." << endl; 
        std::abort(); 
    }

    // Parse filenames
    if (argc - optind != 1)
    {
        cerr << argv[0] << ": missing (or too many) input files! Expecting one filename." << endl;
        print_usage(argv[0]);
        return 1;
    }
    string index = string(argv[optind++]);

    /**
     * Load indexes
     */
    if (verbose)
    {
        cerr << "Loading forward index..." << endl;
    }
    FILE *fp = fopen((index + TextCollection::FMINDEX_EXTENSION).c_str(), "rb");
    if (!fp)
    {
        cerr << argv[0] << ": unable to read input file " << index + TextCollection::FMINDEX_EXTENSION << endl;
        exit(1); 
    }
    TCImplementation* tcOrig = new TCImplementation(fp, 0); // TextCollection::Load(fileOrig, 3);
    fclose(fp);

    if (verbose)
    {
        cerr << "Loading reverse index..." << endl;
    } 
    fp = fopen((index + TextCollection::REVERSE_EXTENSION + TextCollection::FMINDEX_EXTENSION).c_str(), "rb");
    if (!fp)
    {
        cerr << argv[0] << ": unable to read input file " << index + TextCollection::REVERSE_EXTENSION + TextCollection::FMINDEX_EXTENSION << endl;
        exit(1); 
    }
    TCImplementation* tcCompl = new TCImplementation(fp, 0); //TextCollection::Load(fileCompl, 3);
    fclose(fp);
    fp = 0;

    if (tcOrig->isColorCoded() && !color)
    {
        cerr << argv[0] << ": index " << index << " is color coded, please use option -c, --color." << endl;
        return 1;
    }
    if (!tcOrig->isColorCoded() && color)
    {
        cerr <<  argv[0] << ": option -c, --color was given, but the index " << index
             << " is not color coded!" << endl 
             << "Please rebuild the index using -c, --color option, see README for more information." << endl;
        return 1;
    }

    Tools::StartTimer();
    cerr << std::fixed;
    cerr.precision(3);

    unsigned numberOfLines = nReads + skipReads;
    if (nReads == 0 || numberOfLines > tcOrig->GetNumberOfTexts())
        numberOfLines = tcOrig->GetNumberOfTexts();
    
    if (verbose)
    { // Debug print
      	// Find the largest psize we are allowed to choose
        unsigned klimit = threshold%maxerrors == 0 ? threshold/maxerrors : threshold/maxerrors + 1;
        unsigned psize = threshold%(klimit+1) == 0 ? threshold / (klimit + 1) : threshold / (klimit + 1) + 1;
        for (unsigned i = threshold + 1; i < 5000; ++i)
        {
            unsigned tmp = i%maxerrors == 0 ? i/maxerrors : i/maxerrors + 1;
            if (psize > (i%(tmp+1) == 0 ? i / (tmp + 1) : i / (tmp + 1) + 1))
                psize = i%(tmp+1) == 0 ? i / (tmp + 1) : i / (tmp + 1) + 1;
        }
        
        cerr << "parameters: ";

        if (mismatch)
            cerr << "mismatches only";
        else
            cerr << "edit distance";

        if (suffixfilters)
            cerr << " with suffix filters, t = " << threshold << ", e = " << maxerrors << " (" << 1.0/(double)maxerrors << ")"
                 << ", skip = " << skipReads << ", upto = " << numberOfLines << ", psize = " << psize << endl;
        else
            cerr << " with simple backtracking, t = " << threshold << ", k = " << maxerrors
                 << ", skip = " << skipReads << ", upto = " << numberOfLines << endl;
    }

#ifdef PARALLEL_SUPPORT
    if (parallel != 0)
    {
        if (verbose && parallel == 1) 
            cerr << "Using only one core (default setting, see option -P, --parallel)" << endl;
        if (verbose && parallel != 1)
            cerr << "Using " << parallel << " cores." << endl;
        omp_set_num_threads(parallel);
    }
    else
        if (verbose) cerr << "Using all " << omp_get_max_threads() << " cores available." << endl;
#pragma omp parallel 
#endif    
    {
    SXSI::OverlapQuery *oq = 0;
    if (suffixfilters)
    {
        if (mismatch)
            oq = new SXSI::SuffixFilterOverlap(threshold, maxerrors);
        else
            oq = new SXSI::SuffixFilterOverlapEd(threshold, maxerrors);
    }
    else
    {
        if (mismatch)
            oq = new SXSI::MismatchOverlap(threshold, maxerrors);
        else
            oq = new SXSI::EditDistanceOverlap(threshold, maxerrors);
    }

    if (oq == 0)
    {
        cerr << "error: assert failure!" << endl;
        exit(1);
    }
    
#ifdef PARALLEL_SUPPORT
#pragma omp for schedule (dynamic)
#endif
        for(unsigned i = skipReads; i < numberOfLines; ++i)
        {
            uchar* data;
            SXSI::OverlapQuery::full_dist_result fr;
            SXSI::OverlapQuery::full_dist_result::iterator it;
      
            unsigned idA, idB;
            int OHA, OHB, OLA, OLB;
            unsigned lenData, lenA, lenB;
            char O;

	    data = (uchar *) tcOrig->GetText(i);
	    //cerr << "Pattern: " << data << endl;

            lenData = strlen((char*)data);
            //Normal orientation
            O='N';
	    oq->setInsideAlignments(true);
            fr = oq->SuffixPrefixMatch(tcOrig, data);
	    
	    std::sort(fr.begin(), fr.end(), SXSI::OverlapQuery::FullDistResultComp());
	    
	    // Have to output the results of this idA in subsequent manner	    	
#ifdef PARALLEL_SUPPORT
#pragma omp critical (STDOUT)
#endif
            for (it = fr.begin(); it != fr.end(); ++it)
            {
                if(i<=(*it).idB)
                {
                    idA = i;
                    idB = (*it).idB;
                    OLA = (*it).OLA;
                    OLB = (*it).OLB;
    		    
                    lenA = lenData;
                    lenB = strlen((char *)tcOrig->GetText(idB));
                    OHA = lenA-OLA-(*it).offsetB; //((*it).second).first;
                    OHB = lenB-OLB-(*it).offsetB; //((*it).second).first; // OK

                    //idA idB O (N/I) OHA OHB OLA OLB
                }
                else
                {
                    idA = (*it).idB;
                    idB = i;    				
                    OLA = (*it).OLB;
                    OLB = (*it).OLA;

                    lenA = strlen((char *)tcOrig->GetText(idA));
                    lenB = lenData;

                    OHA = -(lenB-OLB-(*it).offsetB); //((*it).second).first);
                    OHB = -(lenA-OLA-(*it).offsetB); //((*it).second).first);
                }

                if(!((idA==idB)&&(OHA==0||OHB==0)))
                {
                    printf("%d\t%d\t%c\t%d\t%d\t%d\t%d\t%u\n", idA,idB,O ,OHA,OHB,OLA,OLB, (*it).k);
                }                
  				
            }
            		
            //Inverted orientation
            O='I';
            fr = oq->SuffixPrefixMatch(tcCompl, data);
	    
	    std::sort(fr.begin(), fr.end(), SXSI::OverlapQuery::FullDistResultComp());

#ifdef PARALLEL_SUPPORT
#pragma omp critical (STDOUT)
#endif
            for (it = fr.begin(); it != fr.end(); ++it)
            {
                if(i<=(*it).idB)
                {
                    idA = i;
                    idB = (*it).idB;
                    OLA = (*it).OLA; 
                    OLB = (*it).OLB; 
    				
                    lenA = lenData;
                    lenB =strlen((char *)tcCompl->GetText(idB)); 
    				
                    OHA = (lenA-OLA)-(*it).offsetB; //((*it).second).first;
                    OHB = (lenB-OLB)-(*it).offsetB; //((*it).second).first;
                    
                    if(!((idA==idB)&&(OHA==0||OHB==0)))
                    {
                        printf("%d\t%d\t%c\t%d\t%d\t%d\t%d\t%u\n", idA,idB,O ,OHA,OHB,OLA,OLB, (*it).k);
                    }
                    //idA idB O (N/I) OHA OHB OLA OLB
                }
    			
                if((i>(*it).idB)&&((*it).offsetB>0))
                {
                    idA = (*it).idB;
                    idB = i;
                    OLA = (*it).OLB; 
                    OLB = (*it).OLA; 
                    
                    lenA = strlen((char *)tcOrig->GetText(idA));
                    lenB = lenData;
    				
                    OHA = (lenA-OLA-(*it).offsetB); //((*it).second).first); // originally -()
                    OHB = (lenB-OLB-(*it).offsetB); //((*it).second).first); // originally -()

                    if(!((idA==idB)&&(OHA==0||OHB==0)))
                    {
                        printf("%d\t%d\t%c\t%d\t%d\t%d\t%d\t%u\n", idA,idB,O ,OHA,OHB,OLA,OLB, (*it).k);
                    }
                }
            }
            //also I orientation
	    // Here we do not output reverse->normal alignments for non-overlaps (since it is equal to normal->reverse)
	    oq->setInsideAlignments(false); 
            data = (uchar *) tcCompl->GetText(i);
            fr = oq->SuffixPrefixMatch(tcOrig, data);

	    std::sort(fr.begin(), fr.end(), SXSI::OverlapQuery::FullDistResultComp());
	    	
#ifdef PARALLEL_SUPPORT
#pragma omp critical (STDOUT)
#endif
            for (it = fr.begin(); it != fr.end(); ++it)
            {
                OLA = (*it).OLA;
                OLB = (*it).OLB;
                if((i<=(*it).idB)&&((*it).offsetB==0)){
                    idA = i;
                    idB = (*it).idB;
    				
                    lenA = lenData;
                    lenB = strlen((char *)tcCompl->GetText(idB)); 
                    OHA = -(lenB-OLB);
                    OHB = -(lenA-OLA);

                    //idA idB O (N/I) OHA OHB OLA OLB
                    if(!((idA==idB)&&(OHA==0||OHB==0)))
                    {
                        printf("%d\t%d\t%c\t%d\t%d\t%d\t%d\t%u\n", idA,idB,O ,OHA,OHB,OLA,OLB, (*it).k);
                    }
                }
            }

            if (verbose && i-skipReads > 0 && i % 10000 == 0)  
            {
#ifdef PARALLEL_SUPPORT
#pragma omp critical (PRINTSTDERR)
#endif
                fprintf(stderr, "Compared %d strings, time %.3f s/pattern, total %.1f s (%.2f h)\n", i-skipReads, 
                        Tools::GetTime()/(double)(i-skipReads), 
                        Tools::GetTime(), Tools::GetTime()/3600.0);
            }
        }
            
        delete oq;

    } // End parallel

    if (verbose)
        fprintf(stderr, "Finished %d strings, time %.3f s/pattern, total %.1f s (%.2f h)\n", 
                (numberOfLines-skipReads),
                Tools::GetTime()/(double)(numberOfLines-skipReads), 
                Tools::GetTime(), Tools::GetTime()/3600.0);

    delete tcOrig;
    delete tcCompl;
}
