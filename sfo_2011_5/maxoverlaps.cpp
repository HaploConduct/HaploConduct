// Returns the overlaps with max length-2*k
// Usage: ./maxoverlaps < overlapfile > outputfile
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
using namespace std;

int main(int argc, char* argv[]) 
{
    if (argc != 1)
    {
        printf("usage: %s < overlapfile > outputfile\n", argv[0]);
        return 1;
    }

    int cA = -1;
    int cB = -1;

    // Init buffers
    int length[8];
    char *row[8];
    // Reset buffers
    for (unsigned i = 0; i < 8; ++i)
    {
        length[i] = 0;
        row[i] = 0;
    }
    
    while (!feof(stdin))
    {
        char line[2000];
        if (fgets(line, 1999, stdin) == 0) 
	  break;

	{
	  unsigned l = strlen(line);
	  unsigned tabs = 0;
	  for (unsigned i = 0; i < l; ++i)
	    if (line[i] == '\t')
	      ++ tabs;
	  if (tabs != 7)
	    {
	      cerr << "maxoverlap error: malformed row: " << line ;
	      exit(1);
	    } 
	}
	
        int A, B, OHA, OHB, OLA, OLB, K;
        char O;
        if (sscanf(line, "%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d", &A, &B, &O ,&OHA, &OHB, &OLA, &OLB, &K) != 8)
	{
	  cerr << "maxoverlap error: scanf failed on row: " << line ;
	  exit(1);
	}

        unsigned i = 0; // Default: OHA >= 0 && OHB >= 0
        if (OHA < 0 && OHB < 0) i = 2;
        if (OHA > 0 && OHB < 0) i = 4;
        if (OHA < 0 && OHB > 0) i = 6;
        if (O == 'I')
            i ++;

        if (A == cA && B == cB)
        {
	    if (OLA - (2*K) > length[i])
            {
	        length[i] = OLA - (2*K);
                delete [] row[i];
                row[i] = new char[2000];
                strncpy(row[i], line, 1999);
		row[i][1999] = 0;
            }
        }
        else
        {
            // Flush buffers
            for (unsigned j = 0; j < 8; ++j)
                if (row[j] != 0)
                {
                    printf("%s", row[j]);
                    delete [] row[j];
                    row[j] = 0;                        
                }

            // Reset buffers
            for (unsigned j = 0; j < 8; ++j)
                length[j] = 0;

            cA = A;
            cB = B;
            length[i] = OLA - (2*K);
            row[i] = new char[2000];
            strncpy(row[i], line, 1999);
	    row[i][1999] = 0;
        }
    }

    // Flush buffers
    for (unsigned j = 0; j < 8; ++j)
        if (row[j] != 0)
        {
            printf("%s", row[j]);
            delete [] row[j];
        }


    return 0;
}
