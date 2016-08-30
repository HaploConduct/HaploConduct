#ifndef _SXSI_EditDistance_h_
#define _SXSI_EditDistance_h_

/**
 * TODO
 * [ ] clean up code
 */

#include <utility>
#include <cstring>
#include "Tools.h"

class EditDistance 
{
public:
    uchar const *pat;
    unsigned  patlen;
    unsigned maxlen;
    unsigned **dp;
    unsigned p;

    bool end()
    {
        if (p < maxlen)
            return false;
        return true;
    }

    EditDistance(uchar const *pattern, unsigned len, unsigned maxk) {
        pat = pattern;
        patlen = len;
        maxlen = len + maxk;
        p = 0;
        dp = new unsigned*[patlen + 1];
        for (unsigned i = 0; i <= patlen; i++) {
            dp[i] = new unsigned[maxlen+1];
            dp[i][0] = i;
        }
    }

    ~EditDistance() {
        for (unsigned i = 0; i <= patlen; i++)
            delete [] dp[i];
        delete [] dp;
    }

    std::pair<unsigned, unsigned> prefixDist(uchar const *cand, unsigned candlen) {
        unsigned mindist = patlen+candlen;
        unsigned minlen = patlen+candlen;
        unsigned *col = new unsigned[patlen+1];
        unsigned *ncol = new unsigned[patlen+1];
        for (unsigned i = 0; i <= patlen; i++)
            col[i] = i;
        for (unsigned j = 1; j <= candlen; j++) {
            ncol[0] = j;
            for (unsigned i = 1; i <= patlen; i++) {
                unsigned match = (cand[j-1] == pat[i-1])? col[i - 1] : (col[i - 1] + 1);
                ncol[i] = std::min(col[i] + 1, ncol[i - 1] + 1);
                ncol[i] = std::min(ncol[i],match);
            }
            if (ncol[patlen]<mindist)
            {
                mindist = ncol[patlen];
                minlen = j;
            }
            for (unsigned i = 0; i <= patlen; i++)
                col[i] = ncol[i];
        }
        delete [] col;
        delete [] ncol;
        return std::make_pair(mindist,minlen);
    }


    unsigned pushChar(uchar c) {
        unsigned mindist = patlen;
        p++;
        if (p>maxlen)
            std::cout << "out of bounds in DP computation, p = " << p << ", max = " << maxlen << std::endl;
        dp[0][p] = p;
        for (unsigned i = 1; i <= patlen; i++) {
            unsigned match = (c == pat[i-1])? dp[i - 1][p-1] : (dp[i - 1][p-1] + 1);
            dp[i][p] = std::min(dp[i][p-1] + 1, dp[i - 1][p] + 1);
            dp[i][p] = std::min(dp[i][p],match);
            if (dp[i][p]<mindist)
                mindist = dp[i][p];
        }
        return mindist;
    }

    void popChar() {
        --p;
        return;
    }

    unsigned dist() {
        return dp[patlen][p];
    }

    unsigned getSuffixLength()
    {
        unsigned mindist = dp[patlen][p];
        unsigned minlen = patlen;
        for (unsigned i = 0; i < patlen; ++i)
            if (dp[i][p] < mindist)
            {
                mindist = dp[i][p];
                minlen = i;
            }
        return minlen;
    }

};


class MyersEditDistanceIncremental
{
  public:

    MyersEditDistanceIncremental(const uchar *pattern, unsigned _length, unsigned _max_errors) :
      length(_length), max_errors(_max_errors), min_row(~0u)
    {
      this->max_length = this->length + this->max_errors;
      this->blocks = (this->length + W - 1) / W;
      this->last_block_length = this->length % W;
      if(this->last_block_length == 0) { this->last_block_length = W; }
      this->last_block_mask = ~0lu >> (this->blocks * W - this->length);

      this->initializeCharVectors(pattern);
      this->initializePlusMinusVectors();
    }

    ~MyersEditDistanceIncremental()
    {
      for(unsigned c = 0; c < 256; c++)
      {
        delete[] this->char_vectors[c]; this->char_vectors[c] = 0;
      }
      delete[] this->char_vectors; this->char_vectors = 0;

      for(unsigned i = 0; i <= this->max_length; i++)
      {
        delete[] this->plus_vectors[i]; this->plus_vectors[i] = 0;
        delete[] this->minus_vectors[i]; this->minus_vectors[i] = 0;
      }
      delete[] this->plus_vectors; this->plus_vectors = 0;
      delete[] this->minus_vectors; this->minus_vectors = 0;

      delete[] this->score; this->score = 0;
    }

    unsigned pushChar(unsigned c)
    {
      if(c >= 256)
      {
        std::cerr << "MyersEditDistance: Invalid character value " << c << "!" << std::endl;
        return ~0u;
      }
      if(this->pos >= this->max_length)
      {
        std::cerr << "MyersEditDistance: Cannot pushChar(c) past limit!" << std::endl;
        std::cerr << "MyersEditDistance: Current text position " << this-> pos << ", maximum " << this->max_length << std::endl;
        return ~0u;
      }

      this->pos++;
      ulong c_vec, diagonal, h_plus = 0, h_minus = 0, plus, minus, temp;
      ulong d_overflow = 0, hp_overflow = 1, hm_overflow = 0;

      /*
        Update rules from http://www.cs.uta.fi/~helmu/pubs/phd.pdf (Figure 25).
      */
      for(unsigned i = 0; i < this->blocks; i++)
      {
        c_vec = this->char_vectors[c][i];
        plus = this->plus_vectors[this->pos - 1][i];
        minus = this->minus_vectors[this->pos - 1][i];

        diagonal = c_vec & plus;
        temp = diagonal + plus + d_overflow;
        d_overflow = ((temp < diagonal) || (temp == diagonal && d_overflow == 1) ? 1 : 0);
        diagonal = (temp ^ plus) | c_vec | minus;

        h_plus = minus | ~(diagonal | plus);
        h_minus = diagonal & plus;

        this->plus_vectors[this->pos][i] =
          (h_minus << 1) | hm_overflow | ~(diagonal | (h_plus << 1) | hp_overflow);
        this->minus_vectors[this->pos][i] = diagonal & ((h_plus << 1) | hp_overflow);

        hp_overflow = (h_plus & HIGH_BIT ? 1 : 0);
        hm_overflow = (h_minus & HIGH_BIT ? 1 : 0);
      }
      this->plus_vectors[this->pos][this->blocks - 1] &= this->last_block_mask;
      this->minus_vectors[this->pos][this->blocks - 1] &= this->last_block_mask;

      this->score[this->pos] = this->score[this->pos - 1];
      if(h_plus & (1lu << (this->length % W - 1)))
      {
        this->score[this->pos]++;
      }
      else if(h_minus & (1lu << (this->length % W - 1)))
      {
        this->score[this->pos]--;
      }

      // Find the minimal value in the current column.
      int minimum = this->score[this->pos], cur = minimum;
      min_row = this->length;
      for(unsigned i = this->blocks; i > 0;)
      {
        i--;
        plus = this->plus_vectors[this->pos][i];
        minus = this->minus_vectors[this->pos][i];

        for(unsigned j = W; j > 0;)
        {
          j -= CHAR_BIT;
          unsigned a = (plus >> j) & 0xFF, b = (minus >> j) & 0xFF;
          if (minimum > cur + myers_optimum[256 * a + b])
          {
              minimum = cur + myers_optimum[256 * a + b];
              min_row = i*W + j + myers_best_pos[256 * a + b];
          }
          cur += myers_total[256 * a + b];
        }
      }
      min_row_value = minimum;
      return minimum;
    }

    unsigned dist() {
        return this->score[this->pos];
    }
    
    void popChar()
    {
      if(this->pos == 0)
      {
        std::cerr << "MyersEditDistance: Cannot popChar() from text position 0!" << std::endl;
        return;
      }

      this->pos--;
    }

    inline unsigned getSuffixLength()
    {
        return min_row;
    }
    inline unsigned getMinimum()
    {
        return min_row_value;
    }

/*
    // FIXME implement
    std::pair<unsigned, unsigned> prefixDist(uchar const *text, unsigned n)
*/

    bool end()
    {
        if (pos < max_length)
            return false;
        return true;
    }

static void initMyersFourRussians()
{
//  myers_total = new signed char[256 * 256];
//  myers_optimum = new signed char[256 * 256];
//   myers_best_pos = new uchar[256 * 256];

  unsigned bp = 8;
  for(unsigned i = 0; i < 256; i++)
  {
    for(unsigned j = 0; j < 256; j++)
    {
      int cur = 0, best = 0;
      for(unsigned k = CHAR_BIT; k > 0;)
      {
        k--;
        if(i & (1 << k)) { cur--; }
        if(j & (1 << k)) { cur++; }
        if(cur < best)
        {
          best = cur;
          bp = k;
        }
      }

      myers_total[256 * i + j] = cur;
      myers_optimum[256 * i + j] = best;
      myers_best_pos[256 * i + j] = bp;
    }
  }
}

    unsigned length;
private:
    static signed char myers_total[];
    static signed char myers_optimum[];
    static uchar       myers_best_pos[];


    const static ulong HIGH_BIT = 1lu << (W - 1);

    ulong **char_vectors;
    unsigned max_errors;
    unsigned max_length;
    unsigned min_row;
    unsigned min_row_value;

    unsigned blocks;            // Number of blocks required for the pattern
    unsigned last_block_length; // Number of bits in the final block
    ulong last_block_mask;

    ulong **plus_vectors;
    ulong **minus_vectors;
    unsigned* score;

    unsigned pos;               // We have matched pos characters of text


    void initializeCharVectors(uchar const *pattern)
    {
      this->char_vectors = new ulong*[256];
      for(unsigned c = 0; c < 256; c++)
      {
        this->char_vectors[c] = new ulong[this->blocks];
        std::memset(this->char_vectors[c], 0, this->blocks * sizeof(ulong));
      }

      for(unsigned i = 0; i < this->length; i++)
      {
        this->char_vectors[pattern[i]][i / W] |= (1lu << (i % W));
      }
    }

    void initializePlusMinusVectors()
    {
      this->plus_vectors = new ulong*[this->max_length + 1];
      this->minus_vectors = new ulong*[this->max_length + 1];
      this->score = new unsigned[this->max_length + 1];

      for(unsigned i = 0; i <= this->max_length; i++)
      {
        this->plus_vectors[i] = new ulong[this->blocks];
        this->minus_vectors[i] = new ulong[this->blocks];
      }

      this->initialStep();
    }

    void initialStep()
    {
      for(unsigned i = 0; i < this->blocks; i++)
      {
        this->plus_vectors[0][i] = ~0lu;
        this->minus_vectors[0][i] = 0;
      }
      this->plus_vectors[0][this->blocks - 1] = this->last_block_mask;

      this->score[0] = this->length;
      this->pos = 0;
    }



};




    class MyersEditDistance
    {
    public:
        ulong **Peq;
        bool S[256];
        unsigned k;
        unsigned m;
        unsigned blockCount;
        unsigned block;
        unsigned lastBlockLength;
        ulong *Pv;
        ulong *Mv;
        ulong *Score;
        ulong twopowmpos;
        ulong twopowwpos;
 
        void init(uchar const *P)
        {
            ulong i;
            ulong j;
            ulong l;
            for(i = 0;i< 256;i++) {
                S[i] = false;
                for (l=0;l <= ((m-1) / W); l++)
                    Peq[l][i] = 0;
            }
            for (l=0;l <= ((m-1) / W); l++)
                for (j = 1; j <= W; j++) {
                    if (W*l+j > m)
                        break;
                    Peq[l][P[W*l+j-1]] = Peq[l][P[W*l+j-1]] | (1lu << (j-1));
                    if (j+W*l <= k+1)
                        S[(unsigned)P[W*l+j-1]] = true;
                }
            for (l=0; l < blockCount; l++) {
                Mv[l] = 0lu;
                Score[l] = (l+1lu)*W;
                Pv[l] = ~0lu;
            }
            Mv[blockCount] = 0;
            Score[blockCount] = m;
            Pv[blockCount] = (~0lu) >> ((blockCount+1lu)*W - m);
        }

        MyersEditDistance(unsigned _m, unsigned maxk) 
            : Peq(0), k(maxk), m(_m)
        {
            ulong l;
            Peq = new ulong *[(m-1) / W + 1]; // (unsigned **)malloc((((m-1) / w)+1)*sizeof(unsigned*));
            for (l=0;l <= ((m-1) / W); l++)
                Peq[l] = new ulong[256]; //(unsigned *)malloc(256*sizeof(unsigned));

            blockCount= ((m-1) / W);
            block = (k-1) / W;
            lastBlockLength = m % W;
            if (lastBlockLength == 0)
                lastBlockLength = W;
            Pv = new ulong[blockCount+1];  // (unsigned *)malloc((blockCount+1)*sizeof(unsigned));
            Mv = new ulong[blockCount+1];  // (unsigned *)malloc((blockCount+1)*sizeof(unsigned));
            Score = new ulong[blockCount+1];//(unsigned *)malloc((blockCount+1)*sizeof(unsigned));
            twopowmpos = 1lu << (lastBlockLength-1);
            twopowwpos = 1lu << (W-1);
        }

        ~MyersEditDistance()
        {
            delete [] Pv; Pv = 0;
            delete [] Mv; Mv = 0;
            delete [] Score; Score = 0;
            unsigned l;
            for (l=0;l <= ((m-1) / W); l++)
            {
                delete [] Peq[l];
                Peq[l] = 0;
            }
            delete [] Peq; Peq = 0;
        }

        std::pair<unsigned, unsigned> prefixDist(uchar const *text, unsigned n)
        {
            unsigned mindist = m+n;
            unsigned minlen = m+n;
            //            int count=0;
            //int blockAverage=0;
            ulong Eq, Xv, Xh, Ph, Mh;
            ulong Phendbit;
            ulong Mhendbit;
            ulong temp;
            ulong j;
            ulong l;
            ulong overk=1;
            block = (k-1) / W;
            for(j = 1; j <= n; j++) {
                Phendbit = 1lu; //0;
                Mhendbit = 0lu;
                for (l=0; l <= blockCount; l++) {
                    Eq = Peq[l][(unsigned char)text[j-1]];
                    Xv = Eq | Mv[l];
                    Xh = ((((Eq | Mhendbit)& Pv[l]) + Pv[l]) ^ Pv[l]) | (Eq | Mhendbit);
                    Ph = Mv[l] | ~ (Xh | Pv[l]);
                    Mh = Pv[l] & Xh;
                    temp = l < blockCount ? twopowwpos : twopowmpos;
                    if (Ph & temp)
                        Score[l] += 1;
                    else if (Mh & temp)
                        Score[l] -= 1;
                    temp = (Ph >> (W-1));
                    Ph <<= 1;
                    Ph |= Phendbit;
                    Phendbit = temp;
                    temp = (Mh >> (W-1));
                    Mh <<= 1;
                    Mh |= Mhendbit;
                    Mhendbit = temp;
                    Pv[l] = (Mh) | ~ (Xv | Ph);
                    Mv[l] = Ph & Xv;
                    if (block == l+1 &&
                        ((block == blockCount &&
                          (Score[l] > k+lastBlockLength ||
                           (Score[l] > k && overk != 0 && j - overk >= lastBlockLength))) ||
                         (block != blockCount &&
                          (Score[l] > k+W ||
                           (Score[l] > k && overk != 0 && j - overk >= W)))))  {
                        // lohkon Score kasvanut niin isoksi, ettei seuraavaa kannata laskea
                        overk = 0;
                        block = l;
                        break;
                    }
                    else if (block == l+1 && overk == 0 && Score[l] > k)
                        // Talletetaan ensimm‰inen diagonaali jolla Score > k. N‰in tiedet‰‰n
                        // jatkossa koska seuraavan lohkon laskeminen on turhaa.
                        overk = j;
                    else if (block == l+1 && Score[l] <= k)
                        // Diagonaalilla Score <= k ==> Nollataan muuttuja overk.
                        overk = 0;
                    else if (block == l && block != blockCount && Score[l] <= k) {
                        // Score <= k ==> jatketaan seuraavaan lohkoon. Koska seuraavaa lohkoa
                        // ei ole laskettu edellisell‰ sarakkeella, pit‰‰ sen arvot p‰ivitt‰‰.
                        Score[l+1] = k+1lu;
                        Pv[l+1] = 0lu;
                        Mv[l+1] = 0lu;
                        overk = 0lu;
                        block = l+1;
                    }
                    else if (block == l)
                        // Seuraavaan lohkoon ei kannata siirty‰, kun Score > k.
                        break;
                }
//                blockAverage += block+1;
                // if (block == blockCount && Score[blockCount] <= k) {
                if (block == blockCount && Score[blockCount] <= mindist)
                {
                    mindist = Score[blockCount];
                    minlen = j;
                
                    //printf("Min %u at len %u\n", mindist, minlen);
                //    count++;
                }
            }
//            return count;
            return std::make_pair(mindist,minlen);
            /*ShowMessage("Esiintymi‰ lˆytyi " + IntToStr(laskuri) + "kpl. " +
              "Aikaa meni " + IntToStr(msec) + "msec");
              ShowMessage("Keskim‰‰r‰inen lohkojen m‰‰r‰ " +
              FormatFloat("0.00",(float)blockAverage / n));*/
        }

};

#endif
