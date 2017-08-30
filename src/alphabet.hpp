#ifndef ALPHABET_HPP_
#define ALPHABET_HPP_

#include <map>
#include <string>
#include <vector>

namespace alphabet{
    extern std::map<unsigned char, unsigned char> AlphaBit;
    extern std::map<unsigned char, unsigned char> BitAlpha;

    extern std::map<unsigned char, unsigned char> AlphaReduct;
    extern std::map<unsigned char, unsigned char> AlphaBitReduct;
    extern std::map<unsigned char, unsigned char> BitReduct;

    extern std::map<std::string, unsigned char> CodonProtein;
    extern std::map<unsigned short,unsigned char> CodonProteinBit;

    extern std::vector<unsigned char> ProteinAlphabet;
    extern std::vector<unsigned char> DnaAlphabet;
    extern std::vector<unsigned char> ProteinUniqAlphabet;
    extern std::vector<unsigned char> PureDna;
    extern std::vector<unsigned char> ProteinBitAlphabet;
    extern std::vector<unsigned char> ProteinUniqBitAlphabet;
    extern std::vector<unsigned char> DnaBitAlphabet;
    extern std::vector<unsigned char> PureDnaBit;

    extern unsigned BitLength;

    std::map<unsigned char, unsigned char> & AlphaMap(int BitLeng);
    std::map<unsigned char, unsigned char> & BitMap(int BitLeng);
    unsigned char DnaInvert(unsigned char Letter);
};
#endif
