#include "alphabet.hpp"

namespace alphabet{
    unsigned BitLength = 5;

    std::map<unsigned char, unsigned char> AlphaBit = {{'A',1},{'B',2},{'C',3}
        ,{'D',4},{'E',5},{'F',6},{'G',7},{'H',8},{'I',9},{'J',10},{'K',11},{'L',12}
        ,{'M',13},{'N',14},{'O',15},{'P',16},{'Q',17},{'R',18},{'S',19},{'T',20}
        ,{'U',21},{'V',22},{'W',23},{'X',24},{'Y',25},{'Z',26},{'*',27}};
    std::map<unsigned char, unsigned char> BitAlpha = {{1,'A'},{2,'B'},{3,'C'}
        ,{4,'D'},{5,'E'},{6,'F'},{7,'G'},{8,'H'},{9,'I'},{10,'J'},{11,'K'},{12,'L'}
        ,{13,'M'},{14,'N'},{15,'O'},{16,'P'},{17,'Q'},{18,'R'},{19,'S'},{20,'T'}
        ,{21,'U'},{22,'V'},{23,'W'},{24,'X'},{25,'Y'},{26,'Z'},{27,'*'}};

    std::map<unsigned char, unsigned char> AlphaReduct = {{'A','S'},{'B','K'},{'C','C'}
        ,{'D','K'},{'E','K'},{'F','F'},{'G','G'},{'H','H'},{'I','I'},{'J','I'},{'K','K'}
        ,{'L','I'},{'M','M'},{'N','K'},{'O','O'},{'P','P'},{'Q','K'},{'R','K'},{'S','S'}
        ,{'T','S'},{'U','U'},{'V','I'},{'W','W'},{'X','X'},{'Y','Y'},{'Z','K'},{'*','*'}};
    std::map<unsigned char, unsigned char> AlphaBitReduct = {{'A',19},{'B',11}
        ,{'C',3},{'D',11},{'E',11},{'F',6},{'G',7},{'H',8},{'I',9},{'J',9}
        ,{'K',11},{'L',9},{'M',13},{'N',11},{'O',15},{'P',16},{'Q',11},{'R',11}
        ,{'S',19},{'T',19},{'U',21},{'V',9},{'W',23},{'X',24},{'Y',25},{'Z',11}
        ,{'*',27}};
    std::map<unsigned char, unsigned char> BitReduct = {{1,19},{2,11},{3,3},{4,11}
        ,{5,11},{6,6},{7,7},{8,8},{9,9},{10,9},{11,11},{12,9},{13,13},{14,11}
        ,{15,15},{16,16},{17,11},{18,11},{19,19},{20,19},{21,21},{22,9},{23,23}
        ,{24,24},{25,25},{26,11},{27,27}};

    std::map<std::string, unsigned char> CodonProtein = {{"AAA",'K'},{"AAC",'N'}
        ,{"AAG",'K'},{"AAT",'N'},{"ACA",'T'},{"ACC",'T'},{"ACG",'T'},{"ACT",'T'}
        ,{"AGA",'R'},{"AGC",'S'},{"AGG",'R'},{"AGT",'S'},{"ATA",'I'},{"ATC",'I'}
        ,{"ATG",'M'},{"ATT",'I'},{"CAA",'Q'},{"CAC",'H'},{"CAG",'Q'},{"CAT",'H'}
        ,{"CCA",'P'},{"CCC",'P'},{"CCG",'P'},{"CCT",'P'},{"CGA",'R'},{"CGC",'R'}
        ,{"CGG",'R'},{"CGT",'R'},{"CTA",'L'},{"CTC",'L'},{"CTG",'L'},{"CTT",'L'}
        ,{"GAA",'E'},{"GAC",'D'},{"GAG",'E'},{"GAT",'D'},{"GCA",'A'},{"GCC",'A'}
        ,{"GCG",'A'},{"GCT",'A'},{"GGA",'G'},{"GGC",'G'},{"GGG",'G'},{"GGT",'G'}
        ,{"GTA",'V'},{"GTC",'V'},{"GTG",'V'},{"GTT",'V'},{"TAA",'*'},{"TAC",'Y'}
        ,{"TAG",'*'},{"TAT",'Y'},{"TCA",'S'},{"TCC",'S'},{"TCG",'S'},{"TCT",'S'}
        ,{"TGA",'*'},{"TGC",'C'},{"TGG",'W'},{"TGT",'C'},{"TTA",'L'},{"TTC",'F'}
        ,{"TTG",'L'},{"TTT",'F'}};
    std::map<unsigned short,unsigned char> CodonProteinBit = {{273,11},{275,14}
        ,{277,11},{284,14},{305,20},{307,20},{309,20},{316,20},{337,18},{339,19}
        ,{341,18},{348,19},{449,9},{451,9},{453,13},{460,9},{785,17},{787,8}
        ,{789,17},{796,8},{817,16},{819,16},{821,16},{828,16},{849,18},{851,18}
        ,{853,18},{860,18},{961,12},{963,12},{965,12},{972,12},{1297,5},{1299,4}
        ,{1301,5},{1308,4},{1329,1},{1331,1},{1333,1},{1340,1},{1361,7},{1363,7}
        ,{1365,7},{1372,7},{1473,22},{1475,22},{1477,22},{1484,22},{3089,27}
        ,{3091,25},{3093,27},{3100,25},{3121,19},{3123,19},{3125,19},{3132,19}
        ,{3153,27},{3155,3},{3157,23},{3164,3},{3265,12},{3267,6},{3269,12},{3276,6}};

    std::vector<unsigned char> ProteinAlphabet = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','*'};
    std::vector<unsigned char> DnaAlphabet = {'A','B','C','D','G','H','K','M','N','R','S','T','V','W','Y'};
    std::vector<unsigned char> ProteinUniqAlphabet = {'E','F','I','J','L','O','P','Q','U','X','Z','*'};
    std::vector<unsigned char> PureDna = {'A','C','G','N','T'};
    std::vector<unsigned char> ProteinBitAlphabet = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27};
    std::vector<unsigned char> DnaBitAlphabet = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    std::vector<unsigned char> PureDnaBit = {1,3,5,9,12};
    unsigned char DnaInvert(unsigned char Letter){
        if(Letter < 16){
            switch(Letter){
                case 1:
                    Letter = 12;
                    break;
                case 3:
                    Letter = 5;
                    break;
                case 5:
                    Letter = 3;
                    break;
                case 12:
                    Letter = 1;
                    break;
            }
        }
        else{
            switch(Letter){
                case 'A':
                    Letter = 'T';
                    break;
                case 'C':
                    Letter = 'G';
                    break;
                case 'G':
                    Letter = 'C';
                    break;
                case 'T':
                    Letter = 'A';
                    break;
            }
        }
        return Letter;
    }
};
