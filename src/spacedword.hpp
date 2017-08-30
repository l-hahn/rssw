#ifndef SPACEDWORD_HPP_
#define SPACEDWORD_HPP_

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include "alphabet.hpp"
#include "pattern.hpp"

class spacedword{
    public:
        spacedword();
        spacedword(std::vector<char> & SeqVec, pattern & Pat, signed short SeqPos, signed SeqNo = -1, signed FamID = -1);
        spacedword(std::vector<unsigned char> & SeqVec, pattern & Pat, signed short SeqPos, signed SeqNo = -1, signed FamID = -1);
        spacedword(int64_t Spw, signed short SeqPos, signed SeqNo = -1, signed FamID = -1);
        

        void push_back(unsigned char Letter);
        std::string to_string() const;

        signed short position() const;
        int64_t bits() const;
        int64_t bits_comp() const;
        unsigned size() const;
        signed sequence() const;
        signed family() const;
        
        void set_position(signed short Pos);
        void set_sequence(signed No);
        void set_family(signed FamID);


        unsigned char operator[](unsigned Idx) const;
        
        bool operator==(const spacedword & Spw) const;
        bool operator!=(const spacedword & Spw) const;
        bool operator<(const spacedword & Spw) const;
        bool operator>(const spacedword & Spw) const;
        
    private:
        template<typename T>
        void create_word(T & Vec, pattern & Pat);

        int64_t _BitSpacedWord;
        signed short _SequencePos;
        signed _SequenceNo;
        signed _FamilyID;
};
#endif
