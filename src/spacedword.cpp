#include "spacedword.hpp"
spacedword::spacedword(): _BitSpacedWord(0), _SequencePos(-1), _SequenceNo(-1), _FamilyID(-1){}
spacedword::spacedword(std::vector<char> & SeqVec, pattern & Pat, signed short SeqPos, signed SeqNo, signed FamID): _SequencePos(SeqPos), _SequenceNo(SeqNo), _FamilyID(FamID){
    create_word(SeqVec, Pat);
}
spacedword::spacedword(std::vector<unsigned char> & SeqVec, pattern & Pat, signed short SeqPos, signed SeqNo, signed FamID): _SequencePos(SeqPos), _SequenceNo(SeqNo), _FamilyID(FamID){
    create_word(SeqVec, Pat);
}
spacedword::spacedword(int64_t Spw, signed short SeqPos, signed SeqNo, signed FamID): _BitSpacedWord(Spw), _SequencePos(SeqPos), _SequenceNo(SeqNo), _FamilyID(FamID){

}


template<typename T>
void spacedword::create_word(T & Vec, pattern & Pat){
    auto SeqIt = Vec.begin() + _SequencePos;
    _BitSpacedWord = 0;
    for(auto & PatIt : Pat.match_pos()){
        push_back(*(SeqIt+PatIt));
    }
}

void spacedword::push_back(unsigned char Letter){
    _BitSpacedWord <<= alphabet::BitLength;
    if(Letter > alphabet::ProteinBitAlphabet[alphabet::ProteinBitAlphabet.size()-1]){
        _BitSpacedWord |= alphabet::AlphaBit[Letter];
    }
    else{
        _BitSpacedWord |= Letter;
    }
}

std::string spacedword::to_string() const{
    std::string AlphaWord = "";
    int64_t BitWord = _BitSpacedWord;
    char BitMask = (1 << alphabet::BitLength)-1; 
    while(BitWord > 0 && alphabet::BitLength > 0){
        AlphaWord.push_back(alphabet::BitAlpha[BitWord & BitMask]);
        BitWord >>= alphabet::BitLength;
    }
    std::reverse(AlphaWord.begin(),AlphaWord.end());
    return AlphaWord;
}

unsigned spacedword::size() const{
    return (std::log2(_BitSpacedWord)+1)/alphabet::BitLength;
}
signed short spacedword::position() const{
    return _SequencePos;
}
int64_t spacedword::bits() const{
    return _BitSpacedWord;
}
int64_t spacedword::bits_comp() const{
    int64_t RevCompBits=0, TmpBits = _BitSpacedWord;
    unsigned char BitMask = (1 << alphabet::BitLength)-1;
    while(TmpBits > 0){
        RevCompBits <<= alphabet::BitLength;
        RevCompBits |= alphabet::DnaInvert(TmpBits & BitMask);
        TmpBits >>= alphabet::BitLength;
    }
    while(RevCompBits > 0){
        TmpBits <<= alphabet::BitLength;
        TmpBits |= (RevCompBits & BitMask);
        RevCompBits >>= alphabet::BitLength;
    }
    RevCompBits = TmpBits;
    return RevCompBits;
}
signed spacedword::sequence() const{
    return _SequenceNo;
}
signed spacedword::family() const{
    return _FamilyID;
}

void spacedword::set_position(signed short Pos){
    _SequencePos = Pos;
}
void spacedword::set_sequence(signed No){
    _SequenceNo = No;
}
void spacedword::set_family(signed FamID){
    _FamilyID = FamID;
}


unsigned char spacedword::operator[](unsigned Idx) const{
    char BitMask = (1 << alphabet::BitLength)-1;
    int64_t BitLetter = _BitSpacedWord;
    for(unsigned i = 0; i < size()-Idx; i++){
        BitLetter >>= alphabet::BitLength;
    }
    return alphabet::BitAlpha[BitLetter & BitMask];
}

bool spacedword::operator==(const spacedword & Spw) const{
    return _BitSpacedWord == Spw._BitSpacedWord;
}
bool spacedword::operator!=(const spacedword & Spw) const{
    return _BitSpacedWord != Spw._BitSpacedWord;
}
bool spacedword::operator<(const spacedword & Spw) const{
    return _BitSpacedWord < Spw._BitSpacedWord;
}
bool spacedword::operator>(const spacedword & Spw) const{
    return _BitSpacedWord > Spw._BitSpacedWord;
}
