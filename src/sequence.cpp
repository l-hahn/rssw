#include "sequence.hpp"

sequence::sequence(){
}

sequence::sequence(std::string & Seq, std::string SeqName, bool AlphaReduct): _SequenceName(SeqName){
    set_sequence(Seq, AlphaReduct);
}

sequence::sequence(std::string && Seq, std::string SeqName, bool AlphaReduct): _SequenceName(SeqName){
    set_sequence(Seq, AlphaReduct);
}

sequence::sequence(std::vector<char> & Seq, std::string SeqName, bool AlphaReduct): _SequenceName(SeqName){
    set_sequence(Seq, AlphaReduct);
}

sequence::sequence(std::vector<char> && Seq, std::string SeqName, bool AlphaReduct): _SequenceName(SeqName){
    set_sequence(Seq, AlphaReduct);
}

sequence::sequence(std::vector<unsigned char> & Seq, std::string SeqName, bool AlphaReduct): _SequenceName(SeqName){
    set_sequence(Seq, AlphaReduct);
}

sequence::sequence(std::vector<unsigned char> && Seq, std::string SeqName, bool AlphaReduct): _SequenceName(SeqName){
    set_sequence(Seq, AlphaReduct);
}

std::vector<spacedword> sequence::spaced_words(pattern & Pat, signed SeqNo, signed short FamID){
    std::vector<spacedword> Spw;
    for(signed short i = 0; i < (signed short)_Sequence.size() - (signed short)Pat.length() + 1; i++){
            Spw.push_back(spacedword(_Sequence, Pat, i, SeqNo, FamID));
    }
    return Spw;
}

void sequence::push_back(unsigned char Letter){
    if(std::isalpha(Letter)){
        _Sequence.push_back(std::toupper(Letter));
    }
    else if(Letter <= alphabet::ProteinBitAlphabet[alphabet::ProteinBitAlphabet.size()-1]){
        _Sequence.push_back(Letter);
    }
    else if(Letter == '*'){
        _Sequence.push_back(Letter);
    }
}

void sequence::clear(){
    _Sequence.clear();
    _SequenceName = "";   
}

void sequence::to_bit(){
    for(auto & Letter : _Sequence){
        if(Letter >= alphabet::ProteinBitAlphabet[alphabet::ProteinBitAlphabet.size()-1]){
            Letter = alphabet::AlphaBit[Letter];
        }
    }
}

std::string sequence::to_string() const{
    if(is_bit()){
        std::string StringSeq;
        for(auto & Letter : _Sequence){
            StringSeq.push_back(alphabet::BitAlpha[Letter]);
        }
        return StringSeq;
    }
    return std::string(_Sequence.begin(), _Sequence.end());
}

bool sequence::is_bit() const{
    if(_Sequence[0] <= alphabet::ProteinBitAlphabet[alphabet::ProteinBitAlphabet.size()-1]){
        return true;
    }
    return false;
}

bool sequence::check_protein(){
    size_t Ctr = 0;
    std::vector<unsigned char> AlphaVec;
    if(is_bit()){    
        AlphaVec = alphabet::PureDnaBit;
    }
    else{
        AlphaVec = alphabet::PureDna;
    }
    for(auto & Letter : _Sequence){
        if(std::find(AlphaVec.begin(),AlphaVec.end(), Letter) != AlphaVec.end()){
            Ctr++;
        }
    }
    if(Ctr < ((_Sequence.size()/4)*3)){
        _IsProtein = true;
        return true;
    }
    _IsProtein = false;
    return false;
}

bool sequence::is_dna() const{
    return _IsProtein == false;
}

bool sequence::is_protein() const{
    return _IsProtein;
}

std::string sequence::name() const{
    return _SequenceName;
}


unsigned char sequence::operator[](size_t Idx) const{
    return _Sequence[Idx];
}

bool sequence::operator<(const sequence & Seq) const{
    return _Sequence.size() < Seq._Sequence.size();
}

bool sequence::operator>(const sequence & Seq) const{
    return _Sequence.size() < Seq._Sequence.size();
}

sequence::iterator sequence::begin(){
    return _Sequence.begin();
}

sequence::iterator sequence::end(){
    return _Sequence.end();
}


size_t sequence::size() const{
    return _Sequence.size();
}

void sequence::set_name(std::string & Name){
    _SequenceName = Name;
}
