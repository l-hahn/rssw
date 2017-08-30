#include "representspw.hpp"

rep_spw::rep_spw(int64_t BitSpacedWord):_BitSpacedWord(BitSpacedWord){	
}

rep_spw::rep_spw(int64_t BitSpacedWord, unsigned Family, unsigned Position, double FamilyScore):_BitSpacedWord(BitSpacedWord){
	push_back(Family, Position, FamilyScore);	
}

rep_spw::rep_spw(int64_t BitSpacedWord, family_score & FamilyScore):_BitSpacedWord(BitSpacedWord){
	push_back(FamilyScore);	
}


int64_t rep_spw::bits() const{
	return _BitSpacedWord;
}

size_t rep_spw::size() const{
	return _Family.size();
}

void rep_spw::clear(){
	_Family.clear();
}

void rep_spw::push_back(unsigned Family, unsigned Position, double FamilyScore){
	_Family.push_back(family_score(Family, Position, FamilyScore));
}

void rep_spw::push_back(family_score & FamilyScore){
	_Family.push_back(FamilyScore);
}


family_score rep_spw::operator[](size_t Idx) const{
	return _Family[Idx];
}

family_score & rep_spw::operator[](size_t Idx){
	return _Family[Idx];
}

bool rep_spw::operator==(const rep_spw & Spw) const{
	return _BitSpacedWord == Spw._BitSpacedWord;
}
bool rep_spw::operator!=(const rep_spw & Spw) const{
	return _BitSpacedWord != Spw._BitSpacedWord;
}

bool rep_spw::operator<(const rep_spw & Spw) const{
	return _BitSpacedWord < Spw._BitSpacedWord;
}

bool rep_spw::operator>(const rep_spw & Spw) const{
	return _BitSpacedWord > Spw._BitSpacedWord;
}


rep_spw::iterator rep_spw::begin(){
	return _Family.begin();
}

rep_spw::iterator rep_spw::end(){
	return _Family.end();
}