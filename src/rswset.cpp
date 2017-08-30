#include "rswset.hpp"

rsw_set::rsw_set(){
}

rsw_set::rsw_set(pattern & Pattern): _Pattern(Pattern){
}
rsw_set::rsw_set(pattern && Pattern): _Pattern(Pattern){
}

void rsw_set::push_back(rep_spw && Spw){
    push_back(Spw);
}
void rsw_set::push_back(rep_spw & Spw){
    _RepresentWords.push_back(Spw);
}
void rsw_set::merge(rsw_set & RswSet){
    if(RswSet._Pattern != _Pattern){
        return;
    }
	std::sort(_RepresentWords.begin(), _RepresentWords.end());
	std::sort(RswSet._RepresentWords.begin(), RswSet._RepresentWords.end());
	auto IterSrc = RswSet._RepresentWords.begin();
    auto IterSink = _RepresentWords.begin();

    while(IterSink != _RepresentWords.end() && IterSrc != RswSet._RepresentWords.end()){
        if(*IterSink < *IterSrc){
            IterSink++;
        }
        else if(*IterSink == *IterSrc){
            for(auto & FamScore : *IterSrc){
                IterSink->push_back(FamScore);
            }
            IterSrc++;
        }
        else{
            IterSrc++;
        }
    }
    
    rsw_set Add(_Pattern);
    auto & Sink = Add._RepresentWords;
    auto & Src1 = RswSet._RepresentWords;
    auto & Src2 = _RepresentWords;
    std::set_difference(Src1.begin(), Src1.end(), Src2.begin(), Src2.end(), std::back_inserter(Sink));
    Src2.insert(Src2.end(), Sink.begin(), Sink.end());
    std::sort(Src2.begin(), Src2.end());
}
void rsw_set::clear(){
	_RepresentWords.clear();
}
void rsw_set::set_pattern(pattern && Pattern){
	set_pattern(Pattern);	
}
void rsw_set::set_pattern(pattern & Pattern){
	_Pattern = Pattern;
	_RepresentWords.clear();
}
pattern rsw_set::pat() const{
    return _Pattern;
}
size_t rsw_set::size() const{
	return _RepresentWords.size();
}
rep_spw & rsw_set::operator[](size_t Idx){
    return _RepresentWords[Idx];
}
bool rsw_set::operator==(const rsw_set & Rsw) const{
	return _Pattern == Rsw.pat();
}

rsw_set::iterator rsw_set::begin(){
    return _RepresentWords.begin();
}

rsw_set::iterator rsw_set::end(){
    return _RepresentWords.end();
}
void rsw_set::insert(iterator FillEnd, iterator InputBegin, iterator InputEnd){
	_RepresentWords.insert(FillEnd, InputBegin, InputEnd);
}