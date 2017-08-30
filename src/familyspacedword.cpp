#include "familyspacedword.hpp"

family_spacedword::family_spacedword(std::vector<spacedword>::iterator Begin, std::vector<spacedword>::iterator End, double Coverage): _Begin(Begin), _End(End), _Coverage(Coverage){
}

void family_spacedword::set_begin(std::vector<spacedword>::iterator Begin){
    _Begin = Begin;
}
void family_spacedword::set_end(std::vector<spacedword>::iterator End){
    _End = End;
}
void family_spacedword::set_coverage(double Coverage){
    _Coverage = Coverage;
}


std::vector<spacedword>::iterator family_spacedword::begin() const{
    return _Begin;
}
std::vector<spacedword>::iterator family_spacedword::end() const{
    return _End;
}
double family_spacedword::coverage() const{
    return _Coverage;
}