#include "familyscore.hpp"

family_score::family_score(unsigned ID, unsigned Position, double Score): _ID(ID), _Position(Position), _Score(Score){
}

unsigned family_score::id() const{
    return _ID;
}

unsigned family_score::position() const{
    return _Position;
}

double family_score::score() const{
    return _Score;
}

void family_score::set_id(unsigned ID){
    _ID = ID;
}
void family_score::set_position(unsigned Pos){
    _Position = Pos;
}
void family_score::set_score(double Score){
    _Score = Score;
}

bool family_score::operator<(const family_score & F) const{
    return _Position < F._Position;
}

bool family_score::operator>(const family_score & F) const{
    return _Position > F._Position;
}

bool family_score::operator==(const family_score & F) const{
    return _ID == F._ID;
}