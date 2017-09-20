#include "spacedworddb.hpp"
spacedword_db::spacedword_db(bool IsProtein, unsigned WordPerSeq): _WordPerSeq(WordPerSeq), _TypeProtein(IsProtein){
}

spacedword_db::spacedword_db(std::string & FileName){
    from_file(FileName);
}

spacedword_db::spacedword_db(rsw_set & SpwFam, bool IsProtein, unsigned WordPerSeq): _WordPerSeq(WordPerSeq), _TypeProtein(IsProtein){
    push_back(SpwFam);
}

void spacedword_db::push_back(rsw_set & SpwFam){
    _SpacedWordsList.push_back(SpwFam);
    std::sort(_SpacedWordsList[_SpacedWordsList.size()-1].begin(), _SpacedWordsList[_SpacedWordsList.size()-1].end());
}


void spacedword_db::merge(spacedword_db & DBAddition){
    // std::vector<std::string> & NameVec = DBAddition._FamilyName;
    // if(_FamilyName.size() < NameVec.size()){
    //     _FamilyName.resize(NameVec.size(),"");
    // }
    // for(unsigned i = 0; i < NameVec.size(); i++){
    //     if(NameVec[i].size() != 0){
    //         _FamilyName[i] = NameVec[i];
    //     }
    // }
    std::unordered_map<std::string,unsigned> & AddFamNameID = DBAddition._FamilyNameID;
    for(auto & Entry : AddFamNameID){
        auto Search = _FamilyNameID.find(Entry.first);
        if(Search != _FamilyNameID.end()){
            _FamilyNameID.insert(Entry);
            _IDFamilyName.insert(std::make_pair(Entry.second,Entry.first));
        }
    }


    for(auto & DBAdd : DBAddition){
        unsigned PatPos;
        auto PatIt = std::find_if(_SpacedWordsList.begin(), _SpacedWordsList.end(), [DBAdd](const rsw_set & SwFam){
            return SwFam == DBAdd;
        });
        if(PatIt != _SpacedWordsList.end()){
            PatPos = std::distance(_SpacedWordsList.begin(), PatIt);
        }
        else{
            PatPos = _SpacedWordsList.size();
            _SpacedWordsList.resize(_SpacedWordsList.size()+1);
        }
        DBAdd.merge(_SpacedWordsList[PatPos]);
    }
}


void spacedword_db::add_family_id(std::string && FamID, unsigned ID){
    add_family_id(FamID, ID);
}
void spacedword_db::add_family_id(std::string & FamID, unsigned ID){
    // if(ID >= _FamilyName.size()){
    //     _FamilyName.resize(ID+1,"");
    // }
    // _FamilyName[ID] = FamID;
    _FamilyNameID.insert(std::make_pair(FamID,ID));
    _IDFamilyName.insert(std::make_pair(ID,FamID));
}

// std::vector<std::string> spacedword_db::get_name_id() const{
//     return _FamilyName;
// }

// unsigned spacedword_db::new_id(){
//     _FamilyName.push_back("");
//     return _FamilyName.size()-1;
// }

std::string spacedword_db::id_fam_name(unsigned ID){
    //return _FamilyName[ID];
    return _IDFamilyName[ID];
}


signed spacedword_db::fam_name_id(std::string && Str){
    return fam_name_id(Str);
}
signed spacedword_db::fam_name_id(std::string & Str){
    // auto Pos = std::find(_FamilyName.begin(), _FamilyName.end(), Str);
    // if(Pos != _FamilyName.end()){
    //     return std::distance(_FamilyName.begin(), Pos);
    // }
    // return -1;
    return _FamilyNameID[Str];
}


rsw_set & spacedword_db::operator[](size_t Idx){
    return _SpacedWordsList[Idx];
}

size_t spacedword_db::size() const{
    return _SpacedWordsList.size();
}

size_t spacedword_db::families_size() const{
    return _FamilyNameID.size();
}

unsigned spacedword_db::words_per_seq() const{
    return _WordPerSeq;
}

bool spacedword_db::is_protein() const{
    return _TypeProtein == true;
}

bool spacedword_db::is_dna() const{
    return is_protein() == false;
}


spacedword_db::iterator spacedword_db::begin(){
    return _SpacedWordsList.begin();
}

spacedword_db::iterator spacedword_db::end(){
    return _SpacedWordsList.end();
}


void spacedword_db::to_file(std::string & FileName, bool AlphaReduction){
    std::string Type = "Protein";
    if(_TypeProtein == false){
        Type = "DNA";
    }
    std::ofstream Output(FileName+".swdb");
    std::string ReductText = "";
    if(AlphaReduction){
        ReductText = " R";
    }
    Output << "#<" << Type << " Spaced Word Data Base WPS " << _WordPerSeq << ReductText << ">#" << std::endl;
    for(auto & PatSpw : _SpacedWordsList){
        Output << "#" << PatSpw.pat().to_string() << "#" << std::endl;
        for(auto & SpwFam : PatSpw){
            Output << "<" << spacedword(SpwFam.bits(),0).to_string() << ">" << std::endl;
            for(auto & Fam : SpwFam){
                Output << id_fam_name(Fam.id()) << "\t" << Fam.score() << "\t" << Fam.position() << "\t" << Fam.pos_scr() << std::endl;
            }
        }
    }
    Output.close();
}

void spacedword_db::from_file(std::string & FileName){
    std::ifstream Input(FileName);
    std::string Parser;
    std::getline(Input, Parser);
    if(Parser.substr(0,2) == "#<"){
        if(Parser[2] == 'P'){
            _TypeProtein = true;
        }
        else if(Parser[2] == 'D'){
            _TypeProtein = false;
        }
        else{
            std::cerr << "FATAL ERROR!\nCannot parse Sequence Type!" <<std::endl;
            std::cerr <<"Header-Format: #<[DNA|PROTEIN] ... [WPS:<INT>]>#" << std::endl;
            std::exit(-1);
        }
        std::string WPS;
        unsigned StartDiffPos = 3;
        if(Parser[Parser.size()-StartDiffPos] == 'R' || Parser[Parser.size()-StartDiffPos] == 'r'){
            StartDiffPos = 5;
        }
        for(unsigned Pos = Parser.size()-StartDiffPos; Pos >= 0; Pos--){
            if(std::isdigit(Parser[Pos])){
                WPS.push_back(Parser[Pos]);
            }
            else{
                break;
            }
        }
        _WordPerSeq = std::stoi(std::string(WPS.rbegin(), WPS.rend()));
        if(_WordPerSeq == 0){
            std::cerr << "FATAL ERROR!\nCannot parse Words per Sequence(WPS) parameter!" <<std::endl;
            std::cerr <<"Header-Format: #<[DNA|PROTEIN] ... [WPS:<INT>]>#" << std::endl;
            std::exit(-1);
        }
        rsw_set TmpSpwFam;
        spacedword Spw;
        while(!Input.eof()){
            Input >> Parser;
            if(Parser.size() != 0){
                switch(Parser[0]){
                    case '<':
                        Parser = Parser.substr(1,Parser.size()-2);
                        Spw = spacedword();
                        for(char c : Parser){
                            Spw.push_back(c);
                        }
                        TmpSpwFam.push_back(rep_spw(Spw.bits()));
                        break;
                    case '#':
                        if(TmpSpwFam.size() != 0){
                            _SpacedWordsList.push_back(TmpSpwFam);
                            TmpSpwFam.clear();
                        }
                        Parser = Parser.substr(1,Parser.size()-2);
                        TmpSpwFam.set_pattern(pattern(Parser));
                        break;
                    default:
                        // signed ID = fam_name_id(Parser);
                        // if(ID < 0){
                        //     ID = new_id();
                        //     add_family_id(Parser, ID);
                        // }
                        auto search = _FamilyNameID.find(Parser);
                        unsigned ID = _FamilyNameID.size();
                        if(search != _FamilyNameID.end()){
                            ID = _FamilyNameID[Parser];
                        }
                        else{
                            add_family_id(Parser, ID);
                        }
                        double Score, PosRelScr;
                        unsigned Position;
                        Input >> Score >> Position >> PosRelScr;
                        TmpSpwFam[TmpSpwFam.size()-1].push_back(ID, Position, Score, PosRelScr);
                }
                Parser.clear();
            }
        }
        if(TmpSpwFam.size() != 0){
            _SpacedWordsList.push_back(TmpSpwFam);
            TmpSpwFam.clear();
        }
    }
    else{
        std::cerr << "FATAL ERROR!\nCannot parse File!" <<std::endl;
        std::cerr << "Header-Format: #<[DNA|PROTEIN] ... [WPS:<INT>]>#" << std::endl;
        std::cerr << "Data-Format:\n#[Pattern]#\n<[SpacedWord]>\n[FamID]\t[Score]\t[Order Number]\n..." << std::endl;
        std::exit(-1);
    }
    Input.close();
}