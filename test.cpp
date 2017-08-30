#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "src/sequencefamily.hpp"
#include "src/spacedword.hpp"
#include "src/pattern.hpp"
#include "src/sequence.hpp"

int main(){
    std::string FamilyFile = "PF07690_seed.txt";
    std::string StringPat = "100001101010111";
    sequence_family SeqFam(FamilyFile,true);
    pattern Pat(StringPat);
    std::vector<spacedword> SpwVec;
    unsigned Ctr = 0;
    for(auto & Seq : SeqFam){
        std::vector<spacedword> Spws = Seq.spaced_words(Pat,Ctr++);
        SpwVec.insert(SpwVec.end(),Spws.begin(),Spws.end());
    }

    std::cout << "Sort" << std::endl;
    std::sort(SpwVec.begin(),SpwVec.end());
    auto Start = SpwVec.begin();
    std::cout << "Set Size" << std::endl;
    for(auto SpwIt = SpwVec.begin(); SpwIt != SpwVec.end(); SpwIt++){
        if(*SpwIt != *Start || (SpwIt+1) == SpwVec.end()){
            for(auto Iter = Start; Iter != SpwIt; Iter++){
                Iter->set_family(std::distance(Start,SpwIt));
            }
            Start = SpwIt;
        }
    }
    std::cout << "Sort size" << std::endl;
    std::sort(SpwVec.begin(),SpwVec.end(),[](spacedword & SpwA, spacedword & SpwB){
        if(SpwA.family() == SpwB.family()){
            return SpwA < SpwB;
        }
        return SpwA.family() > SpwB.family();
    });
    auto CutOff = SpwVec.begin();
    while(CutOff->family() > 1){
        CutOff++;
    }
    SpwVec.resize(std::distance(SpwVec.begin(),CutOff));


    std::unordered_map<uint64_t,size_t> SpwID;
    Ctr = 0;
    Start = SpwVec.begin();
    for(auto SpwIt = SpwVec.begin(); SpwIt != SpwVec.end(); SpwIt++){
        if(*SpwIt != *Start || (SpwIt+1) == SpwVec.end()){
            SpwID.insert(std::make_pair(Start->bits(),Ctr++));
            Start = SpwIt;
        }
    }

    std::vector< std::vector< std::pair<int,int> > > SpwOrdMat(SpwID.size(), std::vector< std::pair<int,int> >(SpwID.size(), std::pair<int,int>(0,0)));
    std::cout << "MatSize: " << SpwID.size() << "x" << SpwID.size() << std::endl;
    std::sort(SpwVec.begin(),SpwVec.end(),[](spacedword & SpwA, spacedword & SpwB){
        if(SpwA.sequence() == SpwB.sequence()){
            if(SpwA.position() == SpwB.position()){
                return SpwA < SpwB;
            }
            return SpwA.position() < SpwB.position();
        }
        return SpwA.sequence() < SpwB.sequence();
    });

    Start = SpwVec.begin();
    for(auto SpwItS = SpwVec.begin(); SpwItS != SpwVec.end()-1; SpwItS++){
        for(auto SpwItE = SpwItS+1; SpwItE != SpwVec.end(); SpwItE++){
            SpwOrdMat[SpwID[SpwItS->bits()]][SpwID[SpwItE->bits()]].first++;
            SpwOrdMat[SpwID[SpwItS->bits()]][SpwID[SpwItE->bits()]].second++;
            SpwOrdMat[SpwID[SpwItE->bits()]][SpwID[SpwItS->bits()]].first--;
            SpwOrdMat[SpwID[SpwItE->bits()]][SpwID[SpwItS->bits()]].second++;
        }
    
    }

    std::sort(SpwOrdMat.begin(), SpwOrdMat.end(),[](std::vector< std::pair<int,int> > & RowA, std::vector< std::pair<int,int> > & RowB){
        return
        std::accumulate(RowA.begin(),RowA.end(),0,[](int A, std::pair<int,int> & B){
            if(B.first > 0){
                return A + 1;
            }
            return A - 1;
        })
        <
        std::accumulate(RowB.begin(),RowB.end(),0,[](int A, std::pair<int,int> & B){
            if(B.first > 0){
                return A + 1;
            }
            return A - 1;
        });
    });

    for(unsigned i = 0; i < SpwOrdMat.size(); i++){
        for(unsigned j = 0; j < SpwOrdMat.size(); j++){
            //std::cout << "[" << SpwOrdMat[i][j].first << "," << SpwOrdMat[i][j].second << "] ";
            std::cout << "[" << (SpwOrdMat[i][j].first/(double)SpwOrdMat[i][j].second) << "] ";
        }
        std::cout << std::endl;
    }
}
