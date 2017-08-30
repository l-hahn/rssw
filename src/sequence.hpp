#ifndef SEQUENCE_HPP_
#define SEQUENCE_HPP_

#include <string>
#include <vector>
#include "alphabet.hpp"
#include "pattern.hpp"
#include "spacedword.hpp"

class sequence{
    friend class sequence_family;
    
    public:
        sequence();
        sequence(std::string & Seq, std::string SeqName = "", bool AlphaReduct = false);
        sequence(std::string && Seq, std::string SeqName = "", bool AlphaReduct = false);
        sequence(std::vector<char> & Seq, std::string SeqName = "", bool AlphaReduct = false);
        sequence(std::vector<char> && Seq, std::string SeqName = "", bool AlphaReduct = false);
        sequence(std::vector<unsigned char> & Seq, std::string SeqName = "", bool AlphaReduct = false);
        sequence(std::vector<unsigned char> && Seq, std::string SeqName = "", bool AlphaReduct = false);

        std::vector<spacedword> spaced_words(pattern & Pat, signed SeqNo = -1, signed short FamID = -1);
        void push_back(unsigned char Letter);
        void clear();
        void to_bit();
        std::string to_string() const;
        bool check_protein();
        bool is_bit() const;
        bool is_protein() const;
        bool is_dna() const;


        std::string name() const;
        size_t size() const;
        void set_name(std::string & Name);
        void set_name(std::string && Name);

        template<typename T>
        void set_sequence(T & Seq, bool AlphaReduct = false);

        unsigned char operator[](size_t Idx) const;
        bool operator<(const sequence & Seq) const;
        bool operator>(const sequence & Seq) const;
        typedef std::vector<unsigned char>::iterator iterator;
        iterator begin();
        iterator end();

    private:
        std::vector<unsigned char> _Sequence;
        std::string _SequenceName;
        bool _IsProtein;
};

template<typename T>
void sequence::set_sequence(T & Seq, bool AlphaReduct){
    if(AlphaReduct){
        for(auto & Letter : Seq){
            if(std::isalpha(Letter)){
                _Sequence.push_back(alphabet::AlphaReduct[std::toupper(Letter)]);
            }
            else if(Letter <= alphabet::ProteinBitAlphabet[alphabet::ProteinBitAlphabet.size()-1]){
                _Sequence.push_back(alphabet::AlphaBitReduct[Letter]);
            }
            else if(Letter == '*'){
                _Sequence.push_back(Letter);
            }
        }    
        _IsProtein = check_protein();
    }
    else{
        for(auto & Letter : Seq){
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
        _IsProtein = check_protein();
    }
}

#endif