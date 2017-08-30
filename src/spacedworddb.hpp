#ifndef SPACEDWORDDB_HPP_
#define SPACEDWORDDB_HPP_

#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <unordered_map>
#include <string>
#include <sstream>
#include <sys/stat.h> 
#include <vector>
#include "alphabet.hpp"
#include "rswset.hpp"
#include "spacedword.hpp"
#include "pattern.hpp"

class spacedword_db{
    public:
        spacedword_db(bool IsProtein, unsigned WordPerSeq = 1);
        spacedword_db(std::string & FileName);
        spacedword_db(rsw_set & SpwFam, bool IsProtein, unsigned WordPerSec = 1);

        void push_back(rsw_set & SpwFam);
        void to_file(std::string & FileName, bool AlphaReduction = false);
        void from_file(std::string & FileName);
        void merge(spacedword_db & DBAddition);

        void add_family_id(std::string && FamID, unsigned ID);
        void add_family_id(std::string & FamID, unsigned ID);
        std::string id_fam_name(unsigned ID);
        signed fam_name_id(std::string && Str);
        signed fam_name_id(std::string & Str);

        size_t size() const;
        size_t families_size() const;
        bool is_protein() const;
        bool is_dna() const;
        unsigned words_per_seq() const;

        rsw_set & operator[](size_t Idx);
        typedef std::vector<rsw_set>::iterator iterator;
        iterator begin();
        iterator end();

    private:
        //std::vector<std::string> get_name_id() const;
        //unsigned new_id();

        std::vector<rsw_set> _SpacedWordsList;
        //std::vector<std::string> _FamilyName;
        std::unordered_map<std::string,unsigned> _FamilyNameID;
        std::unordered_map<unsigned,std::string> _IDFamilyName;
        unsigned _WordPerSeq;
        bool _TypeProtein;
};
#endif