#ifndef RSWSET_HPP_
#define RSWSET_HPP_

#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <vector>
#include "representspw.hpp"
#include "pattern.hpp"
#include "spacedword.hpp"

class rsw_set{
    public:
        rsw_set();
        rsw_set(pattern & Pattern);
        rsw_set(pattern && Pattern);

        void push_back(rep_spw && Spw);
        void push_back(rep_spw & Spw);
        void merge(rsw_set & RswSet);
        void clear();
        void set_pattern(pattern && Pattern);
        void set_pattern(pattern & Pattern);
        pattern pat() const;
        size_t size() const;
        bool operator==(const rsw_set & Rsw) const;
        rep_spw & operator[](size_t Idx);

        typedef std::vector<rep_spw>::iterator iterator;
        iterator begin();
        iterator end();
        void insert(iterator FillEnd, iterator InputBegin, iterator InputEnd);

    private:
        pattern _Pattern;
        std::vector<rep_spw> _RepresentWords;
};
#endif