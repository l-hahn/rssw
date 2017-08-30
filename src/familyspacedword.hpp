#ifndef FAMILYSPACEDWORD_HPP_
#define FAMILYSPACEDWORD_HPP_

#include <iterator>
#include <vector>
#include "spacedword.hpp"


class family_spacedword{
    public:
        family_spacedword(std::vector<spacedword>::iterator Begin, std::vector<spacedword>::iterator End, double Coverage);

        void set_begin(std::vector<spacedword>::iterator Begin);
        void set_end(std::vector<spacedword>::iterator End);
        void set_coverage(double Coverage);

        std::vector<spacedword>::iterator begin() const;
        std::vector<spacedword>::iterator end() const;
        double coverage() const;

    private:
        std::vector<spacedword>::iterator _Begin;
        std::vector<spacedword>::iterator _End;
        double _Coverage;
};
#endif