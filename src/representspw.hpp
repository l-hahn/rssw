#ifndef REPRESENTSPW_HPP_
#define REPRESENTSPW_HPP_

#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <vector>
#include "familyscore.hpp"

class rep_spw{
	public:
		rep_spw(int64_t BitSpacedWord);
		rep_spw(int64_t BitSpacedWord, unsigned Family, unsigned Position, double FamilyScore);
		rep_spw(int64_t BitSpacedWord, family_score & FamilyScore);

		int64_t bits() const;
		size_t size() const;
		void clear();
		void push_back(unsigned Family, unsigned Position, double Score);
		void push_back(family_score & FamilyScore);

		family_score operator[](size_t Idx) const;
		family_score & operator[](size_t Idx);
		bool operator==(const rep_spw & Spw) const;
		bool operator!=(const rep_spw & Spw) const;
		bool operator<(const rep_spw & Spw) const;
		bool operator>(const rep_spw & Spw) const;

		typedef std::vector<family_score>::iterator iterator;
		iterator begin();
		iterator end();

	private:
		int64_t _BitSpacedWord;
		std::vector<family_score> _Family;
};
#endif