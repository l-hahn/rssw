#ifndef CREATEDATABASE_HPP_
#define CREATEDATABASE_HPP_

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <omp.h>
#include <utility>
#include <string>
#include <vector>
#include "alphabet.hpp"
#include "dboptions.hpp"
#include "familyspacedword.hpp"
#include "sequencefamily.hpp"
#include "spacedword.hpp"
#include "spacedworddb.hpp"
#include "pattern.hpp"
#include "patternset.hpp"

void create_database();
void get_seqcov(std::vector<sequence_family> & Families, std::vector<double> & SequenceCoverage);
void get_overlap(spacedword_db & SwDataBase);

void multi_files(std::vector<sequence_family> & FamilyVec, spacedword_db & FamiliesSpacedWords);

void create_rssw(std::vector<sequence_family> & FamilyVec, spacedword_db & FamiliesSpacedWords);
void create_words(std::vector< std::vector<spacedword> > & FamilySpacedWords, std::vector<sequence_family> & FamilyVec, pattern & Pat);
void sort_hash(std::vector< std::vector<spacedword> > & SpacedHashVector, std::vector< std::vector<spacedword> > & FamilySpacedWords, std::vector<int64_t> HashKeyVec, char MaskBlocks, int64_t MaxVal);
void extract_unique_spw(std::vector< std::vector<spacedword> > & SpacedHashVector, std::vector<int64_t> HashKeyVec);
void family_unique(std::vector< std::vector<spacedword> > & SpacedHashVector, std::vector< std::vector<spacedword> > & FamilySpacedWords, std::vector<int64_t> HashKeyVec);
void create_not_unique_spw(std::vector< std::vector<spacedword> > & FamilySpacedWords);
void extract_family_rssw(std::vector< std::vector<spacedword> > & FamilySpacedWords, std::vector<sequence_family> & FamilyVec, rsw_set & RSSW_SPW);

std::vector< std::pair<size_t, double> > spw_order(std::vector<spacedword> & SpwVec)

void clean_line();
void init_omp();

#endif
