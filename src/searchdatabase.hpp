#ifndef SEARCHDATABASE_HPP_
#define SEARCHDATABASE_HPP_

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <omp.h>
#include <vector>
#include "alphabet.hpp"
#include "dsoptions.hpp"
#include "familyscore.hpp"
#include "pattern.hpp"
#include "sequence.hpp"
#include "sequencefamily.hpp"
#include "rswset.hpp"
#include "spacedhit.hpp"
#include "spacedword.hpp"
#include "spacedworddb.hpp"

void search_database();
void read_family(std::vector<sequence_family> & Families);
void read_database(spacedword_db & SwDataBase, std::vector<sequence_family> & Families);
void multi_files(std::vector<sequence_family> & Families, spacedword_db & SwDataBase);
void set_multi_file_name();
void matching_db_seq(std::vector<sequence_family> & Families, spacedword_db & SwDataBase);

std::vector<spacedhit> family_hit(sequence_family & Family, unsigned FamID, spacedword_db & SwDataBase);
std::vector<spacedhit> pattern_hit(sequence_family & Family, unsigned FamID, rsw_set & DbSw);
void classify(spacedword_db & SwDataBase, std::vector<sequence_family> & Families, std::vector<spacedhit> & DataBaseHit);

void init_omp();

void evaluate_sens_spec(spacedword_db & SwDataBase, std::vector<sequence_family> & Families, std::vector<spacedhit> & DataBaseHit);
#endif