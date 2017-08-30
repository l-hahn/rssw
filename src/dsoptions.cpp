#include "dsoptions.hpp"

namespace dsoptions{
    std::string InSeq;
    std::string InSwDb = "";
    std::string OutDetect = "Detection";
    double HitThreshold = 0;
    signed BlockSize = -1;
    signed WordPerSeq = -1;
    unsigned FamilySize = 0;
    unsigned FamilyOffset = 200;
    unsigned ThreadNumber = 1;
    bool AlphaReduction = false;
    bool Translated = false;
};