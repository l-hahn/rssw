#include "createdatabase.hpp"

std::vector<int> FamilySize;

void create_database(){
    patternset PatSet;
    if(dboptions::PatternFile.size() != 0){
        //TODO READ PATTERN FILE
    }
    std::vector<sequence_family> Families = sequence_family::families_file(dboptions::InSeqFam, dboptions::AlphaReduction);
    clean_line();
    if(Families.size() == 1){
        std::cout << Families.size() << " Family was read. Inputtype: ";
    }
    else{
        std::cout << Families.size() << " Families were read. Inputtype: ";
    }
    dboptions::FamilySize = Families.size();
    if(Families[0].is_protein()){
        std::cout << "Protein " << std::endl;
    }
    else{
        std::cout << "DNA" << std::endl;
    }

    std::vector<double> SequenceCoverage(Families.size(), -1);
    init_omp();
    spacedword_db FamiliesSpacedWords(Families[0].is_protein(), dboptions::WordPerSeq);
    multi_files(Families, FamiliesSpacedWords);

    // create_rssw(Families, FamiliesSpacedWords);
    // std::cout << "[  ] Writing database to file '" << dboptions::OutSwDb << ".swdb' ..." << std::flush;
    // FamiliesSpacedWords.to_file(dboptions::OutSwDb, dboptions::AlphaReduction);
    // std::cout << "\r[OK] Writing database to file '" << dboptions::OutSwDb << ".swdb' ..." << std::endl;
}


void multi_files(std::vector<sequence_family> & FamilyVec, spacedword_db & FamiliesSpacedWords){
    bool End = false;
    int Try = 3;
    std::string Parse, Opt;
    while(!End && Try > 0){
        if(Try == 3){
            create_rssw(FamilyVec, FamiliesSpacedWords);
            std::string OutFile = dboptions::OutSwDb + "_c"+std::to_string(dboptions::SpwSeqThrshld)
                +"_n"+std::to_string(dboptions::SpwNoiseThrshld)+"_s"+std::to_string(dboptions::WordPerSeq)
                +"_m"+std::to_string(dboptions::PatternNumber)+"_w"+std::to_string(dboptions::PatternWeight)
                +"_d"+std::to_string(dboptions::PatternMinDC);
            std::cout << "[  ] Writing database to file '" << OutFile << ".swdb' ..." << std::flush;
            FamiliesSpacedWords.to_file(OutFile, dboptions::AlphaReduction);
            std::cout << "\r[OK] Writing database to file '" << OutFile << ".swdb' ..." << std::endl;
        }
        std::cout << "Again? [Y/N]: ";
        std::cin >> Parse;
        if(Parse[0] == 'Y' || Parse[0] == 'y'){
            std::cout << "Used parameter set:" << std::endl;
            std::cout << "-c " << dboptions::SpwSeqThrshld << " | -n " << dboptions::SpwNoiseThrshld << " | -s " << dboptions::WordPerSeq;
            std::cout << " | -m " << dboptions::PatternNumber << " | -w " << dboptions::PatternWeight << " | -d " << dboptions::PatternMinDC;
            std::cout << "-" << dboptions::PatternMaxDC << " | -f " << dboptions::FamilyOffset << " | -t " << dboptions::ThreadNumber << std::endl;
            Try = 3;
            std::cout << "New Values? [Y/N]: ";
            std::cin >> Parse;
            if(Parse[0] == 'Y' || Parse[0] == 'y'){
                std::cout << "Use cmd line parameter pair, and 'END END' for stopping new value submission." << std::endl;
                std::cout << "Param: ";
                std::cin >> Parse >> Opt;
                while(Parse != "END"){
                    switch(Parse[1]){
                        case 'C':
                        case 'c':
                            dboptions::SpwSeqThrshld = std::stof(Opt);
                            break;
                        case 'D':
                        case 'd':
                            dboptions::parse_length(Opt.c_str());
                            break;
                        case 'F':
                        case 'f':
                            dboptions::FamilyOffset = std::stoi(Opt);
                            break;
                        case 'M':
                        case 'm':
                            dboptions::PatternNumber = std::stoi(Opt);
                            break;
                        case 'N':
                        case 'n':
                            dboptions::SpwNoiseThrshld = std::stof(Opt);
                            break;
                        case 'S':
                        case 's':
                            dboptions::WordPerSeq = std::stoi(Opt);
                            break;
                        case 'T':
                        case 't':
                            dboptions::ThreadNumber = std::stoi(Opt);
                            init_omp();
                            break;
                        case 'W':
                        case 'w':
                            dboptions::PatternWeight = std::stoi(Opt);
                            break;
                        default:
                            std::cerr << "Cannot parse '" << Parse << " " << Opt << "'!" << std::endl;
                    }
                    std::cout << "Param: ";
                    std::cin >> Parse >> Opt;
                }
            }
            else{
                std::cout << "No changes!" << std::endl;
            }
            std::cout << std::endl;
        }
        else if (Parse[0] == 'N' || Parse[0] == 'N'){
            End = true;
        }
        else{
            std::cerr << "Please submit 'Y' for yes or 'N' for no! Another try" << std::endl;
            Try--;
        }
    }
    if(Try == 0){
        std::cerr << "Three parsing attempts without correct understanding!" << std::endl;
    }
    std::cout << "Exit program!" << std::endl;
}



void create_rssw(std::vector<sequence_family> & FamilyVec, spacedword_db & FamiliesSpacedWords){
    char MaskBlocks = std::min((unsigned)3,dboptions::PatternWeight);
    int64_t MaxVal = alphabet::ProteinBitAlphabet[alphabet::ProteinBitAlphabet.size()-1];
    int64_t Tmp = MaxVal;
    for(int i = 0; i < MaskBlocks-1; i++){
        MaxVal <<= alphabet::BitLength;
        MaxVal |= Tmp;
    }
    std::vector<int64_t> HashKeyVec;
    for(auto Acid : alphabet::ProteinBitAlphabet){
        for(auto Acid2 : alphabet::ProteinBitAlphabet){
            for(auto Acid3 : alphabet::ProteinBitAlphabet){
                int64_t Prefix = Acid;
                Prefix <<= alphabet::BitLength;
                Prefix |= Acid2;
                Prefix <<=alphabet::BitLength;
                Prefix |= Acid3;
                HashKeyVec.push_back(Prefix);
            }
        }
    }
    patternset CoverSet;
    char Counter = 100;
   
    std::vector< std::vector<unsigned> > FamCovSize(FamilyVec.size());
    for(unsigned i = 0; i < FamCovSize.size(); i++){
            FamCovSize[i] = std::vector<unsigned>(FamilyVec[i].size(),0);
            FamiliesSpacedWords.add_family_id(FamilyVec[i].id(), i);
    }

    for(unsigned i = 0; i < dboptions::PatternNumber; i++){
        pattern Pat(dboptions::PatternWeight, dboptions::PatternMaxDC);
        std::cout << "Pattern " << i+1 << "/" << dboptions::PatternNumber << std::endl;
        while(CoverSet.is_uniq(Pat) == false && Counter > 1){
            Pat = pattern(dboptions::PatternWeight, dboptions::PatternMaxDC);
            Counter--;
        }
        if(Counter == 0){
            break;
        }
        rsw_set RSSW_SPW(Pat);
        CoverSet.push_back(Pat);
        std::vector< std::vector<spacedword> > FamilySpacedWords(FamilyVec.size());
        std::vector< std::vector<spacedword> > SpacedHashVector(MaxVal+1);

        create_words(FamilySpacedWords, FamilyVec, Pat);
        sort_hash(SpacedHashVector, FamilySpacedWords, HashKeyVec, MaskBlocks, MaxVal);
        extract_unique_spw(SpacedHashVector, HashKeyVec);
        family_unique(SpacedHashVector, FamilySpacedWords, HashKeyVec);
        create_not_unique_spw(FamilySpacedWords);
        extract_family_rssw(FamilySpacedWords, FamilyVec, RSSW_SPW);
        FamiliesSpacedWords.push_back(RSSW_SPW);
    }
}


void create_words(std::vector< std::vector<spacedword> > & FamilySpacedWords, std::vector<sequence_family> & FamilyVec, pattern & Pat){
    std::cout << "\t[  ] Creating all spaced words" << std::flush;
    FamilySize = std::vector<int>(FamilyVec.size(),0); 
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic) nowait
        for(unsigned i = 0; i < FamilyVec.size(); i++){
            if(FamilyVec[i].size() >= dboptions::FamilyOffset){
                for(unsigned j = 0; j < FamilyVec[i].size(); j++){
                    std::vector<spacedword> FamilyWords = FamilyVec[i][j].spaced_words(Pat, j, i);
                    FamilySpacedWords[i].insert(FamilySpacedWords[i].end(), FamilyWords.begin(), FamilyWords.end());
                    if(FamilyVec[i][j].size() >= Pat.length()){
                        FamilySize[i]++;
                    }
                }
            }
        }
    }
    std::cout << "\r\t[OK] Creating all spaced words" << std::endl;
}

void sort_hash(std::vector< std::vector<spacedword> > & SpacedHashVector, std::vector< std::vector<spacedword> > & FamilySpacedWords, std::vector<int64_t> HashKeyVec, char MaskBlocks, int64_t MaxVal){
    int64_t HashMask = 0;
    for(int i = 0; i < MaskBlocks; i++){
        HashMask <<= alphabet::BitLength;
        HashMask |= (1 << alphabet::BitLength)-1;
    }
    for(unsigned i = 0; i < dboptions::PatternWeight-MaskBlocks; i++){
        HashMask <<= alphabet::BitLength;
    }

    std::cout << "\t[  ] Sorting and hashing Spaced Words from different Families" << std::flush;
    unsigned & NumThrd = dboptions::ThreadNumber;
    unsigned ThrdSize = (size_t)std::ceil((double)HashKeyVec.size()/NumThrd);
    if(ThrdSize < 1){
        ThrdSize = 1;
    }
    #pragma omp parallel
    {
        unsigned ThrdNum = omp_get_thread_num();
        unsigned TPrefMax = ThrdSize*ThrdNum + (ThrdSize - 1);
        if(ThrdSize*ThrdNum < HashKeyVec.size()){
            if(TPrefMax >= HashKeyVec.size()){
                TPrefMax = HashKeyVec.size()-1;
            }
            int64_t HashMin = HashKeyVec[ThrdSize*ThrdNum];
            int64_t HashMax = HashKeyVec[TPrefMax];
            for(unsigned i = 0; i < FamilySpacedWords.size(); i++){
                for(auto & Spw : FamilySpacedWords[i]){
                    int64_t Hash = ((Spw.bits() & HashMask) >> ((dboptions::PatternWeight-MaskBlocks)*alphabet::BitLength));
                    if(HashMin <= Hash && Hash <= HashMax){
                        SpacedHashVector[Hash].push_back(Spw);
                    }
                }
            }
            // for(unsigned i = HashMin; i <= HashMax; i++){
            //     SpacedHashVector[i].shrink_to_fit();
            // }
        }
    }
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for(unsigned i = 0; i < FamilySpacedWords.size(); i++){
            FamilySpacedWords[i].clear();
        }
    }

    std::cout << "\r\t[OK] Sorting and hashing Spaced Words from different Families" << std::endl;
}

void extract_unique_spw(std::vector< std::vector<spacedword> > & SpacedHashVector, std::vector<int64_t> HashKeyVec){
    std::cout << "\t[  ] Extracting unique spaced words" << std::flush;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic) nowait
        for(unsigned i = 0 ; i < HashKeyVec.size(); i++){
            std::vector<spacedword> & HashVec = SpacedHashVector[HashKeyVec[i]], UniqHashVec;
            if(HashVec.size() == 0){
                continue;
            }
            std::sort(HashVec.begin(), HashVec.end(), [](const spacedword & SpwA, const spacedword & SpwB){
                if(SpwA == SpwB){
                    return SpwA.family() < SpwB.family();
                }
                return SpwA < SpwB;
            });
            std::vector< family_spacedword > HashedFamilies;
            auto BucketBegin = HashVec.begin();
            auto BucketEnd = HashVec.begin();
            bool BucketStop = false;
            HashedFamilies.push_back(family_spacedword(BucketBegin, BucketBegin, 0));
            while(BucketBegin != HashVec.end()){
                if(BucketEnd == HashVec.end()){
                    BucketStop = true;
                }
                else if(*BucketBegin != *BucketEnd){
                    BucketStop = true;
                }
                if(BucketStop){
                    int FamID = HashedFamilies[HashedFamilies.size()-1].begin()->family();
                    HashedFamilies[HashedFamilies.size()-1].set_end(BucketEnd);
                    HashedFamilies[HashedFamilies.size()-1].set_coverage(std::distance(HashedFamilies[HashedFamilies.size()-1].begin(), BucketEnd)/(double) FamilySize[FamID]);
                    std::sort(HashedFamilies.begin(), HashedFamilies.end(), [](const family_spacedword & SpwItA, const family_spacedword & SpwItB){
                        return SpwItA.coverage() > SpwItB.coverage();
                    });
                    if(std::distance(BucketBegin,BucketEnd) != 0){
                        if(HashedFamilies.size() > 1){
                            if(HashedFamilies[0].coverage() >= dboptions::SpwSeqThrshld && HashedFamilies[1].coverage() < dboptions::SpwNoiseThrshld){
                                UniqHashVec.insert(UniqHashVec.end(), HashedFamilies[0].begin(), HashedFamilies[0].end());
                            }
                        }
                        else{
                            if(HashedFamilies[0].coverage() >= dboptions::SpwSeqThrshld){
                                UniqHashVec.insert(UniqHashVec.end(), BucketBegin, BucketEnd);
                            }
                        }
                    }
                    BucketBegin = BucketEnd;
                    HashedFamilies.clear();
                    HashedFamilies.push_back(family_spacedword(BucketBegin, BucketBegin, 0));
                    BucketStop = false;
                }
                else{
                    if(HashedFamilies[HashedFamilies.size()-1].begin()->family() != BucketEnd->family()){
                        HashedFamilies[HashedFamilies.size()-1].set_end(BucketEnd);
                        int FamID = HashedFamilies[HashedFamilies.size()-1].begin()->family();
                        HashedFamilies[HashedFamilies.size()-1].set_coverage(std::distance(HashedFamilies[HashedFamilies.size()-1].begin(), BucketEnd)/(double) FamilySize[FamID]);
                        HashedFamilies.push_back(family_spacedword(BucketEnd, BucketEnd, 0));
                    }
                    BucketEnd++;
                }
            }
            std::swap(HashVec, UniqHashVec);
        }
    }
    std::cout << "\r\t[OK] Extracting unique spaced words" << std::endl;
}

void family_unique(std::vector< std::vector<spacedword> > & SpacedHashVector, std::vector< std::vector<spacedword> > & FamilySpacedWords, std::vector<int64_t> HashKeyVec){
    unsigned & NumThrd = dboptions::ThreadNumber;
    unsigned ThrdSize = (size_t)std::ceil((double)FamilySpacedWords.size()/NumThrd);
    if(ThrdSize < 1){
        ThrdSize = 1;
    }
    std::cout << "\t[  ] Returning unique spaced words to families" << std::flush;
    #pragma omp parallel
    {
        int ThrdNum = omp_get_thread_num();
        int FamIDMin = ThrdSize*ThrdNum;
        int FamIDMax = FamIDMin + (ThrdSize - 1);
        if((unsigned)FamIDMax >= FamilySpacedWords.size()){
            FamIDMax = FamilySpacedWords.size()-1;
        }
        for(unsigned i = 0; i < SpacedHashVector.size(); i++){
            for(auto & Spw : SpacedHashVector[i]){
                if(FamIDMin <= Spw.family() && Spw.family() <= FamIDMax){
                    FamilySpacedWords[Spw.family()].push_back(Spw);
                }
            }
        }
        // for(unsigned i = FamIDMin; i <= FamIDMax; i++){
        //     FamilySpacedWords[i].shrink_to_fit();
        // }
    }
    SpacedHashVector.clear();
    std::cout << "\r\t[OK] Returning unique spaced words to families" << std::endl;
}

void create_not_unique_spw(std::vector< std::vector<spacedword> > & FamilySpacedWords){
    std::cout << "\t[  ] Creating not-unique spaced words per family" << std::flush;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic) nowait
        for(unsigned i = 0; i < FamilySpacedWords.size(); i++){
            std::sort(FamilySpacedWords[i].begin(), FamilySpacedWords[i].end(), [](const spacedword & SpwA, const spacedword & SpwB){
                if(SpwA == SpwB){
                    return SpwA.sequence() < SpwB.sequence();
                }
                return SpwA < SpwB;
            });
            std::vector<spacedword> RSWFam;
            auto BucketBegin = FamilySpacedWords[i].begin();
            auto BucketEnd = FamilySpacedWords[i].begin();
            bool BucketStop = false;
            while(BucketBegin != FamilySpacedWords[i].end()){
                if(BucketEnd == FamilySpacedWords[i].end()){
                    BucketStop = true;
                }
                else if(*BucketBegin != *BucketEnd){
                    BucketStop = true;
                }
                if(BucketStop){
                    unsigned BucketSize = std::distance(BucketBegin,BucketEnd);
                    if(BucketSize >= dboptions::MinSeqWord){
                        std::for_each(BucketBegin, BucketEnd, [BucketSize](spacedword & Spw){
                            Spw.set_family(BucketSize);
                        });
                        RSWFam.insert(RSWFam.end(), BucketBegin, BucketEnd);
                    }
                    BucketBegin = BucketEnd;
                    BucketStop = false;
                }
                else{
                    BucketEnd++;
                }
            }
            //RSWFam.shrink_to_fit();
            std::swap(FamilySpacedWords[i],RSWFam);
        }
    }
    std::cout << "\r\t[OK] Creating not-unique spaced words per family" << std::endl;
}

void extract_family_rssw(std::vector< std::vector<spacedword> > & FamilySpacedWords, std::vector<sequence_family> & FamilyVec, rsw_set & RSSW_SPW){
    std::cout << "\t[  ] Extracting representativ spaced words" << std::flush;
    #pragma omp parallel
    {
        rsw_set RswSet(RSSW_SPW.pat());
        #pragma omp for schedule(dynamic) nowait
        for(unsigned i = 0; i < FamilySpacedWords.size(); i++){
            auto & SpacedBuckets = FamilySpacedWords[i];
            rsw_set LocalRswSet(RSSW_SPW.pat()); 
            std::vector<unsigned short> FamSeqCover(FamilyVec[i].size(), 0);
            std::sort(SpacedBuckets.begin(), SpacedBuckets.end(), [](const spacedword & SpwA, const spacedword & SpwB){
                if(SpwA.family() == SpwB.family()){ //Family := BucketSize
                    return SpwA < SpwB;
                }
                return SpwA.family() > SpwB.family();
            });
            double MaxWordCover = 0;
            if(SpacedBuckets.size() > 0){
                MaxWordCover = SpacedBuckets[0].family();
            }

            std::vector< std::pair<size_t, double> > SpwOrder = spw_order(SpacedBuckets);

            if(dboptions::Greedy == true){
                bool NewSeq = false;
                auto BucketBegin = SpacedBuckets.begin();
                unsigned PosCtr = 0;
                for(auto BucketIter = BucketBegin; BucketIter != SpacedBuckets.end(); BucketIter++){
                    if((*BucketIter != *BucketBegin || (BucketIter+1) == SpacedBuckets.end()) && NewSeq == true){
                        LocalRswSet.push_back(rep_spw(BucketBegin->bits(), i, SpwOrder[PosCtr].first, BucketBegin->family()/MaxWordCover,SpwOrder[PosCtr].second)); //Family := BucketSize
                        BucketBegin = BucketIter;
                        NewSeq = false;
                        PosCtr++;
                    }
                    if(FamSeqCover[BucketIter->sequence()] < dboptions::WordPerSeq){
                        NewSeq = true;
                    }
                    FamSeqCover[BucketIter->sequence()]++;
                }
                if(NewSeq && BucketBegin != SpacedBuckets.end()){
                    LocalRswSet.push_back(rep_spw(BucketBegin->bits(), i, SpwOrder[PosCtr].first, BucketBegin->family()/MaxWordCover,SpwOrder[PosCtr].second)); //Family := BucketSize
                }
                RswSet.merge(LocalRswSet);
            }
            else{
                //TODO: NOT GREEDY!!
            }
        }
        #pragma omp critical
        {
            RSSW_SPW.merge(RswSet);
        }
    }
    std::cout << "\r\t[OK] Extracting representativ spaced words" << std::endl;
}

std::vector< std::pair<size_t, double> > spw_order(std::vector<spacedword> & SpwVec){
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
    std::sort(SpwVec.begin(),SpwVec.end(),[](spacedword & SpwA, spacedword & SpwB){
        if(SpwA.sequence() == SpwB.sequence()){
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


    std::vector<size_t> IndexSort(SpwOrdMat.size());
    std::iota(IndexSort.begin(),IndexSort.end(), 0);
    std::sort(IndexSort.begin(),IndexSort.end(),[&](int L, int R){
        double LM = std::accumulate(SpwOrdMat[L].begin(),SpwOrdMat[L].end(),0.0,[](double A, std::pair<int,int> & B){
            if(B.second > 0){
                return A + (B.first/(double)B.second);    
            }
            return A;
        });
        double RM = std::accumulate(SpwOrdMat[R].begin(),SpwOrdMat[R].end(),0.0,[](double A, std::pair<int,int> & B){
            if(B.second > 0){
                return A + (B.first/(double)B.second);    
            }
            return A;
        });
        return LM < RM; 
    });

    std::vector< std::pair<size_t, double> > PosScr(IndexSort.size());

    std::transform(IndexSort.begin(), IndexSort.end(), PosScr.begin(), [&SpwOrdMat](size_t & Idx){
        double Scr = std::accumulate(SpwOrdMat[Idx].begin(),SpwOrdMat[Idx].end(),0.0,[](double A, std::pair<int,int> & B){
            if(B.second > 0){
                return A + (B.first/(double)B.second);    
            }
            return A;
        });
        return std::make_pair(std::move(Idx),(Scr/IndexSort.size()));
    });

    std::sort(SpwVec.begin(), SpwVec.end(), [](const spacedword & SpwA, const spacedword & SpwB){
        if(SpwA.family() == SpwB.family()){ //Family := BucketSize
            return SpwA < SpwB;
        }
        return SpwA.family() > SpwB.family();
    });

    return PosScr;
}


void clean_line(){
    std::vector<char> EmptyVector(80,' ');
    std::cout << "\r" << std::string(EmptyVector.begin(), EmptyVector.end()) << "\r";
}


void init_omp(){
    #ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(dboptions::ThreadNumber);
    #endif    
}


// bool End = false;
    // int Try = 3;
    // std::string Parse, Opt;
    // while(!End && Try > 0){
    //     if(Try == 3){
    //         create_rssw(Families, FamiliesSpacedWords);
    //     }
    //     std::cout << "Again? [Y/N]: ";
    //     std::cin >> Parse;
    //     if(Parse[0] == 'Y' || Parse[0] == 'y'){
    //         std::cout << "Used parameter set:" << std::endl;
    //         std::cout << "-c " << dboptions::SpwSeqThrshld << " | -n " << dboptions::SpwNoiseThrshld << " | -s " << dboptions::WordPerSeq;
    //         std::cout << " | -m " << dboptions::PatternNumber << " | -w " << dboptions::PatternWeight << " | -d " << dboptions::PatternMinDC;
    //         std::cout << "-" << dboptions::PatternMaxDC << " | -f " << dboptions::FamilyOffset << " | -t " << dboptions::ThreadNumber << std::endl;
    //         Try = 3;
    //         std::cout << "New Values? [Y/N]: ";
    //         std::cin >> Parse;
    //         if(Parse[0] == 'Y' || Parse[0] == 'y'){
    //             std::cout << "Use cmd line parameter pair, and 'END END' for stopping new value submission." << std::endl;
    //             std::cout << "Param: ";
    //             std::cin >> Parse >> Opt;
    //             while(Parse != "END"){
    //                 switch(Parse[1]){
    //                     case 'C':
    //                     case 'c':
    //                         dboptions::SpwSeqThrshld = std::stof(Opt);
    //                         break;
    //                     case 'D':
    //                     case 'd':
    //                         dboptions::parse_length(Opt.c_str());
    //                         break;
    //                     case 'F':
    //                     case 'f':
    //                         dboptions::FamilyOffset = std::stoi(Opt);
    //                         break;
    //                     case 'M':
    //                     case 'm':
    //                         dboptions::PatternNumber = std::stoi(Opt);
    //                         break;
    //                     case 'N':
    //                     case 'n':
    //                         dboptions::SpwNoiseThrshld = std::stof(Opt);
    //                         break;
    //                     case 'S':
    //                     case 's':
    //                         dboptions::WordPerSeq = std::stoi(Opt);
    //                         break;
    //                     case 'T':
    //                     case 't':
    //                         dboptions::ThreadNumber = std::stoi(Opt);
    //                         init_omp();
    //                         break;
    //                     case 'W':
    //                     case 'w':
    //                         dboptions::PatternWeight = std::stoi(Opt);
    //                         break;
    //                     default:
    //                         std::cerr << "Cannot parse '" << Parse << " " << Opt << "'!" << std::endl;
    //                 }
    //                 std::cout << "Param: ";
    //                 std::cin >> Parse >> Opt;
    //             }
    //         }
    //         else{
    //             std::cout << "No changes!" << std::endl;
    //         }
    //         std::cout << std::endl;
    //     }
    //     else if (Parse[0] == 'N' || Parse[0] == 'N'){
    //         End = true;
    //     }
    //     else{
    //         std::cerr << "Please submit 'Y' for yes or 'N' for no! Another try" << std::endl;
    //         Try--;
    //     }
    // }
    // if(Try == 0){
    //     std::cerr << "Three parsing attempts without correct understanding!" << std::endl;
    // }
    // std::cout << "Exit program!" << std::endl;