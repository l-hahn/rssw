#include "searchdatabase.hpp"

std::vector<size_t> FamSizeUsed;

void search_database(){
    if(dsoptions::InSwDb.size() == 0 || dsoptions::InSeq.size() == 0){
        std::cerr << "FATAL ERROR!\nINPUT FILES missing!" << std::endl;
        std::exit(-1);
    }
    std::vector<sequence_family> Families;
    read_family(Families);
    
    spacedword_db SwDataBase(Families[0].is_protein());
    read_database(SwDataBase, Families);

    FamSizeUsed = std::vector<size_t>(Families.size(),0);

    init_omp();
    multi_files(Families, SwDataBase);
    //matching_db_seq(Families, SwDataBase);
}

void read_family(std::vector<sequence_family> & Families){
    if(dsoptions::AlphaReduction == false){
        std::string Parser;
        std::ifstream Header(dsoptions::InSwDb);
        std::getline(Header, Parser);
        if(Parser[Parser.size()-3] == 'R' || Parser[Parser.size()-3] == 'r'){
            dsoptions::AlphaReduction = true;
        }
        Header.close();
    }
    Families = sequence_family::families_file(dsoptions::InSeq, dsoptions::AlphaReduction);
    dsoptions::FamilySize = Families.size();
    if(Families.size() < 2){
        if(Families[0].size() < 2){
            std::cout << Families[0].size() << " Sequence was read. Inputtype: ";
        }
        else{
            std::cout << Families[0].size() << " Sequences were read. Inputtype: ";
        }
    }
    else if(Families.size() == 0){
        std::cerr << "FATAL ERROR!\nNo sequence input!" << std::endl;
        std::exit(-1);
    }
    else{
        std::cout << Families.size() << " Families were read. Inputtype: ";
    }
    if(Families[0].is_protein()){
        std::cout << "Protein " << std::endl;
    }
    else{
        std::cout << "DNA" << std::endl;
    }
    std::cout << std::endl;
}

void read_database(spacedword_db & SwDataBase, std::vector<sequence_family> & Families){
    std::cout << "Reading database '" << dsoptions::InSwDb << "' ... " << std::flush;
    SwDataBase = spacedword_db(dsoptions::InSwDb);
    std::cout << "\n--> Done!" << std::endl;
    if(dsoptions::WordPerSeq < 0){
        dsoptions::WordPerSeq = SwDataBase.words_per_seq();
    }
    bool IsProtein = true;
    if(Families[0].is_protein()){
        if(dsoptions::BlockSize == -1){
            dsoptions::BlockSize = 1167; //DNA default BlockSize divided by 3, codon size.
        }
    }
    else{
        IsProtein = false;
        if(dsoptions::BlockSize == -1){
            dsoptions::BlockSize = 3500; //Typical DNA gene size.
        }
    }
    std::cout << "\n\tBlockSize = " << dsoptions::BlockSize << std::endl;
    std::cout << "\tThreshold = " << dsoptions::HitThreshold << std::endl;
    std::cout << "\tWords per Seq = " << dsoptions::WordPerSeq << std::endl << std::endl;

    if(SwDataBase.is_dna() && IsProtein){
        std::cerr << "TYPE ERROR!\nDataBase Type is DNA, Sequences Type is Protein! Cannot convert SpacedWord to Protein!" << std::endl;
        std::exit(-1);
    }
    if(SwDataBase.is_protein() && IsProtein == false){
        std::cout << "DataBase Type is Protein, Sequences Type is DNA, translating sequences to protein ... " << std::flush;
            for(auto & Fam : Families){
                Fam.translate();
            }
        std::cout << "\n--> Done!" << std::endl;
        dsoptions::Translated = true;
    }
}


void multi_files(std::vector<sequence_family> & Families, spacedword_db & SwDataBase){
    set_multi_file_name();
    bool End = false;
    int Try = 3;
    std::string Parse, Opt;
    while(!End && Try > 0){
        if(Try == 3){
            matching_db_seq(Families, SwDataBase);
        }
        std::cout << "Again with new file(/s)? [Y/N]: ";
        std::cin >> Parse;
        if(Parse[0] == 'Y' || Parse[0] == 'y'){
            Try = 3;
            std::cout << "Use '-f <FAMILY-FILE> ' for new Family-File, -d <DB-FILE> for new spw DB-File and 'END END' for stopping new value submission." << std::endl;
            std::cout << "Param: ";
            std::cin >> Parse >> Opt;
            while(Parse != "END"){
                switch(Parse[1]){
                    case 'D':
                    case 'd':
                        dsoptions::InSwDb = std::string(Opt);
                        set_multi_file_name();
                        read_database(SwDataBase, Families);
                        break;
                    case 'F':
                    case 'f':
                        dsoptions::InSeq = std::string(Opt);
                        set_multi_file_name();
                        read_family(Families);
                        break;
                    default:
                        std::cerr << "Cannot parse '" << Parse << " " << Opt << "'!" << std::endl;
                }
                std::cout << "Param: ";
                std::cin >> Parse >> Opt;
            }
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

void set_multi_file_name(){
    auto PosEnd = dsoptions::InSeq.find_last_of(".");
    std::string TmpSeq, TmpDB;
    if(PosEnd != std::string::npos){
        TmpSeq = dsoptions::InSeq.substr(0,PosEnd);
    }
    else{
        TmpSeq = dsoptions::InSeq;
    }
    PosEnd = dsoptions::InSwDb.find_last_of(".");
    auto PosBegin = dsoptions::InSwDb.find_last_of("/");
    auto PosTmp = dsoptions::InSwDb.find_last_of("\\");
    if(PosBegin >= PosEnd){
        PosEnd = dsoptions::InSwDb.size()-1;
    }
    if(PosBegin != std::string::npos && PosTmp != std::string::npos){
        PosBegin = std::max(PosBegin, PosTmp);
    }
    else{
        PosBegin = std::min(PosBegin, PosTmp);
    }
    if(PosBegin != std::string::npos && PosEnd != std::string::npos){
        TmpDB = dsoptions::InSwDb.substr(PosBegin+1, PosEnd-PosBegin-1);
    }
    else if(PosBegin == std::string::npos && PosEnd != std::string::npos){
        TmpDB = dsoptions::InSwDb.substr(0, PosEnd);
    }
    else if(PosBegin != std::string::npos && PosEnd == std::string::npos){
        TmpDB = dsoptions::InSwDb.substr(PosBegin+1, dsoptions::InSwDb.size()-1-PosBegin-1);
    }
    else{
        TmpDB = dsoptions::InSwDb;
    }
    dsoptions::OutDetect = TmpSeq+"-=-"+TmpDB;
}

void matching_db_seq(std::vector<sequence_family> & Families, spacedword_db & SwDataBase){
    std::vector<spacedhit> DataBaseHit;
    std::cout << "Matching spaced words database to input sequences  ... " << std::endl;
    for(unsigned i = 0; i < Families.size(); i++){
        if(Families[i].size() > 1){
            std::cout << "\rProcessing Family " << i+1 << "/" << Families.size() <<  " " << Families[i].id() << std::flush;
            std::vector<spacedhit> SeqHits = family_hit(Families[i], i, SwDataBase);
            DataBaseHit.insert(DataBaseHit.end(), SeqHits.begin(), SeqHits.end());
        }
    }
    std::sort(DataBaseHit.begin(), DataBaseHit.end());
    std::cout << "\r" << std::string(80,' ') << "\r--> Done!" << std::endl;

    classify(SwDataBase, Families, DataBaseHit);

    std::ofstream Output;
    std::cout << "\nWriting results into file '" << dsoptions::OutDetect << ".swds' ... " << std::flush; 
    Output.open(dsoptions::OutDetect+".swds");
    Output << "#FamName\tSeqName\tORF\tMatch DbID\tPosition\tScore" << std::endl;
    for(auto SwHit : DataBaseHit){
        unsigned ORF = SwHit.orf();
        std::string OrfPrint = "+";
        if(ORF%2 == 1){
            OrfPrint[0] = '-';
        }
        OrfPrint.push_back((ORF/2)+1+'0');
        Output << Families[SwHit.fam_id()].id() << "\t" << Families[SwHit.fam_id()][SwHit.seq_num()].name() << "\t" << OrfPrint << "\t" << SwDataBase.id_fam_name(SwHit.db_fam_id()) << "\t" << SwHit.position() << "\t" << SwHit.score() << std::endl;
    }
    Output.close();
    std::cout << "\n--> Done!" << std::endl;

    evaluate_sens_spec(SwDataBase, Families, DataBaseHit);
}


std::vector<spacedhit> family_hit(sequence_family & Family, unsigned FamID, spacedword_db & SwDataBase){
    std::vector<spacedhit> SwMatches;
    for(auto & DbSw : SwDataBase){
        std::vector<spacedhit> pattern_hits = pattern_hit(Family, FamID, DbSw);
        SwMatches.insert(SwMatches.end(), pattern_hits.begin(), pattern_hits.end());
    }

    std::sort(SwMatches.begin(), SwMatches.end(), [](const spacedhit & SwHA, const spacedhit & SwHB){
        if(SwHA.db_fam_id() == SwHB.db_fam_id()){
            if(SwHA.seq_num() == SwHB.seq_num()){
                return SwHA.position() < SwHB.position();
            }
            return SwHA.seq_num() < SwHB.seq_num();
        }
        return SwHA.db_fam_id() < SwHB.db_fam_id();
    });

    std::vector<spacedhit> TrueSwMatches;
    auto BucketStart = SwMatches.begin(), BucketIter = SwMatches.begin();
    bool NotEnd = true, BlockDist = false;
    if(SwMatches.size() != 0){
        while(NotEnd){
            if(BucketIter == SwMatches.end()){
                NotEnd = false;
                BlockDist = true;
            }
            else{
                if(BucketStart->db_fam_id() != BucketIter->db_fam_id()){
                    BlockDist = true;
                }
                else if(BucketStart->seq_num() != BucketIter->seq_num()){
                    BlockDist = true;
                }
                else if((BucketIter->position() - BucketStart->position()) > (unsigned)dsoptions::BlockSize){
                    BlockDist = true;
                }
            }
            if(BlockDist){
                double BlockScore = BucketStart->score();
                unsigned MatchCount = 1;
                for(auto BuckIt = BucketStart + 1; BuckIt < BucketIter; BuckIt++){
                    if(BuckIt->fam_id() > (BuckIt-1)->fam_id()){
                        BlockScore += BuckIt->score();
                        MatchCount++;
                    }
                }
                if(MatchCount >= (unsigned)dsoptions::WordPerSeq && BlockScore >= dsoptions::HitThreshold){
                    unsigned short SeqOrf = 0;
                    size_t SeqID = BucketStart->seq_num();
                    if(dsoptions::Translated){
                        SeqOrf = SeqID % 6;
                        SeqID /= 6;
                    }
                    size_t Pos = ((BucketIter - 1)->position() - BucketStart->position())/2 + BucketStart->position();
                    TrueSwMatches.push_back(spacedhit(FamID, SeqID, SeqOrf, BucketStart->db_fam_id(), Pos , BlockScore));
                }
                BucketStart = BucketIter;
                BlockDist = false;
            }
            if(NotEnd){
                BucketIter++;
            }
        }
        SwMatches = std::move(TrueSwMatches);
    }
    return SwMatches;
}



std::vector<spacedhit> pattern_hit(sequence_family & Family, unsigned FamID, rsw_set & DbSw){
    std::vector<spacedword> SpacedWords;
    pattern Pat = DbSw.pat();
    size_t Ctr = 0;
    for(unsigned i = 0; i < Family.size(); i++){
        std::vector<spacedword> TmpSw;
        TmpSw = Family[i].spaced_words(Pat, i);
        if(TmpSw.size() > 0){
            TmpSw.erase(std::unique(TmpSw.begin(), TmpSw.end()), TmpSw.end());
            SpacedWords.insert(SpacedWords.end(), TmpSw.begin(), TmpSw.end());
            Ctr++;
        }
    }
    FamSizeUsed[FamID] = Ctr;
    if(Ctr < dsoptions::FamilyOffset){
        FamSizeUsed[FamID] = 0;
        return std::vector<spacedhit>();
    }
    std::sort(SpacedWords.begin(), SpacedWords.end(), [](const spacedword & SwA, const spacedword & SwB){
        return SwA.bits() < SwB.bits();
    });

    std::vector<spacedhit> SwMatches;
    unsigned & NumThrd = dsoptions::ThreadNumber;
    unsigned ThrdSize = (size_t)std::ceil((double)DbSw.size()/NumThrd);
    #pragma omp parallel
    {
        std::vector<spacedhit> LocalSwMatches;
        unsigned ThrdNum = omp_get_thread_num();
        unsigned TPrefMax = ThrdSize*ThrdNum + (ThrdSize - 1);
        if(ThrdSize*ThrdNum < DbSw.size()){
            if(TPrefMax >= DbSw.size()){
                TPrefMax = DbSw.size()-1;
            }
            auto DbSpwMin = DbSw.begin()+ThrdSize*ThrdNum;
            auto DbSpwMax = DbSw.begin()+TPrefMax+1;
            auto SpwMin = SpacedWords.begin();
            auto SpwMax = SpacedWords.end();

            while(SpwMin != SpwMax && DbSpwMin != DbSpwMax){
                if(SpwMin->bits() < DbSpwMin->bits()){
                    SpwMin++;
                }
                else if(SpwMin->bits() > DbSpwMin->bits()){
                    DbSpwMin++;
                }
                else{
                    for(auto & FamScore : *DbSpwMin){
                        int64_t Pos = SpwMin->position();
                        if(dsoptions::Translated){
                            Pos *= 3; //Protein -> DNA retranslation for positions, if origin was DNA!
                            Pos += 1;
                        }
                        LocalSwMatches.push_back(spacedhit(FamScore.position(), SpwMin->sequence(), 0, FamScore.id(), Pos, FamScore.score()));
                    }
                    SpwMin++;
                }
            }
            #pragma omp critical
            {
                SwMatches.insert(SwMatches.end(), LocalSwMatches.begin(), LocalSwMatches.end());
            }
        }
    }
    return SwMatches;
}

void classify(spacedword_db & SwDataBase, std::vector<sequence_family> & Families, std::vector<spacedhit> & DataBaseHit){
    std::vector<spacedhit> DataBaseClassify, Tmp;
    auto BucketStart = DataBaseHit.begin(), BucketIter = BucketStart;
    bool NotEnd = true;
    bool Border = false;
    while(NotEnd){
        if(BucketIter == DataBaseHit.end()){
            NotEnd = false;
            Border = true;
        }
        else{
            if(BucketIter->fam_id() != BucketStart->fam_id()){
                Border = true;
            }
            else if(BucketIter->seq_num() != BucketStart->seq_num()){
                Border = true;
            }
        }
        if(Border){
            auto MaxIter = BucketStart;
            for(auto BIter = BucketStart; BIter != BucketIter; BIter++){
                if(MaxIter->score() <= BIter->score()){
                    if(MaxIter->score() < BIter->score()){
                        Tmp.clear();
                    }
                    MaxIter = BIter;
                    Tmp.push_back(*MaxIter);
                }
            }
            DataBaseClassify.insert(DataBaseClassify.end(), Tmp.begin(), Tmp.end());
            Border = false;
            BucketStart = BucketIter;
            Tmp.clear();
        }
        BucketIter++;
    }
    std::swap(DataBaseClassify, DataBaseHit);
}


void init_omp(){
    #ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(dsoptions::ThreadNumber);
    #endif    
}
















void evaluate_sens_spec(spacedword_db & SwDataBase, std::vector<sequence_family> & Families, std::vector<spacedhit> & DataBaseHit){
    std::map<size_t, size_t> DbIdToFamId;
    for(size_t i = 0; i < SwDataBase.families_size(); i++){
        std::string FamName = SwDataBase.id_fam_name(i);
        auto Pos = std::find_if(Families.begin(), Families.end(), [FamName](const sequence_family & Fam){
            return FamName == Fam.id();
        });
        if(Pos != Families.end()){
            DbIdToFamId.insert(std::pair<size_t,size_t>(i,std::distance(Families.begin(), Pos)));
        }
    }

    std::vector<size_t> TP_Count(Families.size(), 0), FP_Count(Families.size(), 0) , FN_Count(Families.size(), 0), TN_Count(Families.size(), 0);
    for(auto & Hit : DataBaseHit){
        if(DbIdToFamId[Hit.db_fam_id()] == Hit.fam_id()){
            TP_Count[Hit.fam_id()]++;
        }
        else{
            FP_Count[DbIdToFamId[Hit.db_fam_id()]]++;
        }
    }
    for(unsigned i = 0; i < Families.size(); i++){
        if(FamSizeUsed[i] > 0){
            FN_Count[i] = FamSizeUsed[i] - TP_Count[i];
            TN_Count[i] = DataBaseHit.size() - TP_Count[i] - FP_Count[i] - FN_Count[i];
        }
    }

    std::ofstream Output(dsoptions::OutDetect+".eval");
    Output << std::fixed << std::setprecision(5) << std::endl; 
    Output << "#PfamID\tSensitivity\tSpecifity\tTP\t\tFP\t\tTN\t\tFN\t\tSeqSize" << std::endl;
    size_t TP = 0, FP = 0, FN = 0, TN = 0;
    for(unsigned i = 0; i < Families.size(); i++){
        if(FamSizeUsed[i] == 0){
            continue;
        }
        Output << Families[i].id();
        if(TP_Count[i] + FN_Count[i] > 0){
            Output << "\t" << TP_Count[i]/((double)TP_Count[i]+FN_Count[i]);
        }
        else{
            Output << "\t- "; 
        }
        if(TP_Count[i] + FP_Count[i] > 0){
            Output << "\t" << TP_Count[i]/((double)TP_Count[i]+FP_Count[i]);
        }
        else{
            Output << "\t- ";
        }
        Output << "\t" << TP_Count[i] << "\t\t" << FP_Count[i] << "\t\t" << TN_Count[i] << "\t\t" << FN_Count[i] << "\t\t" << FamSizeUsed[i] << std::endl;
        TP += TP_Count[i];
        FP += FP_Count[i];
        FN += FN_Count[i];
        TN += TN_Count[i];
    }
    Output << "\n\nTotal\tSensitivity\tSpecifity" << std::endl;
    Output << "Pfam";
    if(TP + FN > 0){
        Output << "\t" << TP/((double)TP+FN);
    }
    else{
        Output << "\t- "; 
    }
    if(TP + FP > 0){
        Output << "\t" << TP/((double)TP+FP);
    }
    else{
        Output << "\t- ";
    }
    Output << std::endl << std::endl;
    Output << "Sens = TP/(TP+FN)" << std::endl;
    Output << "Spec = TP/(TP+FP)" << std::endl;
    Output.close();   
}