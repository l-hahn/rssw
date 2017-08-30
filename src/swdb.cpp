#include <cstdlib>
#include <iostream>
#include <string>
#include "createdatabase.hpp"
#include "dboptions.hpp"

void SecureMessage(std::string && errmsg, std::string && opt);

int main(int argc, char* argv[]){
    bool OutSet = false;
    if(argc < 2){
        SecureMessage("param",std::string(argv[0]));
        return 0;
    }
    else{
        bool Exit = false;
        std::string tmp;
        for(int i = 1; i < argc; i++){
            tmp = argv[i];
            if(tmp == "--version" || tmp == "-v" || tmp == "-V"){
                SecureMessage("version", "");
                Exit = true;
            }
            else if(tmp == "--help" || tmp == "-h" || tmp == "-H"){
                SecureMessage("param", argv[0]);
                Exit = true;
            }
        }
        if(Exit){
            return 0;
        }
        int ctr;
        for(ctr = 1; ctr < argc-1; ctr++){
            std::string tmp = argv[ctr];
            switch(tmp[1]){
                case 'C':
                case 'c':
                    if (ctr < argc - 2) {
                        ctr++;
                        dboptions::SpwSeqThrshld = std::atof(argv[ctr]);
                    }
                    break;
                case 'D':
                case 'd':
                    if (ctr < argc - 2) {
                        ctr++;
                        dboptions::parse_length(argv[ctr]);
                    }
                    break;
                case 'F':
                case 'f':
                    if (ctr < argc - 2) {
                        ctr++;
                        dboptions::FamilyOffset = std::atoi(argv[ctr]);
                    }
                    break;
                case 'M':
                case 'm':
                    if (ctr < argc - 2) {
                        ctr++;
                        dboptions::PatternNumber = std::atoi(argv[ctr]);
                    }
                    break;
                case 'N':
                case 'n':
                    if (ctr < argc - 2) {
                        ctr++;
                        dboptions::SpwNoiseThrshld = std::atof(argv[ctr]);
                    }
                    break;
                case 'R':
                case 'r':
                    dboptions::AlphaReduction = true;
                    break;
                case 'S':
                case 's':
                    if (ctr < argc - 2) {
                        ctr++;
                        dboptions::WordPerSeq = std::atoi(argv[ctr]);
                    }
                    break;
                case 'T':
                case 't':
                    if (ctr < argc - 2) {
                        ctr++;
                        dboptions::ThreadNumber = std::atoi(argv[ctr]);
                    }
                    break;
                case 'W':
                case 'w':
                    if (ctr < argc - 2) {
                        ctr++;
                        dboptions::PatternWeight = std::atoi(argv[ctr]);
                    }
                    break;
                case '-':
                    if(tmp == "--fast"){
                        dboptions::Greedy = true;
                    }
                    else if(tmp == "--pattern"){
                        if(ctr < argc-2){
                            ctr++;
                            dboptions::PatternFile = std::string(argv[ctr]);
                        }
                    }
                    else if(tmp == "--outfile"){
                        OutSet = true;
                        if(ctr < argc-2){
                            ctr++;
                            dboptions::OutSwDb = std::string(argv[ctr]);
                        }
                    }
                    else if(tmp == "--overlap"){
                        dboptions::Overlap = true;
                    }
                    else if(tmp == "--seqcov"){
                        dboptions::SeqCov = true;
                    }
                    break;
                default:
                    SecureMessage("arg",std::string(argv[ctr]));
            }
        }
        if(ctr == argc-1){
            dboptions::InSeqFam = std::string(argv[ctr++]);
            if(OutSet == false){
                auto Pos = dboptions::InSeqFam.find_last_of(".");
                if(Pos != std::string::npos){
                    dboptions::OutSwDb = dboptions::InSeqFam.substr(0,Pos);
                }
                else{
                    dboptions::OutSwDb = dboptions::InSeqFam;
                }
            }
            else{
                auto Pos = dboptions::OutSwDb.find_last_of(".");
                if(Pos != std::string::npos){
                    dboptions::OutSwDb = dboptions::OutSwDb.substr(0,Pos);
                }
            }
            create_database();
        }
        else{
            std::cerr << "FATAL ERROR!\nINPUT FILE missing!" << std::endl;
        }
    }
    return 0;
}

/*===Functions===============================================================*/
void SecureMessage(std::string && errmsg, std::string && opt){
    if(errmsg == "arg"){
        std::cerr << "Argument '" << opt << "' is unknown!" << std::endl;
        return;
    }
    else if (errmsg == "param") {
        std::cout << "Usage:\n\t" << opt << " [options] <SEQUENCE FILE>\n" << std::endl;
        std::cout << "Options:\n" << std::endl;
        std::cout << "\t\t -m <INT>: \t\t Number of patterns, default: n = " << dboptions::PatternNumber << ".\n" << std::endl;
        std::cout << "\t\t -w <INT>: \t\t Pattern weight, default: w = " << dboptions::PatternWeight << ".\n" << std::endl;
        std::cout << "\t\t -d <INT>: \t\t Number of don't care positions, default: d = " << dboptions::PatternMinDC << std::endl;
        std::cout << "\t\t    <INT>-<INT>: \t Min. and max. number of don't care positions.\n" << std::endl;
        std::cout << "\t\t -s <INT>: \t\t Number of minimal required spaced words per sequence, default s = " << dboptions::WordPerSeq << ".\n" << std::endl;
        std::cout << "\t\t -f <INT>: \t\t Number of minimal required sequences for a family, default f = " << dboptions::FamilyOffset << ".\n" << std::endl;
        std::cout << "\t\t -c <FLOAT>: \t\t Minimal sequence occurence value for a spaced word within a family, 0 <= c <= 1, default c = " << dboptions::SpwSeqThrshld << ".\n" << std::endl;
        std::cout << "\t\t -n <FLOAT>: \t\t Maximal sequence occurence value for a spaced word in other families, 0 <= n <= 1, default n = " << dboptions::SpwNoiseThrshld << ".\n" << std::endl;
        std::cout << "\t\t -r: \t\t\t Activate alphabet reduction for protein sequences.\n" << std::endl;
        std::cout << "\t\t -t <INT>: \t\t Number of CPU-Threads to be used, default: t = " << dboptions::ThreadNumber << ".\n" << std::endl;
        std::cout << "\t\t --outfile <FILE>: \t Save Database to <FILE-prefix>.swdb instead of <SEQUENCE FILE-prefix>.swdb.\n" << std::endl;
        std::cout << "\t\t --fast: \t\t Create database in a greedy, faster, way.\n" << std::endl;
        //std::cout << "\t\t --pattern <FILE>: \t Reading pattern from <File> in pattern format with {0,O,o,-,} = don't care and {1,X,x,#,*,} = match, seperated by ','|' '|'.'|';'|'\\n'|'\\t'.\n" << std::endl;
        std::cout << "\t=== Additional Parameters ====" << std::endl;
        std::cout << "\t\t --overlap: \t\t Calculate the overlap among the database, save into <OUTFILE>.ovlp.\n" << std::endl;
        std::cout << "\t\t --seqcov: \t\t Calculate the spaced words sequence coverage for each family, save into <OUTFILE>.sqcv.\n" << std::endl;
        std::cout << "\t\t --version: \t\t Print the program version.\n" << std::endl;
        std::cout << "\t\t --help: \t\t Print this help.\n" << std::endl;
        return;
    }
    else if (errmsg == "version") {
        std::cout << "SWDB, Version 0.40.4 - (c) 2017 Lars Hahn" << std::endl;
        std::cout << "This program is released under GPLv3." << std::endl << std::endl;
        std::cout << "This program calculates for a given database with sequence families a represen-\ntative set of spaced words (RSSW)." << std::endl;
        std::cout << "The RSSW file output can be used with the search tool SWDS for classifying se-\nquences according to the RSSW." << std::endl;
        return;
    }
}