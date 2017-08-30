#include <iostream>
#include <string>
#include "dsoptions.hpp"
#include "searchdatabase.hpp"


void SecureMessage(std::string && errmsg, std::string && opt);

int main(int argc, char* argv[]){
    bool OutSet = false;
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
    if(argc < 3){
        SecureMessage("param",std::string(argv[0]));
        return 0;
    }
    else{
        int ctr;
        for(ctr = 1; ctr < argc-2; ctr++){
            std::string tmp = argv[ctr];
            switch(tmp[1]){
                case 'B':
                case 'b':
                    if (ctr < argc - 3) {
                        ctr++;
                        dsoptions::BlockSize = std::atoi(argv[ctr]);
                    }
                    break;
                case 'F':
                case 'f':
                    if (ctr < argc - 3) {
                        ctr++;
                        dsoptions::FamilyOffset = std::atoi(argv[ctr]);
                    }
                    break;
                case 'H':
                case 'h':
                    if (ctr < argc - 3) {
                        ctr++;
                        dsoptions::HitThreshold = std::atof(argv[ctr]);
                    }
                    break;
                case 'R':
                case 'r':
                    dsoptions::AlphaReduction = true;
                    break;
                case 'S':
                case 's':
                    if (ctr < argc - 3) {
                        ctr++;
                        dsoptions::WordPerSeq = std::atoi(argv[ctr]);
                    }
                    break;
                case 'T':
                case 't':
                    if (ctr < argc - 3) {
                        ctr++;
                        dsoptions::ThreadNumber = std::atoi(argv[ctr]);
                    }
                    break;
                case '-':
                    if(tmp == "--help"){
                        SecureMessage("param",std::string(argv[0]));
                        return 0;
                    }
                    else if(tmp == "--version"){
                        SecureMessage("version","");
                        return 0;
                    }
                    else if(tmp == "--outfile"){
                        if (ctr < argc - 3) {
                            OutSet = true;
                            ctr++;
                            dsoptions::OutDetect = std::string(argv[ctr]);
                        }
                    }
                    break;
                default:
                    SecureMessage("arg",std::string(argv[ctr]));
            }
        }
        if(ctr == argc-2){
            dsoptions::InSwDb = std::string(argv[ctr++]);
            dsoptions::InSeq = std::string(argv[ctr++]);
            if(OutSet == false){
                auto Pos = dsoptions::InSeq.find_last_of(".");
                if(Pos != std::string::npos){
                    dsoptions::OutDetect = dsoptions::InSeq.substr(0,Pos);
                }
                else{
                    dsoptions::OutDetect = dsoptions::InSeq;
                }
            }
            else{
                auto Pos = dsoptions::OutDetect.find_last_of(".");
                if(Pos != std::string::npos){
                    dsoptions::OutDetect = dsoptions::OutDetect.substr(0,Pos);
                }
            }
            search_database();
        }
        else{
            std::cerr << "FATAL ERROR!\nINPUT FILES missing!" << std::endl;
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
        std::cout << "Usage:\n\t" << opt << " [options] <SWDB FILE> <SEQUENCE FILE>\n" << std::endl;
        std::cout << "Options:\n" << std::endl;
        std::cout << "\t\t -h <FLOAT>: \t\t Minimal threshold hit score for detection, default h = " << dsoptions::HitThreshold << ".\n" << std::endl;
        std::cout << "\t\t -b <INT>: \t\t Window size used for detection, default b = 3500 | 1167 (DNA|Protein).\n" << std::endl;
        std::cout << "\t\t -s <INT>: \t\t Number of minimal required spaced words per sequence, default s = " << dsoptions::WordPerSeq << ".\n" << std::endl;
        std::cout << "\t\t -f <INT>: \t\t Number of minimal required sequences for a family, default f = " << dsoptions::FamilyOffset <<".\n" << std::endl;
        std::cout << "\t\t -r: \t\t\t Activate alphabet reduction for protein sequences.\n" << std::endl;
        std::cout << "\t\t -t <INT>: \t\t Number of CPU-Threads to be used, default: t = " << dsoptions::ThreadNumber << ".\n" << std::endl;
        std::cout << "\t\t --outfile <FILE>: \t Save detection to <FILE-prefix>.swds instead of <SEQUENCE FILE-prefix>.swds.\n" << std::endl;
        std::cout << "\t=== Additional Parameters ====" << std::endl;
        std::cout << "\t\t --version: \t\t Print the program version.\n" << std::endl;
        std::cout << "\t\t --help: \t\t Print this help.\n" << std::endl;
        return;
    }
    else if (errmsg == "version") {
        std::cout << "SWDS, Version 0.40.4 - (c) 2017 Lars Hahn" << std::endl;
        std::cout << "This program is released under GPLv3." << std::endl << std::endl;
        std::cout << "This program classifies a set of sequences/sequence-families with a given repre-\nsentative set of spaced words(RSSW), created by the SWDB tool." << std::endl;
        std::cout << "The output file contains for sequences an RSSW, where it was classified to, and \nadditional information." << std::endl;
        return;
    }
}