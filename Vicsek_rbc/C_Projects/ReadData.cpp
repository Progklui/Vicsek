#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

class ReadData {
private:
    double *dataArray;
    int numberLines;

public:
    ReadData(const string &file) {
        fstream newfile;

        //assuming that the file is within the directory "Parameter"
        newfile.open("../Parameter/" + file, ios::in);

        this->numberLines = countLines(newfile);
        
        this->dataArray = new double[numberLines - 1];

        if (newfile.is_open()) {
            string tp;
            for (int i = 0; i < numberLines; i++) {
                getline(newfile, tp);
                if (i == 0) {
                    continue;
                }
                
                std::stringstream temp(tp);
                std::string segment;
                std::vector<std::string> seglist;

                while (std::getline(temp, segment, ';')) {
                    seglist.push_back(segment); //Spit string at ';' character
                }

                dataArray[i - 1] = stod(seglist[1]);
            }
            newfile.close();
        }

    }

    double *getData();
    
    int numberInput();


private:
    int countLines(fstream &stream);
};

int ReadData::countLines(fstream &stream) {
    int counter = 0;
    if (stream.is_open()) {
        string tp;
        while (getline(stream, tp)) {
            counter++;
        }
    }

    stream.clear();
    stream.seekg(0);

    return counter;
}

double *ReadData::getData() {
    return dataArray;
}

int ReadData::numberInput() {
    return numberLines - 1;
}
