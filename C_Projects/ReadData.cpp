#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <sstream>

using namespace std;

class ReadData {
private:
    string *strings;
    long numberOfParticles;
    long length;
    double timeStep;
    double radius;
    double speed;
    double eta;

public:
    ReadData(const string &file) {
        fstream newfile;

        //assuming that the file is within the directory "Parameter"
        newfile.open("Parameter/" + file, ios::in);

        int numberLines = countLines(newfile);

        this->strings = new string[numberLines - 1];

        if (newfile.is_open()) {
            string tp;
            int counter = 0;
            for (int i = 0; i < numberLines; i++) {
                getline(newfile, tp);
                if (i == 0) {
                    continue;
                }
                this->strings[i - 1] = tp;
            }
            newfile.close();
        }

    }


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

int main(void) {

    //just testing some things

    ReadData newData = *new ReadData("input.csv");

    /*std::string s = "scott>=tiger";
    std::string delimiter = ">=";
    std::string token = s.substr(7, s.find(delimiter)); // token is "scott"

    std::cout << token << "\n";*/

    std::stringstream test("this_is_a_test_string");
    std::string segment;
    std::vector<std::string> seglist;

    while (std::getline(test, segment, '_')) {
        seglist.push_back(segment); //Spit string at '_' character
    }

    for (int i = 0; i < seglist.size(); i++) {
        std::cout << seglist[i] << "\n";
    }


    return EXIT_SUCCESS;
}