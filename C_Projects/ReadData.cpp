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

                if (seglist[0] == "NumberParticles") {
                    this->numberOfParticles = stoi(seglist[1]);
                } else if (seglist[0] == "Length") {
                    this->length = stod(seglist[1]);
                } else if (seglist[0] == "TimeStep") {
                    this->timeStep = stod(seglist[1]);
                } else if (seglist[0] == "Radius") {
                    this->radius = stod(seglist[1]);
                } else if (seglist[0] == "Speed") {
                    this->speed = stod(seglist[1]);
                } else if (seglist[0] == "eta") {
                    this->eta = stod(seglist[1]);
                }

                this->strings[i - 1] = tp;
            }
            newfile.close();
        }

    }

    int getNumberParticles();

    double getLength();

    double getTimeStep();

    double getRadius();

    double getSpeed();

    double getEta();


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

int ReadData::getNumberParticles() {
    return numberOfParticles;
}

double ReadData::getLength() {
    return length;
}

double ReadData::getTimeStep() {
    return timeStep;
}

double ReadData::getRadius() {
    return radius;
}

double ReadData::getSpeed() {
    return speed;
}

double ReadData::getEta() {
    return eta;
}