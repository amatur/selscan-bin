#include <string>
#include <fstream>
#include <iostream>

using namespace std;

#ifndef LOGGER_H
#define LOGGER_H

class Logger{

    public:
            static void open( const string & logFile);
            static void close();
            // write message
            static void write( const string & message);
            ofstream log;

    private:
          Logger();
          
          //Logger instance (singleton)
          static Logger instance;
};

#endif