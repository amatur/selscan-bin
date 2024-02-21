#include "logger.hpp"

Logger Logger::instance;

Logger::Logger(){}

void Logger::open( const string& logFile){
   instance.log.open(logFile.c_str());
}
void Logger::close(){
    instance.log.close();

} 
void Logger::write(const string& message){
    ostream& stream =  instance.log ;
    stream << message<< endl;
}