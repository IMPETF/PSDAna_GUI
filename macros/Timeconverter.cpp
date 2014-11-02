#include <iostream>
#include <ctime>
#include <string>

using namespace std;

unsigned long StringtoTC(string s){	//format of 30/09/2014 18:18:18 is  "20140930181818"
	time_t rawtime, basetime;
	struct tm * timeinfo;

	time(&basetime);

	//construct the base time from 01/01/2013 00:00:00
	timeinfo = localtime(&basetime);
	timeinfo->tm_year = 2013 - 1900;
	timeinfo->tm_mon = 0;
	timeinfo->tm_mday = 1;
	timeinfo->tm_hour = 0;
	timeinfo->tm_min = 0;
	timeinfo->tm_sec = 0;
	basetime = mktime(timeinfo);

	//construct the real time from string s
	timeinfo->tm_year = stoi(s.substr(0, 4)) - 1900;
	timeinfo->tm_mon = stoi(s.substr(4, 2)) - 1;
	timeinfo->tm_mday = stoi(s.substr(6, 2));
	timeinfo->tm_hour = stoi(s.substr(8, 2));
	timeinfo->tm_min = stoi(s.substr(10, 2));
	timeinfo->tm_sec = stoi(s.substr(12, 2));
	rawtime = mktime(timeinfo);

	return rawtime-basetime;
}

string TCtostring(unsigned long t){
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[15];

	time(&rawtime);

	//construct the base time from 01/01/2013 00:00:00
	timeinfo = localtime(&rawtime);
	timeinfo->tm_year = 2013 - 1900;
	timeinfo->tm_mon = 0;
	timeinfo->tm_mday = 1;
	timeinfo->tm_hour = 0;
	timeinfo->tm_min = 0;
	timeinfo->tm_sec = 0;
	rawtime = mktime(timeinfo);

	//add seconds to base time 
	rawtime += t;

	//return to struct tm and string
	timeinfo = localtime(&rawtime);
	strftime(buffer, 15, "%Y%m%d%H%M%S", timeinfo);

	return string(buffer);
}

int main()
{
	unsigned long timecode = StringtoTC("20141011123456");
	cout << "The time code of (2014-10-11 12:34:56) is: " << timecode << endl;
	cout << "In hex format is: 0X" << std::hex << timecode << endl;

	cout << "The real time of 0X"<< std::hex << timecode << " is: " << TCtostring(timecode) << endl;

	return 0;
}