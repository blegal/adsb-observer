// Demo:   maps.cpp
// Author: Evan Pezent (evanpezent.com)
// Date:   7/25/2021
#ifndef _GpsPositions_
#define _GpsPositions_

#include <string>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <locale>

//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//

#ifndef _CONV_
#define _CONV_

inline void utcExample()
{
    std::time_t time = std::time({});
    char timeString[std::size("yyyy-mm-dd hh:mm:ss")];
    std::strftime(std::data(timeString), std::size(timeString), "[%T] ", std::gmtime(&time));
    std::cout << timeString;
}

double conv_lon_to_x(const double lon)
{
    return ((lon + 180.0) / 360.0);
}

double conv_lat_to_y(const double lat) {
    double latrad = lat * PI / 180.0;
    return (1.0 - asinh(tan(latrad)) / PI) / 2.0;
}
#endif

//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//

class GpsPositions
{
private:
    
public:
    GpsPositions()
    {
        std::ifstream ifile("maps.gps");
        if( !ifile.is_open() )
        {
            exit( EXIT_FAILURE );
        }

        std::string line, lat, lon;
        while (std::getline(ifile, line))
        {
            std::istringstream iss(line); // string stream
            std::getline(iss, lat, ' ');
            std::getline(iss, lon, ' ');
            float f_lon = conv_lon_to_x( std::stof(lon) );
            float f_lat = conv_lat_to_y( std::stof(lat) );
            pos_x.push_back( f_lon );
            pos_y.push_back( f_lat );
        }
        
        ifile.close();
    }
    
    ~GpsPositions()
    {
        pos_x.clear();
        pos_y.clear();
    }

    std::vector<float> pos_x;
    std::vector<float> pos_y;
};

//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
#endif
