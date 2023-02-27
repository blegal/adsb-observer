// Demo:   maps.cpp
// Author: Evan Pezent (evanpezent.com)
// Date:   7/25/2021
#ifndef _GpsPositionsRefresh_
#define _GpsPositionsRefresh_

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

class GpsPositionsRefresh
{

private:

    std::ifstream ifile;
    int iters = 0;


public:


    GpsPositionsRefresh() : ifile("avions.gps")
    {
        if( !ifile.is_open() )
        {
            exit( EXIT_FAILURE );
        }
    }


    void update()
    {
        iters += 1;             // On simule une mise à jour des données toutes les 1s
        if( iters%60 != 0 )     // si ce n'est pas le cas, on attend ;-)
        {
            return;
        }

        int cnt = 0;
        std::string line, lat, lon, pidx;
        while (std::getline(ifile, line))
        {
            std::istringstream iss(line);   // string stream
            std::getline(iss, lat, ' ');    //
            
            if(lat.size() == 0)              // un espace debute la ligne !
                std::getline(iss, lat, ' '); //

            std::getline(iss, lon,  ' ');    //
            std::getline(iss, pidx, ' ');    //

            const float f_lon = conv_lon_to_x( std::stof(lon) );    // on convertit la lattitude et la longitude
            const float f_lat = conv_lat_to_y( std::stof(lat) );    // dans des coordonnées compatible 2D
            const int p_id    =                std::stoi(pidx) ;

            pos_x.push_back( f_lon );   // On ajoute les données dans la liste
            pos_y.push_back( f_lat );   // On ajoute les données dans la liste
            idx.push_back  ( p_id  );   // On ajoute les données dans la liste
            
            cnt += 1;                   // On a lu 32 données, on passe donc la main
            if( cnt == 16 )             // pour simuler la réception séquentielle des
                break;                  // trames lors d'une reception réelle
        }

        utcExample();
        std::cout <<  "Data were loaded (" <<  pos_x.size() << ") coordonates" << std::endl;
    }
    

    ~GpsPositionsRefresh()
    {
        ifile.close();
        pos_x.clear();
        pos_y.clear();
        idx.clear();
    }


    std::vector<float> pos_x;
    std::vector<float> pos_y;
    std::vector<  int> idx;
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
