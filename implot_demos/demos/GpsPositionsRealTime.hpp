// Demo:   maps.cpp
// Author: Evan Pezent (evanpezent.com)
// Date:   7/25/2021
#ifndef _GpsPositionsRealTime_
#define _GpsPositionsRealTime_

#include <string>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <ctime>
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

inline double conv_lon_to_x(const double lon)
{
    return ((lon + 180.0) / 360.0);
}

inline double conv_lat_to_y(const double lat) {
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

class GpsPositionsRealTime
{
private:
    std::filesystem::file_time_type ftime;

public:

    GpsPositionsRealTime()
    {

    }
    
    void update()
    {

#ifdef __APPLE__
        std::string filename = "/tmp/avions.gps";
//        FILE* f = fopen("/Volumes/RAMDisk/frame.ech.cf32", "w");
#else
        std::string filename = "/tmp/avions.gps";
        //FILE* f = fopen("/tmp/frame.ech.cf32", "w");
#endif

        //
        //  Si il n'y a pas de fichier dans /tmp cela signifie que
        // le recepteur ADSB n'est pas lancé OU qu'il ne detecte
        // pas d'avion, on passe en mode demo
        //

        if( std::filesystem::exists( filename ) == false )
        {
            filename = "../data/avions.gps";

            if( std::filesystem::exists( filename ) == false )
            {
                std::cout << "(EE) The ADSB receiver is not currentlu working..." << std::endl;
                std::cout << "(EE) and their is not demo file available ?!"       << std::endl;
                exit( EXIT_FAILURE );
            }

        }

        //
        // On regarde si le fichier a été mis à jour depuis notre dernier passage
        // si oui, on va charger les valeurs qui sont dedans, sinon on ne fait rien...
        //
        
        //std::cout << "(DD) Testing the update of (" <<  filename << ")" << std::endl;

        auto curTime = std::filesystem::last_write_time( filename );
        if( curTime == ftime )
        {
            return;
        }

        pos_x.clear();
        pos_y.clear();
        idx.clear();

//      std::cout << "(DD) The file was updated (" <<  filename << ")" << std::endl;

        //
        // On va parser le fichier pour extraire les position (lat, long) des avions
        // un couple de valeur par ligne formaté en "%f %f"
        //
        int nLines = 0;
        std::ifstream ifile( filename );
        std::string line, lat, lon, pidx;
        while (std::getline(ifile, line))
        {
            std::istringstream iss(line);   // string stream

            //
            // On recucpere la lattitude
            //
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

            nLines += 1;
        }

        utcExample();
        std::cout <<  "Data were loaded (" <<  nLines << ") coordonates" << std::endl;

        ifile.close();
        ftime = curTime;
    }
    
    ~GpsPositionsRealTime()
    {
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
