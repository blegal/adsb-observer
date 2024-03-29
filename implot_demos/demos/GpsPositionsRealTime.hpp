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
#include <exception>

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
    bool last_fiable_only;

public:

    GpsPositionsRealTime()
    {
        last_fiable_only = false;
    }
    
    void update(const bool fiable_only)
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

            if( last_fiable_only == fiable_only )
                return;
            last_fiable_only = fiable_only;
        }

        pos_x.clear();
        pos_y.clear();
        idx.clear();
        pos_crc.clear();
        pos_oaci.clear();
        pos_type.clear();


//      std::cout << "(DD) The file was updated (" <<  filename << ")" << std::endl;

        //
        // On va parser le fichier pour extraire les position (lat, long) des avions
        // un couple de valeur par ligne formaté en "%f %f"
        //
        int nLines = 0;
        int nInfos = 0;
        std::ifstream ifile( filename );
        std::string line, lat, lon, pidx, crc, oaci, type;
        while (std::getline(ifile, line))
        {
            std::istringstream iss(line);   // string stream

            //
            // On recucpere la lattitude
            //
            std::getline(iss, lat, ' ');    //
            if(lat.size() == 0)              // un espace debute la ligne !
                std::getline(iss, lat, ' '); //

            std::getline(iss, lon,  ' ');    // longitude
            std::getline(iss, crc,  ' ');    // crc
            std::getline(iss, oaci, ' ');    // oaci
            std::getline(iss, type, ' ');    // type
            std::getline(iss, pidx, ' ');    // plane id


            float f_lon  = 0;
            float f_lat  = 0;
            int   i_crc  = 0;
            int   i_oaci = 0;
            int   i_type = 0;
            int   p_id   = 0;

            try
            {
                f_lon  = conv_lon_to_x( std::stof(lon) );    // on convertit la lattitude et la longitude
                f_lat  = conv_lat_to_y( std::stof(lat) );    // dans des coordonnées compatible 2D
                i_crc  =                std::stoi(crc)  ;    // dans des coordonnées compatible 2D
                //i_oaci =                std::stoi(oaci) ;    // dans des coordonnées compatible 2D
                i_type =                std::stoi(type) ;    // dans des coordonnées compatible 2D
                p_id   =                std::stoi(pidx) ;

                if( (fiable_only == false) || (i_crc == 1) )
                {

                    pos_x.push_back   ( f_lon  );   // On ajoute les données dans la liste
                    pos_y.push_back   ( f_lat  );   // On ajoute les données dans la liste
                    pos_crc.push_back ( i_crc  );   // On ajoute les données dans la liste
                    pos_oaci.push_back( i_oaci );   // On ajoute les données dans la liste
                    pos_type.push_back( i_type );   // On ajoute les données dans la liste
                    idx.push_back     ( p_id   );   // On ajoute les données dans la liste
                    nInfos += 1;
                }
                nLines += 1;
            }
            catch(std::exception &err)
            {
                std::cout <<  "(EE) An errors happened diring the file parsing..." << std::endl;
                std::cout <<  "(EE) line " << nLines << " was (" << line << ")"    << std::endl;
                break;
            }
        }

        utcExample();
        std::cout <<  "Data were loaded (" <<  nLines << ") for (" << nInfos << ") coordonates" << std::endl;

        ifile.close();
        ftime = curTime;
    }
    
    ~GpsPositionsRealTime()
    {
        pos_x.clear();
        pos_y.clear();
        idx.clear();
        pos_crc.clear();
        pos_oaci.clear();
        pos_type.clear();
    }

    std::vector<float> pos_x;
    std::vector<float> pos_y;
    std::vector<  int> idx;
    std::vector<  int> pos_crc;
    std::vector<  int> pos_oaci;
    std::vector<float> pos_type;
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
