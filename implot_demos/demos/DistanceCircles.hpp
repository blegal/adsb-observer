#ifndef _DistanceCircles_
#define _DistanceCircles_

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

class DistanceCircles
{

private:

    const int k_points_per = 50;

    double* x_data_050;
    double* y_data_050;
    
    double* x_data_100;
    double* y_data_100;
    
    double* x_data_200;
    double* y_data_200;


public:


    DistanceCircles(const double lon = -0.605208, const double lat = 44.805643)
    {
        //
        // Allocation mémoire
        //

        x_data_050 = new double[k_points_per];
        y_data_050 = new double[k_points_per];

        x_data_100 = new double[k_points_per];
        y_data_100 = new double[k_points_per];

        x_data_200 = new double[k_points_per];
        y_data_200 = new double[k_points_per];


        double h_pos = 0;
        for(int i = 0; i < 100; i += 1)
        {
            if( distance(lat, lon, lat, lon + h_pos) >= 100 )
                break;
            h_pos += 0.05f;
        }

        double v_pos = 0;
        for(int i = 0; i < 100; i += 1)
        {
            if( distance(lat, lon, lat + v_pos, lon) >= 100 )
                break;
            v_pos += 0.05f;
        }

        for (int p = 0; p < k_points_per; p += 1)
        {
            //
            // 50 km
            //
            x_data_050[p] = conv_lon_to_x( lon + 0.5 * h_pos * cos((double)p/(k_points_per-1) * 6.28) );
            y_data_050[p] = conv_lat_to_y( lat + 0.5 * v_pos * sin((double)p/(k_points_per-1) * 6.28) );

            //
            // 100 km
            //
            x_data_100[p] = conv_lon_to_x( lon +       h_pos * cos((double)p/(k_points_per-1) * 6.28) );
            y_data_100[p] = conv_lat_to_y( lat +       v_pos * sin((double)p/(k_points_per-1) * 6.28) );

            //
            // 200 km
            //
            x_data_200[p] = conv_lon_to_x( lon + 2.0 * h_pos * cos((double)p/(k_points_per-1) * 6.28) );
            y_data_200[p] = conv_lat_to_y( lat + 2.0 * v_pos * sin((double)p/(k_points_per-1) * 6.28) );

        }

    }


    void update()
    {
        ImPlot::PlotLine("", x_data_050, y_data_050, k_points_per);
        ImPlot::PlotLine("", x_data_100, y_data_100, k_points_per);
        ImPlot::PlotLine("", x_data_200, y_data_200, k_points_per);
    }
    

    ~DistanceCircles()
    {
        delete[] x_data_050;
        delete[] y_data_050;

        delete[] x_data_100;
        delete[] y_data_100;

        delete[] x_data_200;
        delete[] y_data_200;
    }

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
