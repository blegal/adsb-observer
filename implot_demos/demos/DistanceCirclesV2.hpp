#ifndef _DistanceCirclesV2_
#define _DistanceCirclesV2_

#include <string>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
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

class DistanceCirclesV2
{

private:

    const int k_points_per = 256;

    double* x_data_050;
    double* y_data_050;
    
    double* x_data_100;
    double* y_data_100;
    
    double* x_data_150;
    double* y_data_150;
    
    double* x_data_200;
    double* y_data_200;

    double* x_data_250;
    double* y_data_250;
    
    double* x_data_300;
    double* y_data_300;

    const double lon;
    const double lat;


public:

    double compute_offset_x(const double lat, const double lon, const double dist)
    {
        double h_pos = 0.0;
        for(int i = 0; i < 100000; i += 1)
        {
            if( distance(lat, lon, lat, lon + h_pos) >= dist )
                break;
            h_pos += 0.0001;
        }
        return h_pos;
    }

    double compute_offset_y(const double lat, const double lon, const double dist)
    {
        double v_pos = 0;
        for(int i = 0; i < 100000; i += 1)
        {
            if( distance(lat, lon, lat + v_pos, lon) >= dist )
                break;
            v_pos += 0.0001;
        }
        return v_pos;
    }



    DistanceCirclesV2(const double rayon = 50.0, const double _lon = -0.605208, const double _lat = 44.805643) : lon(_lon), lat(_lat)
    {
        //
        // Allocation mémoire
        //

        x_data_050 = new double[k_points_per];
        y_data_050 = new double[k_points_per];

        x_data_100 = new double[k_points_per];
        y_data_100 = new double[k_points_per];

        x_data_150 = new double[k_points_per];
        y_data_150 = new double[k_points_per];

        x_data_200 = new double[k_points_per];
        y_data_200 = new double[k_points_per];

        x_data_250 = new double[k_points_per];
        y_data_250 = new double[k_points_per];

        x_data_300 = new double[k_points_per];
        y_data_300 = new double[k_points_per];

        regenerate(rayon);
    }

    void regenerate(const double rayon)
    {
        double off_x_050 = compute_offset_x(lat, lon, 1.0 * rayon);
        double off_y_050 = compute_offset_y(lat, lon, 1.0 * rayon);
        
        double off_x_100 = compute_offset_x(lat, lon, 2.0 * rayon);
        double off_y_100 = compute_offset_y(lat, lon, 2.0 * rayon);
        
        double off_x_150 = compute_offset_x(lat, lon, 3.0 * rayon);
        double off_y_150 = compute_offset_y(lat, lon, 3.0 * rayon);
        
        double off_x_200 = compute_offset_x(lat, lon, 4.0 * rayon);
        double off_y_200 = compute_offset_y(lat, lon, 4.0 * rayon);
        
        double off_x_250 = compute_offset_x(lat, lon, 5.0 * rayon);
        double off_y_250 = compute_offset_y(lat, lon, 5.0 * rayon);

        double off_x_300 = compute_offset_x(lat, lon, 6.0 * rayon);
        double off_y_300 = compute_offset_y(lat, lon, 6.0 * rayon);

        for (int p = 0; p < k_points_per; p += 1)
        {
            //  50 km
            x_data_050[p] = conv_lon_to_x( lon + off_x_050 * cos((double)p/(k_points_per-1) * 6.28) );
            y_data_050[p] = conv_lat_to_y( lat + off_y_050 * sin((double)p/(k_points_per-1) * 6.28) );

            // 100 km
            x_data_100[p] = conv_lon_to_x( lon + off_x_100 * cos((double)p/(k_points_per-1) * 6.28) );
            y_data_100[p] = conv_lat_to_y( lat + off_y_100 * sin((double)p/(k_points_per-1) * 6.28) );

            // 150 km
            x_data_150[p] = conv_lon_to_x( lon + off_x_150 * cos((double)p/(k_points_per-1) * 6.28) );
            y_data_150[p] = conv_lat_to_y( lat + off_y_150 * sin((double)p/(k_points_per-1) * 6.28) );

            // 200 km
            x_data_200[p] = conv_lon_to_x( lon + off_x_200 * cos((double)p/(k_points_per-1) * 6.28) );
            y_data_200[p] = conv_lat_to_y( lat + off_y_200 * sin((double)p/(k_points_per-1) * 6.28) );

            // 250 km
            x_data_250[p] = conv_lon_to_x( lon + off_x_250 * cos((double)p/(k_points_per-1) * 6.28) );
            y_data_250[p] = conv_lat_to_y( lat + off_y_250 * sin((double)p/(k_points_per-1) * 6.28) );

            // 250 km
            x_data_300[p] = conv_lon_to_x( lon + off_x_300 * cos((double)p/(k_points_per-1) * 6.28) );
            y_data_300[p] = conv_lat_to_y( lat + off_y_300 * sin((double)p/(k_points_per-1) * 6.28) );
        }
    }

    void update()
    {
        ImPlot::PlotLine("", x_data_050, y_data_050, k_points_per, ImPlotLineFlags_Segments);
        ImPlot::PlotLine("", x_data_100, y_data_100, k_points_per);
        ImPlot::PlotLine("", x_data_150, y_data_150, k_points_per, ImPlotLineFlags_Segments);
        ImPlot::PlotLine("", x_data_200, y_data_200, k_points_per);
        ImPlot::PlotLine("", x_data_250, y_data_250, k_points_per, ImPlotLineFlags_Segments);
        ImPlot::PlotLine("", x_data_300, y_data_300, k_points_per);
    }
    

    ~DistanceCirclesV2()
    {
        delete[] x_data_050;
        delete[] y_data_050;

        delete[] x_data_100;
        delete[] y_data_100;

        delete[] x_data_150;
        delete[] y_data_150;

        delete[] x_data_200;
        delete[] y_data_200;

        delete[] x_data_250;
        delete[] y_data_250;

        delete[] x_data_300;
        delete[] y_data_300;
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
