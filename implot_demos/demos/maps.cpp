// Demo:   maps.cpp
// Author: Evan Pezent (evanpezent.com)
// Date:   7/25/2021

#include <stdio.h>
#include <curl/curl.h>
#include <string>
#include <filesystem>
#include <iostream>
#include <map>
#include <cmath>
#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <stdexcept>
#include <tuple>
#include <atomic>

#include <cstdio>
#include <cstdlib>
#include <fstream>

#include "App.h"
#include "Image.h"

namespace fs = std::filesystem;

// Useful Links and Resources
//
// https://operations.osmfoundation.org/policies/tiles/
// https://wiki.openstreetmap.org/wiki/Tile_servers
// https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames

#define TILE_SERVER   "http://a.tile.openstreetmap.fr/osmfr/"
#define TILE_SIZE     256                                 // the expected size of tiles in pixels, e.g. 256x256px
#define MAX_ZOOM      19                                  // the maximum zoom level provided by the server
#define MAX_THREADS   2                                   // the maximum threads to use for downloading tiles (OSC strictly forbids more than 2)
#define USER_AGENT    "Mozilla/4.0 (compatible; MSIE 8.0; Windows NT 6.0; Trident/4.0)"                        // change this to represent your own app if you extend this code
#define PI 3.14159265359

#if 0
ImVec4 PlaneColor(int type)
{
    switch (value) {
        case  0 : return IM_COL32(  0,   0,   0, 255); // BLACK
        case  1 : return IM_COL32(  0,   0,   0, 255); // BLACK
        case  2 : return IM_COL32(  0,   0,   0, 255); // BLACK
        case  3 : return IM_COL32(  0,   0,   0, 255); // BLACK
        case  5 : return IM_COL32(  0,   0,   0, 255); // BLACK
        case  6 : return IM_COL32(  0,   0,   0, 255); // BLACK
        case  7 : return IM_COL32(  0,   0,   0, 255); // BLACK
        case 10 : return IM_COL32(  0,   0,   0, 255); // BLACK
        case 20 : return IM_COL32(  0,   0,   0, 255); // BLACK
        case 21 : return IM_COL32(  0, 255,   0, 255); // GREEN SurfaceEmergencyVehicle;
        case 23 : return IM_COL32(  0, 255,   0, 255); // GREEN SurfaceServiceVehicle;
        case 24 : return IM_COL32(  0, 255,   0, 255); // GREEN
        case 25 : return IM_COL32(  0, 255,   0, 255); // GREEN
        case 26 : return IM_COL32(  0, 255,   0, 255); // GREEN
        case 27 : return IM_COL32(  0, 255,   0, 255); // GREEN
        case 30 : return IM_COL32(255,   0,   0, 255); // NoInformation;
        case 31 : return IM_COL32(255,   0,   0, 255); // Sailplane;
        case 32 : return IM_COL32(255,   0,   0, 255); // LighterThanAir;
        case 33 : return IM_COL32(255,   0,   0, 255); // Parachutist;
        case 34 : return IM_COL32(255,   0,   0, 255); // UltralightAircraft;
        case 35 : return IM_COL32(255,   0,   0, 255); // Reserved;
        case 36 : return IM_COL32(255,   0,   0, 255); // UnmannedVehicle;
        case 37 : return IM_COL32(255,   0,   0, 255); // SpaceVehicle;
        case 40 : return IM_COL32(255,   0,   0, 255); // NoInformation;
        case 41 : return IM_COL32(128, 128, 255, 255); // BLUE
        case 42 : return IM_COL32( 64,  64, 255, 255); // BLUE
        case 43 : return IM_COL32( 32,  32, 255, 255); // BLUE
        case 44 : return IM_COL32(  0,   0, 255, 255); // BLUE
        case 45 : return IM_COL32(  0,   0, 255, 255); // BLUE
        case 46 : return IM_COL32(  0,   0, 255, 255); // BLUE
        case 47 : return IM_COL32(255, 255, 0, 255); // Rotorcraft;
        default : return IM_COL32(255,   0,   0, 255);
    }
}
#endif

long double toRadians(const long double degree) {
    long double one_deg = (M_PI) / 180;
    return (one_deg * degree);
}

inline long double distance(long double lat1, long double long1,
                     long double lat2, long double long2) {
    lat1  = toRadians(lat1);
    long1 = toRadians(long1);
    lat2  = toRadians(lat2);
    long2 = toRadians(long2);

    // Haversine Formula
    long double dlong = long2 - long1;
    long double dlat = lat2 - lat1;

    long double ans = pow(sin(dlat / 2), 2) +
            cos(lat1) * cos(lat2) *
            pow(sin(dlong / 2), 2);

    ans = 2 * asin(sqrt(ans));

    // Radius of Earth in
    // Kilometers, R = 6371
    // Use R = 3956 for miles
    long double R = 6371;

    // Calculate the result
    ans = ans * R;

    return ans;
}

//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//

int long2tilex(double lon, int z)
 {
    return (int)(floor((lon + 180.0) / 360.0 * (1 << z)));
}

//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//

int lat2tiley(double lat, int z) {
    double latrad = lat * PI/180.0;
    return (int)(floor((1.0 - asinh(tan(latrad)) / PI) / 2.0 * (1 << z)));
}

//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//

double tilex2long(int x, int z) {
    return x / (double)(1 << z) * 360.0 - 180;
}

//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//

double tiley2lat(int y, int z) {
    double n = PI - 2.0 * PI * y / (double)(1 << z);
    return 180.0 / PI * atan(0.5 * (exp(n) - exp(-n)));
}

//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//

#include "GpsPositions.hpp"
#include "GpsPositionsRefresh.hpp"
#include "GpsPositionsRealTime.hpp"
#include "DistanceCircles.hpp"

//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//

struct TileCoord {
    int z; // zoom    [0......20]
    int x; // x index [0...z^2-1]
    int y; // y index [0...z^2-1]
    inline std::string subdir() const { return std::to_string(z) + "/" + std::to_string(x) + "/"; }
    inline std::string dir()    const { return "tiles/" + subdir(); }
    inline std::string file()   const { return std::to_string(y) + ".png"; }
    inline std::string path()   const { return dir() + file(); }
    inline std::string url()    const { return TILE_SERVER + subdir() + file(); }
    inline std::string label()  const { return subdir() + std::to_string(y); }

    std::tuple<ImPlotPoint,ImPlotPoint> bounds() const {
        double n = std::pow(2,z);
        double t = 1.0 / n;
        return {
                   { x*t     , (1+y)*t } ,
                   { (1+x)*t , (y)*t   }
               };
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

bool operator<(const TileCoord& l, const TileCoord& r ) {
    if ( l.z < r.z )  return true;
    if ( l.z > r.z )  return false;
    if ( l.x < r.x )  return true;
    if ( l.x > r.x )  return false;
    if ( l.y < r.y )  return true;
    if ( l.y > r.y )  return false;
                      return false;
}

//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//

enum TileState : int {
    Unavailable = 0, // tile not available
    Loaded,          // tile has been loaded into  memory
    Downloading,     // tile is downloading from server
    OnDisk           // tile is saved to disk, but not loaded into memory
};

typedef Image TileImage;

//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//

struct Tile {
    Tile() : state(TileState::Unavailable) {  }
    Tile(TileState s) : state(s) { }
    TileState state;
    TileImage image;
};

size_t curl_write_cb(void *ptr, size_t size, size_t nmemb, void *userdata) {
    FILE *stream = (FILE *)userdata;
    if (!stream) {
        printf("No stream\n");
        return 0;
    }
    size_t written = fwrite((FILE *)ptr, size, nmemb, stream);
    return written;
}

//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//

class TileManager {
public:

    TileManager() {
        start_workers();
    }
 
    inline ~TileManager() {
        {
            std::unique_lock<std::mutex> lock(m_queue_mutex);
            m_stop = true;
        }
        m_condition.notify_all();
        for(std::thread &worker: m_workers)
            worker.join();
    }

    const std::vector<std::pair<TileCoord, std::shared_ptr<Tile>>>& get_region(ImPlotRect view, ImVec2 pixels) {
        double min_x  = std::clamp(view.X.Min, 0.0, 1.0);
        double min_y  = std::clamp(view.Y.Min, 0.0, 1.0);
        double size_x = std::clamp(view.X.Size(),0.0,1.0);
        double size_y = std::clamp(view.Y.Size(),0.0,1.0);

        double pix_occupied_x = (pixels.x / view.X.Size()) * size_x;
        double pix_occupied_y = (pixels.y / view.Y.Size()) * size_y;

        double units_per_tile_x = view.X.Size() * (TILE_SIZE / pix_occupied_x);
        double units_per_tile_y = view.Y.Size() * (TILE_SIZE / pix_occupied_y);

        int z    = 0;
        double r = 1.0 / pow(2,z);
        while (r > units_per_tile_x && r > units_per_tile_y && z < MAX_ZOOM)
            r = 1.0 / pow(2,++z);
        
        m_region.clear();
        if (!append_region(z, min_x, min_y, size_x, size_y) && z > 0)
        {
            append_region(--z, min_x, min_y, size_x, size_y);
            std::reverse(m_region.begin(),m_region.end());
        }
        return m_region;
    }

    std::shared_ptr<Tile> request_tile(TileCoord coord)
    {
        std::lock_guard<std::mutex> lock(m_tiles_mutex);
        if (m_tiles.count(coord))
            return get_tile(coord);
        else if (fs::exists(coord.path()))
            return load_tile(coord);
        else
            download_tile(coord);
        return nullptr;
    }

    int tiles_loaded() const     { return m_loads;     }
    int tiles_downloaded() const { return m_downloads; }
    int tiles_cached() const     { return m_loads - m_downloads; }
    int tiles_failed() const     { return m_fails; }
    int threads_working() const  { return m_working; }

private:

    bool append_region(int z, double min_x, double min_y, double size_x, double size_y)
    {
//        printf("double min_x = %f and double min_y = %f / %f / %f\n", min_x, min_y, size_x, size_y);
        int k = pow(2,z);
        double r = 1.0 / k;
        int xa = min_x * k;
        int xb = xa + ceil(size_x / r) + 1;
        int ya = min_y * k;
        int yb = ya + ceil(size_y / r) + 1;
        
        xb = std::clamp(xb,0,k);
        yb = std::clamp(yb,0,k);
        
        bool covered = true;
        for (int x = xa; x < xb; ++x)
        {
            for (int y = ya; y < yb; ++y)
            {
                TileCoord coord{z,x,y};
                std::shared_ptr<Tile> tile = request_tile(coord);
                m_region.push_back({coord,tile});
                if (tile == nullptr || tile->state != TileState::Loaded)
                    covered = false;
            }
        }
        return covered;
    }

    void download_tile(TileCoord coord) {
        auto dir = coord.dir();
        fs::create_directories(dir);
        if (fs::exists(dir)) {
            m_tiles[coord] = std::make_shared<Tile>(Downloading);
            {
                std::unique_lock<std::mutex> lock(m_queue_mutex);
                m_queue.emplace(coord);
            }
            m_condition.notify_one();
        }
    }

    std::shared_ptr<Tile> get_tile(TileCoord coord) {
        if (m_tiles[coord]->state == Loaded)
            return m_tiles[coord];
        else if (m_tiles[coord]->state == OnDisk)
            return load_tile(coord);
        return nullptr;
    }

    std::shared_ptr<Tile> load_tile(TileCoord coord) {
        auto path = coord.path();
        if (!m_tiles.count(coord))
            m_tiles[coord] = std::make_shared<Tile>();
        if (m_tiles[coord]->image.LoadFromFile(path.c_str())) {
            m_tiles[coord]->state = TileState::Loaded;
            m_loads++;
            return m_tiles[coord];
        }
        m_fails++;
        printf("TileManager[00]: Failed to load \"%s\"\n", path.c_str());
        if (!fs::remove(path))
            printf("TileManager[00]: Failed to remove \"%s\"\n", path.c_str());
        printf("TileManager[00]: Removed \"%s\"\n",path.c_str());
        m_tiles.erase(coord);
        return nullptr;
    }

    void start_workers() {
        for(int thrd = 1; thrd < MAX_THREADS+1; ++thrd) {
            m_workers.emplace_back(
                [this, thrd] {
                    printf("TileManager[%02d]: Thread started\n",thrd);
                    CURL* curl = curl_easy_init();
                    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curl_write_cb);
                    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1);
                    curl_easy_setopt(curl, CURLOPT_USERAGENT, USER_AGENT);
                    for(;;)
                    {
                        TileCoord coord;
                        {
                            std::unique_lock<std::mutex> lock(m_queue_mutex);
                            m_condition.wait(lock,
                                [this]{ return m_stop || !m_queue.empty(); });
                            if(m_stop && m_queue.empty()) {
                                curl_easy_cleanup(curl);
                                printf("TileManager[%02d]: Thread terminated\n",thrd);
                                return;
                            }
                            coord = std::move(m_queue.front());
                            m_queue.pop();
                        }
                        m_working++;
                        bool success = true;
                        auto dir  = coord.dir();
                        auto path = coord.path();
                        auto url  = coord.url();
                        FILE *fp = fopen(coord.path().c_str(), "wb");
                        if (fp) {
                            curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
                            curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
                            CURLcode cc = curl_easy_perform(curl);
                            fclose(fp);
                            if (cc != CURLE_OK) {
                                printf("TileManager[%02d]: Failed to download: \"%s\"\n", thrd, url.c_str());
                                long rc = 0;
                                curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &rc);
                                if (!((rc == 200 || rc == 201) && rc != CURLE_ABORTED_BY_CALLBACK))
                                    printf("TileManager[%02d]: Response code: %d\n",thrd,rc);
                                fs::remove(coord.path());
                                success = false;
                            }
                        }
                        else {
                            printf("TileManager[%02d]: Failed to open or create file \"%s\"\n",thrd, path.c_str());
                            success = false;
                        }
                        if (success) {
                            m_downloads++;
                            std::lock_guard<std::mutex> lock(m_tiles_mutex);
                            m_tiles[coord]->state = OnDisk;
                        }
                        else {
                            m_fails++;
                            std::lock_guard<std::mutex> lock(m_tiles_mutex);
                            m_tiles.erase(coord);
                        }
                        m_working--;
                    }
                }
            );
        }
    }

    std::atomic<int> m_loads     = 0;
    std::atomic<int> m_downloads = 0;
    std::atomic<int> m_fails     = 0;
    std::atomic<int> m_working   = 0;
    std::map<TileCoord,std::shared_ptr<Tile>> m_tiles;
    std::mutex m_tiles_mutex;
    std::vector<std::pair<TileCoord, std::shared_ptr<Tile>>> m_region;
    std::vector<std::thread> m_workers;
    std::queue<TileCoord> m_queue;
    std::mutex m_queue_mutex;
    std::condition_variable m_condition;
    bool m_stop = false;
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
static bool debug             = false;

static bool only_reliable_gps = false;

struct ImMaps : public App {
    using App::App;

    void Update() override {

        if( soucoupe == nullptr )
        {
            soucoupe = new Image();
            soucoupe->LoadFromFile("soucoupe.png");
        }

        //
        // This code highlights the tile prefetching process for debugging purpose
        //
        static int renders = 0;
        if (ImGui::IsKeyPressed(ImGui::GetKeyIndex(ImGuiKey_A)))
            debug = !debug;

        //
        // this option enables to remove CRC brute forced points
        //
        if (ImGui::IsKeyPressed(ImGui::GetKeyIndex(ImGuiKey_B)))
            only_reliable_gps = !only_reliable_gps;

        ImGui::SetNextWindowPos({0,0});
        ImGui::SetNextWindowSize(GetWindowSize());
        ImGui::Begin("Map",0,ImGuiWindowFlags_NoResize|ImGuiWindowFlags_NoTitleBar);
        if (debug) {
            int wk = mngr.threads_working();
            int dl = mngr.tiles_downloaded();
            int ld = mngr.tiles_loaded();
            int ca = mngr.tiles_cached();
            int fa = mngr.tiles_failed();
            ImGui::Text("FPS: %.2f    Working: %d    Downloads: %d    Loads: %d    Caches: %d    Fails: %d    Renders: %d", ImGui::GetIO().Framerate, wk, dl, ld, ca, fa, renders);
        }

        //if (ImGui::BeginTabBar("Plots")) {
        //    if (ImGui::BeginTabItem("Magnitude")) {

        ImPlotAxisFlags ax_flags = ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoTickLabels | ImPlotAxisFlags_NoGridLines| ImPlotAxisFlags_Foreground;
        if (ImPlot::BeginPlot("##Map",ImVec2(-1,-1),ImPlotFlags_Equal|ImPlotFlags_NoMouseText)) {
            ImPlot::SetupAxes(NULL,NULL,ax_flags,ax_flags|ImPlotAxisFlags_Invert);
            ImPlot::SetupAxesLimits(0,1,0,1);

            auto pos     = ImPlot::GetPlotPos();
            auto size    = ImPlot::GetPlotSize();
            auto limits  = ImPlot::GetPlotLimits();
            auto& region = mngr.get_region(limits,size);
            renders = 0;

            if (debug) {
                float ys[] = {1,1};
                ImPlot::SetNextFillStyle({1,0,0,1},0.5f);
                ImPlot::PlotShaded("##Bounds",ys,2);
            }

            for (auto& pair : region)
            {
                TileCoord coord            = pair.first;    // les coordonées
                std::shared_ptr<Tile> tile = pair.second;   // la tuile
                auto [bmin,bmax] = coord.bounds();
                if (tile != nullptr) {
                    auto col = debug ? ((coord.x % 2 == 0 && coord.y % 2 != 0) || (coord.x % 2 != 0 && coord.y % 2 == 0))? ImVec4(1,0,1,1) : ImVec4(1,1,0,1) : ImVec4(1,1,1,1);
                    ImPlot::PlotImage("##Tiles",(void*)(intptr_t)tile->image.ID,bmin,bmax,{0,0},{1,1},col);
                }
                if (debug)
                    ImPlot::PlotText(coord.label().c_str(),(bmin.x+bmax.x)/2,(bmin.y+bmax.y)/2);
                renders++;
            }



            {
                ImPlot::PushStyleColor(ImPlotCol_InlayText, IM_COL32(255, 0, 0, 255));

                gps.update( only_reliable_gps );

#if 1
                //
                // On indique que nous souhaitons afficher des lignes noir pour
                // représenter les trajectoires des avions ainsi que des marqueurs
                // de type crois (+)
                //
                ImVec4 black = ImVec4( 0, 0, 0, 1);
                ImPlot::PushStyleColor    (ImPlotCol_Line, black);
                ImPlot::SetNextMarkerStyle(ImPlotMarker_Plus);

                //
                // On parcours l'ensemble des coordonnées à notre disposition
                //
                      int start  = 0;
                const int length = gps.pos_x.size();
                for(int l = 1; l < length; l += 1)
                {
                    if( gps.idx[l] == gps.idx[start] )
                    {
                        // l'avions est toujours le même donc on a rien à faire ...
                    }
                    else
                    {
                        //
                        // On passe à l'avion suivant donc on doit lancer le tracé
                        // de la séquence de points
                        //
                        float* ptr_x = gps.pos_x.data() + start;
                        float* ptr_y = gps.pos_y.data() + start;
                        int    size  = l - start;
                        ImPlot::PlotLine("", ptr_x, ptr_y, size);
                        ImPlot::PushStyleColor(ImPlotCol_Line, black);
                        ImPlot::PlotText("X", gps.pos_x[l-1], gps.pos_y[l-1]);
                        start = l;
                        ImPlot::PushStyleColor(ImPlotCol_Line, black);
                        ImPlot::SetNextMarkerStyle(ImPlotMarker_Plus);
                    }
                }
                //
                // On affiche le dernier avion qui se trouvait en fin de liste
                //
                int    size  = length - start;
                float* ptr_x = gps.pos_x.data() + start;
                float* ptr_y = gps.pos_y.data() + start;
                ImPlot::PlotLine("", ptr_x, ptr_y, size);
                ImPlot::PushStyleColor(ImPlotCol_Line, black);
                ImPlot::PlotText("X", gps.pos_x[length-1], gps.pos_y[length-1]);

//              ImPlotPoint bounds_min;
//                bounds_min.x = gps.pos_x[length-1] - 0.01;
//                bounds_min.y = gps.pos_y[length-1] - 0.01;

//                ImPlotPoint bounds_max;
//                bounds_max.x = gps.pos_x[length-1] + 0.01;
//                bounds_max.y = gps.pos_y[length-1] + 0.01;

//                ImPlot::PlotImage("##Tiles",(void*)(intptr_t)soucoupe, bounds_min, bounds_max, {0,0},{1,1}, ImVec4(1,1,1,1));

#else
                const int length = gps.pos_x.size();
                for(int l = 0; l < length; l += 1)
                {
                    if( l != (length-1) )
                    {
                        if( gps.idx[l] == gps.idx[l+1] )
                            ImPlot::PlotText(".", gps.pos_x[l], gps.pos_y[l]);
                        else
                            ImPlot::PlotText("X", gps.pos_x[l], gps.pos_y[l]);
                    }
                    else
                    {
                        ImPlot::PlotText("X", gps.pos_x[l], gps.pos_y[l]);
                    }
                }
#endif
#if 1
                //
                // Affichage de la position de l'ENSEIRB
                //

                double lon = -0.605208;
                double lat = 44.805643;
                double pX  = conv_lon_to_x( lon );
                double pY  = conv_lat_to_y( lat );
                ImPlot::PlotText("O", pX, pY);


                dist.update();

#if 0                //
                //
                // grille horizontale de 100 km
                //
                double h_pos = 0;
                for(int i = 0; i < 100; i++)
                {
                    if( distance(lat, lon, lat, lon + h_pos) >= 100 )
                        break;
                    h_pos += 0.1f;
                }
//                std::cout << distance(lat, lon, lat, lon + h_pos) << std::endl;
                //
                //
//                static double h_vals[5] = { conv_lon_to_x(lon - h_pos), conv_lon_to_x(lon - h_pos/2), conv_lon_to_x(lon), conv_lon_to_x(lon + h_pos/2), conv_lon_to_x(lon + h_pos) };
//                ImPlot::PlotInfLines("", h_vals, 5);
                //
                //
                // grille verticale de 100 km
                //
                double v_pos = 0;
                for(int i = 0; i < 100; i++)
                {
                    if( distance(lat, lon, lat + v_pos, lon) >= 100 )
                        break;
                    v_pos += 0.1f;
                }
//                std::cout << distance(lat, lon, lat + v_pos, lon) << std::endl;
                //
                //
//                static double v_vals[5] = { conv_lat_to_y(lat - v_pos), conv_lat_to_y(lat - v_pos/2), conv_lat_to_y(lat), conv_lat_to_y(lat + v_pos/2), conv_lat_to_y(lat + v_pos) };
//                ImPlot::PlotInfLines("", v_vals, 5, ImPlotInfLinesFlags_Horizontal);
#endif
                //
                //
                // Affichage de la position du musée de GreenWich
                //
                lon = 0.0;
                lat = 51.476852;
                pX  = conv_lon_to_x( lon );
                pY  = conv_lat_to_y( lat );
                ImPlot::PlotText("X", pX, pY);

#if 0
                lon = -0.605208;
                lat = 44.805643;
                static const int k_points_per = 50;
                static double x_data_050[k_points_per];
                static double y_data_050[k_points_per];
                static double x_data_100[k_points_per];
                static double y_data_100[k_points_per];
                static double x_data_200[k_points_per];
                static double y_data_200[k_points_per];
                for (int p = 0; p < k_points_per; p+=1)
                {
                    x_data_050[p] = conv_lon_to_x( lon + 0.5 * h_pos * cos((double)p/(k_points_per-1) * 6.28) );
                    y_data_050[p] = conv_lat_to_y( lat + 0.5 * v_pos * sin((double)p/(k_points_per-1) * 6.28) );

                    x_data_100[p] = conv_lon_to_x( lon +       h_pos * cos((double)p/(k_points_per-1) * 6.28) );
                    y_data_100[p] = conv_lat_to_y( lat +       v_pos * sin((double)p/(k_points_per-1) * 6.28) );

                    x_data_200[p] = conv_lon_to_x( lon + 2.0 * h_pos * cos((double)p/(k_points_per-1) * 6.28) );
                    y_data_200[p] = conv_lat_to_y( lat + 2.0 * v_pos * sin((double)p/(k_points_per-1) * 6.28) );

                }
                ImPlot::PlotLine("", x_data_050, y_data_050, k_points_per);
                ImPlot::PlotLine("", x_data_100, y_data_100, k_points_per);
                ImPlot::PlotLine("", x_data_200, y_data_200, k_points_per);

#endif
#endif
            }

            //
            // On rajoute le petit logo à la fin de la fenetre
            //

            ImPlot::PushPlotClipRect();
            static const char* label = ICON_FA_COPYRIGHT " OpenStreetMap Contributors";
            auto label_size = ImGui::CalcTextSize(label);
            auto label_off  = ImPlot::GetStyle().MousePosPadding;
            ImPlot::GetPlotDrawList()->AddText({pos.x + label_off.x, pos.y+size.y-label_size.y-label_off.y},IM_COL32_BLACK,label);
            ImPlot::PopPlotClipRect();
            ImPlot::EndPlot();
        }
        ImGui::End();
    }

    Image* soucoupe = nullptr;

    TileManager          mngr;
#if 0
    GpsPositionsRefresh  gps;
#else
    GpsPositionsRealTime gps;
    DistanceCircles dist;
#endif
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

int main(int argc, char const *argv[])
{
    ImMaps app("ImMaps",960,540,argc,argv);
    app.Run();
}
