#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <stdlib.h>
#include <SDL.h>
#include <time.h>
#include <string>
#include <math.h>
#include <constants.h>
#include <diagnostics.h>
#include <thread>

using namespace std;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int screenSize[3] = { 720, 720 , 0 };


double dt = 1e-10; 
double t = 0;

double X = 0.01;
double Y = 0.01;

int k = 0;
int j = 0;
int maxk = 25;

int n = 20000;//Max number of beam particles
int ndt = 1;//number generated per time step
int dtn = 1;//number of time step per generation


double EB = 1000;//eV
double epsilon = 0;
double VB = sqrt(2 * abs(EB) / mh * ev);
double VX = VB * cos(epsilon);
double VY = 0;
double cosb = 1;
double sinb = 0;
double XB = 0;
double YB = 0;
double RB = Y / 2;
double delta = pi / 8;

double T = 10000000;
double Tsp = 40000;

double sputter_y[15] = { 0.007258, 0.121991, 0.279897, 0.439335, 0.591515, 0.734769, 0.869293, 0.995785, 1.115029, 1.227765, 1.334657, 1.436285, 1.533153, 1.625701, 1.71431 };


double dE = 100;
double dteta = pi * 0.5 * 0.025;

bool const UI = true;
bool diagnostics = true;
bool quit = false;
const bool collision = true;
bool control = false;
int xMouse, yMouse;
int xWindow, yWindow;

unsigned seed = chrono::steady_clock::now().time_since_epoch().count();
default_random_engine rnd(seed);


class Boundary
{

public:

    double x1;
    double x2;
    double y1;
    double y2;
    double alpha = 0;
    double sina = 0;
    double cosa = 0;
    double sin2a = 0;
    double cos2a = 0;
    int A;
    float conductivity = 1;
    double T;
    double diposition[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double sputtering[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    double temperature[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    double coords[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    double dx = 0;
    double dy = 0;
    double L;
    double dL = 0;
    double Lx;
    double Ly;
    double Lx_1 = 0;

    Boundary(double X1, double Y1, double X2, double Y2, int a, double t);
    ~Boundary();
    void Sputter(double E, double X, double Y, double phi);
    void TermalConduction();


private:
    double M_1 = 0;

};

Boundary::Boundary(double X1, double Y1, double X2, double Y2, int a, double t)
{
    x1 = X1;
    x2 = X2;
    y1 = Y1;
    y2 = Y2;
    A = a;
    T = t;
    //alpha = acos((x2 - x1) / sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)));
    alpha = atan((y2 - y1) / (x2 - x1));
    //alpha = abs(alpha) - (sgn(alpha) - 1) * pi / 4;
    //cout << alpha << endl;
    sina = sin(alpha);
    cosa = cos(alpha);
    sin2a = sin(2 * alpha);
    cos2a = cos(2 * alpha);

    dx = x1 - x2;
    dy = y1 - y2;
    Lx = abs(dx);
    Ly = abs(dy);
    L = sqrt(Lx * Lx + Ly * Ly);
    dL = L / 10;
    Lx_1 = 1 / Lx;
    M_1 = 1 / (A * mh);
}

Boundary::~Boundary()
{
}

vector<Boundary> bound;//= { Boundary(0 ,Y ,X,0,-1) };//{ Boundary(X/2,Y/2,X*0.6,Y*0.4,-1), Boundary(X * 0.6,Y * 0.4,X*0.6,Y*0.35,-1) ,Boundary(X * 0.6,Y * 0.35,X*0.55,Y*0.3,-1) ,Boundary(X * 0.55,Y * 0.3,X*0.3,Y*0.275,-1) ,Boundary(X * 0.3,Y * 0.275,X*0.2,Y*0.2,-1) };

class Species
{
public:

    double x;
    double y;
    double vx;
    double vy;
    int vxn = sgn(vx);
    int vyn = sgn(vy);

    double x1 = 0;
    double y1 = 0;
    double XC = 0;
    double YC = 0;
    
    Species(bool ref, double x1, double x2, double v1, double v2, double m);
    ~Species();

    double V() {
        double v = sqrt(vx * vx + vy * vy);
        return v;
    }

    double V2() {
        double v = vx * vx + vy * vy;
        return v;
    }

    int push(double dt);

private:
    double M;
    
    double d_1 = 0;
    double d1 = 0;
    double d2 = 0;
    bool reflect = 1;
    double v = 0;
    double phi = 0;
    double cosphi = 0;
    double sinphi = 0;
    double cos2g = 0;
    double sin2g = 0;
    double cosg = 0;
    double sing = 0;

};

vector<Species> ions;
vector<Species> ionsS;

void generateParticlesVector(vector<Species>& part, int n, double X1, double X2, double Y1, double Y2, double Vx, double Vy, int a, double T = 0, double phi = 0) {

    double x, y;
        //sigma, 
        //vx, vy;
        //dvx, dvy;

    //sigma = sqrt(kb * T / (a * mh));
    
    //default_random_engine e(seed);

    for (int i = 0; i < n; i++)
    {
        /*if (T > 0)
        {
            
            normal_distribution<double> distribution(0, sigma);
            dvx = distribution(rnd);
            dvy = distribution(rnd);// hernya
        }*/
        //else { dvx = dvy = 0; }

        x = X1 + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (X2 - X1)));
        y = Y1 + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (Y2 - Y1)));
        //vx = Vx;// +dvx;
        //vy = Vy;// +sin(epsilon - phi / 2 + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (phi)))) * VB;
        Species p = Species(1, x, y, VX, VY, a*mh);
        part.push_back(p);
        beam_born++;
    }
}

void Boundary::Sputter(double E, double X, double Y, double phi) {

    double meanV, vx, vy, gama, v, s;
    int J, n;
    double teta = abs(abs(alpha - phi) - pi/2);

    //cout << teta << endl;

    J = (int)round(E / dE);
    if (J > 14) J = 14;//?????
    s = (sputter_y[J]) * abs(cos(teta));
    n = (int)ceil(s);
    int a = (int)(s/n * 100);
    int b = rand() % 100;
  
    if (b < a) {

        //int g = (X - x1) * Lx_1 * 10;
        //if (g > 9) g = 9; else if (g == -2147483648) g = 0;
        //cout << g << endl;
        //sputtering[g]+=n;

        meanV = sqrt(2 * kb * T * M_1);
        normal_distribution<double> angleDist(pi*0.5 * sgn(phi - alpha), 0.6);
        gama = alpha - angleDist(rnd);//alpha - static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (pi))) * sgn(phi - alpha); //*sgn(teta)
        //cout << gama << endl;
        chi_squared_distribution<double> distribution(3);

        for (int i = 0; i < n; i++)
        {
            v = distribution(rnd)*meanV;
            vx = v * cos(gama);
            vy = v * sin(gama);
            
            Species p = Species(0, X, Y, vx, vy, A*mh);
            ionsS.push_back(p);
            sputter_born++;
        }
        
    }
}

void Boundary::TermalConduction() {

    for (int i = 1; i < 9; i++)
    {
        temperature[i] = temperature[i] + dt / (dL * dL) * (temperature[i - 1] - 2 * temperature[i] + temperature[i + 1]);
    }

}

Species::Species(bool ref, double x1, double x2, double v1, double v2, double m)
{
    M = m;
    x = x1;
    y = x2;
    vx = v1;
    vy = v2;
    reflect = ref;
}

Species::~Species()
{
    if (reflect)
    {
        beam_dead++;
    }
}

int Species::push(double dt)  
{
    x1 = x + vx * dt;
    y1 = y + vy * dt;

    //vxn = sgn(vx);
    vyn = sgn(vy);

    for (int i = 0; i < bound.size(); i++)
    {
        d_1 = 1/((x - x1) * (bound[i].dy) - (y - y1) * (bound[i].dx));
        d1 = ((x - bound[i].x1) * (bound[i].dy) - (y - bound[i].y1) * (bound[i].dx)) * d_1;
        d2 = ((x - bound[i].x1) * (y - y1) - (y - bound[i].y1) * (x - x1)) * d_1;

        if (d1 >= 0 && d1 <= 1 && d2 >= 0 && d2 <= 1)
        {
            XC = x + (d1 * (x1 - x));
            YC = y + (d2 * (y1 - y));

            if (reflect) {
 
                v = V();
                
                cosphi = vx / v;
                sinphi = vy / v;
                phi = acos(cosphi) * vyn;//acos(cosphi) * (vyn + 1) * 0.5 + (2 * pi - acos(cosphi)) * -(vyn - 1) * 0.5;

                cosg = bound[i].cos2a * cosphi + bound[i].sin2a * sinphi;
                sing = bound[i].sin2a * cosphi - bound[i].cos2a * sinphi;

                cos2g = bound[i].cos2a * (cosphi * cosphi - sinphi * sinphi) + 2 * bound[i].sin2a * sinphi * cosphi;
                sin2g = bound[i].sin2a * (cosphi * cosphi - sinphi * sinphi) - 2 * bound[i].cos2a * sinphi * cosphi;

                vx = v * cosg;
                vy = v * sing;

                x = XC + ((x1 - XC) * cos2g - (y1 - YC) * sin2g) * 2;//((x1 - XC) * cos2g - (y1 - YC) * sin2g);//0.001 * X * (sinphi * bound[i].cosa - cosphi * bound[i].sina) * sgn(bound[i].alpha);  //Hernya                
                y = YC + ((x1 - XC) * sin2g + (y1 - YC) * cos2g) * 2;//((x1 - XC) * sin2g + (y1 - YC) * cos2g);//0.001 * Y * (sinphi * bound[i].cosa - cosphi * bound[i].sina);  

                bound[i].Sputter(V2() * mh * 0.5 * ev_1, x, y, phi);

                return i + 1;

            }
            else {

                int g = (XC - bound[i].x1) * bound[i].Lx_1 * 10;
                if (g > 9) g = 9; else if (g == -2147483648) g = 0;
                //cout << g << endl;
                bound[i].diposition[g]++;
                return -(i + 1);
            }
        }
        else
        {
            continue;
        }
    }
    x = x1;
    y = y1;
    if (y1 < 0)
    {
       y = -y1;
        vy = -vy;
    }
    else y = y1;
    return 0;
}

void pushParticles(vector<Species> &part, int type) {

    for (int j = 0; j < part.size(); j++)
    {
        if (part[j].push(dt) < 0) { 
            swap(part[j], part.back());
            part.pop_back();
            sputter_absorbed++;//temporary
        }
        else if (part[j].x<0 || part[j].x > X || part[j].y > Y) { 
            swap(part[j], part.back());
            part.pop_back();
            if (type)
            {
                sputter_away++;
            }
        }
    }
}

void Diagnostics(string path) {

    ofstream myfile;

    /*
    for (size_t i = 0; i < ions.size(); i++)
    {
        if (ions[i].V() > maxV) maxV = ions[i].V();
    }

    for (size_t i = 0; i < nv; i++)
    {
        V[i] = 0;
    }

    for (size_t i = 0; i < ions.size(); i++)
    {
        V[(int)round(ions[i].V() / maxV * nv)] += 1;
    }
    */

    myfile.open("diag/diag.txt", ios_base::app);
    myfile << t << '\t' 
        << ions.size() << '\t' 
        << beam_born << '\t' 
        << ionsS.size() << '\t' 
        << sputter_born  << '\t' 
        << sputter_absorbed << '\t' 
        << sputter_away << endl;
    myfile.close();

    myfile.open("diag/diagV.txt");

    float dv = maxV / nv;

    for (size_t i = 1; i <= nv; i++)
    {
        myfile << dv * i << '\t' << V[i - 1] << endl;
    }
    myfile << maxV << endl;
    myfile.close();

    for (size_t j = 0; j < bound.size(); j++)
    {
        myfile.open(format("diag/boundary{}.txt", j+1));
        for (size_t i = 0; i < 10; i++)
        {
            myfile << bound[j].diposition[i] << '\t'
                   << bound[j].sputtering[i] << endl;
        }
        myfile.close();
    }
    
}

void drawParticles(SDL_Renderer *rend, vector<Species> part, unsigned char r, unsigned char g, unsigned char b) {

    SDL_SetRenderDrawColor(rend, r, g, b, 0xFF);

    for (int i = 0; i < part.size(); i++)
    {
        SDL_RenderDrawPoint(rend, (int)(part[i].x / X * (screenSize[0] - screenSize[2]) + screenSize[2]), (int)(screenSize[1] - part[i].y / Y * screenSize[1]));
    }
}

void drawBoundary(SDL_Renderer* rend, vector<Boundary> bound, unsigned char  r, unsigned char g, unsigned char b) {

    SDL_SetRenderDrawColor(rend, r, g, b, 0xFF);

    for (size_t i = 0; i < bound.size(); i++)
    {
        SDL_RenderDrawLine(rend, (int)(bound[i].x1 / X * (screenSize[0] - screenSize[2]) + screenSize[2]), (int)(screenSize[1] - bound[i].y1 / Y * screenSize[1]), (int)(bound[i].x2 / X * (screenSize[0] - screenSize[2]) + screenSize[2]), (int)(screenSize[1] - bound[i].y2 / Y * screenSize[1]));
    }
}

void input(string path) {

    
    double param[6] = {0,0,0,0,0,0};
    ifstream file;
    string text;
    file.open("input.txt");
    if (file.is_open())
    {
        while (true)
        {
            file >> text;
            if (text == "Domain")
            {
                while (text != "}")
                {
                    file >> text;
                    if (text == "X")
                    {
                        file >> text;
                        file >> text;
                        X = stod(text);
                    }
                    else if (text == "R")
                    {
                        file >> text;
                        file >> text;
                        Y = stod(text);
                    }
                    else if (text == "dt")
                    {
                        file >> text;
                        file >> text;
                        dt = stod(text);
                    }
                    else if (text == "T")
                    {
                        file >> text;
                        file >> text;
                        T = stod(text);
                    }
                    else if (text == "EB")
                    {
                        file >> text;
                        file >> text;
                        EB = stod(text);
                    }
                    else if (text == "YB")
                    {
                        file >> text;
                        file >> text;
                        YB = stod(text);
                    }
                    else if (text == "RB")
                    {
                        file >> text;
                        file >> text;
                        RB = stod(text);
                    }
                    else if (text == "dtn")
                    {
                        file >> text;
                        file >> text;
                        dtn = stoi(text);
                    }
                    else if (text == "delta")
                    {
                        file >> text;
                        file >> text;
                        delta = stod(text);
                    }
                    else if (text == "control")
                    {
                        file >> text;
                        file >> text;
                        control = stoi(text);
                    }
                    else if (text == "screenW")
                    {
                        file >> text;
                        file >> text;
                        screenSize[0] = stoi(text);
                    }
                    else if (text == "screenH")
                    {
                        file >> text;
                        file >> text;
                        screenSize[1] = stoi(text);
                    }
                    else if (text == "diagnostics")
                    {
                        file >> text;
                        file >> text;
                        diagnostics = stoi(text);
                    }
                    else if (text == "ndt")
                    {
                        file >> text;
                        file >> text;
                        ndt = stoi(text);
                    }
                }
            }

            if (text == "Boundary")
            {
                while (text != "}")
                {
                    file >> text;
                    if (text == "X1")
                    {
                        file >> text;
                        file >> text;
                        param[0] = stod(text);
                    }
                    else if (text == "R1")
                    {
                        file >> text;
                        file >> text;
                        param[1] = stod(text);
                    }
                    else if (text == "X2")
                    {
                        file >> text;
                        file >> text;
                        param[2] = stod(text);
                    }
                    else if (text == "R2")
                    {
                        file >> text;
                        file >> text;
                        param[3] = stod(text);
                    }
                    else if (text == "A")
                    {
                        file >> text;
                        file >> text;
                        param[4] = stod(text);
                    }
                    else if (text == "T")
                    {
                        file >> text;
                        file >> text;
                        param[5] = stod(text);
                    }
                }
                bound.push_back(Boundary(param[0], param[1], param[2], param[3], param[4], param[5]));

            }
            if (text == "/") break;
            
        }
        file.close();
    }
}



int main(int argv, char** args)
{
    //INIT DOMAIN
    input("input.txt");

    //bound = { Boundary(0, Y, X, Y, 26, Tsp), Boundary(X, Y, X, 0, 26, Tsp) };



    if (diagnostics)
    {
        ofstream ofs;
        ofs.open("diag/diag.txt", ofstream::out | ofstream::trunc);
        ofs.close();
    }

    if (UI) {

        SDL_Init(SDL_INIT_VIDEO);
        SDL_Event event;
        SDL_Window* window = SDL_CreateWindow("PIC", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, screenSize[0], screenSize[1], 0);
        SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);

        while (!quit)
        {
            while (SDL_PollEvent(&event) != 0)
            {
                //User requests quit
                if (event.type == SDL_QUIT)
                {
                    //cout << "quit" << endl;
                    quit = true;
                }

                if (control)
                {
                    if (event.type == SDL_MOUSEMOTION)
                    {
                        SDL_GetGlobalMouseState(&xMouse, &yMouse);
                        SDL_GetWindowPosition(window, &xWindow, &yWindow);
                    }
                    if (event.type == SDL_MOUSEWHEEL)
                    {
                        if (event.wheel.y > 0) // scroll up
                        {
                            epsilon += pi * 0.025;
                        }
                        else if (event.wheel.y < 0) // scroll down
                        {
                            epsilon -= pi * 0.025;
                        }

                        cosb = cos(epsilon);
                        sinb = sin(epsilon);
                    }
                    if (event.type == SDL_KEYDOWN)
                    {
                        if (event.key.keysym.sym == SDLK_UP)
                        {
                            EB += 10;
                        }
                        else
                            if (event.key.keysym.sym == SDLK_DOWN)
                            {
                                EB -= 10;
                            }
                            else
                                if (event.key.keysym.sym == SDLK_RIGHT)
                                {
                                    delta += pi * 0.03;
                                }
                                else
                                    if (event.key.keysym.sym == SDLK_LEFT)
                                    {
                                        delta -= pi * 0.03;
                                    }

                    }
                }
                
                
            }
            
            if (control)
            {
                XB = (xMouse - screenSize[2] - xWindow) * X / (screenSize[0] - screenSize[2]);
                YB = (screenSize[1] - yMouse + yWindow) * Y / screenSize[1];
                VB = sqrt(2 * abs(EB) * mh_1 * ev);
                VX = VB * cosb;
                VY = VB * sinb;
            }

            
            //Physics

            //seed = chrono::steady_clock::now().time_since_epoch().count();
            //cout << "generating" << endl;
            if (ions.size() < n && j++ == dtn) {
                generateParticlesVector(ions, ndt, XB, XB, YB - RB, YB + RB, VX, VY, 130, 0, delta); j = 0;
            }
            //cout << "pushing" << endl;

            pushParticles(ions,0);
            pushParticles(ionsS,1);

            t += dt;

            //Render
            //cout << "rendering" << endl;
            if (++k == maxk)
            {

                SDL_RenderClear(renderer);
                drawBoundary(renderer, bound, 0xFF, 0, 0);
                drawParticles(renderer, ions, 0xFF, 0xFF, 0xFF);
                drawParticles(renderer, ionsS, 0, 0xFF, 0);

                
                SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
                SDL_RenderPresent(renderer);

                //cout << "diagnostics" << endl;
                if (diagnostics && ++diag_count == diag_count_max) {
                    diag_count = 0;

                    Diagnostics("");
          
                }

                k = 0;
            }

        }

        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
        ions.clear();
        ionsS.clear();
        bound.clear();

    }
    //cout << "exit" << endl;
    return 0;
}
