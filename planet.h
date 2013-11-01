#ifndef PLANET_H
#define PLANET_H
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <sstream>

using namespace std;
class Planet
{
public:
        string name;
        double mass;
        double X,Y,Vx,Vy,Time,Fx,Fy,R, dVx, dVy,K,U,E,L, Aphelion, Perihelion,R2;
        double CMX, CMY, totalmass;
    Planet();
    Planet(string name, double mass, double radius, double angle, double Vo);
    double Radius();
    double Velocity();
    double Distance(Planet planet);
    void Values(int n, int i, Planet *System);
    void Coordinates(int n, Planet *System);

};

#endif // PLANET_H
