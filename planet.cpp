#include "planet.h"
Planet::Planet()
{

}

//constructor with parameter
Planet::Planet(string name, double mass, double Ro, double angle, double Vo)
{
    const double PI = 3.141592653589793;
    //const double AU = 1.496e11;

    Time=0;
    angle *= PI/180;
    X= Ro*cos(angle);
    Y= Ro*sin(angle);
    Vx= -Vo*sin(angle);
    Vy= Vo*cos(angle);
    this->mass= mass;
    this->name= name;
    //Vx= dVx;
    //Vy=dVy;
    Aphelion=Ro;
    Perihelion=Ro;

}
//to calculate the radius
double Planet::Radius()
{
    return sqrt(X*X + Y*Y);
}

 //to calculate the velocity
double Planet::Velocity()
{
    return sqrt(Vx*Vx + Vy*Vy);
}

// to calculate the distance between two planets
double Planet::Distance(Planet planet)
{
    return sqrt(pow(X-planet.X,2)+pow(Y-planet.Y,2));
}

void Planet::Values(int n, int i, Planet *System)
{
    const double PI=3.141592653589793;
    const double G= 4*PI*PI;
    Fx=0;
    Fy=0;
    U=0;
    K=0.5*mass*(Vx*Vx + Vy*Vy);
    R=sqrt(X*X+Y*Y);
    for(int j=0;j<n;j++){
        if(j!=i){
        R2=sqrt(pow(System[j].X-X,2)+pow(System[j].Y-Y,2));
        Fx+=System[j].mass*mass*(System[j].X-X)/pow(R2,3);
        Fy+=System[j].mass*mass*(System[j].Y-Y)/pow(R2,3);
        U-=mass*System[j].mass/R2;
        }
    }
    Fx*=G;
    Fy*=G;
    U*=G;
    E=K+U;
    L=(Vy*X-Vx*Y)/(sqrt((Vx*Vx + Vy*Vy))*R);
    if(Aphelion<R){
        Aphelion=R;
    }
    if(Perihelion>R){
        Perihelion=R;
    }
}

void Planet::Coordinates(int n, Planet *System){
    double cmx=0;
    double cmy=0;
    totalmass=0;
    CMX=0;
    CMY=0;
    for(int i=0;i<n;i++){
        cmx+=System[i].mass*System[i].X;
        cmy+=System[i].mass*System[i].Y;
        totalmass+=System[i].mass;
    }
    CMX=cmx/totalmass;
    CMY=cmy/totalmass;
    X-=CMX;
    Y-=CMY;
}

