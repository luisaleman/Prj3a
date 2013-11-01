#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include "planet.h"
#include <sstream>

using namespace std;

int main()
{
    //const double G=1;
    double Tmax,time;
    int n_steps;
    cout << "Set the maiximum time (years): ";
    cin >> Tmax;
    cout << "Number of steps ";
    cin >> n_steps;
    int n_objects=10;                                    //one planet only
    time=0;

    double h= Tmax/(double)n_steps;
    cout << "Step length = " << h << "years"<< endl<< endl;

    Planet Sun= Planet("Sun", 1, 0, 0,0);
    Planet Earth= Planet("Earth",3e-6,1,20,6.2609472); //6.28318    8.8857659  6.2609472
    Planet Jupiter= Planet("Jupiter",9.5e-4,5.2,285,2.7477737); //2.7477737 285
    Planet Mars= Planet("Mars",3.3e-7,1.52,250,5.06194848);
    Planet Venus= Planet("Venus", 2.451e-6,0.72,45, 7.3626049);
    Planet Mercury= Planet("Mercury", 1.2e-6,0.39,350, 10.0641888);
    Planet Saturn= Planet("Saturn", 2.75e-4,9.54,25, 2.0372256);
    Planet Uranus=Planet("Uranus",4.4e-5,19.19,75,1.4317344);
    Planet Neptune= Planet("Neptune", 5.15e-5,30.06,295, 1.1416032);
    Planet Pluto= Planet("Pluto",6.55e-9,39.53,295,0.988128);

    Planet *System= new Planet[n_objects];
    System[0]= Sun;
    System[1]= Earth;
    System[2]=Jupiter;
    System[3]=Mercury;
    System[4]=Venus;
    System[5]=Mars;
    System[6]=Saturn;
    System[7]=Uranus;
    System[8]=Neptune;
    System[9]=Pluto;

    ofstream *myfile= new ofstream[n_objects];
    for(int i=0;i<n_objects;i++){
        string output;
        stringstream str1,str2;
        output=System[i].name;
        str1<<Tmax;
        str2<<n_steps;
        output+=str1.str();
        output+="years";
        output+=str2.str();
        output+="steps.txt";
        myfile[i].open(output.c_str());
    }

    for(int i=0; i<n_objects; i++){
        System[i].Coordinates(n_objects, System);
    }
    double **Ax= new double*[n_objects];
    double **Ay= new double*[n_objects];
    double **Vx= new double*[n_objects];
    double **Vy= new double*[n_objects];
    for(int j=0;j<n_objects;j++){
        Ax[j]= new double[4];
        Ay[j]= new double[4];
        Vx[j]= new double[4];
        Vy[j]= new double[4];
    }
    double *Xpos= new double[n_objects];
    double *Ypos= new double[n_objects];
    double hh;
    double h6=h/6;

    for(int i=0;i<n_objects;i++){
        System[i].Values(n_objects, i, System);
    }

    while(time<Tmax){
        time +=h;
        //Initial values
        for(int j=0;j<n_objects;j++){
            System[j].Values(n_objects,j, System);
            Ax[j][0]=System[j].Fx/System[j].mass;
            Ay[j][0]=System[j].Fy/System[j].mass;
            Vx[j][0]=System[j].Vx;
            Vy[j][0]=System[j].Vy;
            Xpos[j]=System[j].X;
            Ypos[j]=System[j].Y;
        }
        //Steps 1-3 Runge-Kutta
        for(int j=0;j<3;j++){
            for(int i=0; i<n_objects;i++){
                if(j<2) {hh = h/2.0;} else {hh= h;};
                System[i].Y=Ypos[i]+Vy[i][j]*hh;
                System[i].X=Xpos[i]+Vx[i][j]*hh;
                Vx[i][j+1]=Vx[i][0]+hh*Ax[i][j];
                Vy[i][j+1]=Vy[i][0]+hh*Ay[i][j];
            }
            for(int i=0; i<n_objects;i++){
                System[i].Values(n_objects, i, System);
                Ax[i][j+1]=System[i].Fx/System[i].mass;
                Ay[i][j+1]=System[i].Fy/System[i].mass;
            }
        }

        //Step 4 Runge-Kutta
        for(int i=0; i<n_objects;i++){
            System[i].Vx=Vx[i][0]+h6*(Ax[i][0]+2*Ax[i][1]+2*Ax[i][2]+Ax[i][3]);
            System[i].Vy=Vy[i][0]+h6*(Ay[i][0]+2*Ay[i][1]+2*Ay[i][2]+Ay[i][3]);
            System[i].X= Xpos[i]+h6*(Vx[i][0]+2*Vx[i][1]+2*Vx[i][2]+Vx[i][3]);
            System[i].Y= Ypos[i]+h6*(Vy[i][0]+2*Vy[i][1]+2*Vy[i][2]+Vy[i][3]);
        }
        for(int i=0; i<n_objects;i++){
            myfile[i] << setw(20) << time;
            myfile[i] << setw(20) << System[i].X;
            myfile[i] << setw(20) << System[i].Y;
            myfile[i] << setw(20) << System[i].Vx;
            myfile[i] << setw(20) << System[i].Vy;
            myfile[i] << setw(20) << System[i].K;
            myfile[i] << setw(20) << System[i].U;
            myfile[i] << setw(20) << System[i].L;
            myfile[i] << endl;
        }
    }

    for(int i=0; i<n_objects;i++){
            myfile[i].close();
        }

    for(int i=0; i<n_objects;i++){
        cout << System[i].name << endl;
        cout << "Vx: " << System[i].Vx << " Vy "<<System[i].Vy << endl ;
        cout << "Distance to CM : " << System[i].R << endl;
        cout << "Perihelion : "<< System[i].Perihelion << setw(5) << " Aphelion : " << System[i].Aphelion << endl;
        cout << "L : "<< System[i].L << setw(5) << "K : "<< System[i].K<< endl;
        cout << "U : "<< System[i].U<< endl<<endl;
    }

}


