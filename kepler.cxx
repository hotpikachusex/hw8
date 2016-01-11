#include<cmath>
#include<fstream>
#include<iostream>

using namespace std;

void Hamiltonian(double& H, double* p,  double* q);

int main(){
 const double dt = 0.05, tEnd = 20*M_PI, e = 0.6;
 const int N = int(tEnd/dt);
 double p[2], q[2], t, H, Q;
 q[0] = 1 - e;
 q[1] = 0;
 p[0] = 0;
 p[1] = sqrt((1+e)/(1-e));
 Hamiltonian(H, p, q);

 ofstream out("kepler.txt");
 out << 0 << "\t" << q[0] << "\t" << q[1] << "\t" << p[0] << "\t" << p[1] << "\t" << H << endl;

 for(int i = 1; i < N; i++){
 t = i * dt;

 Q = q[0] * q[0] + q[1] * q[1];

 p[0] -= dt * pow(Q,-3./2) * q[0];
 p[1] -= dt * pow(Q,-3./2) * q[1];
 q[0] += dt * p[0];
 q[1] += dt * p[1];

 Hamiltonian(H, p, q);

 out << t << "\t" << q[0] << "\t" << q[1] << "\t" << p[0] << "\t" << p[1] << "\t" << H << endl;
 }



 out.close();
 return 0;
}


void Hamiltonian(double& H, double* p, double* q){
 H = 0.5 * (p[0] * p[0] + p[1] * p[1])-1./(sqrt(q[0] * q[0] + q[1] * q[1]));
}