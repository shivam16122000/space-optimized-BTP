#include<bits/stdc++.h>
using namespace std;

#define V_C 30.0
#define R_O 0.0001
#define PI 3.14
#define A_RADIUS 0.9975
#define PARTICLE_DENSITY 1000.0
#define AIR_DENSITY 1.225
#define AIR_VISCOSITY 0.0000181
#define TIME_STEP 0.000001
#define GRAV_CONSTANT 9.8

// helper function defination start
double generateRandomNumber(double x, double y) {
    std::mt19937 gen(time(NULL));
    std::uniform_real_distribution<double> dis(x, y);
    double m = dis(gen);
    return m;
}

double getNorm(vector<double>&v){
    double sum = 0.0;
    for(int i=0;i<v.size();i++){
        sum = sum + (v[i]*v[i]);
    }
    return sqrt(sum);
}

void Log3dVector(string &s,double &x,double &y,double &z){
    cout<<s<<" is "<<x<<" "<<y<<" "<<z<<" "<<endl;
}

void logndVector(string &s,vector<double>& v){
    cout<<s<<" is ";
    for(auto i:v){
        cout<<" "<<i<<" ";
    }
    cout<<endl;
}

vector<double> getRandomInitialVelocityVectorInConeParameter(vector<double> & N, vector<double> & A, double Vc){
   double Nx = N[0];
   double Ny = N[1];
   double Nz = N[2];
   
   double Ax = A[0];
   double Ay = A[1];
   double Az = A[2];

   double ETAx = generateRandomNumber(-1.0,1.0)*((rand()%10)/10.0); 
   double ETAz = generateRandomNumber(-1.0,1.0)*((rand()%10)/10.0);
   double ETAy = generateRandomNumber(0.0,1.0)*((rand()%10)/10.0);
   
   double Ni_x = (Nx + (Ax*ETAx));
   double Ni_y = (Ny + (Ay*ETAy));
   double Ni_z = (Nz + (Az*ETAz));

   vector<double> Ni = {Ni_x,Ni_y,Ni_z};
   
   double NormOfNi = getNorm(Ni);
   
   double ni_x = Ni_x/NormOfNi;
   double ni_y = Ni_y/NormOfNi;
   double ni_z = Ni_z/NormOfNi;
   
   double Vx = Vc*ni_x;
   double Vy = Vc*ni_y;
   double Vz = Vc*ni_z;

   vector<double> VelocityVector = {Vx,Vy,Vz};

   return VelocityVector;
   
}

vector<double> generateRandomPosition(){
    double etaX = generateRandomNumber(-0.0002,0.0002)*((rand()%10)/10.0);
    double etaY = generateRandomNumber(-0.0002,0.0002)*((rand()%10)/10.0);
    double etaZ = generateRandomNumber(-0.0002,0.0002)*((rand()%10)/10.0);

    double RoX = 0.0;
    double RoY = 0.0;
    double RoZ = 2.0;

    double Rx = RoX + etaX;
    double Ry = RoY + etaY;
    double Rz = RoZ + etaZ;

    vector<double>Pos = {Rx,Ry,Rz};
    return Pos;
}

double getRadius(double A,double Ro){
    double eta = generateRandomNumber(-1.0,1.0)*((rand()%10)/10.0);
    double r = Ro*(1+(A*eta));
    return r;
    // return Ro;
}

double getMass(double r,double particleDensity,double pi){
    double mass = (particleDensity)*(4.0/3.0)*(pi)*(r)*(r)*(r);
    return mass;
}

vector<double> getGravForce(double mass,double gravConstant){
       double force = -mass*gravConstant;
       vector<double> GF = {0.0,0.0,force};
       return GF;
}

vector<double>getDragForce(double airDensity,double pi, vector<double> & Vf, vector<double>& Vi, double radius){
   double area = pi*radius*radius;

   double relVx = Vf[0] - Vi[0];
   double relVy = Vf[1] - Vi[1];
   double relVz = Vf[2] - Vi[2];

   vector<double> relV = {relVx,relVy,relVz};
   double relVelocityNorm = getNorm(relV);

   double viscosityCoeff = 0.000018;

   double reynoldNumber = (2.0*radius*airDensity*relVelocityNorm)/viscosityCoeff;

   double Cd;
   if(reynoldNumber<1){
       Cd = 24.0/reynoldNumber;
   }else if(reynoldNumber<400){
       Cd = 24.0/pow(reynoldNumber,0.646);
   }else if(reynoldNumber<300000){
       Cd = 0.5;
   }else if(reynoldNumber<2000000){
       Cd = 0.000366*pow(reynoldNumber,0.4275);
   }else{
       Cd = 0.18;
   }

   double constant = 0.5*airDensity*Cd*area*relVelocityNorm;

   double Fx = constant*relVx;
   double Fy = constant*relVy;
   double Fz = constant*relVz;

   vector<double> DF = {Fx,Fy,Fz};
   return DF;

}

// double getVerticalTerminalVelocity(double airDensity,double dragConstant,double pi, double radius, double mass, double gravConstant){
//          double Vt2 = (2.0*mass*gravConstant)/(airDensity*dragConstant*pi*radius*radius);
//          return -sqrt(Vt2);
// }

vector<double> getNewVelocity(vector<double> &prevVelocity,double timeStep,double mass, vector<double>& gravForce,vector<double>& dragForce){
        double prevVx = prevVelocity[0];
        double prevVy = prevVelocity[1];
        double prevVz = prevVelocity[2];

        double gravForce_x = gravForce[0];
        double gravForce_y = gravForce[1];
        double gravForce_z = gravForce[2];

        double dragForce_x = dragForce[0];
        double dragForce_y = dragForce[1];
        double dragForce_z = dragForce[2];

        double totalForce_x = gravForce_x + dragForce_x;
        double totalForce_y = gravForce_y + dragForce_y;
        double totalForce_z = gravForce_z + dragForce_z;

        double newVx = (prevVx +  ((timeStep/mass)*totalForce_x));
        double newVy = (prevVy +  ((timeStep/mass)*totalForce_y));
        double newVz = (prevVz +  ((timeStep/mass)*totalForce_z));

        vector<double>newVelocity = {newVx,newVy,newVz};
        return newVelocity;
        
}



vector<double> getNewPosition(vector<double> &prevPosition,double timeStep,vector<double> &prevVelocity){
        double prevVx = prevVelocity[0];
        double prevVy = prevVelocity[1];
        double prevVz = prevVelocity[2];

        double prevPositionX = prevPosition[0];
        double prevPositionY = prevPosition[1];
        double prevPositionZ = prevPosition[2];

        double newPosX = (prevPositionX + (timeStep*prevVx));
        double newPosY = (prevPositionY + (timeStep*prevVy));
        double newPosZ = (prevPositionZ + (timeStep*prevVz));

        vector<double> newPosition = {newPosX,newPosY,newPosZ};
        return newPosition;
}

// helper function defination over


// main program
int main(){
       int noOfParticles = 50;
       vector<double>particles_Radius(noOfParticles);
       vector<double>particles_Mass(noOfParticles);
       vector<vector<double>>particles_Curr_Velocity(noOfParticles);
       vector<vector<double>>particles_Curr_Position(noOfParticles);
       vector<vector<double>>particles_Curr_DragForce(noOfParticles);
       vector<vector<double>>particles_GravForce(noOfParticles);
    
       vector<double> Nc = { 0.0,1.0,0.0 };
       vector<double> Ac = { 1.0,0.5,1.0 };
       vector<double> Vf = {0.0,0.0,0.0};

       for(int i=0;i<noOfParticles;i++){
           particles_Radius[i] = getRadius(A_RADIUS,R_O);
           particles_Mass[i] = getMass(particles_Radius[i],PARTICLE_DENSITY,PI);
           particles_Curr_Velocity[i] = getRandomInitialVelocityVectorInConeParameter(Nc,Ac,V_C); 
           particles_Curr_Position[i] = generateRandomPosition();
           particles_Curr_DragForce[i] = getDragForce(AIR_DENSITY,PI,Vf,particles_Curr_Velocity[i],particles_Radius[i]);
           particles_GravForce[i] = getGravForce(particles_Mass[i],GRAV_CONSTANT);
       }
       // { CODE FOR PRINTING DATA IN XYZ FILE };
       ofstream MyFile("simulation.xyz");
       MyFile <<noOfParticles<<endl;
       MyFile<<"current time is "<<0<<" sec"<<endl;
       for(int k=0;k<noOfParticles;k++){
           MyFile<<(k+1)<<" "<<particles_Curr_Position[k][0]<<" "<<particles_Curr_Position[k][1]<<" "<<particles_Curr_Position[k][2]<<" "<<particles_Curr_Velocity[k][0]<<" "<<particles_Curr_Velocity[k][1]<<" "<<particles_Curr_Velocity[k][2]<<endl;
       }

       for(int i=1;i<4000000;i++){
           for(int j=0;j<noOfParticles;j++){
               vector<double>prevVelocity = particles_Curr_Velocity[j];
               particles_Curr_Velocity[j] = getNewVelocity(particles_Curr_Velocity[j],TIME_STEP,particles_Mass[j],particles_GravForce[j],particles_Curr_DragForce[j]);
               particles_Curr_Position[j] = getNewPosition(particles_Curr_Position[j],TIME_STEP,prevVelocity);
               if(particles_Curr_Position[j][2]<=0){
                   particles_Curr_Position[j][2] = 0.0;
                   particles_Curr_Velocity[j] = {0.0,0.0,0.0};
               }
               particles_Curr_DragForce[j] = getDragForce(AIR_DENSITY,PI,Vf,particles_Curr_Velocity[j],particles_Radius[j]);   
           }
           // { CODE FOR PRINTING DATA IN XYZ FILE };
           if(i%1000==0){
                MyFile <<noOfParticles<<endl;
                MyFile<<"current time is "<<i*0.000001<<" sec"<<endl;
                for(int k=0;k<noOfParticles;k++){
                    MyFile<<(k+1)<<" "<<particles_Curr_Position[k][0]<<" "<<particles_Curr_Position[k][1]<<" "<<particles_Curr_Position[k][2]<<" "<<particles_Curr_Velocity[k][0]<<" "<<particles_Curr_Velocity[k][1]<<" "<<particles_Curr_Velocity[k][2]<<endl;
                }
           }

       }

       MyFile.close();

    return 0;
}