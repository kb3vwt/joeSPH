#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include "includes/INIReader.h"

//////////////////////////////////////////////////////////////////////////
/*Classes & Structs                                                     */
//////////////////////////////////////////////////////////////////////////
class SimulationParams
{
public:
  string ini_name;
  int nParticles, nSteps, SizeX, SizeY, SizeZ;
  float InteractionDistance, ParticleMass, BulkModulus, Viscosity, RefMassDensity, Gravity;

  SimulationParams (string fname,bool announce)
  {
    INIReader INIRDR(fname);
    ini_name = fname;
    InteractionDistance = INIRDR.GetReal("ParticleInfo","InteractionDistance",-1);
    ParticleMass = INIRDR.GetReal("ParticleInfo","ParticleMass",-1);
    BulkModulus = INIRDR.GetReal("ParticleInfo","BulkModulus",-1);
    Viscosity = INIRDR.GetReal("ParticleInfo","Viscosity",-1);
    RefMassDensity = INIRDR.GetReal("ParticleInfo","RefMassDensity",-1);
    Gravity = INIRDR.GetReal("WorldInfo","Gravity",-1);
    nParticles = INIRDR.GetInteger("WorldInfo","nParticles",-1);
    nSteps = INIRDR.GetInteger("WorldInfo","nSteps",-1);
    SizeX = INIRDR.GetInteger("WorldInfo","SizeX",-1);
    SizeY = INIRDR.GetInteger("WorldInfo","SizeY",-1);
    SizeZ = INIRDR.GetInteger("WorldInfo","SizeZ",-1);

    if (announce == true)
    {
      std::cout << "Loaded in " + ini_name + " Parameters.\n"
                << "--------------------------\n"
                << "Particle Info:\n"
                << "InteractionDistance: " << InteractionDistance << "\n"
                << "ParticleMass: " << ParticleMass << "\n"
                << "BulkModulus: " << BulkModulus << "\n"
                << "Viscosity: " << Viscosity << "\n"
                << "RefMassDensity: " << RefMassDensity << "\n"
                << "--------------------------\n"
                << "World Info:\n"
                << "nParticles: " << nParticles << "\n"
                << "nSteps: " << nSteps << "\n"
                << "SizeX: " << SizeX << "\n"
                << "SizeY: " << SizeY << "\n"
                << "SizeZ: " << SizeZ << "\n";
    }
  }
};

struct Particle
{
  double pos[3];
  double vel[3];
  double acc[3];
  double force[3];
  double rho;
  double grav[3] = {0,0,-10.0};
};



//////////////////////////////////////////////////////////////////////////
/*Functions                                                             */
//////////////////////////////////////////////////////////////////////////
void RectArrangeParticle(SimulationParams Params,Particle *Particles)
{
  //Region to Uniformly Disperse Particles into: SizeX/4 to 3SizeX/4,
  //                                             SizeY/4 to 3SizeY/4,
  //                                             SizeZ/2 to SizeZ
  int nParticles = Params.nParticles;
  int SizeX = Params.SizeX;
  int SizeY = Params.SizeY;
  int SizeZ = Params.SizeZ;

  //Uniform Spacing (along each axis):
  double SpacingXYZ = cbrt((SizeX*SizeY*SizeZ/8.0)/nParticles);
  std::cout << SpacingXYZ << "\n";//Linear spacing between particles in XYZ directions
  //int ParticleXYZ = SizeX / SpacingXYZ; //Length of side of initial volume in terms of particles

  double xPos = SizeX/4.0;
  double yPos = SizeY/4.0;
  double zPos = SizeZ/2.0;
  for(int i = 0; i<nParticles;i++)
  {
    Particles[i].pos[0] = xPos;
    Particles[i].pos[1] = yPos;
    Particles[i].pos[2] = zPos;

    Particles[i].vel[0] = 0.0;
    Particles[i].vel[1] = 0.0;
    Particles[i].vel[2] = 0.0;

    Particles[i].acc[0] = 0.0;
    Particles[i].acc[1] = 0.0;
    Particles[i].acc[2] = 0.0;

    Particles[i].force[0] = 0.0;
    Particles[i].force[1] = 0.0;
    Particles[i].force[2] = 0.0;

    xPos += SpacingXYZ;
    if (xPos > 3*SizeX/4)
    {
      xPos = SizeX/4;
      yPos += SpacingXYZ;
    }
    if (yPos > 3*SizeY/4)
    {
      yPos = SizeY/4;
      zPos += SpacingXYZ;
    }
  }
}

void Calc_Density(SimulationParams Params,Particle *Particles)
{
  int nParticles = Params.nParticles;
  double m = Params.ParticleMass;
  double h = Params.InteractionDistance;
  for(int i = 0; i<=nParticles;i++) //For particle i
  {
    double summation = 0;
    for(int j = 0; j<=nParticles;j++)
    {
      if(i!=j)
      {
        double distance = sqrt(pow(Particles[i].pos[0]-Particles[j].pos[0],2.0)-pow(Particles[i].pos[1]-Particles[j].pos[1],2.0)-pow(Particles[i].pos[2]-Particles[j].pos[2],2.0));
        if(distance <= h);
        {
          summation += pow(((h*h) - (Particles[i].pos[0]*Particles[i].pos[0] + Particles[i].pos[1]*Particles[i].pos[2] + Particles[i].pos[2])),3.0);
        }
      }
    }
    Particles[i].rho = ((4.0*m) / (3.1415*pow(h,2.0))) * summation;
    //std::cout << Params.nParticles << "after rho_" << i << "\n";
  }
}

void Calc_Forces(SimulationParams Params,Particle *Particles)
{
  int nParticles = Params.nParticles;
  double m = Params.ParticleMass;
  double k = Params.BulkModulus;
  double mu = Params.Viscosity;
  double h = Params.InteractionDistance;

  for(int i = 0; i<nParticles;i++)
  {
    double summation[3] = {0,0,0};
    for(int j = 0; i<nParticles;i++)
    {
      if(i!=j)
      {
        double distance = sqrt(pow(Particles[i].pos[0]-Particles[j].pos[0],2.0)-pow(Particles[i].pos[1]-Particles[j].pos[1],2.0)-pow(Particles[i].pos[2]-Particles[j].pos[2],2.0));
        if(distance <= h);
        {
          double rho_0 = Params.RefMassDensity;
          double rho_i = Particles[i].rho;
          double rho_j = Particles[j].rho;
          double r_ij[3] = {Particles[i].pos[0]-Particles[j].pos[0],Particles[i].pos[1]-Particles[j].pos[1],Particles[i].pos[2]-Particles[j].pos[2]};
          double v_ij[3] = {Particles[i].vel[0]-Particles[j].vel[0],Particles[i].vel[1]-Particles[j].vel[1],Particles[i].vel[2]-Particles[j].vel[2]};
          double q_ij = distance/h;

          summation[0] += (m / (3.1415 * pow(h,4)*rho_j))*(1.0 - q_ij)*((15.0*k*(rho_i+rho_j+2.0*rho_0)*(1.0-q_ij)*r_ij[0]/q_ij)-40.0*mu*v_ij[0]);
          summation[1] += (m / (3.1415 * pow(h,4)*rho_j))*(1.0 - q_ij)*((15.0*k*(rho_i+rho_j+2.0*rho_0)*(1.0-q_ij)*r_ij[1]/q_ij)-40.0*mu*v_ij[1]);
          summation[2] += (m / (3.1415 * pow(h,4)*rho_j))*(1.0 - q_ij)*((15.0*k*(rho_i+rho_j+2.0*rho_0)*(1.0-q_ij)*r_ij[2]/q_ij)-40.0*mu*v_ij[2]);
        }
      }
    }
    Particles[i].force[0] = summation[0];
    Particles[i].force[1] = summation[1];
    Particles[i].force[2] = summation[2];
    std::cout << Params.nParticles << "after f_" << i << "\n";
  }
}

void Calc_Acc(SimulationParams Params,Particle *Particles)
{
  int nParticles = Params.nParticles;

  for(int i = 0; i<nParticles;i++)
  {
    Particles[i].acc[0] = Particles[i].force[0] / Particles[i].rho +Particles[i].grav[0];
    Particles[i].acc[1] = Particles[i].force[1] / Particles[i].rho +Particles[i].grav[1];
    Particles[i].acc[2] = Particles[i].force[2] / Particles[i].rho +Particles[i].grav[2];
  }
}

void OutputCSV(SimulationParams Params,Particle *Particles,int stepNumber)
{
  std::ofstream file;
  std::stringstream filename;
  filename << "step" << std::setw(5) << std::setfill('0') << stepNumber << ".csv";
  file.open(filename.str());
  file << "particleID,Density,xPos,yPos,zPos,xVel,yVel,zVel,xAcc,yAcc,zAcc\n";
  for(int i = 0;i<=Params.nParticles;i++)
  {
    file << i << "," << Particles[i].rho << "," << Particles[i].pos[0] << "," << Particles[i].pos[1] << "," << Particles[i].pos[2] << "," << Particles[i].vel[0] << "," << Particles[i].vel[1] << "," << Particles[i].vel[2] << "," << Particles[i].acc[0] << "," << Particles[i].acc[1] << "," << Particles[i].acc[2] << "\n";
  }
  //file.close();
}

//////////////////////////////////////////////////////////////////////////
/*Main Loop                                                             */
//////////////////////////////////////////////////////////////////////////
int main()
{
  //Initialization
  SimulationParams simParams ("sphcfg.ini",true); //Load INI file
  std::cout << simParams.nParticles << "after ini\n";
  Particle ParticleEnsemble[simParams.nParticles]; //Create Particles
  std::cout << simParams.nParticles << "after particle create\n";
  RectArrangeParticle(simParams,ParticleEnsemble); //Set initial Conditions
  std::cout << simParams.nParticles << "after arrange\n";
  Calc_Density(simParams,ParticleEnsemble);
  std::cout << simParams.nParticles << "after density\n";
  //simParams.nParticles = 10000;
  Calc_Forces(simParams,ParticleEnsemble);
  std::cout << simParams.nParticles << "after forces\n";
  Calc_Acc(simParams,ParticleEnsemble);
  std::cout << simParams.nParticles << "after acc\n";

  OutputCSV(simParams,ParticleEnsemble,0);
  std::cout << simParams.nParticles << "after csv\n";
  return 0;
}
