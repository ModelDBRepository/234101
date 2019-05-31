// Simulation software used in Danner SM, Shevtsova NA, Frigon A, Rybak IA.
// Long propriospinal neurons and gait expression in quadrupeds. eLife. submitted
// and Danner SM, Wilshin SD, Shevtsova NA, Rybak IA. Central control of interlimb
// coordination and speed-dependent gait expression in quadrupeds. J Physiol. 2016;
// 594(23):6947-6967.
//
// Network.h
//
#ifndef __iCPG__Network__
#define __iCPG__Network__

#include <stdio.h>
#include <list>
#include <fstream>
#include <string>
#include <map>
#include <chrono>
#include <random>
#include <netcdf>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/algorithm/string.hpp>

#include <ctime>
#include <algorithm>

#include <math.h>

using namespace boost::numeric::odeint;
using namespace boost::numeric::ublas;

// default paramter values
const double c_gBarLeak = 2.8;
const double c_gBarNaP = 5.0;
const double c_ELeak = -60.0;
const double c_wELeak = 0.0;
const double c_ENa = 50.0;
const double c_ESynE = -10.0;
const double c_ESynI = -75.0;
const double c_Cmem = 10.0;
const double c_mV12 = -40;
const double c_mk = -6.5;
const double c_hV12 = -47;
const double c_hk = 5.4;
const double c_htau = 110.0;
const double c_hTauK = 50.0;
const double c_hTauV12 = c_hV12;
const double c_hTau0 = 50.0;
const double c_Vmin = -50.0;
const double c_Vmax = 0.0;
const double c_tauOutput = 5.;
const double c_sigmaNoise = 0.005;
const double c_tauNoise = 10.0;

// check if string num can be cast to type T
template<typename T> bool isValid(const std::string& num) {
    bool flag = true;
    try {
        boost::lexical_cast<T>(num);
    }
    catch (boost::bad_lexical_cast &e) {
        flag = false;
    }
    return flag;
}

// vector type
typedef vector<double,std::vector<double>> myvec;

// struct representing synaptic connection
struct connection{
    int from;
    double *weight; //offset
};

// sparse connection matrix
typedef vector<std::list<connection>> connection_matrix;

// OrnsteinUhlenbeck noise
class OrnsteinUhlenbeck{
private:
    double mu=0.0;
    double sigma=1.0;
    double tau=1.0;
    int N=1;
    std::vector<myvec> data;
public:
    OrnsteinUhlenbeck(){};
    OrnsteinUhlenbeck(int n, double m, double s, double t){N=n;mu=m;sigma=s;tau=t;};
    void calculateUntil(double t);
    double get(int i,double t);
};

// sturcture describing a drive-connection
struct drive{
    int to;         // index of neuron
    double *weight; // pointer to connection weight
    double *offset; // pointer to offset
    drive(int toc, double *weightc, double *offsetc){
        to=toc;
        weight=weightc;
        offset=offsetc;
    }
};

// Observer used to record the state of the simulation with odeint
struct Observer {
    std::list<myvec*> *state_rec;
    
    Observer(){}
    void set_state_rec(std::list<myvec*>* sr){state_rec=sr;};
    void operator()(const myvec&, const double);
};

// Core Network class
class Network{
public:
    
    std::map<int,std::string> names;
    myvec initial;
    myvec gBarLeak;
    myvec gBarNaP;
    
    myvec ELeak;
    myvec wELeak;
    
    myvec ENa;
    myvec ESynE;
    myvec ESynI;
    
    myvec Cmem;
    
    myvec mk;
    myvec mV12;
    
    myvec hk;
    myvec hV12;
    myvec htau;
    myvec hTauK;
    myvec hTauV12;
    myvec hTau0;
    
    myvec Vmax;
    myvec Vmin;
    
    myvec tauOutput;
    
    myvec sigmaNoise;
    myvec tauNoise;
    
    int N_NaP;
    int N_norm;
    
    int iNaP=0;
    int iNorm=0;
    
    int simDirection=1;
    
    range rNaP;
    range rNorm;
    
    connection_matrix connE;
    connection_matrix connI;
    
    std::list<drive> driveE;
    std::list<drive> driveI;
    
    double last_t=-1;
    
    
    double simDuration = 100000;
    double settingPeriod = 10000;
    double scalingFactor = 1.0;
    double sf=1.0;
    
    std::map<std::string,double*> variableMap;
    
    bool stepwise = false;
    int nSteps = 50;
    int stepDur = 1000;
    
    OrnsteinUhlenbeck randomGen;
    int N_Neurons;
    double alphamin = 0.0;
    double alphamax = 0.93;
    double alphaSet=0.0;
    
    Network();
    Network(int N_NP, int N_normal);
    Network(std::string /*filename*/);
    
    void initialize(int /*N_NP*/, int /*N_normal*/);
    
    static myvec create_scalarV(int N, double k);
    static void set_para(myvec &to,double value,int start,int end);
    static void assign_para(myvec &to,myvec from);
    
    myvec genInitialCond() const;
    
    void setConnE(const int from,const int to,double *w);
    void setConnI(const int from,const int to,double *w);
    
    void setDriveE(const int to,  double *weight,  double *offset);
    void setDriveI(const int to,  double *weight,  double *offset);
    
    std::string getName(const int /*neuronID*/) const;
    
    friend std::ostream& operator<<(std::ostream&,const Network&);
    
    //integration step
    void step(const myvec &, myvec &, double);
    
    double calcAlpha(double /*t*/);
    void setAlpha(double as){alphaSet=as;stepwise=false;};
    double getAlpha(){return alphaSet;};
    
    bool updateVariable(std::string, double);
};

// struct supplied to the odeint solver
struct ode_system
{
    Network *net;
    void operator()( const myvec &x , myvec &dxdt , double t){
        (*net).step(x,dxdt,t);
    }
};


// class handling Network simulation
class Simulator{
public:
    std::list<myvec*> state_rec;
    Network* net;
    Observer* observer;
    ode_system sys;
    myvec current_state;
    runge_kutta_fehlberg78< myvec > rkf78;
    double t =0.0;
    
    Simulator(Network *);
    void save_list_to_text(std::string );
    void save_list_to_cdf(std::string );
    void run();
    void run(double /*dur*/);
    void run_dt(double dt);
    
};

#endif /* defined(__iCPG__Network__) */
