// Simulation software used in Danner SM, Shevtsova NA, Frigon A, Rybak IA.
// Long propriospinal neurons and gait expression in quadrupeds. eLife. submitted
// and Danner SM, Wilshin SD, Shevtsova NA, Rybak IA. Central control of interlimb
// coordination and speed-dependent gait expression in quadrupeds. J Physiol. 2016;
// 594(23):6947-6967.
//
// Network.cpp
//
#include "Network.h"
#include <iostream>

// Geline compatible with Unix and Windows
std::istream& safeGetline(std::istream& is, std::string& t)
{
    //from http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
    t.clear();
    
    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();
    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
            case '\n':
                return is;
            case '\r':
                if(sb->sgetc() == '\n')
                    sb->sbumpc();
                return is;
            case EOF:
                if(t.empty())
                    is.setstate(std::ios::eofbit);
                return is;
            default:
                t += (char)c;
        }
    }
}

// definition of sech function
inline double sech(double z){return 2/(exp(z)+exp(-z));};
// returns value of d if positive, otherwise 0.0
inline double pos(double d){return std::signbit(d)?0.0:d;};

// print a vector to osteram
std::ostream& operator<<(std::ostream &strm, const myvec &re){
    auto prec= strm.precision();
    strm.precision(20);
    for(int i=0;i<re.size();++i){
        strm << re[i] << "\t";
    }
    strm.precision(prec);
    return strm;
}

// Tests if all elements of myvec v are the same
bool allEqual(myvec v){
    int diff=0;
    for (auto it = v.begin();it!=v.end();++it){
        if(v[0]!=(*it)){
            diff++;
        }
    }
    if (diff==0){
        return true;
    }else{
        return false;
    }
}

// Network constructor without arguments (initialize needs to be called separately)
Network::Network(){
    names = std::map<int,std::string>();
}

// Network constructor with number of NaP neurons (N_NP) and number of regular neurons (N_normal)
Network::Network(int N_NP, int N_normal){
    //Network();
    initialize(N_NP,N_normal);
};

// Initializes Network class with number of NaP neurons (N_NP) and number of regular neurons (N_normal).
void Network::initialize(int N_NP, int N_normal){
    N_NaP=N_NP;
    N_norm=N_normal;
    rNaP= range(N_norm,N_norm+N_NaP-1);
    rNorm= range(0,N_norm-1);
    
    N_Neurons =N_NaP+N_norm;
    
    ELeak = create_scalarV(N_NaP+N_norm,c_ELeak);
    wELeak = create_scalarV(N_NaP+N_norm,c_wELeak);
    ENa = create_scalarV(N_NaP+N_norm,c_ENa);
    ESynE = create_scalarV(N_NaP+N_norm,c_ESynE);
    ESynI = create_scalarV(N_NaP+N_norm,c_ESynI);
    
    
    Cmem = create_scalarV(N_NaP+N_norm,c_Cmem);
    
    mk = create_scalarV(N_NaP+N_norm,c_mk);
    mV12 = create_scalarV(N_NaP+N_norm,c_mV12);
    
    hk = create_scalarV(N_NaP+N_norm,c_hk);
    hV12 = create_scalarV(N_NaP+N_norm,c_hV12);
    htau = create_scalarV(N_NaP+N_norm,c_htau);
    hTauK = create_scalarV(N_NaP+N_norm,c_hTauK);
    hTauV12 = create_scalarV(N_NaP+N_norm,c_hTauV12);
    hTau0 = create_scalarV(N_NaP+N_norm,c_hTau0);
    
    gBarLeak = create_scalarV(N_NaP+N_norm,c_gBarLeak);
    gBarNaP = create_scalarV(N_NaP+N_norm,c_gBarNaP);
    
    Vmax = create_scalarV(N_NaP+N_norm,c_Vmax);
    Vmin = create_scalarV(N_NaP+N_norm,c_Vmin);
    
    tauOutput = create_scalarV(N_Neurons,c_tauOutput);
    
    sigmaNoise = create_scalarV(N_Neurons,c_sigmaNoise);
    tauNoise = create_scalarV(N_Neurons,c_tauNoise);
    
    connE = connection_matrix(N_NaP+N_norm);
    connI = connection_matrix(N_NaP+N_norm);
    
    driveE = std::list<drive>();
    driveI = std::list<drive>();
    
}

// creates myvec of size N with all entries set to k
myvec Network::create_scalarV(int N, double k){
    myvec ret = myvec(N);
    std::fill(ret.begin(),ret.end(),k);
    return ret;
}

// Returns initial conditions if set, otherwise creates and returns dummy IC
myvec Network::genInitialCond() const{
    myvec ret;
    if(initial.size()==N_norm+N_NaP*2){
        ret= initial;
    }else{
        ret = myvec(N_norm+N_NaP*2);
        for (int i=0;i<N_NaP+N_norm;i++){
            ret[i]=ELeak[i];
        }
        for( int i=0;i<N_NaP;++i){
            ret[N_Neurons+i]=(int((i+1)/2)%2)*0.3+0.3;
        }
    }
    return ret;
}

// create a excitatory synaptic connection from neuron from to neuron to with weight w
void Network::setConnE(const int from,const int to,double* w){
    connection con;
    con.from=from;
    con.weight=w;
    connE(to).push_back(con);
}

// create a inhibitory synaptic connection from neuron from to neuron to with weight w
void Network::setConnI(const int from,const int to,double* w){
    connection con;
    con.from=from;
    con.weight=w;
    connI(to).push_back(con);
}

// set excitatory drive to neuron to with weight weight and offset offset
void Network::setDriveE(const int to,  double *weight,  double *offset){
    driveE.push_back(drive(to,weight,offset));
}

// set inhibitors drive to neuron to with weight weight and offset offset
void Network::setDriveI(const int to,  double *weight,  double *offset){
    driveI.push_back(drive(to,weight,offset));
}

// return name of neuron
std::string Network::getName(const int neuronID) const{
    return (*names.find(neuronID)).second;
}

// print network configuration to stream
std::ostream& operator<<(std::ostream& stream, const Network& net){
    stream << "N_NaP " << net.N_NaP << std::endl;
    stream << "N_Normal " << net.N_norm << std::endl;
    stream << std::endl;
    
    stream << "simDuration " << net.simDuration << std::endl;
    //stream << "scalingFactor " << net.scalingFactor << std::endl;
    
    stream << std::endl;
    
    for (int i=0;i<net.N_Neurons;++i){
        stream << "neuron " << i << ": "<<  net.getName(i) << std::endl;
    }
    stream << std::endl;
    
    for(auto it = net.variableMap.begin();it!=net.variableMap.end();++it){
        stream << "variable " << it->first << " " << std::to_string(*(it->second)) << std::endl;
    }
    
    for (int to = 0;to<net.N_Neurons;++to){
        for (auto it=net.connE(to).begin();it!=net.connE(to).end();++it){
            if(*it->weight!=0.0){
                stream << "connectionE " <<  net.getName(it->from) << " -> " << net.getName(to) << " : " << *it->weight  << std::endl;
            }
        }
        for (auto it=net.connI(to).begin();it!=net.connI(to).end();++it){
            if(*it->weight!=0.0){
                stream << "connectionI " << net.getName(it->from) << " -o " << net.getName(to) << " : -" << *it->weight << std::endl;
            }
        }
    }
    stream << std::endl ;
    for (auto it= net.driveE.begin();it!=net.driveE.end();++it){
        if (*it->weight!=0||*it->offset!=0.0){
            stream << "driveE " <<  *it->weight << " * t + " << *it->offset << " -> "  << net.getName(it->to) << std::endl;
        }
    }
    
    for (auto it= net.driveI.begin();it!=net.driveI.end();++it){
        if (*it->weight!=0||*it->offset!=0.0){
            stream << "driveI " << *it->weight << " * t + " << *it->offset << " -o -"  << net.getName(it->to) << std::endl;
        }
    }
    
    stream << std::endl;
    stream << "gLeak \t" << (allEqual(net.gBarLeak) ? myvec(1,net.gBarLeak[0]): net.gBarLeak ) << std::endl;
    stream << "gBarNaP\t " << (allEqual(net.gBarNaP) ? myvec(1,net.gBarNaP[0]): net.gBarNaP ) << std::endl;
    stream << "Eleak\t " << (allEqual(net.ELeak) ? myvec(1,net.ELeak[0]): net.ELeak ) << std::endl;
    stream << "wEleak\t " << (allEqual(net.wELeak) ? myvec(1,net.wELeak[0]): net.wELeak ) << std::endl;
    stream << "ENa\t " << (allEqual(net.ENa) ? myvec(1,net.ENa[0]): net.ENa ) << std::endl;
    stream << "ESynE\t " << (allEqual(net.ESynE) ? myvec(1,net.ESynE[0]): net.ESynE ) << std::endl;
    stream << "ESynI\t " << (allEqual(net.ESynI) ? myvec(1,net.ESynI[0]): net.ESynI ) << std::endl;
    stream << "Cmem\t " << (allEqual(net.Cmem) ? myvec(1,net.Cmem[0]): net.Cmem ) << std::endl;
    stream << "mk\t " << (allEqual(net.mk) ? myvec(1,net.mk[0]): net.mk ) << std::endl;
    stream << "mV12\t " << (allEqual(net.mV12) ? myvec(1,net.mV12[0]): net.mV12 ) << std::endl;
    
    stream << "hk\t " << (allEqual(net.hk) ? myvec(1,net.hk[0]): net.hk ) << std::endl;
    stream << "hV12\t " << (allEqual(net.hV12) ? myvec(1,net.hV12[0]): net.hV12 ) << std::endl;
    
    stream << "htau\t " << (allEqual(net.htau) ? myvec(1,net.htau[0]): net.htau ) << std::endl;
    stream << "hTauK\t " << (allEqual(net.hTauK) ? myvec(1,net.hTauK[0]): net.hTauK ) << std::endl;
    stream << "hTauV12\t " << (allEqual(net.hTauV12) ? myvec(1,net.hTauV12[0]): net.hTauV12 ) << std::endl;
    
    stream << "hTau0\t " << (allEqual(net.hTau0) ? myvec(1,net.hTau0[0]): net.hTau0 ) << std::endl;
    
    stream << "Vmax\t " << (allEqual(net.Vmax) ? myvec(1,net.Vmax[0]): net.Vmax ) << std::endl;
    stream << "Vmin\t " << (allEqual(net.Vmin) ? myvec(1,net.Vmin[0]): net.Vmin ) << std::endl;
    stream << "sigmaNoise\t " << (allEqual(net.sigmaNoise) ? myvec(1,net.sigmaNoise[0]): net.sigmaNoise ) << std::endl;
    stream << "tauNoise\t " << (allEqual(net.tauNoise) ? myvec(1,net.tauNoise[0]): net.tauNoise ) << std::endl;
    
    stream << "initialConditions\t" << net.genInitialCond() << std::endl;
    
    return stream;
}

// create mask of myvec to be used when reading configuration file
void Network::set_para(myvec &to,double value,int start,int end){
    if((start<=end)&&(end<to.size())){
        for (int i = 0;i<to.size();i++){
            if (i>=start&&i<=end){
                to[i]=value;
            }else{
                to[i]=-1234567890;
            }
        }
    }else{
        std::cout << "error in index for parameter" << start << " " << end <<std::endl;
        return;
    }
}
// assigne myvec mask from to myvec to
void Network::assign_para(myvec &to,myvec from){
    if(to.size()==from.size()){
        for(int i=0;i<to.size();i++){
            if(from[i]!=-1234567890){
                to[i]=from[i];
            }
        }
    }
}

// Network constructor from configuration file with path/name filename
Network::Network(std::string filename){
    std::ifstream myfile (filename);
    std::string line;
    
    N_Neurons=-1;
    if (myfile.is_open())
    {
        int NNaP=-1;
        int NNorm=-1;
        for (int i = 0;i<=1;++i){
            std::string::size_type sz;
            getline(myfile,line);
            std::vector<std::string> strs;
            boost::split(strs, line, boost::is_any_of("\t "));
            
            if (strs[0]=="N_NaP"){
                NNaP = std::stoi(strs[1],&sz);
                std::cout << "NNaP = " << NNaP << std::endl;
            }else if(strs[0]=="N_Normal"){
                NNorm = std::stoi(strs[1],&sz);
                std::cout << "NNorm = " << NNorm << std::endl;
            }
        }
        if(NNaP==-1||NNorm==-1){
            std::cout << "N_NaP and/or N_Normal not defined in the first two lines" << std::endl;
            return;
        }
        initialize(NNaP,NNorm);
        int no_neuron = -1;
        while ( safeGetline (myfile,line) )
        {
            
            std::string::size_type sz;
            std::vector<std::string> strs;
            boost::split(strs, line, boost::is_any_of("\t "));
            
            std::vector<std::string>::iterator i = strs.begin();
            while(i != strs.end())
            {
                if((*i)==std::string("")){
                    strs.erase(i);
                }else{
                    ++i;
                }
            }
            if(strs.size()<2){
                continue;
            }
            if (strs[0]=="neuron"){
                //int no = std::stoi(strs[1],&sz);
                if(++no_neuron<N_Neurons){
                    std::string neuronname;
                    if(strs.size()==2){
                        neuronname=strs[1];
                    }else{
                        neuronname=strs[2];
                    }
                    names.emplace(std::make_pair(no_neuron,neuronname));
                    std::cout << "adding neuron " << neuronname << " nr "  << no_neuron   << std::endl;
                }else{
                    std::cout << "max nr of neurons exceeded with " << no_neuron << "neurons"   << std::endl;
                }
            }else if (strs[0]=="variable"){
                double* value = new double;
                *value = std::stod(strs[2],&sz);
                variableMap[strs[1]]=value;
                std::cout << "new variable " << strs[1] << " = " << *value << std::endl;
            }else if (strs[0]=="connectionE"||strs[0]=="connectionI"){
                int from = -1;
                int to = -1;
                double *weight = nullptr;
                for (auto it=names.begin(); it!=names.end(); ++it){
                    if(it->second==strs[1]){
                        from = it->first;
                    }
                    if(it->second==strs[3]){
                        to = it->first;
                    }
                }
                if(isValid<double>(strs[5])){
                    weight = new double;
                    *weight = std::stod(strs[5],&sz);
                }else if(variableMap.find(strs[5])!=variableMap.end()){
                    auto it = variableMap.find(strs[5]);
                    weight = it->second;
                }
                if(from==-1||to==-1){
                    std::cout << "at line " << line << " neuron names were not recognized "
                    << strs[1] << ":" << from << " " << strs[3] << ":" << to << " " <<  std::endl;
                }else{
                    if (strs[0]=="connectionE"){
                        setConnE(from,to,weight);
                        std::cout << "adding excitatory connection from " << from << strs[1] << " to " << to << strs[3] << " w " << *weight  << std::endl;
                    }
                    if (strs[0]=="connectionI"){
                        setConnI(from,to,weight);
                        std::cout << "adding inhibitory connection from " << from << strs[1] << " to " << to << strs[3] << " w " << *weight << std::endl;
                    }
                }
            }else if(strs[0]=="driveE"||strs[0]=="driveI"){
                double *weight = nullptr;
                double *offset = nullptr;
                int to=-1;
                for (auto it=names.begin(); it!=names.end(); ++it){
                    if(it->second==strs[7]){
                        to = it->first;
                    }
                }
                if(isValid<double>(strs[1])){
                    weight = new double;
                    *weight = std::stod(strs[1],&sz);
                }else if(variableMap.find(strs[1])!=variableMap.end()){
                    auto it = variableMap.find(strs[1]);
                    weight = it->second;
                }
                if(isValid<double>(strs[5])){
                    offset = new double;
                    *offset = std::stod(strs[5],&sz);
                }else if(variableMap.find(strs[5])!=variableMap.end()){
                    auto it = variableMap.find(strs[5]);
                    offset = it->second;
                }
                if(strs[0]=="driveE"){
                    setDriveE(to, weight, offset);
                    std::cout << "adding excitatory drive to " << to << strs[7] << " weight " << *weight << " offset " << *offset << std::endl;
                }else if (strs[0]=="driveI"){
                    setDriveI(to, weight, offset);
                    std::cout << "adding inhibitory drive to " << to << strs[7] << " weight " << *weight << " offset " << *offset << std::endl;
                }
            }else if (strs[0]=="simDuration") {
                simDuration=std::stod(strs[1],&sz);
                std::cout << "setting simDuration to " << simDuration << std::endl;
            }else if (strs[0]=="scalingFactor"){
                sf=std::stod(strs[1],&sz);
                std::cout << "setting scalingFactor to " << sf << std::endl;
            }else if (strs[0]=="stepwise"){
                stepwise=true;
                nSteps = std::stoi(strs[1]);
                stepDur = std::stoi(strs[2]);
                
                if (strs.size()>=5){
                    alphamin = std::stod(strs[3]);
                    alphamax = std::stod(strs[4]);
                }else{
                    alphamax = std::stod(strs[3]);
                }
                std::cout << "turning on stepwise calculation. alpha from " << alphamin << " to " << alphamax << " in " << nSteps << " steps of "<< stepDur << " ms" <<  std::endl;
            }else if(strs[0][0]=='/'){
                // do nothing
                
            }else if (strs[0]=="feedbackAflex"||strs[0]=="feedbackAext"||strs[0]=="feedbackB"){
            }else if(strs[0]=="initialConditions"){
                initial =myvec(N_norm+N_NaP*2);
                for (int j = 1;j<=N_norm+N_NaP*2;++j){
                    if (j>=strs.size()){
                        std::cout<< "not all initial conditions specified" << std::endl;
                        break;
                    }
                    initial[j-1]=std::stod(strs[j],&sz);
                }
            }else{
                myvec para;
                
                if (strs.size()==2){
                    para = create_scalarV(N_NaP+N_norm,std::stod(strs[1],&sz));
                }else if (strs.size()==3){
                    std::vector<std::string> substrs;
                    boost::split(substrs, strs[1], boost::is_any_of("[]:"));
                    auto iter = substrs.begin();
                    while(iter != substrs.end())
                    {
                        if((*iter)==std::string("")){
                            substrs.erase(iter);
                        }else{
                            ++iter;
                        }
                    }
                    if(substrs.size()==2){
                        para = myvec(N_Neurons);
                        set_para(para,std::stod(strs[2],&sz),std::stod(substrs[0],&sz),std::stod(substrs[1],&sz));
                        std::cout<< "update " << strs[0] << " for neurons " << std::stod(substrs[0],&sz) << " : "
                        << std::stod(substrs[1],&sz) << " to " << std::stod(strs[2],&sz) << std::endl;
                    }else{
                        std::cout<< "syntax error in specification for " << strs[0] << std::endl;
                        break;
                    }
                }else if (strs.size()==4){
                    if(strs[1]=="neurons"){
                        para = myvec(N_Neurons,-1234567890);
                        std::string to;
                        for (auto it=names.begin(); it!=names.end(); ++it){
                            if(it->second.find(strs[2])!=std::string::npos){
                                para[it->first]=std::stod(strs[3],&sz);
                                to+=std::to_string(it->first)+", ";
                            }
                        }
                        std::cout << "update " << strs[0] << " for neurons " << to << "(*" << strs[2] << "*) to "
                        << std::stod(strs[3],&sz) << std::endl;
                    }else{
                        std::cout<< "syntax error in specification for " << strs[0] << std::endl;
                    }
                }else if (strs.size()>=N_NaP+N_norm+1){
                    para = myvec(N_Neurons);
                    for (int j = 1;j<=N_NaP+N_norm;j++){
                        if (j>=strs.size()){
                            std::cout<< "not enough values specified for parameter " << strs[0] << std::endl;
                            break;
                        }else{
                            para(j-1)=std::stod(strs[j],&sz);
                        }
                    }
                }
                if (strs[0]=="gLeak"){
                    assign_para(gBarLeak,para);
                    std::cout << "set gLeak to " << (allEqual(gBarLeak) ? myvec(1,gBarLeak[0]): gBarLeak) << std::endl;
                }else if(strs[0] == "gBarNaP"){
                    assign_para(gBarNaP,para);
                    std::cout << "set gBarNaP to " << (allEqual(gBarNaP) ? myvec(1,gBarNaP[0]): gBarNaP) << std::endl;
                }else if(strs[0] == "Eleak"){
                    assign_para(ELeak, para);
                    std::cout << "set Eleak to " << (allEqual(ELeak) ? myvec(1,ELeak[0]): ELeak) << std::endl;
                }else if(strs[0] == "wEleak"){
                    assign_para(wELeak,para);
                    std::cout << "set wEleak to " << (allEqual(wELeak) ? myvec(1,wELeak[0]): wELeak) << std::endl;
                }else if(strs[0] == "ENa"){
                    assign_para(ENa,para);
                    std::cout << "set ENa to " << (allEqual(ENa) ? myvec(1,ENa[0]): ENa) << std::endl;
                }else if(strs[0] == "ESynE"){
                    assign_para(ESynE,para);
                    std::cout << "set ESynE to " << (allEqual(ESynE) ? myvec(1,ESynE[0]): ESynE) << std::endl;
                }else if(strs[0] == "ESynI"){
                    assign_para(ESynI,para);
                    std::cout << "set ESynI to " << (allEqual(ESynI) ? myvec(1,ESynI[0]): ESynI) << std::endl;
                }else if(strs[0] == "Cmem"){
                    assign_para(Cmem,para);
                    std::cout << "set Cmem to " << (allEqual(Cmem) ? myvec(1,Cmem[0]): Cmem) << std::endl;
                }else if(strs[0] == "mk"){
                    assign_para(mk,para);
                    std::cout << "set mk to " << (allEqual(mk) ? myvec(1,mk[0]): mk) << std::endl;
                }else if(strs[0] == "mV12"){
                    assign_para(mV12,para);
                    std::cout << "set mV12 to " << (allEqual(mV12) ? myvec(1,mV12[0]): mV12) << std::endl;
                }else if(strs[0] == "hk"){
                    assign_para(hk,para);
                    std::cout << "set hk to " << (allEqual(hk) ? myvec(1,hk[0]): hk) << std::endl;
                }else if(strs[0] == "hV12"){
                    assign_para(hV12,para);
                    std::cout << "set hV12 to " << (allEqual(hV12) ? myvec(1,hV12[0]): hV12) << std::endl;
                }else if(strs[0] == "htau"){
                    assign_para(htau,para);
                    std::cout << "set htau to " << (allEqual(htau) ? myvec(1,htau[0]): htau) << std::endl;
                }else if(strs[0] == "hTauK"){
                    assign_para(hTauK,para);
                    std::cout << "set hTauK to " << (allEqual(hTauK) ? myvec(1,hTauK[0]): hTauK) << std::endl;
                }else if(strs[0] == "hTau0"){
                    assign_para(hTau0,para);
                    std::cout << "set hTau0 to " << (allEqual(hTau0) ? myvec(1,hTau0[0]): hTau0) << std::endl;
                }else if(strs[0] == "Vmax"){
                    assign_para(Vmax,para);
                    std::cout << "set Vmax to " << (allEqual(Vmax) ? myvec(1,Vmax[0]): Vmax) << std::endl;
                }else if(strs[0] == "Vmin"){
                    assign_para(Vmin,para);
                    std::cout << "set Vmin to " << (allEqual(Vmin) ? myvec(1,Vmin[0]): Vmin) << std::endl;
                }else if(strs[0] == "hTauV12"){
                    assign_para(hTauV12,para);
                    std::cout << "set hTauV12 to " << (allEqual(hTauV12) ? myvec(1,hTauV12[0]): hTauV12) << std::endl;
                }else if(strs[0] == "sigmaNoise"){
                    assign_para(sigmaNoise,para);
                    std::cout << "set sigmaNoise to " << (allEqual(sigmaNoise) ? myvec(1,sigmaNoise[0]): sigmaNoise) << std::endl;
                }else if(strs[0] == "tauNoise"){
                    assign_para(tauNoise,para);
                    std::cout << "set tauNoise to " << (allEqual(tauNoise) ? myvec(1,tauNoise[0]): tauNoise) << std::endl;
                }
            }
        }
        while(no_neuron<N_Neurons-1){
            names.emplace(std::make_pair(++no_neuron,"not initialized"));
        }
        myfile.close();
        
        randomGen = OrnsteinUhlenbeck(NNaP+NNorm,0.0,1.0,tauNoise[0]);
    }
}

//integration step
void Network::step(const myvec &x, myvec &dxdt, double t){
    double alpha = calcAlpha(t);
    myvec transV = myvec(N_Neurons);
    
    for (int i = 0;i<N_Neurons;++i){
        transV(i) = std::min(1.0,std::max(0.0,(x[i]-Vmin[i])/(Vmax[i]-Vmin[i])));
    }
    for (int i = 0;i<N_Neurons;++i){
        if(wELeak[i]==0){
            dxdt[i]=(-gBarLeak[i]*(x[i]-ELeak[i]));
        }else if(wELeak[i]>0){
            dxdt[i]=(-gBarLeak[i]*(x[i]-(ELeak[i]+alpha*wELeak[i]*sf*1e5)));
        }else{
            dxdt[i]=(-gBarLeak[i]*(x[i]-(ELeak[i]*(1.0-alpha))));
        }
        for (auto it=connE(i).begin();it!=connE(i).end();++it){
            dxdt[i]-=pos(*it->weight)*transV(it->from)*(x[i]-ESynE[i]);
            
        }
        for (auto it=connI(i).begin();it!=connI(i).end();++it){
            dxdt[i]-=pos(*it->weight)*transV(it->from)*(x[i]-ESynI[i]);
        }
        dxdt[i]+=-(randomGen.get(i,t)*sigmaNoise[i]);
        dxdt[i]/=Cmem[i];
        if (i < N_NaP){
            double mpInf = 1./(1.+exp((x[i]-mV12[i])/mk[i]));
            double hpInf = 1./(1.+exp((x[i]-hV12[i])/hk[i]));
            double tau_inf = hTau0[i]+(htau[i]-hTau0[i])/(cosh((x[i]-hTauV12[i])/hTauK[i]));
            
            dxdt[i]-=gBarNaP[i]*x[N_NaP+N_norm+i]*mpInf*(x[i]-ENa[i])/Cmem[i];
            dxdt[N_NaP+N_norm+i]=(hpInf-x[N_NaP+N_norm+i])/tau_inf;
        }
    }
    
    for (auto it= driveE.begin();it!=driveE.end();++it){
        dxdt[it->to]-=((alpha*(*it->weight)*1e5)+(*it->offset))*(x[it->to]-ESynE[it->to])/Cmem[it->to];
    }
    
    for (auto it= driveI.begin();it!=driveI.end();++it){
        dxdt[it->to]-=((alpha*(*it->weight)*1e5)+(*it->offset))*(x[it->to]-ESynI[it->to])/Cmem[it->to];
    }
    
}

// get alpha depending on time t
double Network::calcAlpha(double t){
    if(stepwise){
        double tbar = t;
        double m=(nSteps+.5)*stepDur;
        if(tbar>m){
            tbar=2*m-tbar;
        }
        return alphamin+std::floor(tbar/(double)stepDur)*(alphamax-alphamin)/(double)nSteps;
    }else{
        return alphaSet;
    }
}

// update variable name to value
bool Network::updateVariable(std::string name, double value){
    for(auto it = variableMap.begin();it!=variableMap.end();++it){
        if(name==it->first){
            (*it->second)=value;
            return true;
        }
    }
    return false;
}

// function is called every timestep by odeint to record the current state of the network
void Observer::operator()( const myvec &x , const double t )
{
    myvec* y = new myvec(x.size()+1);
    (*y)(0)=t;
    for(int i = 1; i < y->size();++i){
        (*y)(i)=x(i-1);
    }
    (*state_rec).insert((*state_rec).end(),y);
}

// Simulator constructor
Simulator::Simulator(Network *network){
    net=network;
    state_rec=std::list<myvec*>();
    observer = new Observer();
    observer->set_state_rec(&state_rec);
    sys = ode_system();
    sys.net = net;
    current_state=net->genInitialCond();
}

// output simulation results to text
void Simulator::save_list_to_text(std::string postfix=""){
    std::ofstream myfile;
    myfile.open ("./results/example"+postfix+".txt",std::ios::out);
    myfile << "time\t";
    for (int i = 0; i<net->N_Neurons;i++){
        myfile << net->names.at(i) << "\t";
    }
    for (int i = 0; i<net->N_NaP;i++){
        myfile << "h_" << net->names.at(i) << "\t";
    }
    myfile << std::endl;
    for(auto ci = state_rec.begin(); ci!= state_rec.end(); ++ci){
        myfile << **ci << std::endl;
    }
    
    myfile.close();
    std::cout << "file written" << std::endl;
}

// save simulation results to CDF4 file with filename ofliename
void Simulator::save_list_to_cdf(std::string ofilename){
    int NY = (int)(*(state_rec.begin()))->size();
    int NX = (int)state_rec.size();
    std::cout << NX << " " << NY << std::endl;
    float *data_out = new float[NX*NY];
    
    int i=0;
    for(std::list<myvec*>::const_iterator ci = state_rec.begin(); ci!= state_rec.end(); ++ci){
        myvec re =  **ci;
        for(int j=0;j<re.size();++j){
            data_out[i*NY+j]=(float)re[j];
        }
        state_rec.erase(ci);
        ++i;
    }
    try{
        netCDF::NcFile dataFile(ofilename, netCDF::NcFile::replace);
        netCDF::NcDim xDim = dataFile.addDim("x", NX);
        netCDF::NcDim yDim = dataFile.addDim("y", NY);
        std::vector<netCDF::NcDim> dims;
        dims.push_back(xDim);
        dims.push_back(yDim);
        netCDF::NcVar data = dataFile.addVar("data", netCDF::ncFloat, dims);
        data.putVar(data_out);
    }catch(netCDF::exceptions::NcException& e){
        e.what();
        return ;
    }
    delete [] data_out;
}

// run simulation for simDuration specified in configuration file
void Simulator::run(){
    run(net->simDuration);
}

// run simulation for dur ms
void Simulator::run(double dur){
    net->settingPeriod = dur*0.125;
    typedef runge_kutta_fehlberg78< myvec > error_stepper_type;
    typedef runge_kutta_dopri5<myvec > error_stepper_type2;
    runge_kutta4<myvec> rk;
    auto controlled_stepper = make_controlled(1e-5,1e-5,error_stepper_type());
    myvec x = net->genInitialCond();
    
    std::cout << x << std::endl;
    myvec dxdt = myvec(x.size());
    
    if(net->stepwise){
        net->randomGen.calculateUntil((double)((net->nSteps*2+1)*net->stepDur));
        integrate_const( rk , sys , x , 0.0 , 1.0, 0.0001 );
        //std::cout << x << std::endl;
        integrate_adaptive( controlled_stepper, sys, x , 0.0 , (double)(net->nSteps*net->stepDur*2), 0.00001,*observer);
    }
    delete observer;
}

// run single simulation step with delta t of dt.
void Simulator::run_dt(double dt){
    if(t==0.0){
        (*observer)(current_state,t);
    }
    rkf78.do_step(sys,current_state,t,dt);
    t+=dt;
    (*observer)(current_state,t);
}


// precalculate Noise for t ms
void OrnsteinUhlenbeck::calculateUntil(double t){
    data=std::vector<myvec>((int)t+2);
    data[0]=myvec(N);
    std::default_random_engine generator;
    typedef std::chrono::high_resolution_clock myclock;
    unsigned seed= (unsigned)myclock::now().time_since_epoch().count();
    generator.seed(seed);
    std::normal_distribution<double> distribution(0.0,1.0);
    for(int i = 1; i<=t+1;++i){
        myvec x=data[i-1];
        data[i]=myvec(N);
        for(int j=0;j<N;++j){
            (data[i])[j]=x[j]+(mu-x[j])/tau+sigma*(std::sqrt((2./tau)))*distribution(generator);
        }
    }
}

// get noise for neuron i at time t
double OrnsteinUhlenbeck::get(int i, double t){
    int time = std::floor(t);
    if (time < 0 || time > data.size()){
        return 0.0;
    }
    return data[time][i]+(data[time+1][i]-data[time][i])*(t-std::floor(t));
}
