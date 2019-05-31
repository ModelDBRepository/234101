// Simulation software used in Danner SM, Shevtsova NA, Frigon A, Rybak IA.
// Long propriospinal neurons and gait expression in quadrupeds. eLife. submitted
// and Danner SM, Wilshin SD, Shevtsova NA, Rybak IA. Central control of interlimb
// coordination and speed-dependent gait expression in quadrupeds. J Physiol. 2016;
// 594(23):6947-6967.
//
// usage: executable -f config_file [-o output_file] [-u name value] [-a alpha]
//                                  [-U varname tstep value1 value2 [tstep2] [-V varname]]
// -f config_file: text file specifying the neural network model
// -o output_file: the path where the simulationresults should be written to. The file will
//                  be written in CDF-4 file format.
// -u name value:  updates a variable (name) specified in the config_file to the value given
// -a alpha: sets alpha to a constant value (overrides the configuration file)
// -U varname tstep value1 value2 [tstep2]: sets variable varname to value1 and changes
//                  it to value2 at time tstep, if tstep2 is specified, the variable will be
//                  changed back to value1 at time tstep2
// -V varname: adds an additional variable to the update process of -U (has no effect if -U
//                  is not specified.
#include "Network.h"
#include <fstream>
#include <string>

int main(int argc, char **argv)
{
    // if no output filename is specified with [-o filename], ./results/example.cdf will be used
    std::string ofilename = "./results/example.cdf";
    time_t tstart, tend;
    tstart = time(0);
    Network *network;
    bool varyVars=false;
    std::vector<std::string> varnames;
    double tstep = 0.0;
    double valuefirst = 0.0;
    double valuesecond = 0.0;
    double tstep2 = 999999999.0;
    if (std::string(argv[1])=="-f" && argc >=3){
        network=new Network(std::string(argv[2]));
        if(network->N_Neurons!=-1){
            
            Simulator sim = Simulator(network);
            int ia = 3;
            while(ia<argc){
                if(std::string(argv[ia])=="-u" && argc > ia+2){
                    std::string name = std::string(argv[++ia]);
                    double value = std::stod(std::string(argv[++ia]));
                    if(network->updateVariable(name,value)){
                        std::cout << name << " updated to " << value << std::endl;
                    }else{
                        std::cout << "failed to update " << name << " to " << value << std::endl;
                    }
                }else if(std::string(argv[ia])=="-a" && argc > ia+1){
                    network->alphamin = std::stod(std::string(argv[++ia]));
                    network->alphamax = network->alphamin;
                    std::cout << "set alpha to " << network->alphamin << " for all steps" << std::endl;
                }else if(std::string(argv[ia])=="-o" && argc > ia+1){
                    ofilename = std::string(argv[++ia]);
                }else if(std::string(argv[ia])=="-U" && argc > ia+4){
                    varyVars=true;
                    varnames.push_back(std::string(argv[++ia]));
                    tstep = std::stod(std::string(argv[++ia]));
                    valuefirst = std::stod(std::string(argv[++ia]));
                    valuesecond = std::stod(std::string(argv[++ia]));
                    if(argc>ia+1&&argv[ia+1][0]!='-'){
                        tstep2 = std::stod(std::string(argv[++ia]));
                    }
                }else if(std::string(argv[ia])=="-V" && argc > ia+1){
                    varnames.push_back(std::string(argv[++ia]));
                }
                ++ia;
            }
            
            std::cout << "Starting simulation" << std::endl;
            if(!varyVars){
                // regular simulation run
                sim.run();
            }else{
                // if -U is specified
                double duration = 1000.0+tstep*2.0;
                network->randomGen.calculateUntil(duration);
                double dt=0.1;
                std::cout << " updating ";
                for(auto it =varnames.begin();it!=varnames.end();++it){
                    std::string varname = *it;
                    if(varname=="alpha"){
                        network->setAlpha(valuefirst);
                    }else{
                        network->updateVariable(varname,valuefirst);
                    }
                    std::cout << varname;
                }
                std::cout << std::endl;
                for(double t=0.0;t<duration;t+=dt){
                    sim.run_dt(dt);
                    for(auto it =varnames.begin();it!=varnames.end();++it){
                        std::string varname = *it;
                        if (t>=tstep&&t<tstep2){
                            if(varname=="alpha"){
                                network->setAlpha(valuesecond);
                            }else{
                                network->updateVariable(varname,valuesecond);
                            }
                        }else if(t>=tstep2){
                            if(varname=="alpha"){
                                network->setAlpha(valuefirst);
                            }else{
                                network->updateVariable(varname,valuefirst);
                            }
                        }
                    }
                }
            }
            time_t t_savetxtend= time(0);
            sim.save_list_to_cdf(ofilename);
            tend = time(0);
            std::cout << "Save cdf "<< difftime(tend, t_savetxtend) <<" second(s)."<< std::endl;
            std::ofstream file;
            file.open("./runs/run" + std::to_string(tend) + ".txt");
            file << *network << std::endl;
            file.close();
            std::cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< std::endl;
        }else{
            std::cout << "File "<< argv[2] <<" not found."<< std::endl;
        }
    }else{
        std::cout << "specify filename: icpg -f filename [-u var value] [-a alpha] [-U varname tstep value1 value2 tstep2] [-V varname] [-o outfile_ext]"<< std::endl;
    }
}

