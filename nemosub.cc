/*
 *  nemosub.cc
 *  Nemo2
 *
 *  Created by fred on 30.05.07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <time.h>
#include "src/fileparser.h"
#include "src/basicsimulation.h"
#include "src/tstring.h"
#include "src/Uniform.h"
 #include <gsl/gsl_rng.h>

class SUB_DIRECTIVES : public SimComponent {
  
public:
  SUB_DIRECTIVES () {
    set_paramset("sub_cmd", false, NULL);
    add_parameter("sub_program",STR,1,0,0,0); //sbatch, qsub, bsub, etc
    add_parameter("sub_script_directive",STR,1,0,0,0); //this is the #PBS, #OAR, etc in the shell script (without leading #)
    add_parameter("sub_script_header",STR,0,0,0,0); //default is #!/bin/bash
    add_parameter("sub_script_args",STR,0,0,0,0); //the arguments to the shell script, may hold multiple args
    add_parameter("sub_append_ini_to_script",BOOL,0,0,0,0); //tells if the ini file name is written to the script instead of passed as an argument to the script
    add_parameter("sub_jobname",STR,0,0,0,0);
    add_parameter("sub_queue",STR,0,0,0,0);
    add_parameter("sub_parameters",STR,0,0,0,0); //may be multi-args, will be added one per line
    add_parameter("sub_resources",STR,0,0,0,0);  //a single line, without spaces, will be outputed wite pre-pending '-l'
    add_parameter("sub_input_dir",STR,0,0,0,0); //directory where ini files are stored
    add_parameter("sub_output_dir",STR,0,0,0,0); //directory where job output files are saved
  }
  virtual ~SUB_DIRECTIVES() {}
  virtual bool setParameters (){return true;}
  virtual void loadFileServices ( FileServices* loader ){}
  virtual void loadStatServices ( StatServices* loader ){}
};

gsl_rng * RAND::r = 0;

void remove_sub_params (map< string,string >& rec);

using namespace std;

int main (int argc, char **argv)
{
  RAND::init(time(NULL));
  
// ---- sub script writer: ----  
  ParamManager theSim;
  SUB_DIRECTIVES SUB_Params;
  ParamSet* SUB_Pset = SUB_Params.get_paramset();
  theSim.add_paramset(SUB_Pset); //This is the only paramset added to the param manager
  
  FileParser Reader("");
  
  //read the template file which contains the parameters
  ifstream TEMPLATE("sub_template");
  //create an output stream for the submitting shell script
  ostringstream TMP;
  //the list of paramters to pass to the scheduler in the shell script :
  string SUB_cmd;
  
  if( !(TEMPLATE) ) fatal("could not open sub template!!\n");
  
  //paste the content of the template file into our output stream
  TMP << TEMPLATE.rdbuf();
  
  SUB_cmd = TMP.str();
  
  //remove trailing newline, replace with one space char
  if(SUB_cmd[ SUB_cmd.length() -1 ] == '\n') SUB_cmd[ SUB_cmd.length() -1 ] = ' ';
  
  TEMPLATE.close();
  
  //read the iniput filename, and build the simulation parameters list
  //(no external file expansions yet)
  if (argc == 1)
    theSim.build_records(Reader.getParsedParameters("Nemo2.ini"));
  else
    for (int i = 1; i < argc; ++i)
      theSim.build_records(Reader.getParsedParameters(argv[i]));
  
  list< map< string, string > > sims = theSim.get_simRecords();
  list< map< string, string > >::iterator currentSim = sims.begin();
  
  map< string,string > pMap;
  map< string,string >::iterator param;
  
  string filename, jobname, output_dir, input_dir, program, script_args;
  ofstream FH;
  string codename;

  //a codename is added to the simulation filename
  //codename begins with 2 letters, followed by 3 digits:
  unsigned char cap_nocap = (RAND::RandBool() ? 65 : 97);
  unsigned char char_code[5] = {'\0','\0','\0','\0','\0'};
  char_code[0] = cap_nocap + RAND::Uniform(26);
  char_code[1] = (RAND::RandBool() ? 65 : 97) + RAND::Uniform(26);
  char_code[2] = (RAND::RandBool() ? 65 : 97) + RAND::Uniform(26);
  char_code[3] = (RAND::RandBool() ? 65 : 97) + RAND::Uniform(26);
  
  unsigned int num = RAND::Uniform(800);
  
  //vars used to get multiple args of sub parameters
  Param* multiparam;
  vector<string> options;
  
  while(currentSim != sims.end()) {
     
    //the list of current parameters
    pMap = *currentSim;
    
    //set the parameters value for SUB_DIRECTIVES components ONLY!
    theSim.set_parameters(pMap, true);
    
    //setting input dir where init files will be stored
    if(SUB_Pset->isSet("sub_input_dir")) {

      input_dir = SUB_Pset->getArg("sub_input_dir");

      if(input_dir.size() != 0 && input_dir[input_dir.length() -1] != '/')
        input_dir += "/";

    } else
      input_dir = "jobs/init/";
    
    //setting job file output dir
    if(SUB_Pset->isSet("sub_output_dir")) {

      output_dir = SUB_Pset->getArg("sub_output_dir");

      if(output_dir.size() != 0 && output_dir[output_dir.length() -1] != '/')
        output_dir += "/";

    } else
      output_dir = "jobs/";
    
    
    script_args = "";
    
    if(SUB_Pset->isSet("sub_script_args")) {
      
      multiparam = SUB_Pset->get_param("sub_script_args");
      
      if(multiparam->hasMultipleArgs()) {
        
        options = multiparam->getMultiArgs();
        
        for(unsigned int i = 0; i < options.size(); ++i)
          script_args += options[i] + " ";
        
      } else {
        
        script_args = multiparam->getArg();
        
      }
    }
    
    program = SUB_Pset->getArg("sub_program");
    
    //grab the simulation filename
    filename = pMap["filename"];
    
    if(filename.size() > 15) filename = filename.substr(0, 15);
   
    codename = (const char*)char_code + tstring::dble2str(num);
    
    //the digit part gets incremented:
    num++;
    //recode the filename, will be used as the name of the .ini file of this simulation
    filename += codename;
    
    //set the jobname
    if(SUB_Pset->isSet("sub_jobname"))
      jobname =  pMap["sub_jobname"];
    else 
      jobname = filename;
    
    filename += ".ini";

    filename = input_dir + filename;

    //remove the SUB_Pset parameters from the simulation, they are not handled by Nemo
    remove_sub_params(pMap);
 
//----------------------------------------
//--- writing the parameters to the init file:
    FH.open(filename.c_str(),ios::out | ios::trunc);
    
    if(!FH){
      return error("nemosub: %s: open failed!!!\n",filename.c_str());
//      continue;
    }
    FH<<"#sub_codename = "<<codename<<endl;
    FH<<"#sub_jobname = "<<jobname<<endl;
    
    for(param = pMap.begin(); param != pMap.end(); param++) {
      FH<<param->first<<" "<<param->second<<endl;
    }
    FH.close();
    
//---------------------------------------
//--- writing the job submitting script:
    FH.open("sub000",ios::out | ios::trunc);
    if(!FH)
      fatal("sub000: open failed!!!\n");

    string TAG = SUB_Pset->getArg("sub_script_directive");
    
    string header;
    
    if(SUB_Pset->isSet("sub_script_header"))
      header = SUB_Pset->getArg("sub_script_header");
    else
      header = "#!/bin/bash";
    
    FH<<header<<endl;
    
    //some customizations...
    if(TAG == "PBS") {
      FH<<"#PBS -S /bin/bash\n";//this is the default
    
      FH<<"#PBS -N "<<jobname<<endl;
      
    } else if(TAG == "BSUB") {

      FH<<"#BSUB -J "<<jobname<<endl;
      
      //splitting the stdout and stderr:
      FH<<"#BSUB -eo "<<output_dir<<jobname<<".e%J"<<endl;
      FH<<"#BSUB -oo "<<output_dir<<jobname<<".o%J"<<endl;
      
      
    } else if(TAG == "SBATCH") {
    
      FH<<"#SBATCH -J "<<jobname<<endl;
      FH<<"#SBATCH -e "<<output_dir<<jobname<<".e%j"<<endl;
      FH<<"#SBATCH -o "<<output_dir<<jobname<<".o%j"<<endl;
      
    } else if(TAG == "OAR") {
      
      FH<<"#OAR -n "<<jobname<<endl;
      FH<<"#OAR -E "<<output_dir<<jobname<<".e%jobid%"<<endl;
      FH<<"#OAR -O "<<output_dir<<jobname<<".o%jobid%"<<endl;
      
    } else {
      
      warning("not setting jobname of new job as we can't guess the right option.\n");
    
    }
    
    
    if(SUB_Pset->isSet("sub_queue"))
       FH<<"#"<<TAG<<" -q "<<SUB_Pset->getArg("sub_queue")<<endl;
    
    if(SUB_Pset->isSet("sub_resources"))
       FH<<"#"<<TAG<<" -l "<<SUB_Pset->getArg("sub_resources")<<endl;

    if(SUB_Pset->isSet("sub_parameters")) {
      
      multiparam = SUB_Pset->get_param("sub_parameters");
      
      if(multiparam->hasMultipleArgs()) {
        
        options = multiparam->getMultiArgs();
        
        for(unsigned int i = 0; i < options.size(); ++i)
          FH<<"#"<<TAG<<" "<<options[i]<<endl;
      
      } else {
      
        FH<<"#"<<TAG<<" "<<multiparam->getArg()<<endl;
        
      }
    
    }
    
    FH<<SUB_cmd; //walltime and mem and Nemo cmd from the sub_template file

    if(SUB_Pset->isSet("sub_append_ini_to_script")) {
      
      FH<<" "<<filename<<endl;
      FH<<"rm -f "<<filename<<endl;
      
    }
  
    FH<<endl;
    
    FH.close();

    //launch scheduler
    
    string exec_arg, exec_cmd; //the shell command
    
    if( !SUB_Pset->isSet("sub_append_ini_to_script") )
      exec_arg = "sub000 " + filename + " " + script_args;  //this implies that the first argument is always the ini file name
    else
      exec_arg = "sub000 " + script_args;

    //make the script executable, this is often required by the scheduler
    if(system("chmod 755 sub000") != 0)
      fatal("failed to make script \"sub000\" executable\n");

    if(program == "bsub")
      exec_cmd = program + " < " + exec_arg + " >> " + output_dir + "launched_jobs";
    else if(program == "oarsub")
      exec_cmd = program + " ./" + exec_arg + " >> " + output_dir + "launched_jobs";
    else
      exec_cmd = program + " " + exec_arg + " >> " + output_dir + "launched_jobs";
    
#ifdef _DEBUG_

    cout<<pMap["filename"]<<" "<<jobname<<" "<<filename<<endl;
   
    cout<<exec_cmd<<endl;
    
#else

    string cmd =  "printf '" + pMap["filename"] + " " + jobname + " " + filename + " ' >> " + output_dir + "launched_jobs";

    if(system(cmd.c_str()) != 0)
      warning("failed to write in %slaunch_jobs\n",output_dir.c_str());

    if(system(exec_cmd.c_str()) != 0)
      fatal("could not start job with %s\n",program.c_str());

#endif
    
    currentSim++;
  }

return 0;
}

void remove_sub_params (map< string,string >& rec)
{
  rec.erase("sub_program");
  rec.erase("sub_script_directive");
  rec.erase("sub_script_header");
  rec.erase("sub_script_args");
  rec.erase("sub_append_ini_to_script");
  rec.erase("sub_jobname");
  rec.erase("sub_queue");
  rec.erase("sub_parameters");
  rec.erase("sub_resources");
  rec.erase("sub_input_dir");
  rec.erase("sub_output_dir");
}
