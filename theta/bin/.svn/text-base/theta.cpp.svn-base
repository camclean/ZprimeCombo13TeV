#include "interface/cfg-utils.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"
#include "interface/variables-utils.hpp"
#include "interface/variables.hpp"
#include "interface/main.hpp"
#include "interface/model.hpp"
#include "interface/redirect_stdio.hpp"

#include "interface/database.hpp"

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

#ifdef USE_TIMER
#include <boost/timer/timer.hpp>
#endif

#include <termios.h>

#include <fstream>
#include <iostream>

using namespace std;
using namespace theta;
using namespace theta::utils;

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace btime = boost::posix_time;

class MyProgressListener: public ProgressListener{
public:
    virtual void progress(int done_, int total_, int n_error_){
        done = done_;
        total = total_;
        n_error = n_error_;
        print();
    }
    
    void print(){
        if(!is_tty) return;
        btime::ptime now = btime::microsec_clock::local_time();
        if(now < next_update && done < total) return;
        //move back to beginning of terminal line:
        theta::out << "\033[" << chars_written << "D";
        chars_written = 0;  
        char c[200];
        double error_fraction = 100.0 * n_error / done;
        const char * color_start;
        const char * color_end = "\033[0m";
        if(error_fraction > 30){ //red:
           color_start = "\033[31m";
        }
        else if(error_fraction > 5){ //yellow:
           color_start = "\033[1;33m";
        }
        else{ //green:
           color_start = "\033[32m";
        }
        if(total > 0){
            double progress_fraction = 100.0 * done / total;
            chars_written += snprintf(c, 200, "progress: %6d / %-6d [%5.1f%%]   errors: %s%6d [%5.1f%%]%s", done, total, progress_fraction, color_start, n_error, error_fraction, color_end);
        }
        else{
            chars_written += snprintf(c, 200, "progress: %6d   errors: %s%6d [%5.1f%%]%s", done, color_start, n_error, error_fraction, color_end);
        }
        theta::out << c;
        theta::out.flush();
        next_update = now + btime::milliseconds(50);
    }

    MyProgressListener(): stdout_fd(theta::out_fd), done(0), total(0), n_error(0), is_tty(isatty(stdout_fd)),
      chars_written(0), next_update(btime::microsec_clock::local_time()) {
        if(!is_tty) return;
        //disable terminal echoing; we don't expect any input.
        termios settings;
        if (tcgetattr (stdout_fd, &settings) < 0) return; //ignore error
        settings.c_lflag &= ~ECHO;
        tcsetattr (stdout_fd, TCSANOW, &settings);
    }
    
    ~MyProgressListener(){
        if(!is_tty) return;
        //enable terminal echoing again; don't be evil
        termios settings;
        if (tcgetattr (stdout_fd, &settings) < 0) return; //ignore error
        settings.c_lflag |= ECHO;
        tcsetattr (stdout_fd, TCSANOW, &settings);
        if(chars_written > 0)
            theta::out << endl;
    }
    
private:
    int stdout_fd;
    int done, total, n_error;
    bool is_tty;
    int chars_written;
    btime::ptime next_update;
};



namespace{
    //note: these functions are useful to have compatibility
    // with both V2 and V3 of boost::filesystem
    // as in V2, path::filename returns a string whereas in
    // V3, path::filename returns a path.
    std::string to_string(const std::string & s){
        return s;
    }
    std::string to_string(const fs::path & p){
        return p.string();
    }
}

// read the whole file into a string
std::string read_file(const string & fname){
    try{
        size_t s = fs::file_size(fname);
        string result;
        result.resize(s);
        ifstream in(fname.c_str());
        if(!in){
            throw string(); // will be catched below and converted to exception with proper message
        }
        in.read(&result[0], s);
        if(!in){
            throw string();
        }
        return result;
    }
    catch(...){
        throw ConfigurationException("error reading file '" + fname + "'");
    }
}

boost::shared_ptr<Main> build_main(string cfg_filename, bool nowarn, bool print_time){
    boost::optional<Setting> root;
    boost::shared_ptr<SettingUsageRecorder> rec;
    if(not nowarn) rec.reset(new SettingUsageRecorder());
    boost::shared_ptr<Main> main;
    #ifdef USE_TIMER
    std::auto_ptr<boost::timer::cpu_timer> timer;
    #endif
    try {
        #ifdef USE_TIMER
        if(print_time){
            timer.reset(new boost::timer::cpu_timer());
        }
        #endif
        try {
            //as includes in config files should always be resolved relative to the config file's location:
            string old_path = fs::current_path().string();
            //convert any failure to a FileIOException:
            try{
                 if(fs::path(cfg_filename).has_parent_path()){
                    fs::current_path(fs::path(cfg_filename).parent_path());
                 }
                 cfg_filename = to_string(fs::path(cfg_filename).filename());
            }
            catch(fs::filesystem_error & ex){
                 throw ConfigurationException("error during resolving config file path for filename '" + cfg_filename + "'");
            }
            std::string cfg_string = read_file(cfg_filename);
            root = LibconfigSetting::parse(cfg_string, rec);
            fs::current_path(old_path);
        }catch (Exception & ex) {
            stringstream s;
            s << "Error while parsing configuration file '" << cfg_filename << "': " << ex.what();
            throw ConfigurationException(s.str());
        }
        #ifdef USE_TIMER
        if(print_time){
            theta::out << timer->format(4, "Time to read and parse configuration file:        %w sec real, %t sec CPU") << endl;
        }
        #endif
        theta_assert(root);
        Configuration config(*root);
        config.pm->set("default", boost::shared_ptr<VarIdManager>(new VarIdManager));
        
        //process options:
        Configuration cfg_options(config, config.setting["options"]);
        #ifdef USE_TIMER
        if(print_time){
            timer.reset(new boost::timer::cpu_timer());
        }
        #endif
        PluginLoader::execute(cfg_options);
        #ifdef USE_TIMER
        if(print_time){
            theta::out << timer->format(4, "Time to load plugin .so files:                    %w sec real, %t sec CPU") << endl;
            timer.reset(new boost::timer::cpu_timer());
        }
        #endif
        
        //populate VarIdManager from config:
        apply_vm_settings(config);
        //build run:
        main = PluginManager<Main>::build(Configuration(config, (*root)["main"]));
        #ifdef USE_TIMER
        if(print_time){
            theta::out << timer->format(4, "Time to construct object tree from configuration: %w sec real, %t sec CPU") << endl;
        }
        #endif
    }
    catch (Exception & e) {
        theta::err << "Error: " << e.what() << endl;
    }
    if(not main.get()){
        return main;
    }    
    if(rec){
        vector<string> unused;
        rec->get_unused(unused, *root);
        if (unused.size() > 0) {
            theta::out << "WARNING: following setting paths in the configuration file have not been used: " << endl;
            for (size_t i = 0; i < unused.size(); ++i) {
                theta::out << "  " << (i+1) << ". " << unused[i] << endl;
            }
            theta::out << "Comment out these settings to get rid of this message." << endl;
        }
    }
    return main;
}


int main(int argc, char** argv) {
    redirect_stdio();
    fill_theta_dir(argv);
    if(theta_dir==""){
        theta::err << "WARNING: could not determine theta_dir, leaving empty" << endl;
    }
    
    po::options_description desc("Options");
    desc.add_options()("help,h", "show help message")
    ("quiet,q", "quiet mode (suppress progress message)")
    ("print-time", "print some time information to stdout")
    ("nowarn", "do not warn about unused configuration file statements");

    po::options_description hidden("Hidden options");

    hidden.add_options()
    ("cfg-file", po::value<vector<string> >(), "configuration file");

    po::positional_options_description p;
    p.add("cfg-file", -1);

    po::options_description cmdline_options;
    cmdline_options.add(desc).add(hidden);

    po::variables_map cmdline_vars;
    try{
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), cmdline_vars);
    }
    catch(std::exception & ex){
        theta::err << "Error parsing command line options: " << ex.what() << endl;
        return 1;
    }
    po::notify(cmdline_vars);
    
    if(cmdline_vars.count("help")){
        theta::out << desc << endl;
        return 0;
    }
    
    if(cmdline_vars.count("cfg-file")==0){
        theta::err << "Error: you have to specify a configuration file" << endl;
        return 1;
    }
    
    vector<string> cfg_filenames = cmdline_vars["cfg-file"].as<vector<string> >();
    bool quiet = cmdline_vars.count("quiet");
    bool nowarn = cmdline_vars.count("nowarn");
    bool print_time = cmdline_vars.count("print-time");

    try {
        for(size_t i=0; i<cfg_filenames.size(); ++i){
            if(!quiet and cfg_filenames.size() > 1){
                theta::out << "processing file " << (i+1) << " of " << cfg_filenames.size() << ", " << cfg_filenames[i] << endl;
            }
            boost::shared_ptr<Main> main = build_main(cfg_filenames[i], nowarn, print_time);
            if(!main) return 1;
            if(not quiet){
                boost::shared_ptr<ProgressListener> l(new MyProgressListener());
                main->set_progress_listener(l);
            }

            //install signal handler now, not much earlier. Otherwise, plugin loading in build_run()
            // might change it ...
            install_sigint_handler();
            #ifdef USE_TIMER
            std::auto_ptr<boost::timer::cpu_timer> timer;
            if(print_time){
                timer.reset(new boost::timer::cpu_timer());
            }
            #endif
            main->run();
            #ifdef USE_TIMER
            if(print_time){
                std::string timestr = timer->format(4, "Time to run main:                                 %w sec real, %t sec CPU");
                main.reset(); // to destroy progress listener, which outputs newline
                theta::out << timestr  << endl;
            }
            #endif
            if(stop_execution) break;
        }
    }
    catch(ExitException & ex){
       theta::err << "Exit requested: " << ex.message << endl;
       return 1;
    }
    catch (Exception & ex) {
        theta::err << "An error occurred in Main::run: " << ex.what() << endl;
        return 1;
    }
    catch(logic_error & ex){
        theta::err << "A logic error occurred in Main::run: " << ex.what() << endl;
        return 1;
    }
    catch(exception & ex){
        theta::err << "An unspecified exception occurred in Main::run: " << ex.what() << endl;
        return 1;
    }
    if(theta::stop_execution){
        theta::out << "(exiting on SIGINT)" << endl;
    }
    if(!quiet){
        unsigned long n = atomic_get(&n_nll_eval);
        unsigned long nwd = atomic_get(&n_nll_eval_with_derivative);
        if(n != 0){
            theta::out << "Total number of likelihood evaluations: " << n;
            if(nwd != 0){
                theta::out << "; with derivative: " << nwd;
            }
            theta::out << endl;
        }
    }
    return 0;
}
