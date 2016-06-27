#include "interface/redirect_stdio.hpp"

#include <sys/stat.h>
#include <fcntl.h>


// as default: these are synonyms for stdout, stderr:
int theta::out_fd = 1;
int theta::err_fd = 2;
boost::iostreams::stream<boost::iostreams::file_descriptor_sink> theta::out(boost::iostreams::file_descriptor_sink(1, boost::iostreams::never_close_handle));
boost::iostreams::stream<boost::iostreams::file_descriptor_sink> theta::err(boost::iostreams::file_descriptor_sink(2, boost::iostreams::never_close_handle));

void theta::redirect_stdio(){
    static bool already_replaced = false;
    if(!already_replaced){
        theta::out_fd = dup(1);
        theta::out.close();
        theta::out.open(boost::iostreams::file_descriptor_sink(theta::out_fd, boost::iostreams::never_close_handle));
        close(1);
        open("/dev/null", O_WRONLY);
        theta::err_fd = dup(2);
        theta::err.close();
        theta::err.open(boost::iostreams::file_descriptor_sink(theta::err_fd, boost::iostreams::never_close_handle));
        close(2);
        open("/dev/null", O_WRONLY);
        already_replaced = true;
    }
}


