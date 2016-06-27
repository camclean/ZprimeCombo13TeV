#ifndef REDIRECT_STDIO_HPP
#define REDIRECT_STDIO_HPP

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>

namespace theta{

/** some f*** up libararies (most notably ROOT), sometimes dump non-suppressable info to stdout/stderr
 * to write some messages noone is interested in, such as error messages from the minimizer
 * or error messages while opening files, in short, errors which are checked for anyway and just add noise
 * and confusion.
 *
 * To shut them up, the a call to redirect_stdio will:
 * <ul>
 *   <li>Duplicate the current stdout / stderr to file descriptors used by theta::cout and theta::cerr. The only possibility to write 
 *       to stdout / stderr always is through theta::cerr and theta::cout.</li>
 *   <li>Re-open the original stderr and stdout file descriptors to write to "/dev/null".</li>
 * </ul>
 */
void redirect_stdio();

extern boost::iostreams::stream<boost::iostreams::file_descriptor_sink> out;
extern boost::iostreams::stream<boost::iostreams::file_descriptor_sink> err;

// unfortunately, iostreams does not expose file descriptors (at least not yet in boost 1.42). For
// some magic in the default progress listener, we need the file descriptors, so expose them here:
extern int out_fd;
extern int err_fd;


}

#endif

