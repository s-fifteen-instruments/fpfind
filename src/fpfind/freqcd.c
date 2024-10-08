/* freqcd.c : Part of the quantum key distribution software for correcting
              timestamps emitted by a timestamp unit running at a relative
              frequency offset. Proof-of-concept.

   Copyright (C) 2024 Justin Peh, Xu Zifang, Christian Kurtsiefer,
                      National University of Singapore

   This source code is free software; you can redistribute it and/or
   modify it under the terms of the GNU Public License as published
   by the Free Software Foundation; either version 2 of the License,
   or (at your option) any later version.
   This source code is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   Please refer to the GNU Public License for more details.
   You should have received a copy of the GNU Public License along with
   this source code; if not, write to:
   Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

   --

   Program that receives a frequency offset between two parties calculated
   by fpfind, and performs a software correction of the timestamps emitted
   by readevents. Parameters are optimized such that the correction can
   be performed using purely 64-bit constructs.

   Note that there is no support for Windows - many of the functionality
   used are native to Linux.


   usage:
     freqcd [-i infilename] [-o outfilename] [-X]
            [-t timecorr] [-f freqcorr] [-F freqfilename] [-d] [-u]


   DATA STREAM OPTIONS:
     -i infilename:   Filename of source events. Can be a file or a socket
                      and has to supply binary data according to the raw data
                      spec from the timestamp unit. If unspecified, data
                      is read from stdin.
     -o outfilename:  Outfile name for timing corrected events, which can
                      either be a file or socket. Output format is the same as
                      with input. If unspecified, data is written to stdout.
     -F freqfilename: Filename of frequency correction values. Needs to be
                      a readable+writeable socket storing newline-delimited
                      frequency offset values (see '-f' option for format).
                      Takes priority over '-f' supplied offsets.
                      If unspecified, the frequency offset will be static.

   ENCODING OPTIONS:
     -X:              Specifies if both the raw input and output data streams
                      are to be read in legacy format, as specified in the
                      timestamp unit.
     -t timecorr:     Timing offset of the timestamps, in units of ps. If
                      unspecified, offset = 0. The ps unit is chosen for
                      legibility; the actual resolution is 1/256ns (~4ps),
                      in line with the high resolution timestamp spec.
     -f freqcorr:     Frequency offset of the current clock relative to some
                      reference clock, in units of 2^-34 (or 0.6e-10).
                      If unspecified, offset = 0. Maximum absolute value is
                      2097151 (i.e. 2^21-1), see [1] for explanation.
     -d:              Decimal mode. Changes the frequency offset read by '-f'
                      and '-F' to be in units of 0.1 ppb instead, with
                      maximum absolute value of 1220703 (i.e. 2^-13). For
                      human-readability.
     -u:              Update mode. Changes the frequencies read by '-F' to
                      represent frequency differentials instead, i.e. the
                      current frequency offset is shifted by df' rather than
                      replaced with df, given by:

                          df = (1+df)*(1+df')-1
                             ≈ df+df' (to nearest 0.1ppb, if df*df' < 0.05ppb)

   Potential improvements:
     - Consider using epoll() if necessary
     - Merge write procedure with select() call
     - Optimize buffer parameters

   References:
     [1]: Formulae, https://github.com/s-fifteen-instruments/fpfind/issues/6
 */

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>   // open, O_RDONLY, O_NONBLOCK
#include <sys/select.h> // fd_set, usually implicitly declared
#include <unistd.h>  // getopt, select
#include <errno.h>   // select errno
#include <limits.h>  // INT_MAX
#include <math.h>    // pow(..., 10), round(...)
#include <string.h>  // strcmp(...)

/* default definitions */
#define FNAMELENGTH 200            /* length of file name buffers */
#define FNAMEFORMAT "%200s"        /* for sscanf of filenames */
#define FILE_PERMISSONS 0644       /* for all output files */
#define INBUFENTRIES 1024000       /* max. elements in input buffer */
#define OUTBUFENTRIES 1024000      /* max. elements in output buffer */
#define FREQBUFSIZE 100            /* max. length of frequency correction values */
#define RETRYREADWAIT 500000       /* sleep time in usec after an empty read */
#define FCORR_ARESBITS -34         /* absolute resolution of correction, in power of 2 */
#define FCORR_AMAXBITS -13         /* absolute maximum allowed correction, in power of 2 */
#define FCORR_DEFAULT 0            /* frequency correction in units of 2^FCORR_ARESBITS */
#define FCORR_ARESDECIMAL -10      /* absolute resolution of correction, in power of 10, specifically for '-d' mode */
#define FCORR_OFLOWRESBITS 64      /* overflow correction resolution, in power of 2 */
//#define __DEBUG__                /* global debug flag, uncomment to disable, send -v flag for more debug msgs */

// Formatting colors
#define COLOR_ERROR "\x1b[31m"
#define COLOR_RESET "\x1b[0m"

// struct defined in non-legacy format
typedef struct rawevent {
    unsigned int low;
    unsigned int high;
} rawevent;

typedef long long ll;
typedef unsigned long long ull;
#ifdef __SIZEOF_INT128__
typedef unsigned __int128 u128;
typedef __int128 i128;
#else
#warning "128-bit integers unsupported in current compiler: \
slight undercompensation of frequency will occur from rounding errors."
#endif

/* error handling */
char *errormessage[] = {
    "No error",
    "Generic error", /* 1 */
    "Error reading in infilename",
    "Error opening input stream source",
    "Error parsing freq correction as integer",
    "Freq correction value out of range", /* 5 */
    "Unable to allocate memory to inbuffer",
    "Unable to allocate memory to outbuffer",
    "Error reading in outfilename",
    "Error opening output stream",
    "Error reading in freqfilename", /* 10 */
    "Error opening freq correction stream source",
    "Freq correction value not newline terminated",
    "Unable to allocate memory to freqbuffer",
    "Error parsing time correction as integer",
    "Duplicate filename specified for reading and writing", /* 15 */
};
void wmsg(char *message) {
    fprintf(stderr, COLOR_ERROR "%s\n" COLOR_RESET, message);
};
void wmsg_s(size_t n, char *fmt, char *v) {  // dirty fix for overloading signatures
    char *buffer = (char *)malloc(n*sizeof(char));
    snprintf(buffer, n, fmt, v); wmsg(buffer); free(buffer);
};
void wmsg_i(size_t n, char *fmt, int v) {
    char *buffer = (char *)malloc(n*sizeof(char));
    snprintf(buffer, n, fmt, v); wmsg(buffer); free(buffer);
};
int emsg(int code) {
    wmsg(errormessage[code]);
    return code;
};

/* for reading freqcorr values */
int readint(char *buff) {
    char *end;
    long value;
    value = strtol(buff, &end, 10);
    // string is a decimal, zero-terminated, and fits in an int
    if ((end != buff) && (*end == 0) && (abs(value) < INT_MAX))
        return (int)value;
    return INT_MAX;
}

void usage() {
    fprintf(stderr, "\
Usage: freqcd [-i infilename] [-o outfilename] [-X]\n\
              [-t timecorr] [-f freqcorr] [-F freqfilename] [-d] [-u]\n\
Performs frequency correction of timestamps emitted by qcrypto's readevents.\n\
\n\
Supports up to a timing resolution of 1/256ns, but should work with 1/8ns and\n\
2ns timing resolutions as well. An optional static timing correction can be\n\
applied. If '-u' is not supplied, reads from '-F' replace the current freq.\n\
\n\
Data stream options:\n\
  -i infilename    File/socket name for source events. Defaults to stdin.\n\
  -o outfilename   File/socket name for corrected events. Defaults to stdout.\n\
  -F freqfilename  File/socket name of frequency correction values.\n\
\n\
Encoding options:\n\
  -X               Use legacy timestamp format.\n\
  -t timecorr      Timing offset, in units of ps (resolution of 1/256 ns).\n\
  -f freqcorr      Frequency offset, in units of 2^-34 (range: 0-2097151).\n\
  -d               Use units of 0.1ppb for '-f'/'-F'   (range: 0-1220703).\n\
  -u               Modify current frequency relative to reads from '-F'.\n\
\n\
Shows this help with '-h' option. More descriptive documentation:\n\
<https://github.com/s-fifteen-instruments/fpfind/blob/main/src/fpfind/freqcd.c>\n\
");
}


int main(int argc, char *argv[]) {

    /* other constants */
    const int INBUFSIZE = INBUFENTRIES * sizeof(struct rawevent);
    const int OUTBUFSIZE = OUTBUFENTRIES * sizeof(struct rawevent);
    const int FCORR_MAX = 1 << (FCORR_AMAXBITS - FCORR_ARESBITS);
    const int FCORR_TBITS1 = -FCORR_AMAXBITS - 1;  // bit truncations when correcting timestamp
    const int FCORR_TBITS2 = (FCORR_AMAXBITS - FCORR_ARESBITS) + 1;
    const double FCORR_BTO1 = pow(2, FCORR_ARESBITS);

    // Conversion factor for decimal mode, i.e. 1e-10/2^-34 = ~1.7179869
    // Ratio inverted due to negative signs in exponents.
    const double FCORR_DTOB = ((long long)1 << -FCORR_ARESBITS) / pow(10, -FCORR_ARESDECIMAL);

    /* parse options */
    int fcorr = FCORR_DEFAULT;  // frequency correction value
    ll tcorr = 0;               // time correction value
    char infilename[FNAMELENGTH] = {};  // store filename
    char outfilename[FNAMELENGTH] = {};  // store filename
    char freqfilename[FNAMELENGTH] = {};  // store filename
    int islegacy = 0;  // mark if format is legacy
    int isdecimal = 0;  // mark if frequency input is in decimal units
    int isupdate = 0;  // mark if relative frequency is used
    int isdebugverbose = 0;  // mark if need to show more verbose debug
    int opt;  // for getopt options
    opterr = 0;  // be quiet when no options supplied
    while ((opt = getopt(argc, argv, "i:o:F:f:t:xXduhv")) != EOF) {
        switch (opt) {
        case 'i':
            if (sscanf(optarg, FNAMEFORMAT, infilename) != 1) return -emsg(2);
            infilename[FNAMELENGTH-1] = 0;  // security termination
            break;
        case 'o':
            if (sscanf(optarg, FNAMEFORMAT, outfilename) != 1) return -emsg(8);
            outfilename[FNAMELENGTH-1] = 0;  // security termination
            break;
        case 'F':
            if (sscanf(optarg, FNAMEFORMAT, freqfilename) != 1) return -emsg(10);
            freqfilename[FNAMELENGTH-1] = 0;  // security termination
            break;
        case 'f':
            if (sscanf(optarg, "%d", &fcorr) != 1) return -emsg(4);
            break;
        case 't':
            if (sscanf(optarg, "%lld", &tcorr) != 1) return -emsg(14);
            tcorr = (ll) tcorr * 0.256;  // convert ps -> 1/256ns units
            break;
        case 'x':  // retained for legacy purposes, use '-X' instead
            wmsg("'-x' has been deprecated and will be removed not before fpfind:v1.2024.15 - use '-X' instead");
        case 'X':
            islegacy = 1;
            break;
        case 'd':
            isdecimal = 1;
            break;
        case 'u':
            isupdate = 1;
            break;
        case 'h':
            usage(); exit(1);
            break;
        case 'v':  // leave undocumented, for internal use only
            isdebugverbose = 1;
            break;
        }
    }

    /* check specified frequency correction is within-bounds */
    // needs to be done after argument reading for position-agnostic '-f'
    if (isdecimal) fcorr = (int) fcorr * FCORR_DTOB;
    if (abs(fcorr) >= FCORR_MAX) return -emsg(5);

    /* set input and output handler */
    int inhandle = 0;  // stdin by default
    if (infilename[0]) {
        inhandle = open(infilename, O_RDONLY | O_NONBLOCK);
        if (inhandle == -1) return -emsg(3);
    }

    int freqhandle = 0;  // null by default (not stdin)
    if (freqfilename[0]) {
        freqhandle = open(freqfilename, O_RDONLY | O_NONBLOCK);
        if (freqhandle == -1) return -emsg(11);
    }

    int outhandle = 1;  // stdout by default
    if (outfilename[0]) {
        if ((infilename[0]) && (strcmp(outfilename, infilename) == 0)) return -emsg(15);
        if ((freqfilename[0]) && (strcmp(outfilename, freqfilename) == 0)) return -emsg(15);
        outhandle = open(outfilename, O_WRONLY | O_CREAT | O_TRUNC, FILE_PERMISSONS);
        if (outhandle == -1) return -emsg(9);
    }

    /* initialize input and output buffers */
    struct rawevent *inbuffer;
    inbuffer = (struct rawevent *)malloc(INBUFSIZE);
    if (!inbuffer) return -emsg(6);
    struct rawevent *eventptr;  // pointer to position within inbuffer
    int eventnum = 0;  // number of available rawevents for processing
    char *inbufferbytes = (char *)inbuffer;  // lower level byte buffer
    char *inbufferbytes_next;  // pointer to next byte write destination

    struct rawevent *outbuffer;
    outbuffer = (struct rawevent *)malloc(OUTBUFSIZE);
    if (!outbuffer) return -emsg(7);
    int outevents = 0;

    char *freqbuffer;
    int freqbytesread = 0;
    int freqbytesread_next = 0;
    int freqbytespartial = 0;  // size of partial freqcorr value remaining
    char *freqbuffer_next;
    if (freqhandle) {  // allocate memory only if needed
        freqbuffer = (char *)malloc(FREQBUFSIZE);
        if (!freqbuffer) return -emsg(13);
        freqbuffer_next = freqbuffer;  // pointer to next char write destination
    }

    /* parameters for select call */
    fd_set rfds;
    struct timeval tv;
    int retval;

    /* inbuffer reading variables */
    int i, j;
    int inbytesread = 0;
    int inbytesread_next = 0;
    int inbytespartial;  // size of partial rawevent remaining in inbufferbyte

    /* timestamp variables */
    ull tsref = 0;  // reference timestamp to scale by frequency correction,
                    // noting subsequent initializations should zero 'tsref'
    int isset_tsref = 0;    // initialization marker for tsref
    ull ts, tsmeas, tsdiff; // timestamp
    ll tscorr;              // timestamp correction
    ll tsoverflowcorr = tcorr;  // timestamp overflow corrections
                                // overloading with the desired timing corr
#ifdef __SIZEOF_INT128__
    u128 _tsdiff;           // 128-bit equivalents to mitigate rounding error
    i128 _tscorr;
    i128 _tsoverflowcorr = (i128)tcorr << FCORR_OFLOWRESBITS;
#endif
    unsigned int high;      // high word in timestamp
    unsigned int low;       // low word in timestamp
    unsigned int _swp;      // temporary swap variable, support 'legacy' option

    /* main loop */
    while (1) {

        /* discard previously processed rawevents and
           retain partial rawevent left in buffer */
        // TODO: Fix potential bug when 'continue' is called, which
        //       clears the inputbuffer without writing to output stream.
        inbytespartial = inbytesread % sizeof(struct rawevent);
        for (i = inbytesread - inbytespartial, j = 0; j < inbytespartial; i++, j++) {
            inbufferbytes[j] = inbufferbytes[i];
        }
        inbufferbytes_next = &inbufferbytes[inbytespartial];

        /* wait for data on inhandle and freqhandle */
        // TODO: Consider whether to use poll/epoll mechanisms, if frequent
        //       pipe recreation is a concern (high fd).
        FD_ZERO(&rfds);
        FD_SET(inhandle, &rfds);
        if (freqhandle) FD_SET(freqhandle, &rfds);
        tv.tv_sec = 0;
        tv.tv_usec = RETRYREADWAIT;
        retval = select(FD_SETSIZE, &rfds, NULL, NULL, &tv);
        if (retval == -1) {
            wmsg_i(50, "Error %d on select.", errno);
            break;  // graceful close
        }

        if (FD_ISSET(inhandle, &rfds)) {

            /* read data from inhandle */
            // TODO: Highlight corresponding bug in chopper.c. Note that assigning
            //       to inbytesread directly can potentially corrupt events.
            inbytesread_next = read(inhandle, inbufferbytes_next, INBUFSIZE - inbytespartial);
            if (inbytesread_next == 0) {
                wmsg("Input stream closed.");
                break;  // no bytes read (i.e. EOF)
                        // TODO: Check if this should be continue instead,
                        //       when running ad-infinitum
            }
            if (inbytesread_next == -1) {
                wmsg_i(50, "Error %d on input read.", errno);
                break;  // graceful close
            }

            /* concatenate new data */
            inbytesread = inbytespartial + inbytesread_next;
            eventnum = inbytesread / sizeof(struct rawevent);
            eventptr = inbuffer;

            /* micro-optimization to initialize reference timestamp */
            if ((!isset_tsref) && (eventnum > 0)) {
                low = eventptr->low;
                high = eventptr->high;

                // Shift burden of swapping if 'legacy' format is used
                // TODO: Consider a more efficient implementation.
                if (islegacy) {
                    _swp = low;
                    low = high;
                    high = _swp;
                }
                tsref = ((ull)high << 22) | (low >> 10);
                isset_tsref = 1;  // we are done initializing
            }

            /* digest events */
            for (i = 0; i < eventnum; i++) {

                /* extract timestamp value */
                // Assumed 4ps timestamps used
                low = eventptr->low;
                high = eventptr->high;
                if (islegacy) {
                    _swp = low;
                    low = high;
                    high = _swp;
                }
                tsmeas = ((ull)high << 22) | (low >> 10);

                /* calculate timestamp correction */
                tsdiff = (tsmeas - tsref) & 0x3fffffffffffff;  // truncate to 54-bit LSB per timestamp spec
                tscorr = ((ll)(tsdiff >> FCORR_TBITS1) * fcorr) >> FCORR_TBITS2;
                ts = tsmeas + tscorr + tsoverflowcorr;

                /* write corrected timestamp to output buffer */
                eventptr->high = ts >> 22;
                eventptr->low = (ts << 10) | (low & 0x3ff);
#ifdef __DEBUG__
                if (isdebugverbose) {
                    fprintf(stderr, "[debug] Raw event - %08x %08x\n", high, low);
                    fprintf(stderr, "[debug] |   t_i: %014llx (%020llu)\n", tsmeas, tsmeas);
                    fprintf(stderr, "[debug] |  t'_i: %014llx (%020llu)\n", ts, ts);
                    fprintf(stderr, "[debug] +---------- %08x %08x\n", eventptr->high, eventptr->low);
                    fflush(stderr);
                }
#endif
                if (islegacy) {
                    _swp = eventptr->low;
                    eventptr->low = eventptr->high;
                    eventptr->high = _swp;
                }
                outbuffer[outevents++] = *eventptr;
                eventptr++;
            }

            /* accumulate timestamp corrections across batches */
#ifdef __SIZEOF_INT128__
            _tsdiff = ((u128)tsdiff) << FCORR_OFLOWRESBITS;
            _tscorr = ((i128)(_tsdiff >> FCORR_TBITS1) * fcorr) >> FCORR_TBITS2;
            _tsoverflowcorr += _tscorr;  // accumulates using higher resolution units
            tsoverflowcorr = (ll)(_tsoverflowcorr >> FCORR_OFLOWRESBITS);
#else
            tsoverflowcorr += tscorr;
#endif

            /* update reference timestamp to keep within 20 hour overflow condition */
            tsref = tsmeas;
        }

        // Read frequency correction values
        // - Note partial buffer must be maintained, since there is no other
        //   check to verify integrity of values broken up between separate reads.
        // - Note also this falls after reading the input buffer, but there is no
        //   particular reason why this order is chosen.
        if (freqhandle && FD_ISSET(freqhandle, &rfds)) {
            freqbytesread_next = read(freqhandle, freqbuffer_next, FREQBUFSIZE - freqbytespartial - 1);

            // File/pipe closed -> proceed without any further fcorr updates
            if (freqbytesread_next == 0) {
                freqhandle = 0;
                wmsg_s(FNAMELENGTH+25, "File/pipe '%s' closed.", freqfilename);
                // no break here
            }
            if (freqbytesread_next == -1) {
                wmsg_i(32, "Error %d on freqhandle read.", errno);
                break;
            }

            /* concatenate new data */
            freqbytesread = freqbytespartial + freqbytesread_next;

            /* search for valid fcorr values */
            int next_num_idx = 0, fcorr_tmp;
            for (i = 0; i < freqbytesread; i++) {
                if (freqbuffer[i] == '\n') {
                    freqbuffer[i] = 0;  // zero-terminate for readint
                    fcorr_tmp = readint(&freqbuffer[next_num_idx]);
                    if (isdecimal) fcorr_tmp = (int) fcorr_tmp * FCORR_DTOB;
                    if (isupdate) fcorr_tmp = (int)round(((1+fcorr*FCORR_BTO1) * (1+fcorr_tmp*FCORR_BTO1) - 1)/FCORR_BTO1);
                    if (abs(fcorr_tmp) < FCORR_MAX) {
                        fcorr = fcorr_tmp;
#ifdef __DEBUG__
                        if (isdecimal) {
                            fprintf(stderr, "[debug] 'fcorr' updated to '%.3f' ppb.\n", fcorr / FCORR_DTOB / 10);
                        } else {
                            fprintf(stderr, "[debug] 'fcorr' updated to '%d' x 2^-34.\n", fcorr);
                        }
                        fflush(stderr);
#endif
                    }
                    next_num_idx = i+1;  // continue reading
                }
            }

            /* clear parsed numbers from freqbuffer */
            if (next_num_idx > 0) {
                freqbytespartial = freqbytesread - next_num_idx;
                for (i = freqbytesread - freqbytespartial, j = 0; j < freqbytespartial; i++, j++) {
                    freqbuffer[j] = freqbuffer[i];
                }
                freqbuffer_next = &freqbuffer[freqbytespartial];
            }
        }

        // TODO: Shift this back to select call to write only when pipe available,
        //   and increase buffer size of output pipe in case write unavailable.
        //   By same measure, do not flush only when output buffer is full.
        /* write out events */
        retval = write(outhandle, outbuffer, eventnum * sizeof(struct rawevent));
#ifdef __DEBUG__
        if (isdebugverbose) {
            for (i = 0; i < eventnum; i++) {
                fprintf(stderr, "[debug] Verify: %08x %08x\n", outbuffer[i].high, outbuffer[i].low);
            }
        }
#endif
        if (retval != eventnum * sizeof(struct rawevent)) {
            wmsg_i(30, "Error %d on write.", errno);
            break;  // graceful close
        }
        outevents = 0;  // clear outbuffer only after successful write
        eventnum = 0;  // clear events to avoid rewriting:
                       // occurs when 'freqhandle' available for reading, but
                       // 'inhandle' has no more events
    }

    /* free buffers */
    free(inbuffer);
    free(outbuffer);
    if (freqfilename[0]) free(freqbuffer);
    return 0;
}
