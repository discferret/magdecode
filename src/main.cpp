#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cctype>
#include <limits>

// cstdint would be better, but it's in C++0x; gcc's 0x support is still beta.
#include <stdint.h>

extern "C" {
#include "hexdump.h"
}

using namespace std;

// CRC-16/CCITT implementation
// Based on code from http://www.sanity-free.org/133/crc_16_ccitt_in_csharp.html
class CRC16 {
	private:
		static const unsigned int poly = 4129;
		unsigned int table[256];
		unsigned int crcval, initval;
	public:
		CRC16(unsigned int initval = 0xFFFF) {
			this->initval = this->crcval = initval;

			unsigned int temp, a;
			for (size_t i=0; i<(sizeof(this->table)/sizeof(this->table[0])); ++i) {
				temp = 0;

				a = (unsigned int) i << 8;
				
				for (int j=0; j < 8; ++j) {
					if(((temp ^ a) & 0x8000) != 0) {
						temp = (unsigned int)((temp << 1) ^ poly);
					} else {
						temp <<= 1;
					}
					a <<= 1;
				}

				table[i] = temp;
			}
		}

		// calculate crc without updating internal state
		unsigned int calculate(void *buf, size_t len) {
			unsigned int crc = this->crcval;

			for (size_t i=0; i<len; i++) {
				char x = ((unsigned char *)buf)[i];
				crc = ((crc << 8) ^ this->table[((crc >> 8) ^ x) & 0xff]) & 0xffff;
			}

			return crc;
		}

		// calculate crc and update internal state afterwards
		// used for partial crcs
		unsigned int update(void *buf, size_t len) {
			this->crcval = this->calculate(buf, len);
			return this->crcval;
		}

		// reset crc to initialisation default
		void reset(void) {
			this->crcval = this->initval;
		}

		// reset crc to arbitrary value and change init default accordingly
		void reset(unsigned int initval) {
			this->crcval = this->initval = initval;
		}
};

// Load a track image from a file
// TODO: rewrite this using fstream
size_t LoadTrackImage(const char *filename, unsigned int *buffer)
{
	FILE *fp;
	size_t filelen;
	char *buf;
	numeric_limits<unsigned int> int_limits;

	// open the file
	fp=fopen(filename, "rb");
	if (!fp) return -1;

	// find out how big the file is
	fseek(fp, 0, SEEK_END);
	filelen = ftell(fp);
	fseek(fp, 0, SEEK_SET);

	// read entire file into RAM
	if ((buf = new char[filelen]) == NULL) return -1;
	if (fread(buf, 1, filelen, fp) != filelen) return -1;

	// iterate over loop, build timing values
	// This is a little more involved than just copying from one buffer to
	// another. The DA100 will store a '0x00' byte every time the delay
	// counter wraps around (from 127 to 0). Thus, we maintain a carry count
	// (incremented by 128 for each zero byte) which is added onto the next
	// non-zero byte. We also track instances where carrying the value would
	// cause an integer overflow in the data storage buffer -- this causes
	// the function to return -1.
	unsigned int carry = 0;
	size_t count = 0;
	for (size_t s = 0; s < filelen; s++) {
		// check for counter overflow from DA
		if ((buf[s] & 127) == 0) {
			// counter overflow, increase the carry register
			if ((((unsigned long)carry) + 128) > int_limits.max()) {
				// detected data overflow
				delete[] buf;
				return -1;
			} else {
				// otherwise keep carrying
				carry += 128;
			}
		} else {
			// no overflow, just a normal count. store the count and clear
			// the carry counter.
			buffer[count++] = carry + (buf[s] & 127);
			carry = 0;
		}
	}

	// if carry more than zero, dump it into the buffer
	if (carry > 0) {
		buffer[count++] = carry;
	}

	// free the file buffer
	delete[] buf;

	fclose(fp);

	return count;
}

// Decode 16 MFM bits into a single byte
unsigned char decodeMFM(vector<bool> bits, size_t startpos)
{
	uint8_t buf;
	int cnt = 1;

	// Get 8 bits of data, discard every other bit
	for (vector<bool>::iterator i = bits.begin()+startpos; i<bits.begin()+startpos+16; i++) {
		if ((cnt % 2) == 0) {
			buf = (buf << 1) + (*i ? 1 : 0);
		}
		cnt++;
	}

	return buf;
}

// main fnc
int main(int argc, char **argv)
{
	unsigned int buf[128*1024];
	size_t buflen;
	size_t maxval = 0;
	size_t minval = ((size_t)-1);

	const double CLK_FRQ = 100e6;
	const double CLK_TM = 1.0 / CLK_FRQ;

	if (argc < 2) {
		cout << "syntax: " << argv[0] << " filename\n";
		return -1;
	}

	// limit scope of 'x'
	{
		ssize_t x = LoadTrackImage(argv[1], buf);
		if (x < 1) {
			cout << "error reading input file \"" << argv[1] << "\", code " << x << "\n";
			return -1;
		}
		buflen = x;
		printf("buflen = %lu\n", (unsigned long)buflen);
	}

	// calculate RPM and data rate
	unsigned long buftm = 0;
	for (size_t i=0; i<buflen; i++) {
		buftm += buf[i];
		if (maxval < buf[i]) maxval = buf[i];
		if (minval > buf[i]) minval = buf[i];
	}
	printf("total timing val = %lu\n", buftm);
	printf("time(secs) = %f\n", ((float)buftm) * CLK_TM);
	printf("time(ms) = %f\n", ((float)buftm) * CLK_TM * 1000);
	printf("est rpm = %f\n", 60 * (1/(((float)buftm) * CLK_TM)));
	printf("\n");
	printf("maxval = %lu\nminval = %lu\nspan   = %lu\n", maxval, minval, maxval - minval);
	printf("\n");

	// allocate memory for a histogram
	// note that we only allocate memory for 'bins' inside the min-max range
	// calculated above. this saves a bit of memory.
	unsigned int *histogram = new unsigned int[maxval - minval + 1];
	for (size_t i=0; i<=(maxval-minval); i++) histogram[i] = 0;

	// generate the histogram
	for (size_t i=0; i<buflen; i++) {
		histogram[buf[i] - minval]++;
	}

	// Calculate the mean value of the histogram
	float mean = 0;
	for (size_t i=0; i<=(maxval-minval); i++) {
		mean += histogram[i];
	}
	mean /= ((maxval - minval) + 1);
	printf("mean = %0.5f\n", mean);

	// Calculate the standard deviation
	float sd = 0, osqr = 0;
	for (size_t i=0; i<=(maxval - minval); i++) {
		// calculate deviation from mean
		sd = histogram[i] - mean;
		// square the deviation
		sd = sd * sd;
		// add to the variance accumulator
		osqr += sd;
	}

	// Calculate the mean variance
	osqr /= ((float)(maxval - minval + 1));

	printf("o^2 (mean variance) = %0.2f\n", osqr);

	// Take the square root of the variance to get the standard deviation
	printf("standard deviation = %0.4f\n", sqrtf(osqr));

#if 0
	// DEBUG: print out the histogram
	printf("histpoint,time,count\n");
	for (size_t i=0; i<=(maxval - minval); i++) {
		// calculate time in microseconds
		// one count is 1/(40e6) seconds; mul by 1e6 to get usecs
		float t = ((i+minval)*(1.0/40.0e6)) * 1.0e6;
		printf("%d,%0.2f,%d\n", i+minval, t, histogram[i]);
	}
#endif

	// Apply a moving average to the histogram data
	// This eliminates small spikes in the data set, and makes the peak detection
	// a bit more reliable.

	do {
		// Start by making a copy of the original histogram
		unsigned int *histogram_copy = new unsigned int[maxval - minval + 1];
		for (size_t z=0; z<(maxval-minval+1); z++) histogram_copy[z] = histogram[z];

		// Number of previous values to take into account
		const int HIST_PKDET_AVG_BINS = 16;

		// Apply moving average to histogram data
		for (size_t x=1; x<(maxval-minval+1); x++) {
			unsigned long accum = 0;
			unsigned int n = 0;
			for (ssize_t y=(x-HIST_PKDET_AVG_BINS); y<(ssize_t)x; y++) {
				if (y>=0) {
					accum += histogram_copy[y];
					n++;
				}
			}
			if (n>0) accum /= n; else accum = 0;
			histogram[x] = accum;
		}

		delete[] histogram_copy;
	} while (0);

	// Clip histogram noise. Need to find outliers, then eliminate them.
	// For now, just throw away anything with a count less than 1/10th of the mean.
	for (size_t i=0; i<=(maxval-minval); i++) {
		if (histogram[i] < round(mean * 0.1)) histogram[i] = 0;
	}

#if 0
	// DEBUG: print out the clipped histogram
	printf("histpoint\ttime\tcount\n");
	for (size_t i=0; i<=(maxval - minval); i++) {
		// calculate time in microseconds
		// one count is 1/(40e6) seconds; mul by 1e6 to get usecs
		float t = ((i+minval)*(1.0/40.0e6)) * 1.0e6;
		printf("%d\t%0.2f\t%d\n", i+minval, t, histogram[i]);
	}
#endif

	// Peaks are found by measuring the slope of the current and previous histogram
	// points. To put it another way:
	//           (y1-y2)          change in y
	//  slope = ---------  or  = -------------
	//           (x1-x2)          change in x
	//
	// X is the "time" axis (i+minval, or just i), Y is the "count" axis
	// (histogram[i]).
	//
	// Last histogram value
	size_t lasthist = histogram[0];
	// Last delta (dY/dX)
	ssize_t lastdelta = 0;
	// Positions of detected peaks
	size_t peaks[32];
	// Number of detected peaks
	size_t numpeaks = 0;
	for (size_t i=1; i<=(maxval - minval); i++) {
		// break if max number of peaks has been exceeded
		if (numpeaks > (sizeof(peaks)/sizeof(peaks[0]))) {
			printf("WARNING: Maximum number of peaks exceeded!\n");
			break;
		}

		// calculate delta
		ssize_t delta = histogram[i] - lasthist;

		// if this delta is negative and last delta was positive then we've
		// found a peak. Note that the peak is actually in the last histogram
		// 'bin'; we're on a negative slope here, thus the 'last' point was
		// the highest point (peak) as far as we're concerned.
		if ((delta <= 0) && (lastdelta > 0)) {
			peaks[numpeaks++] = i-1;
		}

		// Update last-histogram-point and last-delta
		lasthist = histogram[i];
		lastdelta = delta;
	}

	// List all found peaks
	if (numpeaks > 0) {
		printf("%lu peaks found:\n", numpeaks);
		for (size_t i=0; i<numpeaks; i++) {
			printf("\tpeak #%lu: %3lu\n", i+1, peaks[i]+minval);
		}
	} else {
		printf("No peaks found.\n");
	}

	// Loop over the peaks and see if we have a cluster which are close to the
	// 1t:1.5t:2t profile of MFM.
	//
	// TODO: implement this

	// Basically, we:
	//   - Start at 1.5t and work backwards (and forwards) up to 1t.
	//   - Store the position of the closest peak to 1.5t
	//   - Repeat for 2t
	//
	// The decoder algorithm doesn't need to know where the peaks are, but we
	// can use peak-position testing as part of the Confidence Value
	// calculation. That is, "how confident are we that this really is MFM data?"
	//
	// Note that the MFM CODEC ideally needs to know where the 1t peak is during
	// initialisation. Otherwise we need to know the data rate and sample rate,
	// and set it to a sane default.
	//
	// Also note that the MFM CODEC can handle FM too. So a 2-peak (1t:2t)
	// histogram is a valid FM histogram, and a 3-peak (1t:1.5t:2t) histogram
	// is a valid MFM histogram.

	// Bit-vector to store MFM stream
	vector<bool> mfmbits;

	// Data decoder begins here.
#if 0
	// Current 1t reference. Assumed to be the position of the 1st peak.
	float t = peaks[0] + minval;

	// Iterate over input stream and decode
	const float change_frac = 0.05;
	for (size_t i=0; i<buflen; i++) {
		// Calculate error values for this timing value vs. 1t, 1.5t and 2.0t
		float error_t10 = fabs((buf[i] / 1.0) - t);
		float error_t15 = fabs((buf[i] / 1.5) - t);
		float error_t20 = fabs((buf[i] / 2.0) - t);
		float t_mult;

		// Figure out which error is the lowest
		if ((error_t10 < error_t15) && (error_t10 < error_t20)) {
			// t1.0 is the lowest. "01" sequence.
			mfmbits.push_back(false);
			mfmbits.push_back(true);
			t_mult = 1.0;
		} else if ((error_t15 < error_t10) && (error_t15 < error_t20)) {
			// t1.5 is the lowest. "001" sequence.
			mfmbits.push_back(false);
			mfmbits.push_back(false);
			mfmbits.push_back(true);
			t_mult = 1.5;
		} else {
			// t2.0 is the lowest. "0001" sequence.
			mfmbits.push_back(false);
			mfmbits.push_back(false);
			mfmbits.push_back(false);
			mfmbits.push_back(true);
			t_mult = 2.0;
		}

		// Update reference T
		t = ((1.0 - change_frac) * t) + (change_frac * ((float)buf[i] / t_mult));
	}
#else

#ifdef VCD
	FILE *vcd = fopen("values.vcd", "wt");
	fprintf(vcd, "$version DiscFerret Analyser D2/DPLL 0.1 $end\n"
			"$timescale 1 ns $end\n"
			"$var reg 1 * clock $end\n"
			"$var reg 1 ' pll_clock $end\n"
			"$var reg 1 ! rd_data $end\n"
			"$var reg 1 %% rd_data_latched $end\n"
			"$var reg 1 ^ shaped_data $end\n"
			"$var reg 8 & pjl_shifter $end\n"
			"$var reg 1 ( data_window $end\n"
			"$upscope $end\n"
			"$enddefinitions $end\n"
			"$dumpall\n"
			"0*\n"
			"0'\n"
			"0!\n"
			"0%%\n"
			"0^\n"
			"b00000000 &\n"
			"0(\n"
			"$end\n"
		   );
#endif
	/**
	 * This is a software implementation of Jim Thompson's Phase-Jerked Loop
	 * design, available from AnalogInnovations.com as the PDF file
	 * "FloppyDataExtractor.pdf".
	 *
	 * This consists of:
	 *   A data synchroniser which forces RD_DATA to be active for 2 clock cycles.
	 *   A counter which increments constantly while the PLL is running, and is 
	 *     reset to zero when a data bit is detected.
	 *   A flip-flop which changes state when the counter reaches half its maximum
	 *     value
	 *
	 * The actual values of NSECS_PER_ACQ and NSECS_PER_PLLCK don't really matter.
	 * What actually matters is the ratio between the two -- if you have a 40MHz
	 * acquisition clock and a PLL clock of 16MHz (data rate 500kbps), then the
	 * starting values will be NSECS_PER_ACQ=25 and NSECS_PER_PLLCK=62.5. Problem
	 * is, 62.5 isn't an integer multiple, so we might have issues with
	 * short-term clock jitter. So we multiply by two, which gives us
	 * NSECS_PER_ACQ=50 and NSECS_PER_PLLCK=125, and a timestep of 0.5ns. Much
	 * better.
	 *
	 * We can also change the PJL Counter maximum value if it makes the math
	 * work out better. Be careful though -- reducing the value WILL reduce the
	 * number of available phase-shift steps and thus the PLL accuracy.
	 *
	 * Now we know the ticks-per-acqclk and ticks-per-pllclk values, we can
	 * figure out the optimal timer increment --
	 *   TIMER_INCREMENT = gcd(NSECS_PER_ACQ, NSECS_PER_PLLCK)
	 *   (gcd == greatest common divisor)
	 *
	 * Ideally we want a TIMER_INCREMENT greater than 1. If we get an increment
	 * of 1, then the loop has to run at 1x speed and will be SLOW. Try
	 * increasing NSECS_PER_ACQ and NSECS_PER_PLLCK (multiply them by the same
	 * number e.g. 2, 4, 8, ...), then run the gcd again. Problem is, this isn't
	 * going to gain much if anything in speed because you're going to be running
	 * more loops at a faster rate, thus it's a zero-gain :-/
	 */

	// Nanoseconds counters. Increment once per loop or "virtual" nanosecond.
	unsigned long nsecs1 = 0, nsecs2=0;
	// Number of nanoseconds per acq tick -- (1/freq)*1e9. This is valid for 40MHz.
	const unsigned long NSECS_PER_ACQ = CLK_TM; //(1e9 / 100e6);
	// Number of nanoseconds per PLLCK tick -- (1/16e6)*1e9. 16MHz. 
	// This should be the reciprocal of 32 times the data rate in kbps, multiplied
	// by 1e9 to get the time in nanoseconds.
	// That is, (1/(TRANSITIONS_PER_BITCELL * PJL_COUNTER_MAX * DATA_RATE))*1e9
	const unsigned long NSECS_PER_PLLCK = (1e9 / 32e6);
	// Number of clock increments per loop (timing granularity). Best-case value
	// for this is gcd(NSECS_PER_ACQ, NSECS_PER_PLLCK).
	const unsigned long TIMER_INCREMENT = 1;
	// Maximum value of the PJL counter. Determines the granularity of phase changes.
	const unsigned char PJL_COUNTER_MAX = 64;

	// Iterator for data buffer
	size_t i = 0;

	// True if RD_DATA was high in this clock cycle
	bool rd_latch = false;
	// Same but only active for 2Tcy (SHAPED_DATA)
	int shaped_data = 0;

	// Phase Jerked Loop counter
	unsigned char pjl_shifter = 0;

	// data window
	unsigned char data_window = 0;

	do {
#ifdef VCD
		{ static unsigned long long clocks = 0; fprintf(vcd, "#%lld\n", clocks++); }
#endif
		// Increment counters
		nsecs1 += TIMER_INCREMENT;
		nsecs2 += TIMER_INCREMENT;
#ifdef VCD
		if ((nsecs1 % NSECS_PER_ACQ) == 0) { fprintf(vcd, "1*\n"); }
		if ((nsecs1 % NSECS_PER_ACQ) == (NSECS_PER_ACQ/2)) { fprintf(vcd, "0*\n"); }
		if (nsecs1 == (buf[i] * NSECS_PER_ACQ)) { fprintf(vcd, "0!\n"); } else { fprintf(vcd, "1!\n"); }
		if ((nsecs1 % NSECS_PER_PLLCK) == 0) { fprintf(vcd, "1'\n"); }
		if ((nsecs1 % NSECS_PER_PLLCK) == (NSECS_PER_PLLCK/2)) { fprintf(vcd, "0'\n"); }
#endif
		// Loop 1 -- floppy disc read channel
		if (nsecs1 >= (buf[i] * NSECS_PER_ACQ)) {
			// Flux transition. Set the read latches.
			rd_latch = true;
			shaped_data = 2;

			// Update nanoseconds counter for read channel, retain error factor
			nsecs1 -= (buf[i] * NSECS_PER_ACQ);

			// Update buffer position
			i++;
		}
#ifdef VCD
		fprintf(vcd, "%d%%\n", rd_latch);
		fprintf(vcd, "%d^\n", (shaped_data > 0) ? 1 : 0);
#endif
		// Loop 2 -- DPLL channel
		if (nsecs2 >= NSECS_PER_PLLCK) {
			// Update nanoseconds counter for PLL, retain error factor
			nsecs2 -= NSECS_PER_PLLCK;

			// PJL loop
			if (shaped_data > 0) {
				pjl_shifter = 0;
			} else {
				// increment shifter
				pjl_shifter = (pjl_shifter + 1) % PJL_COUNTER_MAX;
			}

			// DWIN detect
			if (pjl_shifter == (PJL_COUNTER_MAX / 2)) {
				// Data window toggle. Latch the current RD_LATCH blob into the output buffer.
				mfmbits.push_back(rd_latch);
				// Clear the data latch ready for the next data window.
				rd_latch = false;
				// Update DWIN
				data_window ^= 0x01;
			}

			// Update shaped-data time counter
			if (shaped_data > 0) shaped_data--;
		}

#ifdef VCD
		fprintf(vcd, "%d(\n", data_window);
		fputc('b', vcd);
		for (int x=0; x<8; x++) if (pjl_shifter & (1<<(7-x))) fputc('1', vcd); else fputc('0', vcd);
		fprintf(vcd, " &\n");
		if (mfmbits.size() == 1000) { printf("Got 1000 bits, terminating run.\n"); break; }
#endif
	} while (i < buflen);

#ifdef VCD
	fclose(vcd);
#endif
#endif
	printf("mfmbits count = %lu\n", mfmbits.size());


	// Now process the MFM bitstream to find the sync markers
	unsigned long bits = 0;
	unsigned int num_idam = 0, num_dam = 0;
	size_t next_data_dump = 0;
	bool chk_data_crc = false;
	for (size_t i=0; i<mfmbits.size(); i++) {
		size_t dump=0;

		// get next bit
		bits = ((bits << 1) + (mfmbits[i] ? 1 : 0)) & 0xffffffff;

		// compare buffer with sync-longword (sync-A1 : FE, aka IDAM)
		if (bits == 0x44895554) {
			// ID Address Mark
			// i+1 because "i" is the last bit of the IDAM marker; we want the
			// first bit of the new data byte (encoded word).
			printf("IDAM at %lu\n", i+1);
			num_idam++;
			dump = 6;
			chk_data_crc = false;

			// decode the IDAM
			unsigned char *idambuf = new unsigned char[6];
			for (size_t x=0; x<6; x++) {
				idambuf[x] = decodeMFM(mfmbits, i+(x*16)+1);
			}

			printf("\tIDAM = Track %2d, Side %d, Sector %2d; sector size ",
					idambuf[0], idambuf[1], idambuf[2]);
			switch (idambuf[3]) {
				case 0x00:
				case 0x01:
				case 0x02:
				case 0x03:
					printf("%d", 1 << (7+idambuf[3]));
					next_data_dump = 1<<(7+idambuf[3]);
					break;
				default:
					printf("unknown (0x%02X)", idambuf[3]);
					next_data_dump = 0;
					break;
			}

			// check the CRC
			CRC16 c = CRC16();
			c.update((char *)"\xA1\xA1\xA1\xFE", 4);
			unsigned int crc = c.update(idambuf, 4);
			printf("; CRC=%04X %s\n", (idambuf[4] << 8) | idambuf[5],
					((unsigned int)((idambuf[4] << 8) | idambuf[5])==crc) ? "(ok)" : "BAD");

			delete idambuf;
		} else if (bits == 0x44895545) {
			// Data Address Mark
			// i+1 because "i" is the last bit of the DAM marker; we want the
			// first bit of the new data byte (encoded word).
			printf("DAM at %lu%s\n", i+1, (next_data_dump == 0) ? " [ERR: no preceding IDAM]" : "");
			num_dam++;
			dump = next_data_dump;
			next_data_dump = 0;
			chk_data_crc = true;
		}

		if (dump > 0) {
			// dump next few bytes of data in hex
			// TODO: use hex_dump() and a char array instead
			// i+1 because "i" is the last bit of the (I)DAM marker; we want
			// the first bit of the new data byte (encoded word).
			//
			unsigned char *buffer = new unsigned char[dump+2];
			for (size_t x=0; x<dump+2; x++) {
				buffer[x] = decodeMFM(mfmbits, i+(x*16)+1);
			}
			hex_dump(buffer, dump);

			if (chk_data_crc) {
				CRC16 c = CRC16();
				c.update((char *)"\xA1\xA1\xA1\xFB", 4);
				unsigned int crc = c.update(buffer, dump);
				printf("\tData record CRC=%04X %s\n", (buffer[dump+0] << 8) | buffer[dump+1],
						((unsigned int)((buffer[dump+0] << 8) | buffer[dump+1])==crc) ? "(ok)" : "BAD");
			}

			delete[] buffer;
		}
	}

	printf("Seen: %d IDAMs, %d DAMs\n", num_idam, num_dam);

	// clean-up
	delete[] histogram;
	return 0;
}
