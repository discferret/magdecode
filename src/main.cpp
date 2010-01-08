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
int main(void)
{
	unsigned int buf[128*1024];
	size_t buflen;
	size_t maxval = 0;
	size_t minval = ((size_t)-1);

	buflen = LoadTrackImage("memdump.bin", buf);
	printf("buflen = %lu\n", (unsigned long)buflen);

	// calculate RPM and data rate
	unsigned long buftm = 0;
	for (size_t i=0; i<buflen; i++) {
		buftm += buf[i];
		if (maxval < buf[i]) maxval = buf[i];
		if (minval > buf[i]) minval = buf[i];
	}
	printf("total timing val = %lu\n", buftm);
	printf("time(secs) = %f\n", ((float)buftm) * 25e-9);
	printf("time(ms) = %f\n", ((float)buftm) * 25e-9 * 1000);
	printf("est rpm = %f\n", 60 * (1/(((float)buftm) * 25e-9)));
	printf("\n");
	printf("maxval = %d\nminval = %d\nspan   = %d\n", maxval, minval, maxval - minval);
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

	// Start by clipping histogram noise. Need to find outliers, then eliminate them.
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
	int numpeaks = 0;
	for (size_t i=1; i<=(maxval - minval); i++) {
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
		printf("%d peaks found:\n", numpeaks);
		for (int i=0; i<numpeaks; i++) {
			printf("\tpeak #%d: %3d\n", i+1, peaks[i]+minval);
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

	// Data decoder begins here.

	// Current 1t reference. Assumed to be the position of the 1st peak.
	float t = peaks[0] + minval;
	// Bit-vector to store MFM stream
	vector<bool> mfmbits;

	// Iterate over input stream and decode
	for (size_t i=0; i<buflen; i++) {
		// Calculate error values for this timing value vs. 1t, 1.5t and 2.0t
		float error_t10 = fabs((buf[i] / 1.0) - t);
		float error_t15 = fabs((buf[i] / 1.5) - t);
		float error_t20 = fabs((buf[i] / 2.0) - t);

		// Figure out which error is the lowest
		if ((error_t10 < error_t15) && (error_t10 < error_t20)) {
			// t1.0 is the lowest. "01" sequence.
			mfmbits.push_back(false);
			mfmbits.push_back(true);

			// Update reference T
			t = (0.95 * t) + (0.05 * (buf[i] / 1.0));
		} else if ((error_t15 < error_t10) && (error_t15 < error_t20)) {
			// t1.5 is the lowest. "001" sequence.
			mfmbits.push_back(false);
			mfmbits.push_back(false);
			mfmbits.push_back(true);

			// Update reference T
			t = (0.95 * t) + (0.05 * (buf[i] / 1.5));
		} else {
			// t2.0 is the lowest. "0001" sequence.
			mfmbits.push_back(false);
			mfmbits.push_back(false);
			mfmbits.push_back(false);
			mfmbits.push_back(true);

			// Update reference T
			t = (0.95 * t) + (0.05 * (buf[i] / 2.0));
		}
	}

	printf("mfmbits count = %d\n", mfmbits.size());


	// Now process the MFM bitstream to find the sync markers
	unsigned long bits = 0;
	unsigned int num_idam = 0, num_dam = 0;
	for (size_t i=0; i<mfmbits.size(); i++) {
		size_t dump=0;

		// get next bit
		bits = (bits << 1) + (mfmbits[i] ? 1 : 0);

		// compare buffer with sync-longword (sync-A1 : FE, aka IDAM)
		if (bits == 0x44895554) {
			// ID Address Mark
			// i+1 because "i" is the last bit of the IDAM marker; we want the
			// first bit of the new data byte (encoded word).
			printf("IDAM at %d\n", i+1);
			num_idam++;
			dump = 5;
		} else if (bits == 0x44895545) {
			// Data Address Mark
			// i+1 because "i" is the last bit of the DAM marker; we want the
			// first bit of the new data byte (encoded word).
			printf("DAM at %d\n", i+1);
			num_dam++;
			dump = 512;
		}

		if (dump > 0) {
			// dump next few bytes of data in hex
			// TODO: use hex_dump() and a char array instead
			// i+1 because "i" is the last bit of the (I)DAM marker; we want
			// the first bit of the new data byte (encoded word).
			//
			char *buffer = new char[dump];
			for (size_t x=0; x<dump; x++) {
				buffer[x] = decodeMFM(mfmbits, i+(x*16)+1);
			}
			hex_dump(buffer, dump);
			delete[] buffer;
		}
	}

	printf("Seen: %d IDAMs, %d DAMs\n", num_idam, num_dam);

	// clean-up
	delete[] histogram;
	return 0;
}
