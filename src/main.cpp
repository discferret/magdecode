#include <cstdio>
#include <cstdlib>
#include <limits>

using namespace std;

size_t LoadTrackImage(char *filename, unsigned int *buffer)
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
	// For now, just throw away anything with a count less than 10 or so.
	// FIXME: need a better algorithm for this! maybe calculate mean value of histo
	// and clip anything less than 1% of mean?
	for (size_t i=0; i<=(maxval-minval); i++) {
		if (histogram[i] < 10) histogram[i] = 0;
	}

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

		// if this delta is  negative and last delta was positive then we've
		// found a peak
		if ((delta <= 0) && (lastdelta > 0)) {
//			printf("peak: %3d  [delta %2d]\n", i+minval, delta);
			peaks[numpeaks++] = i;
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

	// clean-up
	delete[] histogram;
	return 0;
}
