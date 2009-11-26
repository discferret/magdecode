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
	unsigned int carry = 0;
	size_t count = 0;
	for (size_t s = 0; s < filelen; s++) {
		// check for counter overflow from DA
		if ((buf[s] & 127) == 0) {
			// counter overflow, increase the carry register
			if ((((unsigned long)carry) + 128) > int_limits.max()) {
				// detected data overflow
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

	delete buf;
	fclose(fp);

	return count;
}

// main fnc
int main(void)
{
	unsigned int buf[128*1024];
	size_t buflen;

	buflen = LoadTrackImage("memdump.bin", buf);
	printf("buflen = %lu\n", (unsigned long)buflen);

	// calculate RPM and data rate
	unsigned long buftm = 0;
	for (size_t i=0; i<buflen; i++) {
		buftm += buf[i];
	}
	printf("total timing val = %lu\n", buftm);
	printf("time(secs) = %f\n", ((float)buftm) * 25e-9);
	printf("time(ms) = %f\n", ((float)buftm) * 25e-9 * 1000);
	printf("est rpm = %f\n", 60 * (1/(((float)buftm) * 25e-9)));

	// calculating the data rate involves doing a histogram analysis and finding
	// the first peak.
	// TODO: implement this
}
