#ifndef _HEXDUMP_H
#define _HEXDUMP_H

void hex_dump(void *data, int size);
/* dumps size bytes of *data to stdout. Looks like:
 * [0000] 75 6E 6B 6E 6F 77 6E 20   30 FF 00 00 00 00 39 00 unknown 0.....9.
 */

#endif // _HEXDUMP_H
