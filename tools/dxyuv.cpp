/*

YUF Converter
------------

Form:   Abekas YUV  (OpenDX output)
To:     RAW YUV420P (ffmpeg input)

See:    http://ffmpeg.sourceforge.net/

Author: Lisandro Dalcin

*/


#include <cstdlib>
#include <iostream>

// color component type
typedef unsigned char uint8_t;
// stream sizes type
typedef unsigned int  size_t;

using namespace std;

int main(int argc, char *argv[]) {

  if (argc != 2)
    cerr << "usage: " << argv[0] 
	 << " WIDTHxHEIGHTxFRAMES < istream  > ostream" << endl, exit(0);

  // get stream size info
  bool errsz = false;
  char *s = argv[1];

  const size_t width  = strtoul(s,&s,10);  // horizontal pixels
  if (*s) s++; else errsz = true;

  const size_t height = strtoul(s,&s,10);  // vertical pixels
  if (*s) s++; else errsz = true;

  const size_t frames = strtoul(s,&s,10);  // frames
  if ( !width | !height | !frames ) errsz = true;

  if (errsz)
    cerr << "erroneous sizes..."<< endl, exit(0);

  // pixels info
  const size_t pixsz = 2;            // bytes per pixel
  const size_t frmsz = width*height; // pixels per frame

  // allocate input buffer
  // (one frame, format UYVY, sampling 4:2:2)
  uint8_t *UYVY  = new uint8_t[frmsz*pixsz];

  // allocate output buffers
  // (one frame, format I420, sampling 4:1:1)
  uint8_t *Y = new uint8_t[frmsz  ]; 
  uint8_t *U = new uint8_t[frmsz/4]; 
  uint8_t *V = new uint8_t[frmsz/4];

  // loop over frames
  for (int f=0; f<frames; f++) {

    // pointers to output buffers
    uint8_t *pY = Y;
    uint8_t *pU = U;
    uint8_t *pV = V;

    // read frame from input stream
    cin.read((char *)UYVY,frmsz*pixsz);

    // get Y, U, V planes
    for (int i=0; i<frmsz; i++) {
      *pY++ = UYVY[1+i*2];
    }
    for (int j=0; j<height; j+=2)     // make vertical downsample
      for (int i=0; i<width/2; i++) { // already horizontally downsampled
	*pU++ = UYVY[0+i*4+j*width*pixsz];
	*pV++ = UYVY[2+i*4+j*width*pixsz];
      }

    // write Y, U, V planes to output stream
    cout.write((char *)Y,frmsz  );
    cout.write((char *)U,frmsz/4);
    cout.write((char *)V,frmsz/4);
  }

  // deallocate buffers
  delete[] UYVY;
  delete[] Y;
  delete[] U;
  delete[] V;


  return 0;
}
