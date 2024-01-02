#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "systemc.h"

#define VERBOSE 1

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0
#define BOOSTBLURFACTOR 90.0
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SIGMA 0.6
#define TLOW  0.3
#define THIGH 0.8

#define COLS 3840
#define ROWS 2160
#define SIZE COLS*ROWS
#define VIDEONAME "Keck"
#define IMG_IN    "video/" VIDEONAME "%03d.pgm"
#define IMG_OUT   VIDEONAME "%03d_edges.pgm"
#define IMG_NUM   30 /* number of images processed (1 or more) */
#define AVAIL_IMG 30 /* number of different image frames (1 or more) */

/* upper bound for the size of the gaussian kernel
 * SIGMA must be less than 4.0
 * check for 'windowsize' below
 */
#define WINSIZE 21

/***************************************************
 *               IMAGE DEFINITION                  *
 **************************************************/
typedef struct Image_s
{
    unsigned char img[SIZE];

    Image_s(void)
    {
        for (int i=0; i<SIZE; i++)
        {
            img[i] = 0;
        }
    }

    Image_s& operator=(const Image_s& copy)
    {
        for (int i=0; i<SIZE; i++)
        {
            img[i] = copy.img[i];
        }
        return *this;
    }

    operator unsigned char*()
    {
        return img;
    }

    unsigned char& operator[](const int index)
    {
        return img[index];
    }
} IMAGE;

typedef struct SImage_s {
    short int img[SIZE];

    SImage_s(void) {
        for (int i = 0; i < SIZE; i++) {
            img[i] = 0;
        }
    }

    SImage_s& operator=(const SImage_s& copy) {
        for (int i = 0; i < SIZE; i++) {
            img[i] = copy.img[i];
        }
        return *this;
    }

    operator short int*() {
        return img;
    }

    short int& operator[](const int index) {
        return img[index];
    }
} SIMAGE;


/***************************************************
 *               STIMULUS MODULE                   *
 **************************************************/
SC_MODULE(Stimulus)
{
    // PX_I/O: PORT X=MODULENAME INPUT/OUTPUT
    sc_fifo_out<IMAGE> PS_O;
    IMAGE image;
    
    /******************************************************************************
    * Function: read_pgm_image
    * Purpose: This function reads in an image in PGM format. The image can be
    * read in from either a file or from standard input. The image is only read
    * from standard input when infilename = NULL. Because the PGM format includes
    * the number of columns and the number of rows in the image, these are read
    * from the file. Memory to store the image is allocated OUTSIDE this function.
    * The found image size is checked against the expected rows and cols.
    * All comments in the header are discarded in the process of reading the
    * image. Upon failure, this function returns 0, upon sucess it returns 1.
    ******************************************************************************/
    int read_pgm_image(const char *infilename, unsigned char *image, int rows, int cols)
    {
       FILE *fp;
       char buf[71];
       int r, c;

       /***************************************************************************
       * Open the input image file for reading if a filename was given. If no
       * filename was provided, set fp to read from standard input.
       ***************************************************************************/
       if(infilename == NULL) fp = stdin;
       else{
          if((fp = fopen(infilename, "r")) == NULL){
             fprintf(stderr, "Error reading the file %s in read_pgm_image().\n",
                infilename);
             return(0);
          }
       }

       /***************************************************************************
       * Verify that the image is in PGM format, read in the number of columns
       * and rows in the image and scan past all of the header information.
       ***************************************************************************/
       fgets(buf, 70, fp);
       if(strncmp(buf,"P5",2) != 0){
          fprintf(stderr, "The file %s is not in PGM format in ", infilename);
          fprintf(stderr, "read_pgm_image().\n");
          if(fp != stdin) fclose(fp);
          return(0);
       }
       do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */
       sscanf(buf, "%d %d", &c, &r);
       if(c != cols || r != rows){
          fprintf(stderr, "The file %s is not a %d by %d image in ", infilename, cols, rows);
          fprintf(stderr, "read_pgm_image().\n");
          if(fp != stdin) fclose(fp);
          return(0);
       }
       do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */

       /***************************************************************************
       * Read the image from the file.
       ***************************************************************************/
       if((unsigned)rows != fread(image, cols, rows, fp)){
          fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
          if(fp != stdin) fclose(fp);
          return(0);
       }

       if(fp != stdin) fclose(fp);
       return(1);
    }

    void stimulus(void)
    {
        cout << "IN STIM THREAD!" << endl;
        char infilename[70];
        int i=0, n=0;

        // for loop for write
        for (int i = 0; i < IMG_NUM; i++)
        {
            /* code */
            n = i % AVAIL_IMG;
            sprintf(infilename, IMG_IN, n+1);
            /****************************************************************************
            * Read in the image. This read function allocates memory for the image.
            ****************************************************************************/
            if(VERBOSE) printf("Reading the image %s.\n", infilename);
            if(read_pgm_image(infilename, image, ROWS, COLS) == 0){
                fprintf(stderr, "Error reading the input image, %s.\n", infilename);
                exit(1);
            }
            else
                PS_O.write(image);

        }

    }

    // Stimulus Constructor
    SC_CTOR(Stimulus):
    PS_O("PS_O") // Only initialize ports here
    {
        SC_THREAD(stimulus);
        set_stack_size(128*1024*1024);
    };

};

/***************************************************
 *                  DIN MODULE                   *
 **************************************************/
SC_MODULE(DIN)
{
    // PX_I/O: PORT X=MODULENAME INPUT/OUTPUT
    sc_fifo_in<IMAGE> PDI_I;
    sc_fifo_out<IMAGE> PDI_O;
    
    IMAGE image;
    
    void dIN()
    {
        cout << "STARTED DIN THREAD!" << endl;
        while(1){
            PDI_O.write(PDI_I.read());
        }
    }

    SC_CTOR(DIN):
    PDI_I("PDI_I"),
    PDI_O("PDI_O")
    {
        SC_THREAD(dIN);
        set_stack_size(128*1024*1024);
    };
    
}; // END OF DIN MODULE


/***************************************************
 *  GAUSSIANSMOOTH  (GS) MODULE    *
 **************************************************/

SC_MODULE(GS) {

    IMAGE image;
    SIMAGE smoothedim;
    
    sc_fifo_in<IMAGE> img;
    sc_fifo_out<SIMAGE> smtP;
    
    /*******************************************************************************
    * PROCEDURE: gaussian_smooth
    * PURPOSE: Blur an image with a gaussian filter.
    * NAME: Mike Heath
    * DATE: 2/15/96
    *******************************************************************************/
    void gaussian_smooth(unsigned char *image, int rows, int cols, float sigma,
            short int *smoothedim)
    {
        int r, c, rr, cc,     /* Counter variables. */
              windowsize,        /* Dimension of the gaussian kernel. */
              center;            /* Half of the windowsize. */
           float tempim[SIZE]    /* Buffer for separable filter gaussian smoothing. */
                = {0.0},
                 kernel[WINSIZE] /* A one dimensional gaussian kernel. */
                = {0.0},
                 dot,            /* Dot product summing variable. */
                 sum;            /* Sum of the kernel weights variable. */

       /****************************************************************************
       * Create a 1-dimensional gaussian smoothing kernel.
       ****************************************************************************/
       if(VERBOSE) printf("   Computing the gaussian smoothing kernel.\n");
       make_gaussian_kernel(sigma, kernel, &windowsize);
       center = windowsize / 2;

       /****************************************************************************
       * Blur in the x - direction.
       ****************************************************************************/
       if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
       for(r=0;r<rows;r++){
          for(c=0;c<cols;c++){
             dot = 0.0;
             sum = 0.0;
             for(cc=(-center);cc<=center;cc++){
                if(((c+cc) >= 0) && ((c+cc) < cols)){
                   dot += (float)image[r*cols+(c+cc)] * kernel[center+cc];
                   sum += kernel[center+cc];
                }
             }
             tempim[r*cols+c] = dot/sum;
          }
       }

       /****************************************************************************
       * Blur in the y - direction.
       ****************************************************************************/
       if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
       for(c=0;c<cols;c++){
          for(r=0;r<rows;r++){
             sum = 0.0;
             dot = 0.0;
             for(rr=(-center);rr<=center;rr++){
                if(((r+rr) >= 0) && ((r+rr) < rows)){
                   dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
                   sum += kernel[center+rr];
                }
             }
             smoothedim[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
          }
       }
    }
    
    /*******************************************************************************
    * PROCEDURE: make_gaussian_kernel
    * PURPOSE: Create a one dimensional gaussian kernel.
    * NAME: Mike Heath
    * DATE: 2/15/96
    *******************************************************************************/
    void make_gaussian_kernel(float sigma, float *kernel, int *windowsize)
    {
       int i, center;
       float x, fx, sum=0.0;

       *windowsize = 1 + 2 * ceil(2.5 * sigma);
       center = (*windowsize) / 2;

       if(VERBOSE) printf("      The kernel has %d elements.\n", *windowsize);

       for(i=0;i<(*windowsize);i++){
          x = (float)(i - center);
          fx = pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853));
          kernel[i] = fx;
          sum += fx;
       }

       for(i=0;i<(*windowsize);i++) kernel[i] /= sum;

       if(VERBOSE){
          printf("The filter coefficients are:\n");
          for(i=0;i<(*windowsize);i++)
             printf("kernel[%d] = %f\n", i, kernel[i]);
       }
    }
    
    void GaussianSmooth(void) {
        while (1) {
            img.read(image);

            gaussian_smooth(image, ROWS, COLS, SIGMA,smoothedim);
            
            smtP.write(smoothedim);

            
        }
    }

    SC_CTOR(GS) {
        SC_THREAD(GaussianSmooth);
        set_stack_size(128*1024*1024);
    }
    
};



/***************************************************
 *       DERRIVATIVE MODULE              *
 **************************************************/

SC_MODULE(DERRIVE) {

    SIMAGE imageIn;

    SIMAGE imageX;
    SIMAGE imageY;

    sc_fifo_in<SIMAGE> smoothing; // SIMAGE
    
    sc_fifo_out<SIMAGE> dx1;
    sc_fifo_out<SIMAGE> dy1;
    
    sc_fifo_out<SIMAGE> dx2;
    sc_fifo_out<SIMAGE> dy2;
    
    /*******************************************************************************
    * PROCEDURE: derrivative_x_y
    * PURPOSE: Compute the first derivative of the image in both the x any y
    * directions. The differential filters that are used are:
    *
    *                                          -1
    *         dx =  -1 0 +1     and       dy =  0
    *                                          +1
    *
    * NAME: Mike Heath
    * DATE: 2/15/96
    *******************************************************************************/
    void derrivative_x_y(short int *smoothedim, int rows, int cols,
            short int *delta_x, short int *delta_y)
    {
       int r, c, pos;

       /****************************************************************************
       * Compute the x-derivative. Adjust the derivative at the borders to avoid
       * losing pixels.
       ****************************************************************************/
       if(VERBOSE) printf("   Computing the X-direction derivative.\n");
       for(r=0;r<rows;r++){
          pos = r * cols;
          delta_x[pos] = smoothedim[pos+1] - smoothedim[pos];
          pos++;
          for(c=1;c<(cols-1);c++,pos++){
             delta_x[pos] = smoothedim[pos+1] - smoothedim[pos-1];
          }
          delta_x[pos] = smoothedim[pos] - smoothedim[pos-1];
       }

       /****************************************************************************
       * Compute the y-derivative. Adjust the derivative at the borders to avoid
       * losing pixels.
       ****************************************************************************/
       if(VERBOSE) printf("   Computing the Y-direction derivative.\n");
       for(c=0;c<cols;c++){
          pos = c;
          delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos];
          pos += cols;
          for(r=1;r<(rows-1);r++,pos+=cols){
             delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos-cols];
          }
          delta_y[pos] = smoothedim[pos] - smoothedim[pos-cols];
       }
    }
    
    void Derrive(void) {
        // cout << "****** In DERRIVATIVE_X_Y **********" << endl;
        while (1) {
            smoothing.read(imageIn);

            derrivative_x_y(imageIn, ROWS, COLS, imageX, imageY);

            dx1.write(imageX);
            dy1.write(imageY);
            
            dx2.write(imageX);
            dy2.write(imageY);
        }
    }

    SC_CTOR(DERRIVE) {
        SC_THREAD(Derrive);
        set_stack_size(128*1024*1024);
    }
    
    
};

/***************************************************
 *   MAGNITUDE (MAG)  MODULE         *
 **************************************************/
SC_MODULE(MAG) {
    sc_fifo_in<SIMAGE> dx;
    sc_fifo_in<SIMAGE> dy;

    sc_fifo_out<SIMAGE> mag1;
    sc_fifo_out<SIMAGE> mag2;

    SIMAGE imageX;
    SIMAGE imageY;

    SIMAGE imageMag;
    
    /*******************************************************************************
    * PROCEDURE: magnitude_x_y
    * PURPOSE: Compute the magnitude of the gradient. This is the square root of
    * the sum of the squared derivative values.
    * NAME: Mike Heath
    * DATE: 2/15/96
    *******************************************************************************/
    void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols,
            short int *magnitude)
    {
       int r, c, pos, sq1, sq2;

       for(r=0,pos=0;r<rows;r++){
          for(c=0;c<cols;c++,pos++){
             sq1 = (int)delta_x[pos] * (int)delta_x[pos];
             sq2 = (int)delta_y[pos] * (int)delta_y[pos];
             magnitude[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
          }
       }

    }
    
    void Magnitude(void) {
        while (1) {
            
            dy.read(imageY);
            dx.read(imageX);

            magnitude_x_y(imageX, imageY, ROWS, COLS, imageMag);

            mag1.write(imageMag);
            mag2.write(imageMag);
        }
    }
    
    SC_CTOR(MAG){
        SC_THREAD(Magnitude);
        set_stack_size(128*1024*1024);
    }
    
};


/***************************************************
 *  NON MAX SUPPLY (NMS)  MODULE         *
 **************************************************/
SC_MODULE(NMS) {
    sc_fifo_in<SIMAGE> mag;
    sc_fifo_in<SIMAGE> dx;
    sc_fifo_in<SIMAGE> dy;

    sc_fifo_out<IMAGE> nms;

    SIMAGE imageX;
    SIMAGE imageY;
    SIMAGE imageMag;

    IMAGE imageNMS;
    
    /*******************************************************************************
    * PROCEDURE: non_max_supp
    * PURPOSE: This routine applies non-maximal suppression to the magnitude of
    * the gradient image.
    * NAME: Mike Heath
    * DATE: 2/15/96
    *******************************************************************************/
    void non_max_supp(short *mag, short *gradx, short *grady, int nrows, int ncols,
        unsigned char *result)
    {
        int rowcount, colcount,count;
        short *magrowptr,*magptr;
        short *gxrowptr,*gxptr;
        short *gyrowptr,*gyptr,z1,z2;
        short m00; short gx=0; short gy=0;
        float mag1,mag2; float xperp=0;float yperp=0;
        unsigned char *resultrowptr, *resultptr;

       /****************************************************************************
       * Zero the edges of the result image.
       ****************************************************************************/
        for(count=0,resultrowptr=result,resultptr=result+ncols*(nrows-1);
            count<ncols; resultptr++,resultrowptr++,count++){
            *resultrowptr = *resultptr = (unsigned char) 0;
        }

        for(count=0,resultptr=result,resultrowptr=result+ncols-1;
            count<nrows; count++,resultptr+=ncols,resultrowptr+=ncols){
            *resultptr = *resultrowptr = (unsigned char) 0;
        }

       /****************************************************************************
       * Suppress non-maximum points.
       ****************************************************************************/
       for(rowcount=1,magrowptr=mag+ncols+1,gxrowptr=gradx+ncols+1,
          gyrowptr=grady+ncols+1,resultrowptr=result+ncols+1;
          rowcount<=nrows-2;    // bug fix 3/29/17, RD
          rowcount++,magrowptr+=ncols,gyrowptr+=ncols,gxrowptr+=ncols,
          resultrowptr+=ncols){
          for(colcount=1,magptr=magrowptr,gxptr=gxrowptr,gyptr=gyrowptr,
             resultptr=resultrowptr;colcount<=ncols-2;    // bug fix 3/29/17, RD
             colcount++,magptr++,gxptr++,gyptr++,resultptr++){
             m00 = *magptr;
             if(m00 == 0){
                *resultptr = (unsigned char) NOEDGE;
             }
             else{
                xperp = -(gx = *gxptr)/((float)m00);
                yperp = (gy = *gyptr)/((float)m00);
             }

             if(gx >= 0){
                if(gy >= 0){
                        if (gx >= gy)
                        {
                            /* 111 */
                            /* Left point */
                            z1 = *(magptr - 1);
                            z2 = *(magptr - ncols - 1);

                            mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;

                            /* Right point */
                            z1 = *(magptr + 1);
                            z2 = *(magptr + ncols + 1);

                            mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                        }
                        else
                        {
                            /* 110 */
                            /* Left point */
                            z1 = *(magptr - ncols);
                            z2 = *(magptr - ncols - 1);

                            mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

                            /* Right point */
                            z1 = *(magptr + ncols);
                            z2 = *(magptr + ncols + 1);

                            mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp;
                        }
                    }
                    else
                    {
                        if (gx >= -gy)
                        {
                            /* 101 */
                            /* Left point */
                            z1 = *(magptr - 1);
                            z2 = *(magptr + ncols - 1);

                            mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;

                            /* Right point */
                            z1 = *(magptr + 1);
                            z2 = *(magptr - ncols + 1);

                            mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
                        }
                        else
                        {
                            /* 100 */
                            /* Left point */
                            z1 = *(magptr + ncols);
                            z2 = *(magptr + ncols - 1);

                            mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

                            /* Right point */
                            z1 = *(magptr - ncols);
                            z2 = *(magptr - ncols + 1);

                            mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp;
                        }
                    }
                }
                else
                {
                    if ((gy = *gyptr) >= 0)
                    {
                        if (-gx >= gy)
                        {
                            /* 011 */
                            /* Left point */
                            z1 = *(magptr + 1);
                            z2 = *(magptr - ncols + 1);

                            mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

                            /* Right point */
                            z1 = *(magptr - 1);
                            z2 = *(magptr + ncols - 1);

                            mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
                        }
                        else
                        {
                            /* 010 */
                            /* Left point */
                            z1 = *(magptr - ncols);
                            z2 = *(magptr - ncols + 1);

                            mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

                            /* Right point */
                            z1 = *(magptr + ncols);
                            z2 = *(magptr + ncols - 1);

                            mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
                        }
                    }
                    else
                    {
                        if (-gx > -gy)
                        {
                            /* 001 */
                            /* Left point */
                            z1 = *(magptr + 1);
                            z2 = *(magptr + ncols + 1);

                            mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

                            /* Right point */
                            z1 = *(magptr - 1);
                            z2 = *(magptr - ncols - 1);

                            mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
                        }
                        else
                        {
                            /* 000 */
                            /* Left point */
                            z1 = *(magptr + ncols);
                            z2 = *(magptr + ncols + 1);

                            mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

                            /* Right point */
                            z1 = *(magptr - ncols);
                            z2 = *(magptr - ncols - 1);

                            mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
                        }
                    }
                }

                /* Now determine if the current point is a maximum point */

                if ((mag1 > 0.0) || (mag2 > 0.0))
                {
                    *resultptr = (unsigned char) NOEDGE;
                }
                else
                {
                    if (mag2 == 0.0)
                        *resultptr = (unsigned char) NOEDGE;
                    else
                        *resultptr = (unsigned char) POSSIBLE_EDGE;
                }
            }
        }
    }
    
    void NonMaxSupply(void) {
        while(1) {
            mag.read(imageMag);
            dy.read(imageY);
            dx.read(imageX);

            non_max_supp(imageMag, imageX, imageY, ROWS, COLS, imageNMS);
            nms.write(imageNMS);
        }
    }
    
    SC_CTOR(NMS){
        SC_THREAD(NonMaxSupply);
        set_stack_size(128*1024*1024);
    }
};





/***************************************************
 *   APPLY HYSTERESIS (AH) MODULE         *
 **************************************************/
SC_MODULE(AH) {
    sc_fifo_in<SIMAGE> mag;
    sc_fifo_in<IMAGE> nms;

    sc_fifo_out<IMAGE> edgeP;

    SIMAGE imageMAG;
    IMAGE imageNMS;

    IMAGE edge;
    
    /*******************************************************************************
    * PROCEDURE: follow_edges
    * PURPOSE: This procedure edges is a recursive routine that traces edgs along
    * all paths whose magnitude values remain above some specifyable lower
    * threshhold.
    * NAME: Mike Heath
    * DATE: 2/15/96
    *******************************************************************************/
    void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval,
       int cols)
    {
       short *tempmagptr;
       unsigned char *tempmapptr;
       int i;
       int x[8] = {1,1,0,-1,-1,-1,0,1},
           y[8] = {0,1,1,1,0,-1,-1,-1};

       for(i=0;i<8;i++){
          tempmapptr = edgemapptr - y[i]*cols + x[i];
          tempmagptr = edgemagptr - y[i]*cols + x[i];

          if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)){
             *tempmapptr = (unsigned char) EDGE;
             follow_edges(tempmapptr,tempmagptr, lowval, cols);
          }
       }
    }



    /*******************************************************************************
    * PROCEDURE: apply_hysteresis
    * PURPOSE: This routine finds edges that are above some high threshhold or
    * are connected to a high pixel by a path of pixels greater than a low
    * threshold.
    * NAME: Mike Heath
    * DATE: 2/15/96
    *******************************************************************************/
    void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols,
        float tlow, float thigh, unsigned char *edge)
    {
       int r, c, pos, numedges, highcount, lowthreshold, highthreshold, hist[32768];
       short int maximum_mag=0;

       /****************************************************************************
       * Initialize the edge map to possible edges everywhere the non-maximal
       * suppression suggested there could be an edge except for the border. At
       * the border we say there can not be an edge because it makes the
       * follow_edges algorithm more efficient to not worry about tracking an
       * edge off the side of the image.
       ****************************************************************************/
       for(r=0,pos=0;r<rows;r++){
          for(c=0;c<cols;c++,pos++){
         if(nms[pos] == POSSIBLE_EDGE) edge[pos] = POSSIBLE_EDGE;
         else edge[pos] = NOEDGE;
          }
       }

       for(r=0,pos=0;r<rows;r++,pos+=cols){
          edge[pos] = NOEDGE;
          edge[pos+cols-1] = NOEDGE;
       }
       pos = (rows-1) * cols;
       for(c=0;c<cols;c++,pos++){
          edge[c] = NOEDGE;
          edge[pos] = NOEDGE;
       }

       /****************************************************************************
       * Compute the histogram of the magnitude image. Then use the histogram to
       * compute hysteresis thresholds.
       ****************************************************************************/
       for(r=0;r<32768;r++) hist[r] = 0;
       for(r=0,pos=0;r<rows;r++){
          for(c=0;c<cols;c++,pos++){
         if(edge[pos] == POSSIBLE_EDGE) hist[mag[pos]]++;
          }
       }

       /****************************************************************************
       * Compute the number of pixels that passed the nonmaximal suppression.
       ****************************************************************************/
       for(r=1,numedges=0;r<32768;r++){
          if(hist[r] != 0) maximum_mag = r;
          numedges += hist[r];
       }

       highcount = (int)(numedges * thigh + 0.5);

       /****************************************************************************
       * Compute the high threshold value as the (100 * thigh) percentage point
       * in the magnitude of the gradient histogram of all the pixels that passes
       * non-maximal suppression. Then calculate the low threshold as a fraction
       * of the computed high threshold value. John Canny said in his paper
       * "A Computational Approach to Edge Detection" that "The ratio of the
       * high to low threshold in the implementation is in the range two or three
       * to one." That means that in terms of this implementation, we should
       * choose tlow ~= 0.5 or 0.33333.
       ****************************************************************************/
       r = 1;
       numedges = hist[1];
       while((r<(maximum_mag-1)) && (numedges < highcount)){
          r++;
          numedges += hist[r];
       }
       highthreshold = r;
       lowthreshold = (int)(highthreshold * tlow + 0.5);

       if(VERBOSE){
          printf("The input low and high fractions of %f and %f computed to\n",
         tlow, thigh);
          printf("magnitude of the gradient threshold values of: %d %d\n",
         lowthreshold, highthreshold);
       }

       /****************************************************************************
       * This loop looks for pixels above the highthreshold to locate edges and
       * then calls follow_edges to continue the edge.
       ****************************************************************************/
       for(r=0,pos=0;r<rows;r++){
          for(c=0;c<cols;c++,pos++){
         if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)){
                edge[pos] = EDGE;
                follow_edges((edge+pos), (mag+pos), lowthreshold, cols);
         }
          }
       }

       /****************************************************************************
       * Set all the remaining possible edges to non-edges.
       ****************************************************************************/
       for(r=0,pos=0;r<rows;r++){
          for(c=0;c<cols;c++,pos++) if(edge[pos] != EDGE) edge[pos] = NOEDGE;
       }
    }
    
    void ApplyHystersis(void) {
        while (1) {
            mag.read(imageMAG);
            nms.read(imageNMS);

            apply_hysteresis(imageMAG, imageNMS, ROWS, COLS, TLOW, THIGH, edge);

            edgeP.write(edge);
        }
    }
    SC_CTOR(AH) {
        SC_THREAD(ApplyHystersis);
        set_stack_size(128*1024*1024);
    }
    
};


/***************************************************
 *                  CANNY MODULE                   *
 **************************************************/
SC_MODULE(Canny) // DUT
{
    // PX_I/O: PORT X=MODULENAME INPUT/OUTPUT
    sc_fifo_in<IMAGE> PC_I;
    sc_fifo_out<IMAGE> PC_O;
    
    // CHANNELS
    sc_fifo<SIMAGE> q1;
    sc_fifo<SIMAGE> q2;
    sc_fifo<SIMAGE> q3;
    sc_fifo<SIMAGE> q4;
    sc_fifo<SIMAGE> q5;
    sc_fifo<SIMAGE> q6;
    sc_fifo<SIMAGE> q7;
    
    sc_fifo<IMAGE> q8;
    

    IMAGE image;
    IMAGE edge;
    
    // MODULES
    GS gsM;
    DERRIVE derrivM;
    MAG magM;
    NMS nmsM;
    AH ahM;
    
    void before_end_of_elaboration() {
        
        gsM.img.bind(PC_I);
        gsM.smtP.bind(q1);

        derrivM.smoothing.bind(q1);
        
        derrivM.dx1.bind(q2);
        derrivM.dy1.bind(q3);
        derrivM.dx2.bind(q4);
        derrivM.dy2.bind(q5);
        

        magM.dx.bind(q2);
        magM.dy.bind(q3);
        
        magM.mag1.bind(q6);
        magM.mag2.bind(q7);
        

        nmsM.mag.bind(q6);
        nmsM.dx.bind(q4);
        nmsM.dy.bind(q5);
        nmsM.nms.bind(q8);

        
        ahM.mag.bind(q7);
        ahM.nms.bind(q8);
        
        ahM.edgeP.bind(PC_O);

    }
    
    void DUT_CANNY()
    {
        cout << "STARTED Canny THREAD!" << endl;
        /****************************************************************************
        * Perform the edge detection. All of the work takes place here.
        ****************************************************************************/
        while(1){
            image = PC_I.read();
            
            if(VERBOSE){
                printf("Starting Canny edge detection.\n");
//                canny(image, ROWS, COLS, SIGMA, TLOW, THIGH, edge);
                PC_O.write(image);
                PC_O.write(edge);
            }

        }
    }

    SC_CTOR(Canny):
    PC_I("PC_I"),
    PC_O("PC_O"),
    q1("q1"),
    q2("q2"),
    q3("q3"),
    q4("q4"),
    q5("q5"),
    q6("q6"),
    q7("q7"),
    q8("q8"),
    gsM("gsM"),
    derrivM("derrivM"),
    magM("magM"),
    nmsM("nmsM"),
    ahM("ahM")
    {
        SC_THREAD(DUT_CANNY);
        set_stack_size(128*1024*1024);
    };
    // Canny CONSTRUCTOR

//    /*******************************************************************************
//    * PROCEDURE: canny
//    * PURPOSE: To perform canny edge detection.
//    * NAME: Mike Heath
//    * DATE: 2/15/96
//    *******************************************************************************/
//    void canny(unsigned char *image, int rows, int cols, float sigma,
//             float tlow, float thigh, unsigned char *edge)
//    {
//       unsigned char nms[SIZE]    /* Points that are local maximal magnitude. */
//                = {0};
//       short int smoothedim[SIZE] /* The image after gaussian smoothing.      */
//                = {0},
//                 delta_x[SIZE]    /* The first devivative image, x-direction. */
//                = {0},
//                 delta_y[SIZE]    /* The first derivative image, y-direction. */
//                = {0},
//                 magnitude[SIZE]  /* The magnitude of the gadient image.      */
//                = {0};
//
//       /****************************************************************************
//       * Perform gaussian smoothing on the image using the input standard
//       * deviation.
//       ****************************************************************************/
//       if(VERBOSE) printf("Smoothing the image using a gaussian kernel.\n");
//       gaussian_smooth(image, rows, cols, sigma, smoothedim);
//
//       /****************************************************************************
//       * Compute the first derivative in the x and y directions.
//       ****************************************************************************/
//       if(VERBOSE) printf("Computing the X and Y first derivatives.\n");
//       derrivative_x_y(smoothedim, rows, cols, delta_x, delta_y);
//
//       /****************************************************************************
//       * Compute the magnitude of the gradient.
//       ****************************************************************************/
//       if(VERBOSE) printf("Computing the magnitude of the gradient.\n");
//       magnitude_x_y(delta_x, delta_y, rows, cols, magnitude);
//
//       /****************************************************************************
//       * Perform non-maximal suppression.
//       ****************************************************************************/
//       if(VERBOSE) printf("Doing the non-maximal suppression.\n");
//       non_max_supp(magnitude, delta_x, delta_y, rows, cols, nms);
//
//       /****************************************************************************
//       * Use hysteresis to mark the edge pixels.
//       ****************************************************************************/
//       if(VERBOSE) printf("Doing hysteresis thresholding.\n");
//       apply_hysteresis(magnitude, nms, rows, cols, tlow, thigh, edge);
//    }
    
   
}; // END OF CANNY MODULE


/***************************************************
 *                  DOUT MODULE                   *
 **************************************************/
SC_MODULE(DOUT)
{
    // PX_I/O: PORT X=MODULENAME INPUT/OUTPUT
    sc_fifo_in<IMAGE> PDO_I;
    sc_fifo_out<IMAGE> PDO_O;
    
    
    void dOUT()
    {
        cout << "Started DOUT THREAD!" << endl;
        /****************************************************************************
        * Write out the edge image to a file.
        ****************************************************************************/
        while(1){
            PDO_O.write(PDO_I.read());
        }
        
    }

    SC_CTOR(DOUT):
    PDO_I("PDO_I"),
    PDO_O("PDO_O")
    {
        SC_THREAD(dOUT);
        set_stack_size(128*1024*1024);
    };
}; // END OF DOUT MODULE

/***************************************************
 *                MONITOR MODULE                   *
 **************************************************/
SC_MODULE(Monitor)
{
    // PX_I/O: PORT X=MODULENAME INPUT/OUTPUT
    sc_fifo_in<IMAGE> PM_I;
    
    IMAGE image;
    IMAGE edge;
    
    /* Name of the output "edge" image */
    char outfilename[128];
    int i=0, n=0;
    
    /******************************************************************************
    * Function: write_pgm_image
    * Purpose: This function writes an image in PGM format. The file is either
    * written to the file specified by outfilename or to standard output if
    * outfilename = NULL. A comment can be written to the header if coment != NULL.
    ******************************************************************************/
    int write_pgm_image(const char *outfilename, unsigned char *image, int rows,
        int cols, const char *comment, int maxval)
    {
       FILE *fp;

       /***************************************************************************
       * Open the output image file for writing if a filename was given. If no
       * filename was provided, set fp to write to standard output.
       ***************************************************************************/
       if(outfilename == NULL) fp = stdout;
       else{
          if((fp = fopen(outfilename, "w")) == NULL){
             fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
                outfilename);
             return(0);
          }
       }

       /***************************************************************************
       * Write the header information to the PGM file.
       ***************************************************************************/
       fprintf(fp, "P5\n%d %d\n", cols, rows);
       if(comment != NULL)
          if(strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
       fprintf(fp, "%d\n", maxval);

       /***************************************************************************
       * Write the image data to the file.
       ***************************************************************************/
       if((unsigned)rows != fwrite(image, cols, rows, fp)){
          fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
          if(fp != stdout) fclose(fp);
          return(0);
       }

       if(fp != stdout) fclose(fp);
       return(1);
    }
    
    
    void monitor()
    {
        cout << "IN Monitor THREAD!" << endl;
        
        //for loop for read image
        for (int i = 0; i < IMG_NUM; i++)
        {
            /* code */
            image = PM_I.read();
            edge = PM_I.read();
            
            n = i % AVAIL_IMG;
            sprintf(outfilename, IMG_OUT, n+1);
            if(VERBOSE) printf("Writing the edge iname in the file %s.\n", outfilename);
            if(write_pgm_image(outfilename, edge, ROWS, COLS, "", 255) == 0){
               fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
               exit(1);
            }
        }
        
        sc_stop();
    }
    
    // MONITOR CONSTRUCTOR
    SC_CTOR(Monitor):
    PM_I("PM_I")
    {
        SC_THREAD(monitor);
        set_stack_size(128*1024*1024);
    };
    
}; // END OF MONITOR MODULE


/***************************************************
 *               PLATFORM MODULE                   *
 **************************************************/
SC_MODULE(Platform)
{
    // Instantiate modules as member objects
    DIN din;
    Canny canny;
    DOUT dout;
    
    // COMMUNICATION CHANNELS (FIFO)
    sc_fifo<IMAGE> q3;
    sc_fifo<IMAGE> q4;
    
    sc_fifo_in<IMAGE> PP_I;
    sc_fifo_out<IMAGE> PP_O;
    
    void before_end_of_elaboration()
    {
        // Internal binding of channels to submodule ports
        din.PDI_I.bind(PP_I); // Connect PP_I to PDI_I of DIN
        dout.PDO_O.bind(PP_O); // Connect PP_O to PDO_O of DOUT
        
        // BIND Channel "q3"
        din.PDI_O.bind(q3);
        canny.PC_I.bind(q3);
        
        // BIND Channel "q4"
        canny.PC_O.bind(q4);
        dout.PDO_I.bind(q4);
        
    }
    
    
//    void platform(void)
//    {
//
//
//    }
    
    // MONITOR CONSTRUCTOR
    SC_CTOR(Platform):
    din("din"),
    canny("canny"),
    dout("dout"),
    q3("q3"),
    q4("q4"),
    PP_I("PP_I"), // Initialize external ports
    PP_O("PP_O")
    {
//        SC_THREAD(platform);
//        set_stack_size(128*1024*1024);
    };
    
}; // END OF PLATFORM MODULE

/***************************************************
 *                   TOP MODULE                    *
 **************************************************/
SC_MODULE(Top)
{
    // Declare Components
    // stim, monitor, platform
    
    Stimulus stimulus;
    Monitor monitor;
    Platform platform;
    
    // COMMUNICATION CHANNELS (FIFO)
    sc_fifo<IMAGE> q1;
    sc_fifo<IMAGE> q2;
    
    void before_end_of_elaboration()
    {
        // BIND Channel "q1"
        stimulus.PS_O.bind(q1);
        platform.PP_I.bind(q1);

        // BIND Channel "q2"
        platform.PP_O.bind(q2);
        monitor.PM_I.bind(q2);
    }
    
    
    // Top Constructor
    SC_CTOR(Top):
    stimulus("stimulus"),
    platform("platform"),
    monitor("monitor"),
    q1("q1"),
    q2("q2")
    {
//        SC_THREAD(top);
//        set_stack_size(128*1024*1024);
    }
    
};

// END OF TOP MODULE

/***************************************************
 *                    SC_MAIN                      *
 **************************************************/
int sc_main(int argc, char* argv[])
{
    Top top("top"); // Instantiate one top-level module
    
    sc_start(); // End elaboration, run simulation
    
    return 0;
}

