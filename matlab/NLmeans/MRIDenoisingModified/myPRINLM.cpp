
#include <stdio.h>

#include <nifti1_io.h>


#include "math.h"
//#include "mex.h"
#include <stdlib.h>
//#include "matrix.h"
#include <pthread.h>

typedef struct {
    int rows;
    int cols;
    int slices;
    float * in_image;
    float * out_image;
    float * means_image;
    float * ref_image;
    float * pesos;
    float * coilsens_image;
    int ini;
    int fin;
    int radiusB;
    int radiusS;
    float sigma;
    int order;
} myargument;

bool rician;




void* ThreadFunc(void* pArguments)
{
    float *ima, *fima, *medias, *pesos, *ref, *coilsens, sigma, w, d, hhh, hh, t1, alfa, x, beta;
    int ii, jj, kk, ni, nj, nk, i, j, k, ini, fin, rows, cols, slices, v, p, p1, f, order, rc;
    extern bool rician;
    /*extern float *table;*/

    myargument arg;
    arg = *(myargument *) pArguments;

    rows = arg.rows;
    cols = arg.cols;
    slices = arg.slices;
    ini = arg.ini;
    fin = arg.fin;
    ima = arg.in_image;
    fima = arg.out_image;
    medias = arg.means_image;
    ref = arg.ref_image;
    coilsens = arg.coilsens_image;
    pesos = arg.pesos;
    v = arg.radiusB;
    f = arg.radiusS;
    sigma = arg.sigma;
    order = arg.order;

    hh = 2 * sigma * sigma;
    alfa = 0.5;
    hhh = 0.5 * sigma * sigma;
    rc = rows * cols;

    /* filter*/

    if (order == 1)
    {
        for (k = ini; k < fin; k++)
        {
            for (j = 0; j < rows; j++)
            {
                for (i = 0; i < cols; i++)
                {
                    p = k * rc + j * cols + i;

                    if (ima[p] == 0) continue;

                    for (kk = 0; kk <= v; kk++)
                    {
                        nk = k + kk;
                        for (jj = -v; jj <= v; jj++)
                        {
                            nj = j + jj;
                            for (ii = -v; ii <= v; ii++)
                            {
                                ni = i + ii;
                                if (kk == 0 && jj < 0) continue;
                                if (kk == 0 && jj == 0 && ii <= 0) continue;

                                if (ni >= 0 && nj >= 0 && nk >= 0 && ni < cols && nj < rows && nk < slices)
                                {
                                    p1 = nk * rc + nj * cols + ni;

                                    t1 = abs(medias[p] - medias[p1]);

                                    if (t1 > sigma) continue;

                                    w = (ref[p] - ref[p1]);

                                    w = (w * w) / hh - 1;
                                    d = (t1 * t1) / hhh - 1;

                                    w = w < 0 ? 0 : w;
                                    d = d < 0 ? 0 : d;

                                    d = w * 0.5 + d * 0.5;

                                    //Coil sensitivity (B1) correction
                                    beta = (coilsens[p] - coilsens[p1]);
                                    d *= (beta * beta);


                                    if (d <= 0) w = 1.0;
                                    else if (d > 10)w = 0;
                                    else w = exp(-d);

                                    if (rician > 0)
                                    {
                                        fima[p] += w * ima[p1] * ima[p1];
                                        pesos[p] += w;

                                        fima[p1] += w * ima[p] * ima[p];
                                        pesos[p1] += w;
                                    }
                                    else
                                    {
                                        fima[p] += w * ima[p1];
                                        pesos[p] += w;

                                        fima[p1] += w * ima[p];
                                        pesos[p1] += w;
                                    }

                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        for (k = ini; k < fin; k++)
        {
            for (j = rows - 1; j >= 0; j--)
            {
                for (i = cols - 1; i >= 0; i--)
                {
                    p = k * rc + j * cols + i;
                    if (ima[p] == 0) continue;

                    for (kk = 0; kk <= v; kk++)
                    {
                        nk = k + kk;
                        for (jj = -v; jj <= v; jj++)
                        {
                            nj = j + jj;
                            for (ii = -v; ii <= v; ii++)
                            {
                                ni = i + ii;
                                if (kk == 0 && jj < 0) continue;
                                if (kk == 0 && jj == 0 && ii <= 0) continue;

                                if (ni >= 0 && nj >= 0 && nk >= 0 && ni < cols && nj < rows && nk < slices)
                                {
                                    p1 = nk * rc + nj * cols + ni;

                                    t1 = abs(medias[p] - medias[p1]);

                                    if (t1 > sigma) continue;

                                    w = (ref[p] - ref[p1]);

                                    w = (w * w) / hh - 1;
                                    d = (t1 * t1) / hhh - 1;

                                    w = w < 0 ? 0 : w;
                                    d = d < 0 ? 0 : d;

                                    d = w * 0.5 + d * 0.5;
                                    //Coil sensitivity (B1) correction
                                    beta = (coilsens[p] - coilsens[p1]);
                                    d *= (beta * beta);

                                    if (d <= 0) w = 1.0;
                                    else if (d > 10) w = 0;
                                    else w = exp(-d);

                                    if (rician > 0)
                                    {
                                        fima[p] += w * ima[p1] * ima[p1];
                                        pesos[p] += w;

                                        fima[p1] += w * ima[p] * ima[p];
                                        pesos[p1] += w;
                                    }
                                    else
                                    {
                                        fima[p] += w * ima[p1];
                                        pesos[p] += w;

                                        fima[p1] += w * ima[p];
                                        pesos[p1] += w;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    pthread_exit(0);


    return 0;
}



int show_help(void)
{
    printf(
        "myRINLM3D: non-local means with B1 correction\n"
        "\n"
        "    This program is to demonstrate how to read a NIfTI dataset,\n"
        "    set output filenames and write a NIfTI dataset, all via the\n"
        "    standard NIfTI C library.\n"
        "\n"
        "    basic usage: myRINLM3D -input InputImage -output OutputImage -means ODCTImage -coil CoilSensitivityImage -radiusB v -radiusS f -sigma sigma -rician TRUE\n"
        "\n"
        "    options:     -help           : show this help\n"
        "                 -input <filename>  : Input image\n"
        "                 -output <filename> : Output image\n"
        "                 -means <filename>  : ODCT means image\n"
        "                 -coil <filename>   : Coils Sensitivity (B1 correction)\n"
        "                 -radiusB value     : Set batch radius\n"
        "                 -radiusS value     : Set search radius\n"
        "                 -sigma value     : Set sigma value\n"
        "                 -rician bool     : Enable rician (default=true)\n"
        "                 -v      : Verbose output\n"

        "\n");
    return 0;
}

int main(int argc, char * argv[])
{
    nifti_image * nimg = NULL, * mimg = NULL, * coilimg = NULL;
    char        * inputfilename = NULL, * outputfilename = NULL, *meansfilename = NULL, *coilfilename = NULL;
    int           ac;
    float *ima, *fima, *lista, *pesos, *ref, *medias, *coilsens;
    float sigma, vr, R, h, w, average, totalweight, wmax, d, media, var, th, hh, hhh, wm, t1, t2, alfa;
    int inc, i, j, k, ii, jj, kk, ni, nj, nk, radiusB, radiusS, ndim, indice, Nthreads, ini, fin, order;
    const int  *dims;
    myargument *ThreadArgs;
    pthread_t * ThreadList;


    if (argc < 2) return show_help();    /* typing '-help' is sooo much work */

    /* process user options: 4 are valid presently */
    for (ac = 1; ac < argc; ac++) {
        if (! strncmp(argv[ac], "-h", 2)) {
            return show_help();
        }
        else if (! strcmp(argv[ac], "-input")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -input\n");
                return 1;
            }
            inputfilename = argv[ac];  /* no string copy, just pointer assignment */
        }
        else if (! strcmp(argv[ac], "-output")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -output\n");
                return 2;
            }
            outputfilename = argv[ac];
        }
        else if (! strcmp(argv[ac], "-means")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -means\n");
                return 2;
            }
            meansfilename = argv[ac];
        }
        else if (! strcmp(argv[ac], "-coil")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -coil\n");
                return 2;
            }
            coilfilename = argv[ac];
        }
        else if (! strcmp(argv[ac], "-sigma")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -sigma\n");
                return 2;
            }
            sigma = atof(argv[ac]);
        }
        else if (! strcmp(argv[ac], "-radiusB")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -radiusB\n");
                return 2;
            }
            radiusB = atof(argv[ac]);
        }
        else if (! strcmp(argv[ac], "-radiusS")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -radiusS\n");
                return 2;
            }
            radiusS = atof(argv[ac]);
        }
        else if (! strcmp(argv[ac], "-rician")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -rician\n");
                return 2;
            }
            rician = atoi(argv[ac]);
        }
        else if (! strcmp(argv[ac], "-verb")) {
            if (++ac >= argc) {
                fprintf(stderr, "** missing argument for -verb\n");
                return 2;
            }
            nifti_set_debug_level(atoi(argv[ac]));
        }
        else {
            fprintf(stderr, "** invalid option, '%s'\n", argv[ac]);
            return 1;
        }
    }
    h = sigma / 3.0f; /* this is due to the removal of part of the noise*/

    if (!inputfilename) {
        fprintf(stderr, "** missing option '-input'\n");
        return 1;
    }
    if (!outputfilename) {
        fprintf(stderr, "** missing option '-output'\n");
        return 1;
    }
    if (!meansfilename) {
        fprintf(stderr, "** missing option '-means'\n");
        return 1;
    }
    if (!coilfilename) {
        coilfilename = meansfilename;
    }

    /* read input dataset, including data */
    nimg = nifti_image_read(inputfilename, 1);
    if (!nimg) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", inputfilename);
        return 2;
    }
    mimg = nifti_image_read(meansfilename, 1);
    if (!mimg) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", meansfilename);
        return 2;
    }
    coilimg = nifti_image_read(coilfilename, 1);
    if (!coilimg) {
        fprintf(stderr, "** failed to read NIfTI image from '%s'\n", coilfilename);
        return 2;
    }

    ndim = nimg->dim[0];
    dims = &(nimg->dim[1]);
    int npixels = dims[0];
    for (int i = 1; i < ndim; i++) npixels *= dims[i];

    fima = (float*) calloc(npixels, sizeof(float)); //(float*) nimg->data;
    medias = (float*) calloc(npixels, sizeof(float)); //(float*) mimg>data;
    coilsens = (float*) calloc(npixels, sizeof(float)); //(float*) coilimg->data;
    pesos = (float*) calloc(npixels, sizeof(float));
    for (int i = 1; i < ndim; i++) {
        fima[i] = ((float) nimg->data[i]) * (nimg->scl_slope) + nimg->scl_inter;
        medias[i] = ((float)mimg->data[i]) * (mimg->scl_slope) + mimg->scl_inter;
        coilsens[i] = ((float)coilimg->data[i]) * (coilimg->scl_slope) + coilimg->scl_inter;
    }

    /* calculate means*/

    for (k = 0; k < dims[2]; k++)
    {
        for (j = 0; j < dims[1]; j++)
        {
            for (i = 0; i < dims[0]; i++)
            {
                media = ref[k * (dims[0] * dims[1]) + (j * dims[0]) + i];
                for (ii = -1; ii <= 1; ii++)
                {
                    ni = i + ii;
                    if (ni < 0) ni = -ni;
                    if (ni >= dims[0]) ni = 2 * dims[0] - ni - 1;
                    for (jj = -1; jj <= 1; jj++)
                    {
                        nj = j + jj;
                        if (nj < 0) nj = -nj;
                        if (nj >= dims[1]) nj = 2 * dims[1] - nj - 1;
                        for (kk = -1; kk <= 1; kk++)
                        {
                            nk = k + kk;
                            if (nk < 0) nk = -nk;
                            if (nk >= dims[2]) nk = 2 * dims[2] - nk - 1;

                            if (sqrt(ii * ii + jj * jj + kk * kk) > 1)continue;

                            media = media + ref[nk * (dims[0] * dims[1]) + (nj * dims[0]) + ni];
                        }
                    }
                }
                media = media / 8;
                medias[k * (dims[0]*dims[1]) + (j * dims[0]) + i] = media;
            }
        }
    }

    for (k = 0; k < dims[2]*dims[1]*dims[0]; k++)
    {
        if (rician > 0) fima[k] = ima[k] * ima[k];
        else fima[k] = ima[k];
        pesos[k] = 1.0;
    }

    Nthreads = dims[2] < 8 ? dims[2] : 8;
    if (Nthreads < 1) Nthreads = 1;


    /* Reserve room for handles of threads in ThreadList*/
    ThreadList = (pthread_t *) calloc(Nthreads, sizeof(pthread_t));
    ThreadArgs = (myargument*) calloc(Nthreads, sizeof(myargument));

    order = -1;
    for (i = 0; i < Nthreads; i++)
    {
        /* Make Thread Structure*/
        order = -order;
        ini = (i * dims[2]) / Nthreads;
        fin = ((i + 1) * dims[2]) / Nthreads;
        ThreadArgs[i].cols = dims[0];
        ThreadArgs[i].rows = dims[1];
        ThreadArgs[i].slices = dims[2];
        ThreadArgs[i].in_image = ima;
        ThreadArgs[i].out_image = fima;
        ThreadArgs[i].ref_image = ref;
        ThreadArgs[i].means_image = medias;
        ThreadArgs[i].pesos = pesos;
        ThreadArgs[i].ini = ini;
        ThreadArgs[i].fin = fin;
        ThreadArgs[i].radiusB = radiusB;
        ThreadArgs[i].radiusS = radiusS;
        ThreadArgs[i].sigma = h;
        ThreadArgs[i].order = order;
    }

    for (i = 0; i < Nthreads; i++)
    {
        if (pthread_create(&ThreadList[i], NULL, ThreadFunc, &ThreadArgs[i]))
        {
            printf("Threads cannot be created\n");
            exit(1);
        }
    }

    for (i = 0; i < Nthreads; i++)
    {
        pthread_join(ThreadList[i], NULL);
    }

    free(ThreadArgs);
    free(ThreadList);


    sigma = 2 * sigma * sigma;
    for (k = 0; k < dims[2]*dims[1]*dims[0]; k++)
    {
        if (rician > 0)
        {
            vr = (fima[k] / pesos[k]) - sigma;
            if (vr < 0) fima[k] = 0;
            else fima[k] = sqrt(vr);
            if (ima[k] <= 0) fima[k] = 0;
        }
        else fima[k] /= pesos[k];
        if (fima[k] < fima_min) fima_min = fima[k];
        if (fima[k] > fima_max) fima_max = fima[k];
    }

    free(pesos);
    nimg->scl_slope = fima_max - fima_min;
    nimg->scl_inter = fima_min;

    for (int i = 1; i < npixels; i++)
        nimg->data[i] = (signed int)((fima[i] - nimg->scl_inter) / nimg->scl_slope);


    /* assign nifti_image fname/iname pair, based on output filename
       (request to 'check' image and 'set_byte_order' here) */
    if (nifti_set_filenames(nimg, outputfilename, 1, 1)) return 1;

    /* if we get here, write the output dataset */
    nifti_image_write(nimg);

    /* and clean up memory */
    nifti_image_free(nimg); nifti_image_free(mimg); nifti_image_free(coilimg);

    return 0;
}
