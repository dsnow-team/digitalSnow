
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */





// main/pbrt.cpp*
#include "stdafx.h"
#include "api.h"
#include "probes.h"
#include "parser.h"
#include "parallel.h"
#include <cstdlib>


// main program
int main(int argc, char *argv[]) {
    Options options;
    vector<string> filenames;


	//[DGtal on met que un coeur pour avoir tous les calculs séquentiels et les résultats dans un seul fichier
	bool wavelength(false), resPix(false), dimensionX(false),dimensionY(false),dimensionZ(false), ImagePhoton(false);
		

    // Process command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--ncores")) options.nCores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--outfile")) options.imageFile = argv[++i];
        else if (!strcmp(argv[i], "--quick")) options.quickRender = true;
        else if (!strcmp(argv[i], "--quiet")) options.quiet = true;
        else if (!strcmp(argv[i], "--verbose")) options.verbose = true;
        else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
            printf("usage: pbrt  [--image || -i ] file.pbrt \n"
                   "pbrt [--photon || -p] [--wavelength wavelength(nm) || -w wavelength(nm)] [-x dimImageY] [-y dimImageY] [-z dimImageZ] [--resPixel PixelResolution(micrometer) || -r PixelResolution(micrometer)] [ <filenamePhoton.pbrt> ...\n");
           return 0;
        }
	//[DGtal ajout option pour faire de l'absorption]
	else if (!strcmp(argv[i],"--wavelength") || !strcmp(argv[i],"-w")) { options.lOnde=atof(argv[++i]); wavelength=true;}
	else if (!strcmp(argv[i],"-x")) { options.dimx=atoi(argv[++i]); dimensionX=true; }
	else if (!strcmp(argv[i],"-y")) { options.dimy=atoi(argv[++i]); dimensionY=true; }
	else if (!strcmp(argv[i],"-z")) { options.dimz=atoi(argv[++i]); dimensionZ=true; }
	else if ((!strcmp(argv[i],"--resPixel")) || (!strcmp(argv[i],"-r"))) { options.resolPixel=atof(argv[++i]); resPix=true; }
	else if ((!strcmp(argv[i],"--photon")) || (!strcmp(argv[i],"-p"))){ImagePhoton=true; options.photon=true;options.nCores=1;}	
	else if ((!strcmp(argv[i],"--image")) || (!strcmp(argv[i],"-i"))) { ImagePhoton=true; options.photon=false;}

        else {
		filenames.push_back(argv[i]);
		options.filename=argv[i];
		}
    }

	//[DGtal : test arguments]
	if (!ImagePhoton) {printf("usage: pbrt  [--image || -i ] file.pbrt \n"
                   "pbrt [--photon || -p] [--wavelength wavelength(nm) || -w wavelength(nm)] [-x dimImageY] [-y dimImageY] [-z dimImageZ] [--resPixel PixelResolution(micrometer) || -r PixelResolution(micrometer)] [ <filenamePhoton.pbrt> ...\n"); exit(1);}
	else if (options.photon && ((!wavelength) || (!dimensionX) || (!dimensionY) || (!dimensionZ) || (!resPix)))
	{
            printf("usage: pbrt  [--image || -i ] file.pbrt \n"
                   "pbrt [--photon || -p] [--wavelength wavelength(nm) || -w wavelength(nm)] [-x dimImageY] [-y dimImageY] [-z dimImageZ] [--resPixel PixelResolution(micrometer) || -r PixelResolution(micrometer)] [ <filenamePhoton.pbrt> ...\n");
	exit(1);
	}


    // Print welcome banner
    if (!options.quiet) {
        printf("pbrt version %s of %s at %s [Detected %d core(s)]\n",
               PBRT_VERSION, __DATE__, __TIME__, NumSystemCores());
        printf("Copyright (c)1998-2010 Matt Pharr and Greg Humphreys.\n");
        printf("The source code to pbrt (but *not* the book contents) is covered by the GNU GPL.\n");
        printf("See the file COPYING.txt for the conditions of the license.\n");
        fflush(stdout);
    }
    pbrtInit(options);
    // Process scene description
    PBRT_STARTED_PARSING();
    if (filenames.size() == 0) {
        // Parse scene from standard input
        ParseFile("-");
    } else {
        // Parse scene from input files
        for (u_int i = 0; i < filenames.size(); i++)
            if (!ParseFile(filenames[i]))
                Error("Couldn't open scene file \"%s\"", filenames[i].c_str());
    }
    pbrtCleanup();
    return 0;
}


