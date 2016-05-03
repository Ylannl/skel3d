// Copyright (c) 2012, 2013
// Ravi Peters -- r.y.peters@tudelft.nl
// All rights reserved
// 
// This file is part of Surfonoi.
// 
// Surfonoi is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Surfonoi is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Surfonoi.  If not, see <http://www.gnu.org/licenses/>.

#include "CgalProcessor.h"
#include <tclap/CmdLine.h>

// cnpy
#include <cnpy/cnpy.h>

using namespace std;

int main(int argc, const char * argv[])
{
    try {
        
        TCLAP::CmdLine cmd("Triangulate points, generate obj", ' ', "none", false);
        
        // TCLAP::ValueArg<int> sArg("s","smooth","Number of times to smooth",true,2,"int", cmd);
        // TCLAP::SwitchArg uSwitch("u","unsafe","Smooth without attempting to respect bathymetric safety constraint", cmd, false);
        
        
        TCLAP::UnlabeledValueArg<std::string> inputArg( "input", "path to .xyz file", true, "", "input file", cmd);
        TCLAP::UnlabeledValueArg<std::string> outputArg( "ouput", "path to output", true, "", "output file", cmd);
        cmd.parse(argc,argv);

        std::string input_coords_path = inputArg.getValue()+"/coords.npy";
        std::string input_mask_path = inputArg.getValue()+"/decimate_lfs.npy";
        std::string output_error_path = inputArg.getValue()+"/decimate_error.npy";

        cnpy::NpyArray coords_npy = cnpy::npy_load( input_coords_path.c_str() );
        float* coords_carray = reinterpret_cast<float*>(coords_npy.data);

        cnpy::NpyArray mask_npy = cnpy::npy_load( input_mask_path.c_str() );
        bool* mask_carray = reinterpret_cast<bool*>(mask_npy.data);

        int m = coords_npy.shape[0];

        int count=0;
        for (int i=0; i<m; i++) if (mask_carray[i]) count++;

        float *coords_filterd = new float[count*3];

        int cntr = 0;
        for (int i=0; i<m; i++)
            if (mask_carray[i]){
                for (int j=0; j<3; j++)
                    coords_filterd[3*cntr+j] = coords_carray[i*3+j];
                cntr++;
            }

        CgalProcessor cp(coords_filterd, count);

        float *out = new float[m*3];
        cp.metricL2potri(coords_carray, m, out);
        
        {
          const unsigned int shape[] = { m,3 };
          cnpy::npy_save(output_error_path.c_str(), out, shape, 2, "w");
        }

        cp.dumpOBJ(outputArg.getValue().c_str());
        delete[] coords_filterd; coords_filterd = NULL;
        delete[] out; out = NULL;
        coords_npy.destruct();
        mask_npy.destruct();
        
    } catch (TCLAP::ArgException &e)  // catch any exceptions
	{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
    
    return 0;
}

