///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):      Version of a Creative Commons License,
//                     for details contact Boris Jeremic, jeremic@ucdavis.edu
// PROJECT:            Real ESSI Simulator
// PROGRAMMER:         Han Yang
// DATE:               Mar 2016
// UPDATE HISTORY:     Under test, no update has been made yet
//
// LEGACY/DEFUNCT COPYLEFT (C):
//                     Woody's viral GPL-like license (adapted by BJ):
//                     ``This    source  code is Copyrighted in
//                     worldwide for  an  indefinite  period,  and anybody
//                     caught  using it without our permission, will be
//                     mighty good friends of ourn, cause we don't give
//                     a  darn.  Hack it. Compile it. Debug it. Run it.
//                     Yodel  it.  Enjoy it. We wrote it, that's all we
//                     wanted to do.''
//
/////////////////////////////////////////////////////////////////////////////



// Purpose: This file calculate the energy density from HDF5 output file and write the results back into the same file.


#include <hdf5.h>
#include "H5Cpp.h"
#include <string>
#include <math.h>
#include <map>
#include <iostream>
#include <fstream>
#include "Energy_Post_Processing.h"

using namespace std;

Energy_Post_Processing::Energy_Post_Processing()
{
    is_initialized = false;
}

Energy_Post_Processing::~Energy_Post_Processing()
{
    clean_up();
}

void Energy_Post_Processing::clean_up()
{
    Energy_Out << "Trying to clean everything up!!\n";

    H5Fclose(id_file);
    //H5Fclose(id_file_energy);

    Energy_Out << "Files are successfully closed!!\n";

    delete[] index_kinetic_energy_density;
    delete[] kinetic_energy_density;

    delete[] index_strain_energy_density;
    delete[] strain_energy_density_step;
    delete[] strain_energy_density_cumu;

    delete[] plastic_work_density_step;
    delete[] plastic_work_density_cumu;

    delete[] plastic_free_energy_density_step;
    delete[] plastic_free_energy_density_cumu;

    delete[] plastic_dissipation_density_step;
    delete[] plastic_dissipation_density_cumu;

    delete[] kinetic_energy_density_eleavg;
    delete[] strain_energy_density_eleavg;
    delete[] plastic_work_density_eleavg;
    delete[] mechanical_energy_density_eleavg;
    delete[] index_energy_density_eleavg;

    delete[] element_volume;
    delete[] total_mechanical_energy;
    delete[] total_plastic_work;
    delete[] total_kinetic_energy;
    delete[] total_strain_energy;

    delete[] weight_27nodebrick;

    Energy_Out << "Everything is successfully cleaned up!\n";
}

void Energy_Post_Processing::initialize()
{
    //===========================================================================
    // Define material properties
    //===========================================================================

    mass_density = 2000;

    is_linear_KH = true;
    a1 = 1e8;                   // Material constant for linear kinematic hardening

//    AF_ha = 1e4/20000.0;
//    a1 = AF_ha*2.0/3.0;                   // Material constant for AF kinematic hardening

//    is_AF_KH = false;
//    AF_ha = 10000.0;
//    AF_cr = 1e9 / 1470000.0;



    //===========================================================================
    // Define finite element constants
    //===========================================================================

    // Weight of 27NodeBrick at Gauss points
    weight_27nodebrick = new double[27];

    weight_27nodebrick[0] = 125.0/729;
    weight_27nodebrick[2] = 125.0/729;
    weight_27nodebrick[6] = 125.0/729;
    weight_27nodebrick[8] = 125.0/729;
    weight_27nodebrick[18] = 125.0/729;
    weight_27nodebrick[20] = 125.0/729;
    weight_27nodebrick[24] = 125.0/729;
    weight_27nodebrick[26] = 125.0/729;

    weight_27nodebrick[1] = 200.0/729;
    weight_27nodebrick[3] = 200.0/729;
    weight_27nodebrick[5] = 200.0/729;
    weight_27nodebrick[7] = 200.0/729;
    weight_27nodebrick[9] = 200.0/729;
    weight_27nodebrick[11] = 200.0/729;
    weight_27nodebrick[15] = 200.0/729;
    weight_27nodebrick[17] = 200.0/729;
    weight_27nodebrick[19] = 200.0/729;
    weight_27nodebrick[21] = 200.0/729;
    weight_27nodebrick[23] = 200.0/729;
    weight_27nodebrick[25] = 200.0/729;

    weight_27nodebrick[4] = 320.0/729;
    weight_27nodebrick[10] = 320.0/729;
    weight_27nodebrick[12] = 320.0/729;
    weight_27nodebrick[14] = 320.0/729;
    weight_27nodebrick[16] = 320.0/729;
    weight_27nodebrick[22] = 320.0/729;

    weight_27nodebrick[13] = 512.0/729;

    //Energy_Out << "weight_27nodebrick: " << weight_27nodebrick[13] << endl;


    // Local coordinates of 27NodeBrick at Gauss points
    GP_local_coordinate_27nodebrick = new double*[27];
    for (int i = 0; i < 27; i++)
    {
        GP_local_coordinate_27nodebrick[i] = new double[3];
    }

    GP_local_coordinate_27nodebrick[0][0] = 0.77459666924;   GP_local_coordinate_27nodebrick[0][1] = 0.77459666924;   GP_local_coordinate_27nodebrick[0][2] = 0.77459666924;
    GP_local_coordinate_27nodebrick[1][0] = 0.77459666924;   GP_local_coordinate_27nodebrick[1][1] = 0.77459666924;   GP_local_coordinate_27nodebrick[1][2] = 0.;
    GP_local_coordinate_27nodebrick[2][0] = 0.77459666924;   GP_local_coordinate_27nodebrick[2][1] = 0.77459666924;   GP_local_coordinate_27nodebrick[2][2] = -0.77459666924;
    GP_local_coordinate_27nodebrick[3][0] = 0.77459666924;   GP_local_coordinate_27nodebrick[3][1] = 0.;              GP_local_coordinate_27nodebrick[3][2] = 0.77459666924;
    GP_local_coordinate_27nodebrick[4][0] = 0.77459666924;   GP_local_coordinate_27nodebrick[4][1] = 0.;              GP_local_coordinate_27nodebrick[4][2] = 0.;
    GP_local_coordinate_27nodebrick[5][0] = 0.77459666924;   GP_local_coordinate_27nodebrick[5][1] = 0.;              GP_local_coordinate_27nodebrick[5][2] = -0.77459666924;
    GP_local_coordinate_27nodebrick[6][0] = 0.77459666924;   GP_local_coordinate_27nodebrick[6][1] = -0.77459666924;  GP_local_coordinate_27nodebrick[6][2] = 0.77459666924;
    GP_local_coordinate_27nodebrick[7][0] = 0.77459666924;   GP_local_coordinate_27nodebrick[7][1] = -0.77459666924;  GP_local_coordinate_27nodebrick[7][2] = 0.;
    GP_local_coordinate_27nodebrick[8][0] = 0.77459666924;   GP_local_coordinate_27nodebrick[8][1] = -0.77459666924;  GP_local_coordinate_27nodebrick[8][2] = -0.77459666924;
    GP_local_coordinate_27nodebrick[9][0] = 0.;              GP_local_coordinate_27nodebrick[9][1] = 0.77459666924;   GP_local_coordinate_27nodebrick[9][2] = 0.77459666924;
    GP_local_coordinate_27nodebrick[10][0] = 0.;             GP_local_coordinate_27nodebrick[10][1] = 0.77459666924;  GP_local_coordinate_27nodebrick[10][2] = 0.;
    GP_local_coordinate_27nodebrick[11][0] = 0.;             GP_local_coordinate_27nodebrick[11][1] = 0.77459666924;  GP_local_coordinate_27nodebrick[11][2] = -0.77459666924;
    GP_local_coordinate_27nodebrick[12][0] = 0.;             GP_local_coordinate_27nodebrick[12][1] = 0.;             GP_local_coordinate_27nodebrick[12][2] = 0.77459666924;
    GP_local_coordinate_27nodebrick[13][0] = 0.;             GP_local_coordinate_27nodebrick[13][1] = 0.;             GP_local_coordinate_27nodebrick[13][2] = 0.;
    GP_local_coordinate_27nodebrick[14][0] = 0.;             GP_local_coordinate_27nodebrick[14][1] = 0.;             GP_local_coordinate_27nodebrick[14][2] = -0.77459666924;
    GP_local_coordinate_27nodebrick[15][0] = 0.;             GP_local_coordinate_27nodebrick[15][1] = -0.77459666924; GP_local_coordinate_27nodebrick[15][2] = 0.77459666924;
    GP_local_coordinate_27nodebrick[16][0] = 0.;             GP_local_coordinate_27nodebrick[16][1] = -0.77459666924; GP_local_coordinate_27nodebrick[16][2] = 0.;
    GP_local_coordinate_27nodebrick[17][0] = 0.;             GP_local_coordinate_27nodebrick[17][1] = -0.77459666924; GP_local_coordinate_27nodebrick[17][2] = -0.77459666924;
    GP_local_coordinate_27nodebrick[18][0] = -0.77459666924; GP_local_coordinate_27nodebrick[18][1] = 0.77459666924;  GP_local_coordinate_27nodebrick[18][2] = 0.77459666924;
    GP_local_coordinate_27nodebrick[19][0] = -0.77459666924; GP_local_coordinate_27nodebrick[19][1] = 0.77459666924;  GP_local_coordinate_27nodebrick[19][2] = 0.;
    GP_local_coordinate_27nodebrick[20][0] = -0.77459666924; GP_local_coordinate_27nodebrick[20][1] = 0.77459666924;  GP_local_coordinate_27nodebrick[20][2] = -0.77459666924;
    GP_local_coordinate_27nodebrick[21][0] = -0.77459666924; GP_local_coordinate_27nodebrick[21][1] = 0.;             GP_local_coordinate_27nodebrick[21][2] = 0.77459666924;
    GP_local_coordinate_27nodebrick[22][0] = -0.77459666924; GP_local_coordinate_27nodebrick[22][1] = 0.;             GP_local_coordinate_27nodebrick[22][2] = 0.;
    GP_local_coordinate_27nodebrick[23][0] = -0.77459666924; GP_local_coordinate_27nodebrick[23][1] = 0.;             GP_local_coordinate_27nodebrick[23][2] = -0.77459666924;
    GP_local_coordinate_27nodebrick[24][0] = -0.77459666924; GP_local_coordinate_27nodebrick[24][1] = -0.77459666924; GP_local_coordinate_27nodebrick[24][2] = 0.77459666924;
    GP_local_coordinate_27nodebrick[25][0] = -0.77459666924; GP_local_coordinate_27nodebrick[25][1] = -0.77459666924; GP_local_coordinate_27nodebrick[25][2] = 0.;
    GP_local_coordinate_27nodebrick[26][0] = -0.77459666924; GP_local_coordinate_27nodebrick[26][1] = -0.77459666924; GP_local_coordinate_27nodebrick[26][2] = -0.77459666924;


    // Local coordinates of 27NodeBrick at nodes
    node_local_coordinate_27nodebrick = new int*[27];
    for (int i = 0; i < 27; i++)
    {
        node_local_coordinate_27nodebrick[i] = new int[3];
    }

    node_local_coordinate_27nodebrick[0][0] = 1;    node_local_coordinate_27nodebrick[0][1] = 1;    node_local_coordinate_27nodebrick[0][2] = 1;
    node_local_coordinate_27nodebrick[1][0] = -1;   node_local_coordinate_27nodebrick[1][1] = 1;    node_local_coordinate_27nodebrick[1][2] = 1;
    node_local_coordinate_27nodebrick[2][0] = -1;   node_local_coordinate_27nodebrick[2][1] = -1;   node_local_coordinate_27nodebrick[2][2] = 1;
    node_local_coordinate_27nodebrick[3][0] = 1;    node_local_coordinate_27nodebrick[3][1] = -1;   node_local_coordinate_27nodebrick[3][2] = 1;
    node_local_coordinate_27nodebrick[4][0] = 1;    node_local_coordinate_27nodebrick[4][1] = 1;    node_local_coordinate_27nodebrick[4][2] = -1;
    node_local_coordinate_27nodebrick[5][0] = -1;   node_local_coordinate_27nodebrick[5][1] = 1;    node_local_coordinate_27nodebrick[5][2] = -1;
    node_local_coordinate_27nodebrick[6][0] = -1;   node_local_coordinate_27nodebrick[6][1] = -1;   node_local_coordinate_27nodebrick[6][2] = -1;
    node_local_coordinate_27nodebrick[7][0] = 1;    node_local_coordinate_27nodebrick[7][1] = -1;   node_local_coordinate_27nodebrick[7][2] = -1;
    node_local_coordinate_27nodebrick[8][0] = 0;    node_local_coordinate_27nodebrick[8][1] = 1;    node_local_coordinate_27nodebrick[8][2] = 1;
    node_local_coordinate_27nodebrick[9][0] = -1;   node_local_coordinate_27nodebrick[9][1] = 0;    node_local_coordinate_27nodebrick[9][2] = 1;
    node_local_coordinate_27nodebrick[10][0] = 0;   node_local_coordinate_27nodebrick[10][1] = -1;  node_local_coordinate_27nodebrick[10][2] = 1;
    node_local_coordinate_27nodebrick[11][0] = 1;   node_local_coordinate_27nodebrick[11][1] = 0;   node_local_coordinate_27nodebrick[11][2] = 1;
    node_local_coordinate_27nodebrick[12][0] = 0;   node_local_coordinate_27nodebrick[12][1] = 1;   node_local_coordinate_27nodebrick[12][2] = -1;
    node_local_coordinate_27nodebrick[13][0] = -1;  node_local_coordinate_27nodebrick[13][1] = 0;   node_local_coordinate_27nodebrick[13][2] = -1;
    node_local_coordinate_27nodebrick[14][0] = 0;   node_local_coordinate_27nodebrick[14][1] = -1;  node_local_coordinate_27nodebrick[14][2] = -1;
    node_local_coordinate_27nodebrick[15][0] = 1;   node_local_coordinate_27nodebrick[15][1] = 0;   node_local_coordinate_27nodebrick[15][2] = -1;
    node_local_coordinate_27nodebrick[16][0] = 1;   node_local_coordinate_27nodebrick[16][1] = 1;   node_local_coordinate_27nodebrick[16][2] = 0;
    node_local_coordinate_27nodebrick[17][0] = -1;  node_local_coordinate_27nodebrick[17][1] = 1;   node_local_coordinate_27nodebrick[17][2] = 0;
    node_local_coordinate_27nodebrick[18][0] = -1;  node_local_coordinate_27nodebrick[18][1] = -1;  node_local_coordinate_27nodebrick[18][2] = 0;
    node_local_coordinate_27nodebrick[19][0] = 1;   node_local_coordinate_27nodebrick[19][1] = -1;  node_local_coordinate_27nodebrick[19][2] = 0;
    node_local_coordinate_27nodebrick[20][0] = 0;   node_local_coordinate_27nodebrick[20][1] = 0;   node_local_coordinate_27nodebrick[20][2] = 0;
    node_local_coordinate_27nodebrick[21][0] = 0;   node_local_coordinate_27nodebrick[21][1] = 1;   node_local_coordinate_27nodebrick[21][2] = 0;
    node_local_coordinate_27nodebrick[22][0] = -1;  node_local_coordinate_27nodebrick[22][1] = 0;   node_local_coordinate_27nodebrick[22][2] = 0;
    node_local_coordinate_27nodebrick[23][0] = 0;   node_local_coordinate_27nodebrick[23][1] = -1;  node_local_coordinate_27nodebrick[23][2] = 0;
    node_local_coordinate_27nodebrick[24][0] = 1;   node_local_coordinate_27nodebrick[24][1] = 0;   node_local_coordinate_27nodebrick[24][2] = 0;
    node_local_coordinate_27nodebrick[25][0] = 0;   node_local_coordinate_27nodebrick[25][1] = 0;   node_local_coordinate_27nodebrick[25][2] = 1;
    node_local_coordinate_27nodebrick[26][0] = 0;   node_local_coordinate_27nodebrick[26][1] = 0;   node_local_coordinate_27nodebrick[26][2] = -1;



	//===========================================================================
    // Turn off the auto-printing when failure occurs 
    // so that we can handle the errors appropriately
    //===========================================================================

	H5::Exception::dontPrint();



    //===========================================================================
    // Read ESSI output
    //===========================================================================

    readESSIOutput();   



    //===========================================================================
    // Create file, groups and datasets for energy density components
    //===========================================================================

    createEnergyDatasets();



    //===========================================================================
    // Initialize arrays for energy density components
    //===========================================================================

    index_kinetic_energy_density = new int[max_node_tag];
    //index_kinetic_energy_density[0] = -1;

    kinetic_energy_density = new double[number_of_nodes*number_of_timesteps];


    index_strain_energy_density = new int[max_element_tag];
    //index_strain_energy_density[0] = -1;

    strain_energy_density_step = new double[number_of_GPs*number_of_timesteps];
    strain_energy_density_cumu = new double[number_of_GPs*number_of_timesteps];

    plastic_work_density_step = new double[number_of_GPs*number_of_timesteps];
    plastic_work_density_cumu = new double[number_of_GPs*number_of_timesteps];

    plastic_free_energy_density_step = new double[number_of_GPs*number_of_timesteps];
    plastic_free_energy_density_cumu = new double[number_of_GPs*number_of_timesteps];

    plastic_dissipation_density_step = new double[number_of_GPs*number_of_timesteps];
    plastic_dissipation_density_cumu = new double[number_of_GPs*number_of_timesteps];



    kinetic_energy_density_eleavg = new double[number_of_elements*number_of_timesteps];
    strain_energy_density_eleavg = new double[number_of_elements*number_of_timesteps];
    plastic_work_density_eleavg = new double[number_of_elements*number_of_timesteps];
    plastic_free_energy_density_eleavg = new double[number_of_elements*number_of_timesteps];
    plastic_dissipation_density_eleavg = new double[number_of_elements*number_of_timesteps];
    mechanical_energy_density_eleavg = new double[number_of_elements*number_of_timesteps];

    index_energy_density_eleavg = new int[max_element_tag];
    //index_energy_density_eleavg[0] = -1;



    total_mechanical_energy = new double[number_of_timesteps];
    total_plastic_work = new double[number_of_timesteps];
    total_plastic_free_energy = new double[number_of_timesteps];
    total_plastic_dissipation = new double[number_of_timesteps];
    total_kinetic_energy = new double[number_of_timesteps];
    total_strain_energy = new double[number_of_timesteps];



    element_volume = new double[max_element_tag];

    //back_stress = new double[number_of_GPs*6*number_of_timesteps];



//    center_coordinates = new double[3*number_of_elements];



    //===========================================================================
    // Set status to initialized and ready to compute energy density components
    //===========================================================================

    is_initialized = true;

}

void Energy_Post_Processing::computeAndWriteEnergyDensity()
{
    if (not is_initialized)
    {
        initialize();
    }

    

    //===========================================================================
    // Compute kinetic_energy_density
    //===========================================================================

    for (int i = 0; i < max_node_tag; i++)
    {
        computeKEnergyDensity(i);
        //Energy_Out << "Kinetic energy density at node " << i << " is successfully computed.\n";
    }
    Energy_Out << "Kinetic energy density is successfully computed.\n";

    //Temporary time step and node for testing -- will be commented later
    //int temp_timestep = 500;
    //int temp_node = 500;
    //for (int i = 0; i < 100; i++)
    //{
    //    Energy_Out << "At time step " << i << ", the kinetic energy density at node " << temp_node << " is " << kinetic_energy_density[temp_node][i] << "\n";
    //}

    //===========================================================================
    // Write kinetic_energy_density
    //===========================================================================

    int rank = 2;
    hsize_t dataspace_dims[2] = {number_of_nodes, number_of_timesteps};
    hsize_t dataspace_maxdims[2] = {number_of_nodes, number_of_timesteps};
    hid_t id_dataspace = H5Screate_simple(rank, dataspace_dims, dataspace_maxdims);

    hsize_t start[2]; // Start of hyperslab
    hsize_t stride[2]; // Stride of hyperslab
    hsize_t count[2];  // Block count
    hsize_t block[2];  // Block sizes
    start[0]  = 0; start[1]  = 0;
    stride[0] = 1; stride[1] = number_of_timesteps;
    count[0]  = number_of_nodes; count[1]  = 1;
    block[0]  = 1; block[1]  = number_of_timesteps;
    herr_t status = H5Sselect_hyperslab (id_dataspace, H5S_SELECT_SET, start, stride, count, block);


    rank = 1;
    hsize_t memspace_dims[1] = {number_of_nodes*number_of_timesteps};
    hsize_t memspace_maxdims[1] = {number_of_nodes*number_of_timesteps};
    hid_t id_memspace = H5Screate_simple(rank, memspace_dims, memspace_maxdims);

    start[0]  = 0;
    stride[0] = number_of_timesteps;
    count[0]  = number_of_nodes;
    block[0]  = number_of_timesteps;
    status = H5Sselect_hyperslab (id_memspace, H5S_SELECT_SET, start, stride, count, block);

    status = H5Dwrite(id_KEnergy_dataset, H5T_NATIVE_DOUBLE, id_memspace, id_dataspace, H5P_DEFAULT, kinetic_energy_density);
    //hdf5_check_error(status);

    //Close stuff
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);
    H5Dclose(id_KEnergy_dataset);

   
    //===========================================================================
    // Write index_kinetic_energy_density
    //===========================================================================

    rank = 1;
    hsize_t data_dims[1] = {max_node_tag};
    id_dataspace = H5Screate_simple(rank, data_dims, data_dims);
    id_memspace = H5Screate_simple(rank, data_dims, data_dims);

    status = H5Dwrite(id_index_to_KEnergy_dataset, H5T_NATIVE_INT, id_memspace, id_dataspace, H5P_DEFAULT, index_kinetic_energy_density);
    //hdf5_check_error(status);

    //Close stuff
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);
    H5Dclose(id_index_to_KEnergy_dataset);

    Energy_Out << "Kinetic energy density is successfully output.\n";



    //===========================================================================
    // Compute strain_energy_density
    //===========================================================================

    last_valid_element = 0;
    is_first_valid_element = true;

    for (int i = 0; i < max_element_tag; i++)
    {
        computeSEnergyDensity(i);
        //Energy_Out << "Index to strain energy density = " << index_strain_energy_density[i] << endl;
        //Energy_Out << "Strain energy density at node " << i << " is successfully computed.\n";
    }
    Energy_Out << "Strain energy density is successfully computed.\n";


    //===========================================================================
    // Write strain_energy_density
    //===========================================================================

    rank = 2;
    hsize_t strain_dataspace_dims[2] = {number_of_GPs, number_of_timesteps};
    hsize_t strain_dataspace_maxdims[2] = {number_of_GPs, number_of_timesteps};
    hid_t id_strain_dataspace = H5Screate_simple(rank, strain_dataspace_dims, strain_dataspace_maxdims);

    start[0]  = 0; start[1]  = 0;
    stride[0] = 1; stride[1] = number_of_timesteps;
    count[0]  = number_of_GPs; count[1]  = 1;
    block[0]  = 1; block[1]  = number_of_timesteps;
    status = H5Sselect_hyperslab (id_strain_dataspace, H5S_SELECT_SET, start, stride, count, block);


    rank = 1;
    hsize_t strain_memspace_dims[1] = {number_of_GPs*number_of_timesteps};
    hsize_t strain_memspace_maxdims[1] = {number_of_GPs*number_of_timesteps};
    hid_t id_strain_memspace = H5Screate_simple(rank, strain_memspace_dims, strain_memspace_maxdims);

    start[0]  = 0;
    stride[0] = number_of_timesteps;
    count[0]  = number_of_GPs;
    block[0]  = number_of_timesteps;
    status = H5Sselect_hyperslab (id_strain_memspace, H5S_SELECT_SET, start, stride, count, block);


    status = H5Dwrite(id_SEnergy_step_dataset, H5T_NATIVE_DOUBLE, id_strain_memspace, id_strain_dataspace, H5P_DEFAULT, strain_energy_density_step);
    status = H5Dwrite(id_SEnergy_cumu_dataset, H5T_NATIVE_DOUBLE, id_strain_memspace, id_strain_dataspace, H5P_DEFAULT, strain_energy_density_cumu);
    //hdf5_check_error(status);


    //Close stuff
    H5Sclose(id_strain_dataspace);
    H5Sclose(id_strain_memspace);
    H5Dclose(id_SEnergy_step_dataset);
    H5Dclose(id_SEnergy_cumu_dataset);

    Energy_Out << "Strain energy density is successfully output.\n";



    //===========================================================================
    // Write plastic_work_density_step and plastic_work_density_cumu
    //===========================================================================

    rank = 2;
    hsize_t plastic_dataspace_dims[2] = {number_of_GPs, number_of_timesteps};
    hsize_t plastic_dataspace_maxdims[2] = {number_of_GPs, number_of_timesteps};
    hid_t id_plastic_dataspace = H5Screate_simple(rank, plastic_dataspace_dims, plastic_dataspace_maxdims);

    start[0]  = 0; start[1]  = 0;
    stride[0] = 1; stride[1] = number_of_timesteps;
    count[0]  = number_of_GPs; count[1]  = 1;
    block[0]  = 1; block[1]  = number_of_timesteps;
    status = H5Sselect_hyperslab (id_plastic_dataspace, H5S_SELECT_SET, start, stride, count, block);


    rank = 1;
    hsize_t plastic_memspace_dims[1] = {number_of_GPs*number_of_timesteps};
    hsize_t plastic_memspace_maxdims[1] = {number_of_GPs*number_of_timesteps};
    hid_t id_plastic_memspace = H5Screate_simple(rank, plastic_memspace_dims, plastic_memspace_maxdims);

    start[0]  = 0;
    stride[0] = number_of_timesteps;
    count[0]  = number_of_GPs;
    block[0]  = number_of_timesteps;
    status = H5Sselect_hyperslab (id_plastic_memspace, H5S_SELECT_SET, start, stride, count, block);


    status = H5Dwrite(id_PWork_step_dataset, H5T_NATIVE_DOUBLE, id_plastic_memspace, id_plastic_dataspace, H5P_DEFAULT, plastic_work_density_step);
    status = H5Dwrite(id_PWork_cumu_dataset, H5T_NATIVE_DOUBLE, id_plastic_memspace, id_plastic_dataspace, H5P_DEFAULT, plastic_work_density_cumu);
    
    //hdf5_check_error(status);


    //Close stuff
    //H5Sclose(id_plastic_dataspace);
    //H5Sclose(id_plastic_memspace);
    H5Dclose(id_PWork_step_dataset);
    H5Dclose(id_PWork_cumu_dataset);

    Energy_Out << "Plastic work (each time step) is successfully output.\n";
    Energy_Out << "Plastic work (cumulative) is successfully output.\n";



    //===========================================================================
    // Write plastic_free_energy_density_step and plastic_free_energy_density_cumu
    //===========================================================================

    status = H5Dwrite(id_plastic_free_energy_density_step_dataset, H5T_NATIVE_DOUBLE, id_plastic_memspace, id_plastic_dataspace, H5P_DEFAULT, plastic_free_energy_density_step);
    status = H5Dwrite(id_plastic_free_energy_density_cumu_dataset, H5T_NATIVE_DOUBLE, id_plastic_memspace, id_plastic_dataspace, H5P_DEFAULT, plastic_free_energy_density_cumu);
    
    //hdf5_check_error(status);


    //Close stuff
    //H5Sclose(id_plastic_dataspace);
    //H5Sclose(id_plastic_memspace);
    H5Dclose(id_plastic_free_energy_density_step_dataset);
    H5Dclose(id_plastic_free_energy_density_cumu_dataset);

    Energy_Out << "Plastic free energy (each time step) is successfully output.\n";
    Energy_Out << "Plastic free energy (cumulative) is successfully output.\n";



    //===========================================================================
    // Write plastic_dissipation_density_step and plastic_dissipation_density_cumu
    //===========================================================================

    status = H5Dwrite(id_plastic_dissipation_density_step_dataset, H5T_NATIVE_DOUBLE, id_plastic_memspace, id_plastic_dataspace, H5P_DEFAULT, plastic_dissipation_density_step);
    status = H5Dwrite(id_plastic_dissipation_density_cumu_dataset, H5T_NATIVE_DOUBLE, id_plastic_memspace, id_plastic_dataspace, H5P_DEFAULT, plastic_dissipation_density_cumu);
    
    //hdf5_check_error(status);


    //Close stuff
    H5Sclose(id_plastic_dataspace);
    H5Sclose(id_plastic_memspace);
    H5Dclose(id_plastic_dissipation_density_step_dataset);
    H5Dclose(id_plastic_dissipation_density_cumu_dataset);

    Energy_Out << "Plastic dissipation (each time step) is successfully output.\n";
    Energy_Out << "Plastic dissipation (cumulative) is successfully output.\n";



    //===========================================================================
    // Write index_strain_energy_density
    //===========================================================================

    rank = 1;
    data_dims[0] = max_element_tag;
    id_dataspace = H5Screate_simple(rank, data_dims, data_dims);
    id_memspace = H5Screate_simple(rank, data_dims, data_dims);

    status = H5Dwrite(id_index_to_SEnergy_dataset, H5T_NATIVE_INT, id_memspace, id_dataspace, H5P_DEFAULT, index_strain_energy_density);
    //hdf5_check_error(status);

    //Close stuff
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);
    H5Dclose(id_index_to_SEnergy_dataset);

    Energy_Out << "Strain energy density index is successfully output.\n";



    //===========================================================================
    // Compute element average energy density
    //===========================================================================

    for (int i = 0; i < max_element_tag; i++)
    {
        //Energy_Out << "Average energy density at elements " << i+1 << " is being computed.\n";
        if (class_tags[i] == 8001)
        {
            computeElementAverage_8NodeBrick(i);
        }
        else if (class_tags[i] == 8002)
        {
            computeElementAverage_27NodeBrick(i);
        }
        else if (class_tags[i] == -1)
        {
            //Energy_Out << "Void element tag...\n";
        }
        else
        {
            Energy_Error << "Element " << i << ": Element type currently not implemented!\n";
        }
        
        //Energy_Out << "Average energy density at elements " << i+1 << " is successfully computed.\n";
    }
    Energy_Out << "Element average energy density is successfully computed.\n";



    //===========================================================================
    // Write plastic_work_density_eleavg
    //===========================================================================

    rank = 2;
    hsize_t pavg_dataspace_dims[2] = {number_of_elements, number_of_timesteps};
    hsize_t pavg_dataspace_maxdims[2] = {number_of_elements, number_of_timesteps};
    hid_t id_pavg_dataspace = H5Screate_simple(rank, pavg_dataspace_dims, pavg_dataspace_maxdims);

    start[0]  = 0; start[1]  = 0;
    stride[0] = 1; stride[1] = number_of_timesteps;
    count[0]  = number_of_elements; count[1]  = 1;
    block[0]  = 1; block[1]  = number_of_timesteps;
    status = H5Sselect_hyperslab (id_pavg_dataspace, H5S_SELECT_SET, start, stride, count, block);


    rank = 1;
    hsize_t pavg_memspace_dims[1] = {number_of_elements*number_of_timesteps};
    hsize_t pavg_memspace_maxdims[1] = {number_of_elements*number_of_timesteps};
    hid_t id_pavg_memspace = H5Screate_simple(rank, pavg_memspace_dims, pavg_memspace_maxdims);

    start[0]  = 0;
    stride[0] = number_of_timesteps;
    count[0]  = number_of_elements;
    block[0]  = number_of_timesteps;
    status = H5Sselect_hyperslab (id_pavg_memspace, H5S_SELECT_SET, start, stride, count, block);


    status = H5Dwrite(id_PWork_avg_dataset, H5T_NATIVE_DOUBLE, id_pavg_memspace, id_pavg_dataspace, H5P_DEFAULT, plastic_work_density_eleavg);
    //hdf5_check_error(status);


    //Close stuff
    H5Sclose(id_pavg_dataspace);
    H5Sclose(id_pavg_memspace);
    H5Dclose(id_PWork_avg_dataset);

    Energy_Out << "Element average plastic energy density is successfully output.\n";



    //===========================================================================
    // Write mechanical_energy_density_eleavg
    //===========================================================================

    rank = 2;
    hsize_t tavg_dataspace_dims[2] = {number_of_elements, number_of_timesteps};
    hsize_t tavg_dataspace_maxdims[2] = {number_of_elements, number_of_timesteps};
    hid_t id_tavg_dataspace = H5Screate_simple(rank, tavg_dataspace_dims, tavg_dataspace_maxdims);

    start[0]  = 0; start[1]  = 0;
    stride[0] = 1; stride[1] = number_of_timesteps;
    count[0]  = number_of_elements; count[1]  = 1;
    block[0]  = 1; block[1]  = number_of_timesteps;
    status = H5Sselect_hyperslab (id_tavg_dataspace, H5S_SELECT_SET, start, stride, count, block);


    rank = 1;
    hsize_t tavg_memspace_dims[1] = {number_of_elements*number_of_timesteps};
    hsize_t tavg_memspace_maxdims[1] = {number_of_elements*number_of_timesteps};
    hid_t id_tavg_memspace = H5Screate_simple(rank, tavg_memspace_dims, tavg_memspace_maxdims);

    start[0]  = 0;
    stride[0] = number_of_timesteps;
    count[0]  = number_of_elements;
    block[0]  = number_of_timesteps;
    status = H5Sselect_hyperslab (id_tavg_memspace, H5S_SELECT_SET, start, stride, count, block);


    status = H5Dwrite(id_mechanical_energy_avg_dataset, H5T_NATIVE_DOUBLE, id_tavg_memspace, id_tavg_dataspace, H5P_DEFAULT, mechanical_energy_density_eleavg);
    //hdf5_check_error(status);


    //Close stuff
    H5Sclose(id_tavg_dataspace);
    H5Sclose(id_tavg_memspace);
    H5Dclose(id_mechanical_energy_avg_dataset);

    Energy_Out << "Element average total energy density is successfully output.\n";

  

    //===========================================================================
    // Write index_energy_density_eleavg
    //===========================================================================

    rank = 1;
    data_dims[0] = max_element_tag;
    id_dataspace = H5Screate_simple(rank, data_dims, data_dims);
    id_memspace = H5Screate_simple(rank, data_dims, data_dims);

    status = H5Dwrite(id_index_to_Energy_avg_dataset, H5T_NATIVE_INT, id_memspace, id_dataspace, H5P_DEFAULT, index_energy_density_eleavg);
    //hdf5_check_error(status);

    //Close stuff
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);
    H5Dclose(id_index_to_Energy_avg_dataset);

    Energy_Out << "Element average energy density index is successfully output.\n";


/*
    //===========================================================================
    // Compute center_coordinates
    //===========================================================================

    for (int i = 0; i < number_of_elements; i++)
    {
        getCenterofElement(i+1);
        //Energy_Out << "Center coordinates of elements " << i << " is successfully computed.\n";
    }
    Energy_Out << "Center coordinates is successfully computed.\n";//



    //===========================================================================
    // Write center_coordinates
    //===========================================================================

    rank = 2;
    hsize_t center_dataspace_dims[2] = {number_of_elements, 3};
    hsize_t center_dataspace_maxdims[2] = {number_of_elements, 3};
    hid_t id_center_dataspace = H5Screate_simple(rank, center_dataspace_dims, center_dataspace_maxdims);

    start[0]  = 0; start[1]  = 0;
    stride[0] = 1; stride[1] = 3;
    count[0]  = number_of_elements; count[1]  = 1;
    block[0]  = 1; block[1]  = 3;
    status = H5Sselect_hyperslab (id_center_dataspace, H5S_SELECT_SET, start, stride, count, block);


    rank = 1;
    hsize_t center_memspace_dims[1] = {number_of_elements*3};
    hsize_t center_memspace_maxdims[1] = {number_of_elements*3};
    hid_t id_center_memspace = H5Screate_simple(rank, center_memspace_dims, center_memspace_maxdims);

    start[0]  = 0;
    stride[0] = 3;
    count[0]  = number_of_elements;
    block[0]  = 3;
    status = H5Sselect_hyperslab (id_center_memspace, H5S_SELECT_SET, start, stride, count, block);


    status = H5Dwrite(id_center_dataset, H5T_NATIVE_DOUBLE, id_center_memspace, id_center_dataspace, H5P_DEFAULT, center_coordinates);
    //hdf5_check_error(status);


    //Close stuff
    H5Sclose(id_center_dataspace);
    H5Sclose(id_center_memspace);
    H5Dclose(id_center_dataset);
*/   

    //===========================================================================
    // Compute element volume
    //===========================================================================

    for (int i = 0; i < max_element_tag; i++)
    {
        if (class_tags[i] == 8001)
        {
            computeElementVolume_8NodeBrick(i);
        }
        else if (class_tags[i] == 8002)
        {
            computeElementVolume_27NodeBrick(i);
        }
        else if (class_tags[i] == -1)
        {
            element_volume[i] = -1;
        }
        else
        {
            Energy_Error << "Element " << i << ": Element type currently not implemented!\n";
        }
        
        //Energy_Out << "Volume of element " << i+1 << " is successfully computed.\n";
    }
    Energy_Out << "Element volume is successfully computed.\n";



    //===========================================================================
    // Write element_volume
    //===========================================================================

    rank = 1;
    data_dims[0] = max_element_tag;
    id_dataspace = H5Screate_simple(rank, data_dims, data_dims);
    id_memspace = H5Screate_simple(rank, data_dims, data_dims);

    status = H5Dwrite(id_element_volume_dataset, H5T_NATIVE_DOUBLE, id_memspace, id_dataspace, H5P_DEFAULT, element_volume);
    //hdf5_check_error(status);

    //Close stuff
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);
    H5Dclose(id_element_volume_dataset);

    Energy_Out << "Element volume is successfully output.\n";


/*
    //===========================================================================
    // Write back_stress
    //===========================================================================

    rank = 2;
    hsize_t back_stress_dataspace_dims[2] = {number_of_GPs*6, number_of_timesteps};
    hsize_t back_stress_dataspace_maxdims[2] = {number_of_GPs*6, number_of_timesteps};
    hid_t id_back_stress_dataspace = H5Screate_simple(rank, back_stress_dataspace_dims, back_stress_dataspace_maxdims);

    start[0]  = 0; start[1]  = 0;
    stride[0] = 6; stride[1] = number_of_timesteps;
    count[0]  = number_of_GPs; count[1]  = 1;
    block[0]  = 6; block[1]  = number_of_timesteps;
    status = H5Sselect_hyperslab (id_back_stress_dataspace, H5S_SELECT_SET, start, stride, count, block);


    rank = 1;
    hsize_t back_stress_memspace_dims[1] = {number_of_GPs*6*number_of_timesteps};
    hsize_t back_stress_memspace_maxdims[1] = {number_of_GPs*6*number_of_timesteps};
    hid_t id_back_stress_memspace = H5Screate_simple(rank, back_stress_memspace_dims, back_stress_memspace_maxdims);

    start[0]  = 0;
    stride[0] = number_of_timesteps;
    count[0]  = number_of_GPs*6;
    block[0]  = number_of_timesteps;
    status = H5Sselect_hyperslab (id_back_stress_memspace, H5S_SELECT_SET, start, stride, count, block);

    Energy_Out << "number_of_GPs: " << number_of_GPs << endl;
    Energy_Out << "back_stress: " << back_stress[23998] << endl;

    status = H5Dwrite(id_back_stress_dataset, H5T_NATIVE_DOUBLE, id_back_stress_memspace, id_back_stress_dataspace, H5P_DEFAULT, back_stress);
    //hdf5_check_error(status);


    //Close stuff
    H5Sclose(id_back_stress_dataspace);
    H5Sclose(id_back_stress_memspace);
    H5Dclose(id_back_stress_dataset);

    Energy_Out << "Back stress is successfully output.\n";
*/
 

    //===========================================================================
    // Compute total energy of the whole model at each time step
    //===========================================================================

    computeTotalEnergy();
    Energy_Out << "Total energy and total plastic dissipation is successfully computed.\n";



    //===========================================================================
    // Write total_mechanical_energy and total_plastic_work
    //===========================================================================
   
    rank = 1;
    data_dims[0] = number_of_timesteps;
    id_dataspace = H5Screate_simple(rank, data_dims, data_dims);
    id_memspace = H5Screate_simple(rank, data_dims, data_dims);

    status = H5Dwrite(id_total_mechanical_energy_dataset, H5T_NATIVE_DOUBLE, id_memspace, id_dataspace, H5P_DEFAULT, total_mechanical_energy);
    status = H5Dwrite(id_total_plastic_work_dataset, H5T_NATIVE_DOUBLE, id_memspace, id_dataspace, H5P_DEFAULT, total_plastic_work);
    status = H5Dwrite(id_total_plastic_free_energy_dataset, H5T_NATIVE_DOUBLE, id_memspace, id_dataspace, H5P_DEFAULT, total_plastic_free_energy);
    status = H5Dwrite(id_total_plastic_dissipation_dataset, H5T_NATIVE_DOUBLE, id_memspace, id_dataspace, H5P_DEFAULT, total_plastic_dissipation);
    status = H5Dwrite(id_total_kinetic_energy_dataset, H5T_NATIVE_DOUBLE, id_memspace, id_dataspace, H5P_DEFAULT, total_kinetic_energy);
    status = H5Dwrite(id_total_strain_energy_dataset, H5T_NATIVE_DOUBLE, id_memspace, id_dataspace, H5P_DEFAULT, total_strain_energy);
    //hdf5_check_error(status);

    //Close stuff
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);
    H5Dclose(id_total_mechanical_energy_dataset);
    H5Dclose(id_total_plastic_work_dataset);
    H5Dclose(id_total_plastic_free_energy_dataset);
    H5Dclose(id_total_plastic_dissipation_dataset);
    H5Dclose(id_total_kinetic_energy_dataset);
    H5Dclose(id_total_strain_energy_dataset);

    Energy_Out << "Total energy and total plastic dissipation is successfully output.\n";

}



void Energy_Post_Processing::readESSIOutput()
{

    //===========================================================================
    // Open the specified file for read and write
    //===========================================================================

    hid_t file_access_plist = H5Pcreate(H5P_FILE_ACCESS);


    unsigned flags = H5F_ACC_RDWR;
    id_file = H5Fopen(HDF5filename.c_str(), flags, file_access_plist);

    if (id_file < 0)
    {
        Energy_Error << "Unable to open file: " << HDF5filename << endl;
        return;
    }



    //===========================================================================
    // Read "Number_of_Elements" array, get number of elements
    //===========================================================================

    hsize_t data_dims[1];
    hsize_t data_maxdims[1];
    hsize_t rank_one_array = 1;
    hsize_t id_memspace = 0;

    hsize_t id_elements = H5Dopen(id_file, "Number_of_Elements", H5P_DEFAULT);

    if (id_elements < 0)
    {
        Energy_Error << "Could not open \"Number_of_Elements\" array!\n ";
        return;
    }

    //Open dataspace and get number of elements from it
    int * temp_number_of_elements = 0;
    temp_number_of_elements = new int[1];

    hsize_t id_dataspace = H5Dget_space(id_elements);

    data_dims[0] = 1;
    data_maxdims[0] = 1;
    id_memspace  = H5Screate_simple(rank_one_array, data_dims, data_maxdims);       // create dataspace of memory
    H5Dread(id_elements, H5T_NATIVE_INT, id_memspace, id_dataspace, H5P_DEFAULT, temp_number_of_elements);
    number_of_elements = temp_number_of_elements[0];

    if (number_of_elements == -1)
    {
        Energy_Error << "\"Number_of_Elements\" is not right!\n";
    }
    else
    {
        //Energy_Out << "Dataset \"" << HDF5filename << "\" has " << number_of_elements << " total elements.\n";
    }

    //Cleanup
    H5Dclose(id_elements);
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);



    //===========================================================================
    // Read "Number_of_Nodes" array, get number of nodes
    //===========================================================================

    hsize_t id_num_nodes = H5Dopen(id_file, "Number_of_Nodes", H5P_DEFAULT);

    if (id_num_nodes < 0)
    {
        Energy_Error << "Could not open \"Number_of_Nodes\" array!\n ";
        return;
    }

    //Open dataspace and get number of nodes from it
    int * temp_number_of_nodes = 0;
    temp_number_of_nodes = new int[1];

    id_dataspace = H5Dget_space(id_num_nodes);

    id_memspace  = H5Screate_simple(rank_one_array, data_dims, data_maxdims);       // create dataspace of memory
    H5Dread(id_num_nodes, H5T_NATIVE_INT, id_memspace, id_dataspace, H5P_DEFAULT, temp_number_of_nodes);
    number_of_nodes = temp_number_of_nodes[0];

    if (number_of_nodes == -1)
    {
        Energy_Error << "\"Number_of_Nodes\" is not right!\n";
    }
    else
    {
        //Energy_Out << "Dataset \"" << HDF5filename << "\" has " << number_of_nodes << " total nodes.\n";
    }

    //Cleanup
    H5Dclose(id_num_nodes);
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);



    //===========================================================================
    // Read "Time" arrays
    //===========================================================================
    hsize_t id_time = H5Dopen(id_file, "time", H5P_DEFAULT);

    if (id_time < 0)
    {
        Energy_Error << "Could not open \"time\" array!\n ";
        return;
    }

    id_dataspace = H5Dget_space(id_time);

    int ndims;
    ndims =  H5Sget_simple_extent_ndims(id_dataspace);

    if (ndims != 1)
    {
        Energy_Error << "\"time\" array should be a 1-D array.\n";
        return;
    }

    H5Sget_simple_extent_dims(id_dataspace, data_dims, data_maxdims);
    number_of_timesteps = data_dims[0];

    double* all_timesteps = new double[number_of_timesteps];

    id_memspace  = H5Screate_simple(rank_one_array, data_dims, data_maxdims);       // create dataspace of memory

    H5Dread(id_time, H5T_NATIVE_DOUBLE, id_memspace, id_dataspace, H5P_DEFAULT, all_timesteps);

    //Energy_Out << "Dataset \"" << HDF5filename  << "\" has " <<  number_of_timesteps << " time-steps.\n";

    tfinal = all_timesteps[number_of_timesteps - 1];
    //Energy_Out << "The final time-step of dataset \"" << HDF5filename  << "\" is " <<  tfinal << "\n";

    delta_t = all_timesteps[1] - all_timesteps[0];
    //Energy_Out << "The time interval of dataset \"" << HDF5filename  << "\" is " <<  delta_t << "\n";
    
    // delete all_timesteps;
    H5Dclose(id_time);
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);



    //===========================================================================
    // Read "Index_to_Generalized_Displacements" arrays
    //===========================================================================

    hsize_t id_index_GDs = H5Dopen(id_file, "Model/Nodes/Index_to_Generalized_Displacements", H5P_DEFAULT);

    if (id_index_GDs < 0)
    {
        Energy_Error << "Could not open \"Index_to_Generalized_Displacements\" array!\n ";
        return;
    }

    id_dataspace = H5Dget_space(id_index_GDs);

    ndims =  H5Sget_simple_extent_ndims(id_dataspace);

    if (ndims != 1)
    {
        Energy_Error << "\"Index_to_Generalized_Displacements\" array should be a 1-D array.\n";
        return;
    }

    H5Sget_simple_extent_dims(id_dataspace, data_dims, data_maxdims);

    max_node_tag = data_dims[0];

    index_GDs = new int[max_node_tag];

    id_memspace  = H5Screate_simple(rank_one_array, data_dims, data_maxdims);       // create dataspace of memory

    H5Dread(id_index_GDs, H5T_NATIVE_INT, id_memspace, id_dataspace, H5P_DEFAULT, index_GDs);

    //Energy_Out << "Dataset \"" << HDF5filename  << "\" has " <<  max_node_tag << " indices of GDs.\n";
    
    // delete all_timesteps;
    H5Dclose(id_index_GDs);
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);



    //===========================================================================
    // Open "Generalized_Displacements" array for reading
    //===========================================================================

    id_displacements = H5Dopen(id_file, "Model/Nodes/Generalized_Displacements", H5P_DEFAULT);

    if (id_displacements < 0)
    {
        Energy_Error << "Could not open \"Displacements\" array!\n ";
        return;
    }

    //Energy_Out << "\"Generalized Displacements\" is successfully opened!\n";

    id_displacements_dataspace = H5Dget_space(id_displacements);



    //===========================================================================
    // Open "Outputs" array for reading
    //===========================================================================

    id_outputs = H5Dopen(id_file, "Model/Elements/Outputs", H5P_DEFAULT);

    if (id_outputs < 0)
    {
        Energy_Error << "Could not open \"Outputs\" array!\n ";
        return;
    }

    //Energy_Out << "\"Outputs\" is successfully opened!\n";

    id_outputs_dataspace = H5Dget_space(id_outputs);



    //===========================================================================
    // Read "Number_of_Gauss_Points"
    //===========================================================================

    hsize_t id_number_of_GPs_element = H5Dopen(id_file, "Model/Elements/Number_of_Gauss_Points", H5P_DEFAULT);

    if (id_number_of_GPs_element < 0)
    {
        Energy_Error << "Could not open \"Number_of_Gauss_Points\" array!\n ";
        return;
    }

    id_dataspace = H5Dget_space(id_number_of_GPs_element);

    ndims =  H5Sget_simple_extent_ndims(id_dataspace);

    if (ndims != 1)
    {
        Energy_Error << "\"Number_of_Gauss_Points\" array should be a 1-D array.\n";
        return;
    }

    H5Sget_simple_extent_dims(id_dataspace, data_dims, data_maxdims);

    max_element_tag = data_dims[0];

    number_of_GPs_element = new int[max_element_tag];

    id_memspace  = H5Screate_simple(rank_one_array, data_dims, data_maxdims);       // create dataspace of memory

    H5Dread(id_number_of_GPs_element, H5T_NATIVE_INT, id_memspace, id_dataspace, H5P_DEFAULT, number_of_GPs_element);

    //Energy_Out << "Dataset \"" << HDF5filename  << "\" has " <<  max_element_tag << " indices of \"Number_of_Gauss_Points\".\n";

    number_of_GPs = 0;
    for (int i = 0; i < max_element_tag; i++)
    {
        if (number_of_GPs_element[i] > -1)
        {
            number_of_GPs += number_of_GPs_element[i];
        }
    }
    //Energy_Out << number_of_GPs << endl;
    
    //Clean up
    H5Dclose(id_number_of_GPs_element);
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);


    //===========================================================================
    // Read "Number_of_Nodes"
    //===========================================================================

    hsize_t id_number_of_nodes_element = H5Dopen(id_file, "Model/Elements/Number_of_Nodes", H5P_DEFAULT);

    if (id_number_of_nodes_element < 0)
    {
        Energy_Error << "Could not open \"Number_of_Nodes\" array!\n ";
        return;
    }

    id_dataspace = H5Dget_space(id_number_of_nodes_element);

    ndims =  H5Sget_simple_extent_ndims(id_dataspace);

    if (ndims != 1)
    {
        Energy_Error << "\"Number_of_Nodes\" array should be a 1-D array.\n";
        return;
    }

    H5Sget_simple_extent_dims(id_dataspace, data_dims, data_maxdims);

    number_of_nodes_element = new int[max_element_tag];

    id_memspace  = H5Screate_simple(rank_one_array, data_dims, data_maxdims);       // create dataspace of memory

    H5Dread(id_number_of_nodes_element, H5T_NATIVE_INT, id_memspace, id_dataspace, H5P_DEFAULT, number_of_nodes_element);

    //Energy_Out << "Dataset \"" << HDF5filename  << "\" has " <<  max_element_tag << " indices of \"Number_of_Nodes\".\n";
    
    //Clean up
    H5Dclose(id_number_of_nodes_element);
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);



    //===========================================================================
    // Read "Index_to_Outputs"
    //===========================================================================

    hsize_t id_index_to_outputs = H5Dopen(id_file, "Model/Elements/Index_to_Outputs", H5P_DEFAULT);

    if (id_index_to_outputs < 0)
    {
        Energy_Error << "Could not open \"Index_to_Outputs\" array!\n ";
        return;
    }

    id_dataspace = H5Dget_space(id_index_to_outputs);

    ndims =  H5Sget_simple_extent_ndims(id_dataspace);

    if (ndims != 1)
    {
        Energy_Error << "\"Index_to_Outputs\" array should be a 1-D array.\n";
        return;
    }

    H5Sget_simple_extent_dims(id_dataspace, data_dims, data_maxdims);

    index_to_outputs = new int[max_element_tag];

    id_memspace  = H5Screate_simple(rank_one_array, data_dims, data_maxdims);       // create dataspace of memory

    H5Dread(id_index_to_outputs, H5T_NATIVE_INT, id_memspace, id_dataspace, H5P_DEFAULT, index_to_outputs);

    //Energy_Out << "Dataset \"" << HDF5filename  << "\" has " <<  max_element_tag << " indices of \"Index_to_Outputs\".\n";
    
    //Clean up
    H5Dclose(id_index_to_outputs);
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);


    //===========================================================================
    // Read "Index_to_Connectivity"
    //===========================================================================

    hsize_t id_index_to_connectivity = H5Dopen(id_file, "Model/Elements/Index_to_Connectivity", H5P_DEFAULT);

    if (id_index_to_connectivity < 0)
    {
        Energy_Error << "Could not open \"Index_to_Connectivity\" array!\n ";
        return;
    }

    id_dataspace = H5Dget_space(id_index_to_connectivity);

    ndims =  H5Sget_simple_extent_ndims(id_dataspace);

    if (ndims != 1)
    {
        Energy_Error << "\"Index_to_Connectivity\" array should be a 1-D array.\n";
        return;
    }

    H5Sget_simple_extent_dims(id_dataspace, data_dims, data_maxdims);

    index_to_connectivity = new int[max_element_tag];

    id_memspace  = H5Screate_simple(rank_one_array, data_dims, data_maxdims);       // create dataspace of memory

    H5Dread(id_index_to_connectivity, H5T_NATIVE_INT, id_memspace, id_dataspace, H5P_DEFAULT, index_to_connectivity);

    //Energy_Out << "Dataset \"" << HDF5filename  << "\" has " <<  max_element_tag << " indices of \"Index_to_Connectivity\".\n";
    
    //Clean up
    H5Dclose(id_index_to_connectivity);
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);



    //===========================================================================
    // Read "Connectivity"
    //===========================================================================

    hsize_t id_connectivity = H5Dopen(id_file, "Model/Elements/Connectivity", H5P_DEFAULT);

    if (id_connectivity < 0)
    {
        Energy_Error << "Could not open \"Connectivity\" array!\n ";
        return;
    }

    id_dataspace = H5Dget_space(id_connectivity);

    ndims =  H5Sget_simple_extent_ndims(id_dataspace);

    if (ndims != 1)
    {
        Energy_Error << "\"Connectivity\" array should be a 1-D array.\n";
        return;
    }

    H5Sget_simple_extent_dims(id_dataspace, data_dims, data_maxdims);

    int number_of_indices = data_dims[0];

    connectivity = new int[number_of_indices];

    id_memspace  = H5Screate_simple(rank_one_array, data_dims, data_maxdims);       // create dataspace of memory

    H5Dread(id_connectivity, H5T_NATIVE_INT, id_memspace, id_dataspace, H5P_DEFAULT, connectivity);

    //Energy_Out << "Dataset \"" << HDF5filename  << "\" has " <<  number_of_indices << " indices of \"Connectivity\".\n";
    
    //Clean up
    H5Dclose(id_connectivity);
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);



    //===========================================================================
    // Read "Index_to_Coordinates"
    //===========================================================================

    hsize_t id_index_to_node_coordinates = H5Dopen(id_file, "Model/Nodes/Index_to_Coordinates", H5P_DEFAULT);

    if (id_index_to_node_coordinates < 0)
    {
        Energy_Error << "Could not open \"Index_to_Coordinates\" array!\n ";
        return;
    }

    id_dataspace = H5Dget_space(id_index_to_node_coordinates);

    ndims =  H5Sget_simple_extent_ndims(id_dataspace);

    if (ndims != 1)
    {
        Energy_Error << "\"Index_to_Coordinates\" array should be a 1-D array.\n";
        return;
    }

    H5Sget_simple_extent_dims(id_dataspace, data_dims, data_maxdims);

    index_to_node_coordinates = new int[max_node_tag];

    id_memspace  = H5Screate_simple(rank_one_array, data_dims, data_maxdims);       // create dataspace of memory

    H5Dread(id_index_to_node_coordinates, H5T_NATIVE_INT, id_memspace, id_dataspace, H5P_DEFAULT, index_to_node_coordinates);

    //Energy_Out << "Dataset \"" << HDF5filename  << "\" has " <<  max_node_tag << " indices of \"Index_to_Coordinates\".\n";
    
    //Clean up
    H5Dclose(id_index_to_node_coordinates);
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);



    //===========================================================================
    // Read "Coordinates"
    //===========================================================================

    hsize_t id_node_coordinates = H5Dopen(id_file, "Model/Nodes/Coordinates", H5P_DEFAULT);

    if (id_node_coordinates < 0)
    {
        Energy_Error << "Could not open \"Coordinates\" array!\n ";
        return;
    }

    id_dataspace = H5Dget_space(id_node_coordinates);

    ndims =  H5Sget_simple_extent_ndims(id_dataspace);

    if (ndims != 1)
    {
        Energy_Error << "\"Coordinates\" array should be a 1-D array.\n";
        return;
    }

    H5Sget_simple_extent_dims(id_dataspace, data_dims, data_maxdims);

    number_of_indices = data_dims[0];

    node_coordinates = new double[number_of_indices];

    id_memspace  = H5Screate_simple(rank_one_array, data_dims, data_maxdims);       // create dataspace of memory

    H5Dread(id_node_coordinates, H5T_NATIVE_DOUBLE, id_memspace, id_dataspace, H5P_DEFAULT, node_coordinates);

    //Energy_Out << "Dataset \"" << HDF5filename  << "\" has " <<  number_of_indices << " indices of \"Coordinates\".\n";
    
    //Clean up
    H5Dclose(id_node_coordinates);
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);



    //===========================================================================
    // Read "Class_Tags"
    //===========================================================================

    hsize_t id_class_tags = H5Dopen(id_file, "Model/Elements/Class_Tags", H5P_DEFAULT);

    if (id_class_tags < 0)
    {
        Energy_Error << "Could not open \"Class_Tags\" array!\n ";
        return;
    }

    id_dataspace = H5Dget_space(id_class_tags);

    ndims =  H5Sget_simple_extent_ndims(id_dataspace);

    if (ndims != 1)
    {
        Energy_Error << "\"Class_Tags\" array should be a 1-D array.\n";
        return;
    }

    H5Sget_simple_extent_dims(id_dataspace, data_dims, data_maxdims);

    class_tags = new int[max_element_tag];

    id_memspace  = H5Screate_simple(rank_one_array, data_dims, data_maxdims);       // create dataspace of memory

    H5Dread(id_class_tags, H5T_NATIVE_INT, id_memspace, id_dataspace, H5P_DEFAULT, class_tags);

    //Energy_Out << "Dataset \"" << HDF5filename  << "\" has " <<  max_element_tag << " indices of \"Class_Tags\".\n";
    
    //Clean up
    H5Dclose(id_class_tags);
    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);
}



void Energy_Post_Processing::createEnergyDatasets()
{
/*
    // Create energy file
    std::string energy_filename = HDF5filename;

    size_t f = energy_filename.find(".feioutput");
    energy_filename.replace(f, std::string(".feioutput").length(), ".feipost");

    id_file_energy = H5Fcreate(energy_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
*/

    // Set dataset dimensions
    int rank = 2;
    hsize_t dims_KEnergy[rank] = {number_of_nodes, number_of_timesteps};
    hsize_t maxdims_KEnergy[rank] = {number_of_nodes, number_of_timesteps};
    hid_t id_dataspace_KEnergy = H5Screate_simple(rank, dims_KEnergy, maxdims_KEnergy);

    hsize_t dims_SEnergy[rank] = {number_of_GPs, number_of_timesteps};
    hsize_t maxdims_SEnergy[rank] = {number_of_GPs, number_of_timesteps};
    hid_t id_dataspace_SEnergy = H5Screate_simple(rank, dims_SEnergy, maxdims_SEnergy);

    hsize_t dims_PWork[rank] = {number_of_GPs, number_of_timesteps};
    hsize_t maxdims_PWork[rank] = {number_of_GPs, number_of_timesteps};
    hid_t id_dataspace_PWork = H5Screate_simple(rank, dims_PWork, maxdims_PWork);

    hsize_t dims_TEnergy_avg[rank] = {number_of_elements, number_of_timesteps};
    hsize_t maxdims_TEnergy_avg[rank] = {number_of_elements, number_of_timesteps};
    hid_t id_dataspace_TEnergy_avg = H5Screate_simple(rank, dims_TEnergy_avg, maxdims_TEnergy_avg);

    hsize_t dims_PWork_avg[rank] = {number_of_elements, number_of_timesteps};
    hsize_t maxdims_PWork_avg[rank] = {number_of_elements, number_of_timesteps};
    hid_t id_dataspace_PWork_avg = H5Screate_simple(rank, dims_PWork_avg, maxdims_PWork_avg);

    hsize_t dims_center[rank] = {number_of_elements, 3};
    hsize_t maxdims_center[rank] = {number_of_elements, 3};
    hid_t id_dataspace_center = H5Screate_simple(rank, dims_center, maxdims_center);
/*
    hsize_t dims_back_stress[rank] = {number_of_GPs*6, number_of_timesteps};
    hsize_t maxdims_back_stress[rank] = {number_of_GPs*6, number_of_timesteps};
    hid_t id_dataspace_back_stress = H5Screate_simple(rank, dims_back_stress, maxdims_back_stress);
*/

    rank = 1;
    hsize_t dims_index_to_KEnergy[rank] = {max_node_tag};
    hsize_t maxdims_index_to_KEnergy[rank] = {max_node_tag};
    hid_t id_dataspace_index_to_KEnergy = H5Screate_simple(rank, dims_index_to_KEnergy, maxdims_index_to_KEnergy);

    hsize_t dims_index_to_SEnergy[rank] = {max_element_tag};
    hsize_t maxdims_index_to_SEnergy[rank] = {max_element_tag};
    hid_t id_dataspace_index_to_SEnergy = H5Screate_simple(rank, dims_index_to_SEnergy, maxdims_index_to_SEnergy);

    hsize_t dims_index_to_Energy_avg[rank] = {max_element_tag};
    hsize_t maxdims_index_to_Energy_avg[rank] = {max_element_tag};
    hid_t id_dataspace_index_to_Energy_avg = H5Screate_simple(rank, dims_index_to_Energy_avg, maxdims_index_to_Energy_avg);

    hsize_t dims_total_mechanical_energy[rank] = {number_of_timesteps};
    hsize_t maxdims_total_mechanical_energy[rank] = {number_of_timesteps};
    hid_t id_dataspace_total_mechanical_energy = H5Screate_simple(rank, dims_total_mechanical_energy, maxdims_total_mechanical_energy);


    // Create "Energy_Post_Processing" group
    herr_t status = H5Gget_objinfo (id_file, "Energy_Post_Processing", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    The group either does NOT exist or some other error occurred.\n";
        hid_t energy_density_group = H5Gcreate(id_file, "Energy_Post_Processing", H5P_DEFAULT, H5P_DEFAULT , H5P_DEFAULT);
        //Energy_Out << "    The group is successfully created!\n";     
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing\" exists.\n";
    }


    //Create "Energy_Density" group
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    The group either does NOT exist or some other error occurred.\n";
        hid_t energy_density_group = H5Gcreate(id_file, "Energy_Post_Processing/Energy_Density", H5P_DEFAULT, H5P_DEFAULT , H5P_DEFAULT);
        //Energy_Out << "    The group is successfully created!\n";     
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density\" exists.\n";
    }
    
    //Create "Kinetic_Energy_Density" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density/Kinetic_Energy_Density", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density/Kinetic_Energy_Density\" either does NOT exist or some other error occurred.\n";
        id_KEnergy_dataset = H5Dcreate2(id_file, "Energy_Post_Processing/Energy_Density/Kinetic_Energy_Density", H5T_NATIVE_DOUBLE, id_dataspace_KEnergy, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density/Kinetic_Energy_Density\" exists.\n";
        id_KEnergy_dataset = H5Dopen(id_file, "Energy_Post_Processing/Energy_Density/Kinetic_Energy_Density", H5P_DEFAULT);
    }

    //Create "Strain_Energy_Density_Step" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density/Strain_Energy_Density_Step", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density/Strain_Energy_Density\" either does NOT exist or some other error occurred.\n";
        id_SEnergy_step_dataset = H5Dcreate2(id_file, "Energy_Post_Processing/Energy_Density/Strain_Energy_Density_Step", H5T_NATIVE_DOUBLE, id_dataspace_SEnergy, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density/Strain_Energy_Density_Step\" exists.\n";
        id_SEnergy_step_dataset = H5Dopen(id_file, "Energy_Post_Processing/Energy_Density/Strain_Energy_Density_Step", H5P_DEFAULT);
    }

    //Create "Strain_Energy_Density_Cumulative" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density/Strain_Energy_Density_Cumulative", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density/Strain_Energy_Density\" either does NOT exist or some other error occurred.\n";
        id_SEnergy_cumu_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Energy_Density/Strain_Energy_Density_Cumulative", H5T_NATIVE_DOUBLE, id_dataspace_SEnergy, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density/Strain_Energy_Density_Cumulative\" exists.\n";
        id_SEnergy_cumu_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Energy_Density/Strain_Energy_Density_Cumulative", H5P_DEFAULT);
    }

    //Create "Plastic_Work_Density_Step" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density/Plastic_Work_Density_Step", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density/Plastic_Work_Density_Step\" either does NOT exist or some other error occurred.\n";
        id_PWork_step_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Energy_Density/Plastic_Work_Density_Step", H5T_NATIVE_DOUBLE, id_dataspace_PWork, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density/Plastic_Work_Density_Step\" exists.\n";
        id_PWork_step_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Energy_Density/Plastic_Work_Density_Step", H5P_DEFAULT);
    }

    //Create "Plastic_Work_Density_Cumulative" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density/Plastic_Work_Density_Cumulative", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density/Plastic_Work_Density_Cumulative\" either does NOT exist or some other error occurred.\n";
        id_PWork_cumu_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Energy_Density/Plastic_Work_Density_Cumulative", H5T_NATIVE_DOUBLE, id_dataspace_PWork, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density/Plastic_Work_Density_Cumulative\" exists.\n";
        id_PWork_cumu_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Energy_Density/Plastic_Work_Density_Cumulative", H5P_DEFAULT);
    }

    //Create "Plastic_Free_Energy_Density_Step" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density/Plastic_Free_Energy_Density_Step", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density/Plastic_Free_Energy_Density_Step\" either does NOT exist or some other error occurred.\n";
        id_plastic_free_energy_density_step_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Energy_Density/Plastic_Free_Energy_Density_Step", H5T_NATIVE_DOUBLE, id_dataspace_PWork, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density/Plastic_Free_Energy_Density_Step\" exists.\n";
        id_plastic_free_energy_density_step_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Energy_Density/Plastic_Free_Energy_Density_Step", H5P_DEFAULT);
    }

    //Create "Plastic_Free_Energy_Density_Cumulative" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density/Plastic_Free_Energy_Density_Cumulative", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density/Plastic_Free_Energy_Density_Cumulative\" either does NOT exist or some other error occurred.\n";
        id_plastic_free_energy_density_cumu_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Energy_Density/Plastic_Free_Energy_Density_Cumulative", H5T_NATIVE_DOUBLE, id_dataspace_PWork, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density/Plastic_Free_Energy_Density_Cumulative\" exists.\n";
        id_plastic_free_energy_density_cumu_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Energy_Density/Plastic_Free_Energy_Density_Cumulative", H5P_DEFAULT);
    }

    //Create "Plastic_Dissipation_Density_Step" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density/Plastic_Dissipation_Density_Step", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density/Plastic_Dissipation_Density_Step\" either does NOT exist or some other error occurred.\n";
        id_plastic_dissipation_density_step_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Energy_Density/Plastic_Dissipation_Density_Step", H5T_NATIVE_DOUBLE, id_dataspace_PWork, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density/Plastic_Dissipation_Density_Step\" exists.\n";
        id_plastic_dissipation_density_step_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Energy_Density/Plastic_Dissipation_Density_Step", H5P_DEFAULT);
    }

    //Create "Plastic_Dissipation_Density_Cumulative" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density/Plastic_Dissipation_Density_Cumulative", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density/Plastic_Dissipation_Density_Cumulative\" either does NOT exist or some other error occurred.\n";
        id_plastic_dissipation_density_cumu_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Energy_Density/Plastic_Dissipation_Density_Cumulative", H5T_NATIVE_DOUBLE, id_dataspace_PWork, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density/Plastic_Dissipation_Density_Cumulative\" exists.\n";
        id_plastic_dissipation_density_cumu_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Energy_Density/Plastic_Dissipation_Density_Cumulative", H5P_DEFAULT);
    }

    //Create "Index_to_Kinetic_Energy_Density" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density/Index_to_Kinetic_Energy_Density", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density/Index_to_Kinetic_Energy_Density\" either does NOT exist or some other error occurred.\n";
        id_index_to_KEnergy_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Energy_Density/Index_to_Kinetic_Energy_Density", H5T_NATIVE_INT, id_dataspace_index_to_KEnergy, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density/Index_to_Kinetic_Energy_Density\" exists.\n";
        id_index_to_KEnergy_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Energy_Density/Kinetic_Energy_Density", H5P_DEFAULT);
    }

    //Create "Index_to_Strain_Energy_Density" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density/Index_to_Strain_Energy_Density", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density/Index_to_Strain_Energy_Density\" either does NOT exist or some other error occurred.\n";
        id_index_to_SEnergy_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Energy_Density/Index_to_Strain_Energy_Density", H5T_NATIVE_INT, id_dataspace_index_to_SEnergy, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density/Index_to_Strain_Energy_Density\" exists.\n";
        id_index_to_SEnergy_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Energy_Density/Index_to_Strain_Energy_Density", H5P_DEFAULT);
    }



    //Create "Energy_Density_Element_Average" group
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density_Element_Average", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    The group either does NOT exist or some other error occurred.\n";
        hid_t energy_density_element_average_group = H5Gcreate(id_file, "Energy_Post_Processing/Energy_Density_Element_Average", H5P_DEFAULT, H5P_DEFAULT , H5P_DEFAULT);
        //Energy_Out << "    The group is successfully created!\n";    
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density_Element_Average\" exists.\n";
    }

    //Create "Plastic_Work_Density_Element_Average" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Work_Density_Element_Average", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Work_Density_Element_Average\" either does NOT exist or some other error occurred.\n";
        id_PWork_avg_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Work_Density_Element_Average", H5T_NATIVE_DOUBLE, id_dataspace_PWork_avg, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Work_Density_Element_Average\" exists.\n";
        id_PWork_avg_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Work_Density_Element_Average", H5P_DEFAULT);
    }

    //Create "Plastic_Free_Energy_Density_Element_Average" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Free_Energy_Density_Element_Average", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Free_Energy_Density_Element_Average\" either does NOT exist or some other error occurred.\n";
        id_plastic_free_energy_avg_dataset = H5Dcreate2(id_file, "Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Free_Energy_Density_Element_Average", H5T_NATIVE_DOUBLE, id_dataspace_PWork_avg, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Free_Energy_Density_Element_Average\" exists.\n";
        id_plastic_free_energy_avg_dataset = H5Dopen(id_file, "Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Free_Energy_Density_Element_Average", H5P_DEFAULT);
    }

    //Create "Plastic_Dissipation_Density_Element_Average" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Dissipation_Density_Element_Average", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Dissipation_Density_Element_Average\" either does NOT exist or some other error occurred.\n";
        id_plastic_dissipation_avg_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Dissipation_Density_Element_Average", H5T_NATIVE_DOUBLE, id_dataspace_PWork_avg, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Dissipation_Density_Element_Average\" exists.\n";
        id_plastic_dissipation_avg_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Energy_Density_Element_Average/Plastic_Dissipation_Density_Element_Average", H5P_DEFAULT);
    }

    //Create "Energy_Post_Processing/Mechanical_Energy_Density_Element_Average" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density_Element_Average/Mechanical_Energy_Density_Element_Average", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density_Element_Average/Mechanical_Energy_Density_Element_Average\" either does NOT exist or some other error occurred.\n";
        id_mechanical_energy_avg_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Energy_Density_Element_Average/Mechanical_Energy_Density_Element_Average", H5T_NATIVE_DOUBLE, id_dataspace_TEnergy_avg, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density_Element_Average/Mechanical_Energy_Density_Element_Average\" exists.\n";
        id_mechanical_energy_avg_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Energy_Density_Element_Average/Mechanical_Energy_Density_Element_Average", H5P_DEFAULT);
    }

    //Create "Energy_Post_Processing/Index_to_Energy_Density_Element_Average" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density_Element_Average/Index_to_Energy_Density_Element_Average", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Energy_Density_Element_Average/Index_to_Energy_Density_Element_Average\" either does NOT exist or some other error occurred.\n";
        id_index_to_Energy_avg_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Energy_Density_Element_Average/Index_to_Energy_Density_Element_Average", H5T_NATIVE_INT, id_dataspace_index_to_Energy_avg, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density_Element_Average/Index_to_Energy_Density_Element_Average\" exists.\n";
        id_index_to_Energy_avg_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Energy_Density_Element_Average/Index_to_Energy_Density_Element_Average", H5P_DEFAULT);
    }

/*
    //Create "Center_Coordinates" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Energy_Density/Center_Coordinates", 0, 0);
    if (status != 0)
    {
        Energy_Out << "    \"Energy_Post_Processing/Energy_Density/Center_Coordinates\" either does NOT exist or some other error occurred.\n";
        id_center_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Energy_Density/Center_Coordinates", H5T_NATIVE_DOUBLE, id_dataspace_center, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Energy_Density/Center_Coordinates\" exists.\n";
        id_center_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Energy_Density/Center_Coordinates", H5P_DEFAULT);
    }
*/

    //Create "Element_Volume" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Element_Volume", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Element_Volume\" either does NOT exist or some other error occurred.\n";
        id_element_volume_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Element_Volume", H5T_NATIVE_DOUBLE, id_dataspace_index_to_Energy_avg, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Element_Volume\" exists.\n";
        id_element_volume_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Element_Volume", H5P_DEFAULT);
    }

/*
    //Create "Back_Stress" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Back_Stress", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Back_Stress\" either does NOT exist or some other error occurred.\n";
        id_back_stress_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Back_Stress", H5T_NATIVE_DOUBLE, id_dataspace_back_stress, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Back_Stress\" exists.\n";
        id_back_stress_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Back_Stress", H5P_DEFAULT);
    }
*/

    //Create "Total_Energy" group
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Total_Energy", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    The group either does NOT exist or some other error occurred.\n";
        hid_t energy_density_group = H5Gcreate(id_file, "Energy_Post_Processing/Total_Energy", H5P_DEFAULT, H5P_DEFAULT , H5P_DEFAULT);
        //Energy_Out << "    The group is successfully created!\n";       
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Total_Energy\" exists.\n";
    }


    //Create "Total_Mechanical_Energy" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Total_Energy/Total_Mechanical_Energy", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Total_Energy/Total_Mechanical_Energy\" either does NOT exist or some other error occurred.\n";
        id_total_mechanical_energy_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Total_Energy/Total_Mechanical_Energy", H5T_NATIVE_DOUBLE, id_dataspace_total_mechanical_energy, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Total_Energy/Total_Mechanical_Energy\" exists.\n";
        id_total_mechanical_energy_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Total_Energy/Total_Mechanical_Energy", H5P_DEFAULT);
    }


    //Create "Total_Plastic_Work" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Total_Energy/Total_Plastic_Work", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Total_Energy/Total_Plastic_Work\" either does NOT exist or some other error occurred.\n";
        id_total_plastic_work_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Total_Energy/Total_Plastic_Work", H5T_NATIVE_DOUBLE, id_dataspace_total_mechanical_energy, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Total_Energy/Total_Plastic_Work\" exists.\n";
        id_total_plastic_work_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Total_Energy/Total_Plastic_Work", H5P_DEFAULT);
    }


    //Create "Total_Plastic_Free_Energy" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Total_Energy/Total_Plastic_Free_Energy", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Total_Energy/Total_Plastic_Free_Energy\" either does NOT exist or some other error occurred.\n";
        id_total_plastic_free_energy_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Total_Energy/Total_Plastic_Free_Energy", H5T_NATIVE_DOUBLE, id_dataspace_total_mechanical_energy, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Total_Energy/Total_Plastic_Free_Energy\" exists.\n";
        id_total_plastic_free_energy_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Total_Energy/Total_Plastic_Free_Energy", H5P_DEFAULT);
    }


    //Create "Total_Plastic_Dissipation" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Total_Energy/Total_Plastic_Dissipation", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Total_Energy/Total_Plastic_Dissipation\" either does NOT exist or some other error occurred.\n";
        id_total_plastic_dissipation_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Total_Energy/Total_Plastic_Dissipation", H5T_NATIVE_DOUBLE, id_dataspace_total_mechanical_energy, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Total_Energy/Total_Plastic_Dissipation\" exists.\n";
        id_total_plastic_dissipation_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Total_Energy/Total_Plastic_Dissipation", H5P_DEFAULT);
    }


    //Create "Total_Kinetic_Energy" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Total_Energy/Total_Kinetic_Energy", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Total_Energy/Total_Kinetic_Energy\" either does NOT exist or some other error occurred.\n";
        id_total_kinetic_energy_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Total_Energy/Total_Kinetic_Energy", H5T_NATIVE_DOUBLE, id_dataspace_total_mechanical_energy, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Total_Energy/Total_Kinetic_Energy\" exists.\n";
        id_total_kinetic_energy_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Total_Energy/Total_Kinetic_Energy", H5P_DEFAULT);
    }


    //Create "Total_Strain_Energy" dataset
    status = H5Gget_objinfo (id_file, "Energy_Post_Processing/Total_Energy/Total_Strain_Energy", 0, 0);
    if (status != 0)
    {
        //Energy_Out << "    \"Energy_Post_Processing/Total_Energy/Total_Strain_Energy\" either does NOT exist or some other error occurred.\n";
        id_total_strain_energy_dataset = H5Dcreate2(id_file, "/Energy_Post_Processing/Total_Energy/Total_Strain_Energy", H5T_NATIVE_DOUBLE, id_dataspace_total_mechanical_energy, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //Energy_Out << "    The dataset is successfully created!\n";
    }
    else
    {
        Energy_Warning << "    \"Energy_Post_Processing/Total_Energy/Total_Strain_Energy\" exists.\n";
        id_total_strain_energy_dataset = H5Dopen(id_file, "/Energy_Post_Processing/Total_Energy/Total_Strain_Energy", H5P_DEFAULT);
    }
}



void Energy_Post_Processing::computeKEnergyDensity(int node)
{
    if (index_GDs[node] == -1)
    {
        index_kinetic_energy_density[node] = -1;
        return;
    }
    else if (index_GDs[node]%3 == 0)
    {
        index_kinetic_energy_density[node] = index_GDs[node] / 3;
    }
    else
    {
        Energy_Error << "Something is wrong with \"Index_to_Generalized_Displacement\"...";
        return;
    }
    
    
    //Get displacement
    double disp[3][number_of_timesteps];

    int index_disp = index_GDs[node];
    //Energy_Out << "Index_GDs of node " << node << " is " << index_disp << "\n";

    hsize_t start[2]  = {index_disp, 0};
    hsize_t stride[2] = {3, number_of_timesteps};
    hsize_t count[2]  = {1, 1};
    hsize_t block[2]  = {3, number_of_timesteps};

    H5Sselect_hyperslab(id_displacements_dataspace, H5S_SELECT_SET, start, stride, count, block);

    hsize_t data_dims[2] = {3, number_of_timesteps};
    hsize_t data_maxdims[2] = {3, number_of_timesteps};

    hid_t id_memspace_disp = H5Screate_simple(2, data_dims, data_maxdims);

    H5Dread(id_displacements, H5T_NATIVE_DOUBLE, id_memspace_disp, id_displacements_dataspace, H5P_DEFAULT, disp);
    //Energy_Out << "GDs of node " << node << " is successfully read.\n";

    //Temporary time step for testing -- will be commented later
    //int temp_timestep = 300;
    //Energy_Out << "At time step " << temp_timestep << ", the GDs of node " << node << " are (" << disp[0][temp_timestep] << ", " << disp[1][temp_timestep] << ", " << disp[2][temp_timestep] << ")\n";
 

    //Calculate velocity
    double vel[3][number_of_timesteps];

    vel[0][0] = vel[1][0] = vel[2][0] = 0.0;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 1; j < number_of_timesteps; j++)
        {
            vel[i][j] = (disp[i][j] - disp[i][j-1]) / delta_t;
        }
    }
    //Energy_Out << "At time step " << temp_timestep << ", the velocity components of node " << node << " are (" << vel[0][temp_timestep] << ", " << vel[1][temp_timestep] << ", " << vel[2][temp_timestep] << ")\n";
  

    //Calculate kinetic energy density   
    kinetic_energy_density[index_kinetic_energy_density[node]*number_of_timesteps] = 0;

    for (int j = 1; j < number_of_timesteps; j++)
    {
        kinetic_energy_density[index_kinetic_energy_density[node]*number_of_timesteps+j] = 0.5 * mass_density * (vel[0][j]*vel[0][j] + vel[1][j]*vel[1][j] + vel[2][j]*vel[2][j]);
    }

}



void Energy_Post_Processing::computeSEnergyDensity(int element)
{
    if (number_of_GPs_element[element] < 0)
    {
        index_strain_energy_density[element] = -1;
        index_energy_density_eleavg[element] = -1;
        //Energy_Out << element << " " << number_of_GPs_element[element] << " " << index_strain_energy_density[element] << endl;
        return;
    }
    else if (is_first_valid_element == true)
    {
        index_strain_energy_density[element] = 0;
        index_energy_density_eleavg[element] = 0;
        is_first_valid_element = false;
        last_valid_element = element;
    }
    else
    {
        index_strain_energy_density[element] = index_strain_energy_density[last_valid_element] + number_of_GPs_element[element];
        index_energy_density_eleavg[element] = index_energy_density_eleavg[last_valid_element] + 1;
        last_valid_element = element;
    }
    //Energy_Out << element << " " << number_of_GPs_element[element] << " " << index_strain_energy_density[element] << endl;
    //Energy_Out << element << " " << index_strain_energy_density[element] << endl;

    double output_of_element[number_of_GPs_element[element]*SOLID_ELEMENT_GP_OUTPUT_SIZE][number_of_timesteps];

    int index_output_of_element = index_to_outputs[element];  
    //Energy_Out << element << " " << e << " " << index_strain_energy_density[element] << " " << index_output_of_element << endl;

    hsize_t start[2]  = {index_output_of_element, 0};
    hsize_t stride[2] = {number_of_GPs_element[element]*SOLID_ELEMENT_GP_OUTPUT_SIZE, number_of_timesteps};
    hsize_t count[2]  = {1, 1};
    hsize_t block[2]  = {number_of_GPs_element[element]*SOLID_ELEMENT_GP_OUTPUT_SIZE, number_of_timesteps};

    H5Sselect_hyperslab(id_outputs_dataspace, H5S_SELECT_SET, start, stride, count, block);    

    hsize_t data_dims[2] = {number_of_GPs_element[element]*SOLID_ELEMENT_GP_OUTPUT_SIZE, number_of_timesteps};
    hsize_t data_maxdims[2] = {number_of_GPs_element[element]*SOLID_ELEMENT_GP_OUTPUT_SIZE, number_of_timesteps};    

    hid_t id_output_memspace = H5Screate_simple(2, data_dims, data_maxdims);    

    H5Dread(id_outputs, H5T_NATIVE_DOUBLE, id_output_memspace, id_outputs_dataspace, H5P_DEFAULT, output_of_element);

    //Energy_Out << number_of_GPs << endl;
    //Energy_Out << element << " " << e << " " << output_of_element[0][0] << " " << output_of_element[17][0] << endl;


    for (int e = 0; e < number_of_GPs_element[element]; e++)
    { 
/*        // Calculate stress, elastic strain, plastic strain, deviatoric plastic strain
        double stress[6][number_of_timesteps];
        double elastic_strain[6][number_of_timesteps];
        double plastic_strain[6][number_of_timesteps];
        double deviatoric_plastic_strain[6][number_of_timesteps];
        double mean_plastic_strain[number_of_timesteps];

        for (int i = 0; i < number_of_timesteps; i++)
        {
            mean_plastic_strain[i] = 0;

            for (int j = 0; j < 3; j++)
            {
                stress[j][i] = output_of_element[SOLID_ELEMENT_GP_OUTPUT_SIZE*e+12+j][i];
                elastic_strain[j][i] = output_of_element[SOLID_ELEMENT_GP_OUTPUT_SIZE*e+j][i] - output_of_element[SOLID_ELEMENT_GP_OUTPUT_SIZE*e+j+6][i];
                plastic_strain[j][i] = output_of_element[SOLID_ELEMENT_GP_OUTPUT_SIZE*e+j+6][i];
                mean_plastic_strain[i] += plastic_strain[j][i]/3.0;
            }

            for (int j = 0; j < 3; j++)
            {
                deviatoric_plastic_strain[j][i] = plastic_strain[j][i] - mean_plastic_strain[i];
            }

            for (int j = 3; j < 6; j++)
            {
                stress[j][i] = output_of_element[SOLID_ELEMENT_GP_OUTPUT_SIZE*e+j+12][i];
                elastic_strain[j][i] = output_of_element[SOLID_ELEMENT_GP_OUTPUT_SIZE*e+j][i] - output_of_element[SOLID_ELEMENT_GP_OUTPUT_SIZE*e+j+6][i];
                plastic_strain[j][i] = output_of_element[SOLID_ELEMENT_GP_OUTPUT_SIZE*e+j+6][i];
                deviatoric_plastic_strain[j][i] = plastic_strain[j][i];
            }
        }

        Energy_Out << "Plastic strain of GP #: " << e << "\n";
        for (int j = 0; j < 6; j++)
        {
            Energy_Out << plastic_strain[j][100] << "\n";
        }

        Energy_Out << "Mean plastic strain of GP #: " << e << "\n";
        Energy_Out << mean_plastic_strain[100] << "\n";

        Energy_Out << "Deviatoric plastic strain of GP #: " << e << "\n";
        for (int j = 0; j < 6; j++)
        {
            Energy_Out << deviatoric_plastic_strain[j][100] << "\n";
        }


        double equivalent_plastic_strain[number_of_timesteps];

        for (int i = 0; i < number_of_timesteps; i++)
        {
            equivalent_plastic_strain[i] = 0;

            for (int j = 0; j < 3; j++)
            {
                equivalent_plastic_strain[i] += deviatoric_plastic_strain[j][i]*deviatoric_plastic_strain[j][i];
            }

            for (int j = 3; j < 6; j++)
            {
                equivalent_plastic_strain[i] += deviatoric_plastic_strain[j][i]*deviatoric_plastic_strain[j][i]*2;
            }

            equivalent_plastic_strain[i] = sqrt(equivalent_plastic_strain[i]*2.0/3);
        }

        //Energy_Out << "Equivalent plastic strain of GP #: " << e << "\n";
        //Energy_Out << equivalent_plastic_strain[100] << "\n";
*/

        //===========================================================================
        // Get plastic strain, elastic stran, stress, and back stress
        //===========================================================================

        double stress[6][number_of_timesteps];
        double back_stress[6][number_of_timesteps];
        double elastic_strain[6][number_of_timesteps];
        double plastic_strain[6][number_of_timesteps];

        for (int i = 0; i < number_of_timesteps; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                stress[j][i] = output_of_element[SOLID_ELEMENT_GP_OUTPUT_SIZE*e+12+j][i];
                back_stress[j][i] = output_of_element[SOLID_ELEMENT_GP_OUTPUT_SIZE*e+18+j][i];
                elastic_strain[j][i] = output_of_element[SOLID_ELEMENT_GP_OUTPUT_SIZE*e+j][i] - output_of_element[SOLID_ELEMENT_GP_OUTPUT_SIZE*e+j+6][i];
                plastic_strain[j][i] = output_of_element[SOLID_ELEMENT_GP_OUTPUT_SIZE*e+j+6][i];
            }
        }
        

        //Calculate strain energy density

        strain_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps] = 
            elastic_strain[0][0]*stress[0][0]/2 +
            elastic_strain[1][0]*stress[1][0]/2 +
            elastic_strain[2][0]*stress[2][0]/2 +
            elastic_strain[3][0]*stress[3][0] +
            elastic_strain[4][0]*stress[4][0] +
            elastic_strain[5][0]*stress[5][0];

        strain_energy_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps] = strain_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps];

        for (int j = 1; j < number_of_timesteps; j++)
        {
            //Energy_Out << (index_strain_energy_density[element]+e)*number_of_timesteps+j << endl;
            strain_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps+j] = 
                (elastic_strain[0][j]-elastic_strain[0][j-1])*(stress[0][j]+stress[0][j-1])/2 +
                (elastic_strain[1][j]-elastic_strain[1][j-1])*(stress[1][j]+stress[1][j-1])/2 +
                (elastic_strain[2][j]-elastic_strain[2][j-1])*(stress[2][j]+stress[2][j-1])/2 +
                (elastic_strain[3][j]-elastic_strain[3][j-1])*(stress[3][j]+stress[3][j-1]) +
                (elastic_strain[4][j]-elastic_strain[4][j-1])*(stress[4][j]+stress[4][j-1]) +
                (elastic_strain[5][j]-elastic_strain[5][j-1])*(stress[5][j]+stress[5][j-1]);

            strain_energy_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps+j] = 
                strain_energy_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps+j-1] + 
                strain_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps+j];
        }



        //Calculate plastic work density

        plastic_work_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps] = 
            plastic_strain[0][0]*stress[0][0]/2 +
            plastic_strain[1][0]*stress[1][0]/2 +
            plastic_strain[2][0]*stress[2][0]/2 +
            plastic_strain[3][0]*stress[3][0] +
            plastic_strain[4][0]*stress[4][0] +
            plastic_strain[5][0]*stress[5][0];

        plastic_work_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps] = plastic_work_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps];

        for (int j = 1; j < number_of_timesteps; j++)
        {
            plastic_work_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps+j] = 
                (plastic_strain[0][j]-plastic_strain[0][j-1])*(stress[0][j]+stress[0][j-1])/2 +
                (plastic_strain[1][j]-plastic_strain[1][j-1])*(stress[1][j]+stress[1][j-1])/2 +
                (plastic_strain[2][j]-plastic_strain[2][j-1])*(stress[2][j]+stress[2][j-1])/2 +
                (plastic_strain[3][j]-plastic_strain[3][j-1])*(stress[3][j]+stress[3][j-1]) +
                (plastic_strain[4][j]-plastic_strain[4][j-1])*(stress[4][j]+stress[4][j-1]) +
                (plastic_strain[5][j]-plastic_strain[5][j-1])*(stress[5][j]+stress[5][j-1]);

            plastic_work_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps+j] = 
                plastic_work_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps+j-1] + 
                plastic_work_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps+j];
        }

/*
        //===========================================================================
        // Compute plastic free energy
        // This is the old algorithm, not using any more...
        //===========================================================================

        // First step
        plastic_free_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps] = a1*(
            plastic_strain[0][0]*stress[0][0]/2 +
            plastic_strain[1][0]*stress[1][0]/2 +
            plastic_strain[2][0]*stress[2][0]/2 +
            plastic_strain[3][0]*stress[3][0] +
            plastic_strain[4][0]*stress[4][0] +
            plastic_strain[5][0]*stress[5][0]);

        plastic_free_energy_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps] = plastic_free_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps];

        // Following steps
        if (is_linear_KH == true)
        {
            for (int j = 1; j < number_of_timesteps; j++)
            {
                for (int k = 0; k < 6; k++)
                {
                    back_stress[(6*(index_strain_energy_density[element]+e)+k)*number_of_timesteps+j] = a1*plastic_strain[k][j];
                    //Energy_Out << "back_stress at # " << 6*((index_strain_energy_density[element]+e)*number_of_timesteps+j)+k << " is: " << back_stress[6*((index_strain_energy_density[element]+e)*number_of_timesteps+j)+k] << endl;
                }
                //Energy_Out << "back_stress at # " << 6*((index_strain_energy_density[element]+e)*number_of_timesteps+j)+4 << " is: " << back_stress[6*((index_strain_energy_density[element]+e)*number_of_timesteps+j)+4] << endl;
                

                plastic_free_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps+j] = 
                    computePlasticFreeEnergy_Linear_KH((double *)plastic_strain, j);

                plastic_free_energy_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps+j] = 
                    plastic_free_energy_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps+j-1] + 
                    plastic_free_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps+j];
            }
        }
        else if (is_AF_KH == true)
        {
            for (int j = 1; j < number_of_timesteps; j++)
            {
                plastic_free_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps+j] =
                    computePlasticFreeEnergy_AF_KH((double *)equivalent_plastic_strain, (double *)back_stress, (double *)deviatoric_plastic_strain, j);

                plastic_free_energy_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps+j] = 
                    plastic_free_energy_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps+j-1] + 
                    plastic_free_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps+j];
            }
        }
        else
        {
            Energy_Error << "Something is wrong with hardening settings!\n";
            return;
        }
*/        


        //===========================================================================
        // Compute plastic free energy density
        // This is the new algorithm!!! (Aug 2016)
        //===========================================================================

        // First step
        plastic_free_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps] = 1.0/a1*(
            back_stress[0][0]*back_stress[0][0]/2 +
            back_stress[1][0]*back_stress[1][0]/2 +
            back_stress[2][0]*back_stress[2][0]/2 +
            back_stress[3][0]*back_stress[3][0] +
            back_stress[4][0]*back_stress[4][0] +
            back_stress[5][0]*back_stress[5][0]);

        plastic_free_energy_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps] = plastic_free_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps];

        for (int j = 1; j < number_of_timesteps; j++)
        {
            plastic_free_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps+j] = 1.0/a1*(
                (back_stress[0][j]-back_stress[0][j-1])*(back_stress[0][j]+back_stress[0][j-1])/2 +
                (back_stress[1][j]-back_stress[1][j-1])*(back_stress[1][j]+back_stress[1][j-1])/2 +
                (back_stress[2][j]-back_stress[2][j-1])*(back_stress[2][j]+back_stress[2][j-1])/2 +
                (back_stress[3][j]-back_stress[3][j-1])*(back_stress[3][j]+back_stress[3][j-1]) +
                (back_stress[4][j]-back_stress[4][j-1])*(back_stress[4][j]+back_stress[4][j-1]) +
                (back_stress[5][j]-back_stress[5][j-1])*(back_stress[5][j]+back_stress[5][j-1]));

            plastic_free_energy_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps+j] = 
                plastic_free_energy_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps+j-1] + 
                plastic_free_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps+j];
        }



        //===========================================================================
        // Compute plastic dissipation density
        //===========================================================================

        for (int j = 0; j < number_of_timesteps; j++)
        {
            plastic_dissipation_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps+j] = 
                plastic_work_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps+j] -
                plastic_free_energy_density_step[(index_strain_energy_density[element]+e)*number_of_timesteps+j];

            plastic_dissipation_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps+j] = 
                plastic_work_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps+j] -
                plastic_free_energy_density_cumu[(index_strain_energy_density[element]+e)*number_of_timesteps+j];
        }
    }

}



//===========================================================================
// Compute element average of energy density (8NodeBrick)
//===========================================================================
void Energy_Post_Processing::computeElementAverage_8NodeBrick(int element)
{
    //Energy_Out << "Computing average energy density of element " << element << endl;
    for (int j = 0; j < number_of_timesteps; j++)
    {
        int tag_eleavg = index_energy_density_eleavg[element]*number_of_timesteps+j;

        kinetic_energy_density_eleavg[tag_eleavg] = 0;

        for (int i = 0; i < 8; i++)
        {
            kinetic_energy_density_eleavg[tag_eleavg] += kinetic_energy_density[index_kinetic_energy_density[connectivity[index_to_connectivity[element]+i]]*number_of_timesteps+j];
        }

        kinetic_energy_density_eleavg[tag_eleavg] = kinetic_energy_density_eleavg[tag_eleavg] / 8;
        //Energy_Out << "Successfully computed kinetic energy density of element " << element << " at time step " << j << endl;



        strain_energy_density_eleavg[tag_eleavg] = 0;
        
        for (int i = 0; i < 8; i++)
        {
            strain_energy_density_eleavg[tag_eleavg] += strain_energy_density_cumu[(index_strain_energy_density[element]+i)*number_of_timesteps+j];
        }

        strain_energy_density_eleavg[tag_eleavg] = strain_energy_density_eleavg[tag_eleavg] / 8;
        //Energy_Out << "Successfully computed strain energy density of element " << element << " at time step " << j << endl;



        plastic_dissipation_density_eleavg[tag_eleavg] = 0;

        for (int i = 0; i < 8; i++)
        {
            plastic_dissipation_density_eleavg[tag_eleavg] += plastic_dissipation_density_cumu[(index_strain_energy_density[element]+i)*number_of_timesteps+j];
        }

        plastic_dissipation_density_eleavg[tag_eleavg] = plastic_dissipation_density_eleavg[tag_eleavg] / 8;



        plastic_free_energy_density_eleavg[tag_eleavg] = 0;

        for (int i = 0; i < 8; i++)
        {
            plastic_free_energy_density_eleavg[tag_eleavg] += plastic_free_energy_density_cumu[(index_strain_energy_density[element]+i)*number_of_timesteps+j];
        }

        plastic_free_energy_density_eleavg[tag_eleavg] = plastic_free_energy_density_eleavg[tag_eleavg] / 8;



        plastic_work_density_eleavg[tag_eleavg] = 0;

        for (int i = 0; i < 8; i++)
        {
            plastic_work_density_eleavg[tag_eleavg] += plastic_work_density_cumu[(index_strain_energy_density[element]+i)*number_of_timesteps+j];
        }

        plastic_work_density_eleavg[tag_eleavg] = plastic_work_density_eleavg[tag_eleavg] / 8;


        mechanical_energy_density_eleavg[tag_eleavg] = kinetic_energy_density_eleavg[tag_eleavg] + strain_energy_density_eleavg[tag_eleavg] + plastic_free_energy_density_eleavg[tag_eleavg];
    }
}



//===========================================================================
// Compute element average of energy density (27NodeBrick)
//===========================================================================
void Energy_Post_Processing::computeElementAverage_27NodeBrick(int element)
{
    for (int j = 0; j < number_of_timesteps; j++)
    {
        int tag_eleavg = index_energy_density_eleavg[element]*number_of_timesteps+j;

        kinetic_energy_density_eleavg[tag_eleavg] = 0;

        for (int i = 0; i < number_of_nodes_element[element]; i++)
        {
            kinetic_energy_density_eleavg[tag_eleavg] += kinetic_energy_density[index_kinetic_energy_density[connectivity[index_to_connectivity[element]+i]]*number_of_timesteps+j];
        }

        kinetic_energy_density_eleavg[tag_eleavg] = kinetic_energy_density_eleavg[tag_eleavg] / number_of_nodes_element[element];



        strain_energy_density_eleavg[tag_eleavg] = 0;

        for (int i = 0; i < number_of_GPs_element[element]; i++)
        {
            strain_energy_density_eleavg[tag_eleavg] += strain_energy_density_cumu[(index_strain_energy_density[element]+i)*number_of_timesteps+j]*weight_27nodebrick[i];
        }

        strain_energy_density_eleavg[tag_eleavg] = strain_energy_density_eleavg[tag_eleavg] / 8;



        plastic_work_density_eleavg[tag_eleavg] = 0;

        for (int i = 0; i < number_of_GPs_element[element]; i++)
        {
            plastic_work_density_eleavg[tag_eleavg] += plastic_work_density_cumu[(index_strain_energy_density[element]+i)*number_of_timesteps+j]*weight_27nodebrick[i];
        }

        plastic_work_density_eleavg[tag_eleavg] = plastic_work_density_eleavg[tag_eleavg] / 8;



        plastic_free_energy_density_eleavg[tag_eleavg] = 0;

        for (int i = 0; i < number_of_GPs_element[element]; i++)
        {
            plastic_free_energy_density_eleavg[tag_eleavg] += plastic_free_energy_density_cumu[(index_strain_energy_density[element]+i)*number_of_timesteps+j]*weight_27nodebrick[i];
        }

        plastic_free_energy_density_eleavg[tag_eleavg] = plastic_free_energy_density_eleavg[tag_eleavg] / 8;



        plastic_dissipation_density_eleavg[tag_eleavg] = 0;

        for (int i = 0; i < number_of_GPs_element[element]; i++)
        {
            plastic_dissipation_density_eleavg[tag_eleavg] += plastic_dissipation_density_cumu[(index_strain_energy_density[element]+i)*number_of_timesteps+j]*weight_27nodebrick[i];
        }

        plastic_dissipation_density_eleavg[tag_eleavg] = plastic_dissipation_density_eleavg[tag_eleavg] / 8;



        mechanical_energy_density_eleavg[tag_eleavg] = kinetic_energy_density_eleavg[tag_eleavg] + strain_energy_density_eleavg[tag_eleavg];
    }
}


/*
//===========================================================================
// Compute center coordinates of certain element
// Currently not useful...
//===========================================================================
void Energy_Post_Processing::getCenterofElement(int element)
{
    //Simplification of index to elements, NOT GENERAL!!!
    //element = element-1;

    int num_nodes_ele = number_of_nodes_element[element];
    //Energy_Out << "num_nodes_ele: " << num_nodes_ele << "\n";
  

    center_coordinates[3*index_energy_density_eleavg[element]+0] = 0;
    center_coordinates[3*index_energy_density_eleavg[element]+1] = 0;
    center_coordinates[3*index_energy_density_eleavg[element]+2] = 0;

    for (int i = 0; i < num_nodes_ele; i++)
    {
        int node_tag = connectivity[index_to_connectivity[element]+i];
        //Energy_Out << "node_tag: " << node_tag << "\n";

        center_coordinates[3*index_energy_density_eleavg[element]+0] += node_coordinates[index_to_node_coordinates[node_tag]+0];
        center_coordinates[3*index_energy_density_eleavg[element]+1] += node_coordinates[index_to_node_coordinates[node_tag]+1];
        center_coordinates[3*index_energy_density_eleavg[element]+2] += node_coordinates[index_to_node_coordinates[node_tag]+2];  
        //Energy_Out << "index_to_node_coordinates: " << index_to_node_coordinates[node_tag] << "\n";
        //Energy_Out << "node_coordinates: " << node_coordinates[index_to_node_coordinates[node_tag]+0] << " " << node_coordinates[index_to_node_coordinates[node_tag]+1] << " " << node_coordinates[index_to_node_coordinates[node_tag]+2] << "\n";     
        //Energy_Out << "center_coordinates: " << center_coordinates[3*index_energy_density_eleavg[element]+0] << " " << center_coordinates[3*index_energy_density_eleavg[element]+1] << " " << center_coordinates[3*index_energy_density_eleavg[element]+2] << "\n";     
    }
    
    center_coordinates[3*index_energy_density_eleavg[element]+0] = center_coordinates[3*index_energy_density_eleavg[element]+0] / num_nodes_ele;
    center_coordinates[3*index_energy_density_eleavg[element]+1] = center_coordinates[3*index_energy_density_eleavg[element]+1] / num_nodes_ele;
    center_coordinates[3*index_energy_density_eleavg[element]+2] = center_coordinates[3*index_energy_density_eleavg[element]+2] / num_nodes_ele;
    //Energy_Out << "center_coordinates: " << center_coordinates[3*index_energy_density_eleavg[element]+0] << " " << center_coordinates[3*index_energy_density_eleavg[element]+1] << " " << center_coordinates[3*index_energy_density_eleavg[element]+2] << "\n";     

}
*/


//===========================================================================
// Compute energy components of the entire model
//===========================================================================
void Energy_Post_Processing::computeTotalEnergy()
{
    for (int j = 0; j < number_of_timesteps; j++)
    {
        total_mechanical_energy[j] = 0;
        total_plastic_work[j] = 0;
        total_plastic_dissipation[j] = 0;
        total_kinetic_energy[j] = 0;
        total_strain_energy[j] = 0;

        for (int i = 0; i < max_element_tag; i++)
        {
            if (index_energy_density_eleavg[i] != -1)
            {
                total_mechanical_energy[j] += mechanical_energy_density_eleavg[index_energy_density_eleavg[i]*number_of_timesteps+j]*element_volume[i];
                total_plastic_work[j] += plastic_work_density_eleavg[index_energy_density_eleavg[i]*number_of_timesteps+j]*element_volume[i];
                total_plastic_free_energy[j] += plastic_free_energy_density_eleavg[index_energy_density_eleavg[i]*number_of_timesteps+j]*element_volume[i];
                total_plastic_dissipation[j] += plastic_dissipation_density_eleavg[index_energy_density_eleavg[i]*number_of_timesteps+j]*element_volume[i];
                total_kinetic_energy[j] += kinetic_energy_density_eleavg[index_energy_density_eleavg[i]*number_of_timesteps+j]*element_volume[i];
                total_strain_energy[j] += strain_energy_density_eleavg[index_energy_density_eleavg[i]*number_of_timesteps+j]*element_volume[i];
            }
            
        }
    }
}



//===========================================================================
// Compute volume of each element, ONLY FOR 8 NODE BRICK
// Use the algorithm given by J. Grandy (1997), more efficient than using Jacobian
//===========================================================================
void Energy_Post_Processing::computeElementVolume_8NodeBrick(int element)
{
    /*if (number_of_nodes_element[element] < 0)
    {
        //Energy_Out << element << " " << number_of_nodes_element[element] << endl;
        element_volume[element] = -1;
        return;
    }*/
    //Energy_Out << element << " " << number_of_nodes_element[element] << endl;


    // Get nodal coordinates (8 vertices) of an element
    double element_nodal_coordinates[24];

    //Energy_Out << "The nodal coordinates of element " << element << " is:\n";
    for (int i = 0; i < 8; i++)
    {
        int node_tag = connectivity[index_to_connectivity[element]+i];
        for (int j = 0; j < 3; j++)
        {
            element_nodal_coordinates[3*i+j] = node_coordinates[index_to_node_coordinates[node_tag]+j];
            //Energy_Out << element_nodal_coordinates[3*i+j] << " ";
        }
        //Energy_Out << endl;
    }


    // Compute element volume
    double entry_1[3], entry_2[3], entry_3[3];

    for (int i = 0; i < 3; i++)
    {
        entry_1[i] = element_nodal_coordinates[3*7+i] - element_nodal_coordinates[3*2+i] + element_nodal_coordinates[3*4+i] - element_nodal_coordinates[3*1+i];
        entry_2[i] = element_nodal_coordinates[3*7+i] - element_nodal_coordinates[3*0+i];
        entry_3[i] = element_nodal_coordinates[3*3+i] - element_nodal_coordinates[3*1+i];
    }

    double determinant_1 = entry_1[0]*entry_2[1]*entry_3[2] + 
                           entry_1[1]*entry_2[2]*entry_3[0] +
                           entry_1[2]*entry_2[0]*entry_3[1] -
                           entry_1[2]*entry_2[1]*entry_3[0] -
                           entry_1[1]*entry_2[0]*entry_3[2] -
                           entry_1[0]*entry_2[2]*entry_3[1];

    for (int i = 0; i < 3; i++)
    {
        entry_1[i] = element_nodal_coordinates[3*4+i] - element_nodal_coordinates[3*1+i];
        entry_2[i] = element_nodal_coordinates[3*7+i] - element_nodal_coordinates[3*0+i] + element_nodal_coordinates[3*6+i] - element_nodal_coordinates[3*1+i];
        entry_3[i] = element_nodal_coordinates[3*7+i] - element_nodal_coordinates[3*5+i];
    }

    double determinant_2 = entry_1[0]*entry_2[1]*entry_3[2] + 
                           entry_1[1]*entry_2[2]*entry_3[0] +
                           entry_1[2]*entry_2[0]*entry_3[1] -
                           entry_1[2]*entry_2[1]*entry_3[0] -
                           entry_1[1]*entry_2[0]*entry_3[2] -
                           entry_1[0]*entry_2[2]*entry_3[1];

    for (int i = 0; i < 3; i++)
    {
        entry_1[i] = element_nodal_coordinates[3*7+i] - element_nodal_coordinates[3*2+i];
        entry_2[i] = element_nodal_coordinates[3*6+i] - element_nodal_coordinates[3*1+i];
        entry_3[i] = element_nodal_coordinates[3*7+i] - element_nodal_coordinates[3*5+i] + element_nodal_coordinates[3*3+i] - element_nodal_coordinates[3*1+i];
    }

    double determinant_3 = entry_1[0]*entry_2[1]*entry_3[2] + 
                           entry_1[1]*entry_2[2]*entry_3[0] +
                           entry_1[2]*entry_2[0]*entry_3[1] -
                           entry_1[2]*entry_2[1]*entry_3[0] -
                           entry_1[1]*entry_2[0]*entry_3[2] -
                           entry_1[0]*entry_2[2]*entry_3[1];

    element_volume[element] = (determinant_1 + determinant_2 + determinant_3) / 12;

    if (element_volume[element] < 0)
    {
        element_volume[element] = -element_volume[element];
    }
}



//===========================================================================
// Compute volume of each element, ONLY FOR 27 NODE BRICK
// Use Jacobian directly
//===========================================================================
void Energy_Post_Processing::computeElementVolume_27NodeBrick(int element)
{
    /*if (number_of_nodes_element[element] < 0)
    {
        //Energy_Out << element << " " << number_of_nodes_element[element] << endl;
        element_volume[element] = -1;
        return;
    }*/
    //Energy_Out << element << " " << number_of_nodes_element[element] << endl;


    // Get nodal coordinates (27 nodes) of an element
    double node_global_coordinates[27][3];

    //Energy_Out << "The global nodal coordinates of element " << element << " is:\n";
    for (int i = 0; i < 27; i++)
    {
        int node_tag = connectivity[index_to_connectivity[element]+i];
        for (int j = 0; j < 3; j++)
        {
            node_global_coordinates[i][j] = node_coordinates[index_to_node_coordinates[node_tag]+j];
            //Energy_Out << node_global_coordinates[i][j] << " ";
        }
        //Energy_Out << endl;
    }

    
    // Compute Jacobian matrix at each Gauss point
    double J_GP[3][3][27];
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 27; k++)
                J_GP[i][j][k] = 0;
        }
    }

    // Partial derivative of the shape function
    double partial_N;

    for (int GP_tag = 0; GP_tag < 27; GP_tag++)
    {
        for (int column = 0; column < 3; column++)
        {
            for (int node_tag = 0; node_tag < 8; node_tag++)
            {
                partial_N = 1.0/8.0*(1 + 2*GP_local_coordinate_27nodebrick[GP_tag][0]*node_local_coordinate_27nodebrick[node_tag][0])*
                            GP_local_coordinate_27nodebrick[GP_tag][1]*(1 + GP_local_coordinate_27nodebrick[GP_tag][1]*node_local_coordinate_27nodebrick[node_tag][1])*
                            GP_local_coordinate_27nodebrick[GP_tag][2]*(1 + GP_local_coordinate_27nodebrick[GP_tag][2]*node_local_coordinate_27nodebrick[node_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][0]*node_local_coordinate_27nodebrick[node_tag][1]*node_local_coordinate_27nodebrick[node_tag][2];

                J_GP[0][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = 1.0/8.0*(1 + 2*GP_local_coordinate_27nodebrick[GP_tag][1]*node_local_coordinate_27nodebrick[node_tag][1])*
                            GP_local_coordinate_27nodebrick[GP_tag][2]*(1 + GP_local_coordinate_27nodebrick[GP_tag][2]*node_local_coordinate_27nodebrick[node_tag][2])*
                            GP_local_coordinate_27nodebrick[GP_tag][0]*(1 + GP_local_coordinate_27nodebrick[GP_tag][0]*node_local_coordinate_27nodebrick[node_tag][0])*
                            node_local_coordinate_27nodebrick[node_tag][0]*node_local_coordinate_27nodebrick[node_tag][1]*node_local_coordinate_27nodebrick[node_tag][2];

                J_GP[1][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = 1.0/8.0*(1 + 2*GP_local_coordinate_27nodebrick[GP_tag][2]*node_local_coordinate_27nodebrick[node_tag][2])*
                            GP_local_coordinate_27nodebrick[GP_tag][0]*(1 + GP_local_coordinate_27nodebrick[GP_tag][0]*node_local_coordinate_27nodebrick[node_tag][0])*
                            GP_local_coordinate_27nodebrick[GP_tag][1]*(1 + GP_local_coordinate_27nodebrick[GP_tag][1]*node_local_coordinate_27nodebrick[node_tag][1])*
                            node_local_coordinate_27nodebrick[node_tag][0]*node_local_coordinate_27nodebrick[node_tag][1]*node_local_coordinate_27nodebrick[node_tag][2];

                J_GP[2][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;
            }

            for (int node_tag = 8; node_tag < 15; node_tag = node_tag + 2)
            {
                partial_N = 1.0/4.0*(-2*GP_local_coordinate_27nodebrick[GP_tag][0])*
                            GP_local_coordinate_27nodebrick[GP_tag][1]*(1 + GP_local_coordinate_27nodebrick[GP_tag][1]*node_local_coordinate_27nodebrick[node_tag][1])*
                            GP_local_coordinate_27nodebrick[GP_tag][2]*(1 + GP_local_coordinate_27nodebrick[GP_tag][2]*node_local_coordinate_27nodebrick[node_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][1]*node_local_coordinate_27nodebrick[node_tag][2];

                J_GP[0][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = 1.0/4.0*(1 - GP_local_coordinate_27nodebrick[GP_tag][0]*GP_local_coordinate_27nodebrick[GP_tag][0])*
                            (1 + 2*GP_local_coordinate_27nodebrick[GP_tag][1]*node_local_coordinate_27nodebrick[node_tag][1])*
                            GP_local_coordinate_27nodebrick[GP_tag][2]*(1 + GP_local_coordinate_27nodebrick[GP_tag][2]*node_local_coordinate_27nodebrick[node_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][1]*node_local_coordinate_27nodebrick[node_tag][2];

                J_GP[1][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = 1.0/4.0*(1 - GP_local_coordinate_27nodebrick[GP_tag][0]*GP_local_coordinate_27nodebrick[GP_tag][0])*
                            GP_local_coordinate_27nodebrick[GP_tag][1]*(1 + GP_local_coordinate_27nodebrick[GP_tag][1]*node_local_coordinate_27nodebrick[node_tag][1])*
                            (1 + 2*GP_local_coordinate_27nodebrick[GP_tag][2]*node_local_coordinate_27nodebrick[node_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][1]*node_local_coordinate_27nodebrick[node_tag][2];

                J_GP[2][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;
            }

            for (int node_tag = 9; node_tag < 16; node_tag = node_tag + 2)
            {
                partial_N = 1.0/4.0*(1 + 2*GP_local_coordinate_27nodebrick[GP_tag][0]*node_local_coordinate_27nodebrick[node_tag][0])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][1]*GP_local_coordinate_27nodebrick[GP_tag][1])*
                            GP_local_coordinate_27nodebrick[GP_tag][2]*(1 + GP_local_coordinate_27nodebrick[GP_tag][2]*node_local_coordinate_27nodebrick[node_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][0]*node_local_coordinate_27nodebrick[node_tag][2];

                J_GP[0][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = 1.0/4.0*GP_local_coordinate_27nodebrick[GP_tag][0]*(1 + GP_local_coordinate_27nodebrick[GP_tag][0]*node_local_coordinate_27nodebrick[node_tag][0])*
                            (-2*GP_local_coordinate_27nodebrick[GP_tag][1])*
                            GP_local_coordinate_27nodebrick[GP_tag][2]*(1 + GP_local_coordinate_27nodebrick[GP_tag][2]*node_local_coordinate_27nodebrick[node_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][0]*node_local_coordinate_27nodebrick[node_tag][2];

                J_GP[1][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = 1.0/4.0*GP_local_coordinate_27nodebrick[GP_tag][0]*(1 + GP_local_coordinate_27nodebrick[GP_tag][0]*node_local_coordinate_27nodebrick[node_tag][0])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][1]*GP_local_coordinate_27nodebrick[GP_tag][1])*
                            (1 + 2*GP_local_coordinate_27nodebrick[GP_tag][2]*node_local_coordinate_27nodebrick[node_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][0]*node_local_coordinate_27nodebrick[node_tag][2];

                J_GP[2][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;
            }

            for (int node_tag = 16; node_tag < 20; node_tag++)
            {
                partial_N = 1.0/4.0*(1 + 2*GP_local_coordinate_27nodebrick[GP_tag][0]*node_local_coordinate_27nodebrick[node_tag][0])*
                            GP_local_coordinate_27nodebrick[GP_tag][1]*(1 + GP_local_coordinate_27nodebrick[GP_tag][1]*node_local_coordinate_27nodebrick[node_tag][1])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][2]*GP_local_coordinate_27nodebrick[GP_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][0]*node_local_coordinate_27nodebrick[node_tag][1];

                J_GP[0][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = 1.0/4.0*GP_local_coordinate_27nodebrick[GP_tag][0]*(1 + GP_local_coordinate_27nodebrick[GP_tag][0]*node_local_coordinate_27nodebrick[node_tag][0])*
                            (1 + 2*GP_local_coordinate_27nodebrick[GP_tag][1]*node_local_coordinate_27nodebrick[node_tag][1])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][2]*GP_local_coordinate_27nodebrick[GP_tag][2])*                       
                            node_local_coordinate_27nodebrick[node_tag][0]*node_local_coordinate_27nodebrick[node_tag][1];

                J_GP[1][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = 1.0/4.0*GP_local_coordinate_27nodebrick[GP_tag][0]*(1 + GP_local_coordinate_27nodebrick[GP_tag][0]*node_local_coordinate_27nodebrick[node_tag][0])*
                            GP_local_coordinate_27nodebrick[GP_tag][1]*(1 + GP_local_coordinate_27nodebrick[GP_tag][1]*node_local_coordinate_27nodebrick[node_tag][1])*
                            (-2*GP_local_coordinate_27nodebrick[GP_tag][2])*                       
                            node_local_coordinate_27nodebrick[node_tag][0]*node_local_coordinate_27nodebrick[node_tag][1];

                J_GP[2][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;
            }

            for (int node_tag = 20; node_tag < 21; node_tag++)
            {
                partial_N = (-2*GP_local_coordinate_27nodebrick[GP_tag][0])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][1]*GP_local_coordinate_27nodebrick[GP_tag][1])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][2]*GP_local_coordinate_27nodebrick[GP_tag][2]);

                J_GP[0][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = (1 - GP_local_coordinate_27nodebrick[GP_tag][0]*GP_local_coordinate_27nodebrick[GP_tag][0])*
                            (-2*GP_local_coordinate_27nodebrick[GP_tag][1])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][2]*GP_local_coordinate_27nodebrick[GP_tag][2]);

                J_GP[1][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = (1 - GP_local_coordinate_27nodebrick[GP_tag][0]*GP_local_coordinate_27nodebrick[GP_tag][0])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][1]*GP_local_coordinate_27nodebrick[GP_tag][1])*
                            (-2*GP_local_coordinate_27nodebrick[GP_tag][2]);

                J_GP[2][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;
            }

            for (int node_tag = 21; node_tag < 24; node_tag = node_tag + 2)
            {
                partial_N = 1.0/2.0*(-2*GP_local_coordinate_27nodebrick[GP_tag][0])*
                            GP_local_coordinate_27nodebrick[GP_tag][1]*(1 + GP_local_coordinate_27nodebrick[GP_tag][1]*node_local_coordinate_27nodebrick[node_tag][1])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][2]*GP_local_coordinate_27nodebrick[GP_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][1];

                J_GP[0][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = 1.0/2.0*(1 - GP_local_coordinate_27nodebrick[GP_tag][0]*GP_local_coordinate_27nodebrick[GP_tag][0])*
                            (1 + 2*GP_local_coordinate_27nodebrick[GP_tag][1]*node_local_coordinate_27nodebrick[node_tag][1])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][2]*GP_local_coordinate_27nodebrick[GP_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][1];

                J_GP[1][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = 1.0/2.0*(1 - GP_local_coordinate_27nodebrick[GP_tag][0]*GP_local_coordinate_27nodebrick[GP_tag][0])*
                            GP_local_coordinate_27nodebrick[GP_tag][1]*(1 + GP_local_coordinate_27nodebrick[GP_tag][1]*node_local_coordinate_27nodebrick[node_tag][1])*
                            (-2*GP_local_coordinate_27nodebrick[GP_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][1];

                J_GP[2][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;
            }

            for (int node_tag = 22; node_tag < 25; node_tag = node_tag + 2)
            {
                partial_N = 1.0/2.0*(1 + 2*GP_local_coordinate_27nodebrick[GP_tag][0]*node_local_coordinate_27nodebrick[node_tag][0])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][1]*GP_local_coordinate_27nodebrick[GP_tag][1])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][2]*GP_local_coordinate_27nodebrick[GP_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][0];

                J_GP[0][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = 1.0/2.0*GP_local_coordinate_27nodebrick[GP_tag][0]*(1 + GP_local_coordinate_27nodebrick[GP_tag][0]*node_local_coordinate_27nodebrick[node_tag][0])*
                            (-2*GP_local_coordinate_27nodebrick[GP_tag][1])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][2]*GP_local_coordinate_27nodebrick[GP_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][0];

                J_GP[1][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = 1.0/2.0*GP_local_coordinate_27nodebrick[GP_tag][0]*(1 + GP_local_coordinate_27nodebrick[GP_tag][0]*node_local_coordinate_27nodebrick[node_tag][0])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][1]*GP_local_coordinate_27nodebrick[GP_tag][1])*
                            (-2*GP_local_coordinate_27nodebrick[GP_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][0];

                J_GP[2][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;
            }

            for (int node_tag = 25; node_tag < 27; node_tag++)
            {
                partial_N = 1.0/2.0*(-2*GP_local_coordinate_27nodebrick[GP_tag][0])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][1]*GP_local_coordinate_27nodebrick[GP_tag][1])*
                            GP_local_coordinate_27nodebrick[GP_tag][2]*(1 + GP_local_coordinate_27nodebrick[GP_tag][2]*node_local_coordinate_27nodebrick[node_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][2];

                J_GP[0][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = 1.0/2.0*(1 - GP_local_coordinate_27nodebrick[GP_tag][0]*GP_local_coordinate_27nodebrick[GP_tag][0])*
                            (-2*GP_local_coordinate_27nodebrick[GP_tag][1])*
                            GP_local_coordinate_27nodebrick[GP_tag][2]*(1 + GP_local_coordinate_27nodebrick[GP_tag][2]*node_local_coordinate_27nodebrick[node_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][2];

                J_GP[1][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;

                partial_N = 1.0/2.0*(1 - GP_local_coordinate_27nodebrick[GP_tag][0]*GP_local_coordinate_27nodebrick[GP_tag][0])*
                            (1 - GP_local_coordinate_27nodebrick[GP_tag][1]*GP_local_coordinate_27nodebrick[GP_tag][1])*
                            (1 + 2*GP_local_coordinate_27nodebrick[GP_tag][2]*node_local_coordinate_27nodebrick[node_tag][2])*
                            node_local_coordinate_27nodebrick[node_tag][2];

                J_GP[2][column][GP_tag] += node_global_coordinates[node_tag][column]*partial_N;
            }
        }
    }
/*
    for (int print_GP_tag = 0; print_GP_tag < 27; print_GP_tag++)
    {
        Energy_Out << element << " " << print_GP_tag << endl;
        Energy_Out << J_GP[0][0][print_GP_tag] << " " << J_GP[0][1][print_GP_tag] << " " << J_GP[0][2][print_GP_tag] << endl;
        Energy_Out << J_GP[1][0][print_GP_tag] << " " << J_GP[1][1][print_GP_tag] << " " << J_GP[1][2][print_GP_tag] << endl;
        Energy_Out << J_GP[2][0][print_GP_tag] << " " << J_GP[2][1][print_GP_tag] << " " << J_GP[2][2][print_GP_tag] << endl;   
    }
*/   

    // Compute Jacobian matrix of the element
    double J[3][3];
    for (int row = 0; row < 3; row++)
    {
        for (int column = 0; column < 3; column++)
        {
            J[row][column] = 0;
            for (int GP_tag = 0; GP_tag < 27; GP_tag++)
            {
                J[row][column] += J_GP[row][column][GP_tag]*weight_27nodebrick[GP_tag];
            }
            J[row][column] = J[row][column]/8;
        }
    }


    // Compute element volume
    element_volume[element] = 8*(J[0][0]*J[1][1]*J[2][2] + 
                              J[0][1]*J[1][2]*J[2][0] + 
                              J[0][2]*J[1][0]*J[2][1] -
                              J[0][0]*J[1][2]*J[2][1] -
                              J[0][1]*J[1][0]*J[2][2] -
                              J[0][2]*J[1][1]*J[2][0]);

    if (element_volume[element] < 0)
    {
        element_volume[element] = -element_volume[element];
    }

    //Energy_Out << element << " " << element_volume[element] << endl;
}



//===========================================================================
// Compute plastic free energy for linear kinematic hardening
//===========================================================================
double Energy_Post_Processing::computePlasticFreeEnergy_Linear_KH(double *plastic_strain, int step)
{
    //Energy_Out << step << " " << *((ss+6*number_of_timesteps)+step) << " " << *((ss+6*number_of_timesteps)+step-1) << endl;

    double plastic_free_energy_linear_KH = 0;   

    for (int i = 0; i < 3; i++)
    {
        plastic_free_energy_linear_KH += (*((plastic_strain+i*number_of_timesteps)+step) - *((plastic_strain+i*number_of_timesteps)+step-1)) * (*((plastic_strain+i*number_of_timesteps)+step) + *((plastic_strain+i*number_of_timesteps)+step-1))/2;
    }
    
    for (int i = 3; i < 6; i++)
    {
        plastic_free_energy_linear_KH += (*((plastic_strain+i*number_of_timesteps)+step) - *((plastic_strain+i*number_of_timesteps)+step-1)) * (*((plastic_strain+i*number_of_timesteps)+step) + *((plastic_strain+i*number_of_timesteps)+step-1));
    }
    
    plastic_free_energy_linear_KH = a1 * plastic_free_energy_linear_KH;

    //Energy_Out << step << " " << plastic_free_energy_linear_KH << endl;
    return plastic_free_energy_linear_KH;
}



//===========================================================================
// Compute plastic free energy for Armstrong-Frederic kinematic hardening
//===========================================================================
double Energy_Post_Processing::computePlasticFreeEnergy_AF_KH(double *equivalent_plastic_strain, double *back_stress, double *deviatoric_plastic_strain, int step)
{
    double incremental_back_stress[6];

    for (int i = 0; i < 6; i++)
    {
        incremental_back_stress[i] = 0;
    }

    for (int j = 0; j < 6; j++)
    {
        incremental_back_stress[j] = 
            2.0/3*AF_ha * (*((deviatoric_plastic_strain+j*number_of_timesteps)+step) - *((deviatoric_plastic_strain+j*number_of_timesteps)+step-1)) - 
            AF_cr * (*((back_stress+j*number_of_timesteps)+step-1)) * (*(equivalent_plastic_strain+step) - *(equivalent_plastic_strain+step-1));

        *((back_stress+j*number_of_timesteps)+step) = *((back_stress+j*number_of_timesteps)+step-1) + incremental_back_stress[j];
    }


    double plastic_free_energy_AF_KH = 0;

    for (int j = 0; j < 3; j++)
    {
        plastic_free_energy_AF_KH += 1.0/(2.0/3*AF_ha) * incremental_back_stress[j] * (*((back_stress+j*number_of_timesteps)+step) + *((back_stress+j*number_of_timesteps)+step-1))/2;
    }
    
    for (int j = 3; j < 6; j++)
    {
        plastic_free_energy_AF_KH += 1.0/(2.0/3*AF_ha) * incremental_back_stress[j] * (*((back_stress+j*number_of_timesteps)+step) + *((back_stress+j*number_of_timesteps)+step-1));
    }

    //Energy_Out << step << " " << plastic_free_energy_AF_KH << endl;
    return plastic_free_energy_AF_KH;


}