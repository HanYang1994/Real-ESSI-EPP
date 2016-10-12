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

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

#define Energy_Error std::cerr << std::endl << __PRETTY_FUNCTION__ << " - ERROR: \n"
#define Energy_Warning std::cout << std::endl <<  __PRETTY_FUNCTION__ << " - Warning: \n"
#define Energy_Out std::cout << std::endl << __PRETTY_FUNCTION__ << " : \n"

#define NUM_DOF 3
//#define EightNodeBrickLT_OUTPUT_SIZE 8*6*4
#define SOLID_ELEMENT_GP_OUTPUT_SIZE (6*4+1)

class Energy_Post_Processing
{	

public:

    Energy_Post_Processing();
    ~Energy_Post_Processing();


    void initialize();                              //Open output file and read data
    void clean_up();                                //Clean up all data
    void computeAndWriteEnergyDensity();                     //Compute and output energy density components at all time steps


    void readESSIOutput();
    void createEnergyDatasets();


    void computeKEnergyDensity(int node);           //Compute kinetic energy density at each time step at certain node
    void computeSEnergyDensity(int element);        //Compute strain energy density at each time step at certain element (each Gauss Point)

    
    void computeElementAverage_8NodeBrick(int element);
    void computeElementAverage_27NodeBrick(int element);

    
    void computeTotalEnergy();
    void computeElementVolume_8NodeBrick(int element);
    void computeElementVolume_27NodeBrick(int element);


    void getCenterofElement(int element);


    double computePlasticFreeEnergy_Linear_KH(double *plastic_strain, int step);
    double computePlasticFreeEnergy_AF_KH(double *equivalent_plastic_strain, double *back_stress, double *deviatoric_plastic_strain, int step);


    //Material properties
    double mass_density;

    int hardening_type;
    double a1;
    double k1;

    bool is_AF_KH;
    double AF_ha;
    double AF_cr;


    bool is_initialized;
    std::string HDF5filename;
    int number_of_elements;
    int number_of_timesteps;
    int number_of_GPs;
    int number_of_nodes;
    double delta_t, t1, t2, tfinal;

    int max_node_tag;
    int max_element_tag;               


    int* number_of_nodes_element;                    //Number of nodes in each element
    int* index_GDs;                                  //Index to generalized displacements
    int* index_kinetic_energy_density;
    double* kinetic_energy_density;


    int* number_of_GPs_element;                      //Number of Gauss Points in each element
    int* index_to_outputs;
    int last_valid_element;
    bool is_first_valid_element;
    int* index_strain_energy_density;

    double* strain_energy_density_step;
    double* strain_energy_density_cumu;

    double* plastic_work_density_step;
    double* plastic_work_density_cumu;

    double* plastic_free_energy_density_step;
    double* plastic_free_energy_density_cumu;

    double* plastic_dissipation_density_step;
    double* plastic_dissipation_density_cumu;


    int* index_to_connectivity;
    int* connectivity;

    int* index_to_node_coordinates;
    double* node_coordinates;
    double* center_coordinates;

    int* index_energy_density_eleavg;
    double* kinetic_energy_density_eleavg;
    double* strain_energy_density_eleavg;
    double* plastic_work_density_eleavg;
    double* plastic_free_energy_density_eleavg;
    double* plastic_dissipation_density_eleavg;
    double* mechanical_energy_density_eleavg;            // Mechanical Enrgy = Kinetic Energy + Strain Energy
    double* stored_energy_density_eleavg;            // Stored Enrgy = Kinetic Energy + Strain Energy + Plastic Free Energy


    int* class_tags;
    double* element_volume;

    double* back_stress;


    double* total_mechanical_energy;
    double* total_stored_energy;
    double* total_plastic_work;
    double* total_plastic_free_energy;
    double* total_plastic_dissipation;
    double* total_kinetic_energy;
    double* total_strain_energy;
  

    double* weight_27nodebrick;
    double** GP_local_coordinate_27nodebrick;
    int** node_local_coordinate_27nodebrick;

    //HDF5 handles
    hid_t id_file;
    hid_t id_displacements;
    hid_t id_displacements_dataspace;

    hid_t id_outputs;
    hid_t id_outputs_dataspace;


    hid_t id_file_energy;
    hid_t id_KEnergy_dataset;
    hid_t id_index_to_KEnergy_dataset;
    hid_t id_index_to_SEnergy_dataset;
    hid_t id_SEnergy_step_dataset;
    hid_t id_SEnergy_cumu_dataset;
    hid_t id_PWork_step_dataset;
    hid_t id_PWork_cumu_dataset;
    hid_t id_plastic_free_energy_density_step_dataset;
    hid_t id_plastic_free_energy_density_cumu_dataset;
    hid_t id_plastic_dissipation_density_step_dataset;
    hid_t id_plastic_dissipation_density_cumu_dataset;
    

    hid_t id_PWork_avg_dataset;
    hid_t id_plastic_free_energy_avg_dataset;
    hid_t id_plastic_dissipation_avg_dataset;
    hid_t id_mechanical_energy_avg_dataset;
    hid_t id_stored_energy_avg_dataset;
    hid_t id_index_to_Energy_avg_dataset;

    //hid_t id_center_dataset;

    hid_t id_element_volume_dataset;
    //hid_t id_back_stress_dataset;

    hid_t id_total_mechanical_energy_dataset;
    hid_t id_total_stored_energy_dataset;
    hid_t id_total_plastic_work_dataset;
    hid_t id_total_plastic_free_energy_dataset;
    hid_t id_total_plastic_dissipation_dataset;
    hid_t id_total_kinetic_energy_dataset;
    hid_t id_total_strain_energy_dataset;


protected:


private:


};