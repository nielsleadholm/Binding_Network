#include "Spike/Spike.hpp"
#include "UtilityFunctionsLeadholm.hpp"
#include <array>
#include <iostream>
#include <cstring>
#include <string>


// Network with 5x5 neurons in each layer
// If unable to get interesting activity, then can make the most trivial case of literally two 
// parallel networks with no interactivity (as a proxy for winner-take-all connectivity)
// Inhibitory population is the same size as the excitatory population
// Start actually with Gaussian connectivity and SOMO like architecture to see if possible
// Can later use all-to-all connectivity if necessary-


// Things to add:
// Background neurons inputting to all layers to prevent dead neurons following plasticity changes
// *** need to check in the future this isn't cause some odd correlated activity by each 
// background neuron simultaneously activating neurons in multiple layers etc. ***

// The function which will autorun when the executable is created
int main (int argc, char *argv[]){

  /*
      CHOOSE THE COMPONENTS OF YOUR SIMULATION
  */

  // Create an instance of the Model
  SpikingModel* BinaryModel = new SpikingModel();
  /* Explanation of above notation:
    BinaryModel is intiliazed as a pointer to an object of class SpikingModel
    The 'new' operator is essentially the C++ equivalent of 'malloc' allocates memory for the un-named object, and returns the pointer to this object,
    or if it is an array, the first element. The memory allocation performed by new is with 'dynamic storage duration', such that the lifetime of the 
    object isn't limited to the scope in which it was created. This is also known as allocating memory to the 'heap' (as opposed to the stack)
    and as such memory *de*-allocation is critical in order to prevent a memory leak/'garbage' building up
  */

  //Simulation parameters; these can be relatively easily adjusted to observe the affect on the network
  int training_epochs = 5; // Number of epochs to have STDP active
  int display_epochs = 5; // Number of epochs where the each stimulus is presented with STDP inactive
  float exc_inh_weight_ratio = 2.0; //parameter that determines how much stronger inhibitory synapses are than excitatory synapses
  int background_firing_rate = 50; //approximate firing rate of noisy neurons feeding into all layers, and preventing dead neurons
  
  // Initialize core model parameters
  int x_dim = 5;
  int y_dim = 5;
  int num_images = 2; 
  int input_firing_rate = 20; //approximate firing rate of input stimuli; note multiplier used later to generate actual stimuli
  float competitive_connection_prob = 0.75; // Probability parameter that controls how the two competing halves of the network are connected
  float timestep = 0.0001;  // In seconds
  float lower_weight_limit = 0.01;
  float upper_weight_limit = 0.1;
  BinaryModel->SetTimestep(timestep);

  // Choose an input neuron type
  PatternedPoissonInputSpikingNeurons* patterned_poisson_input_neurons = new PatternedPoissonInputSpikingNeurons();
  // Choose the neuron type
  LIFSpikingNeurons* lif_spiking_neurons = new LIFSpikingNeurons();
  // Choose the synapse type
  ConductanceSpikingSynapses * conductance_spiking_synapses = new ConductanceSpikingSynapses();

  // Allocate your chosen components to the simulator
  BinaryModel->input_spiking_neurons = patterned_poisson_input_neurons;
  BinaryModel->spiking_neurons = lif_spiking_neurons;
  BinaryModel->spiking_synapses = conductance_spiking_synapses;

  // *** Allocate chosen plasticity rule
  custom_stdp_plasticity_parameters_struct * Excit_STDP_PARAMS = new custom_stdp_plasticity_parameters_struct;
  Excit_STDP_PARAMS->a_plus = (upper_weight_limit + lower_weight_limit)/2; //Set to the mean of the excitatory weight distribution
  Excit_STDP_PARAMS->a_minus = 1.0f;
  Excit_STDP_PARAMS->tau_plus = 0.01f;
  Excit_STDP_PARAMS->tau_minus = 0.01f;
  Excit_STDP_PARAMS->learning_rate = 0.001f;
  Excit_STDP_PARAMS->a_star = 0; //Excit_STDP_PARAMS->a_plus * Excit_STDP_PARAMS->tau_minus * Inhib_STDP_PARAMS->targetrate;

  CustomSTDPPlasticity * excitatory_stdp = new CustomSTDPPlasticity((SpikingSynapses *) conductance_spiking_synapses, (SpikingNeurons *) lif_spiking_neurons, (SpikingNeurons *) patterned_poisson_input_neurons, (stdp_plasticity_parameters_struct *) Excit_STDP_PARAMS);  
  
  BinaryModel->AddPlasticityRule(excitatory_stdp);

  /*
      ADD ANY ACTIVITY MONITORS
  */
  SpikingActivityMonitor* spike_monitor_main = new SpikingActivityMonitor(lif_spiking_neurons);
  BinaryModel->AddActivityMonitor(spike_monitor_main);
  // Add activity monitor for poisson input neurons
  SpikingActivityMonitor* spike_monitor_input = new SpikingActivityMonitor(patterned_poisson_input_neurons);
  BinaryModel->AddActivityMonitor(spike_monitor_input);


  // SETTING UP INPUT NEURONS

  std::cout << "\n\n.......\nCreating network neuron groups...\n.......\n\n";
  
  // Creating an input neuron parameter structure
  // Initialize a 2D vector to store the neuron group IDs of each excitatory layer, including the input as the 0th layer
  // Note however this vector will not include the background activity neuron group
  std::vector<std::vector<int>> neuron_params_vec;

  // Note the first dimension corresponds to the layer, indexed from 0, corresponding to the input neurons
  // The second dimension corresponds to the 'left' or 'right' side of the network, indexed by 0 and 1 respectively

  patterned_poisson_input_spiking_neuron_parameters_struct* neuron_params_0_0 = new patterned_poisson_input_spiking_neuron_parameters_struct();
  neuron_params_0_0->group_shape[0] = x_dim;    // x-dimension of the input neuron layer
  neuron_params_0_0->group_shape[1] = y_dim;   // y-dimension of the input neuron layer
  neuron_params_vec.push_back(std::vector<int>());
  neuron_params_vec[0].push_back(BinaryModel->AddInputNeuronGroup(neuron_params_0_0));

  patterned_poisson_input_spiking_neuron_parameters_struct* neuron_params_0_1 = new patterned_poisson_input_spiking_neuron_parameters_struct();
  neuron_params_0_1->group_shape[0] = x_dim;    // x-dimension of the input neuron layer
  neuron_params_0_1->group_shape[1] = y_dim;   // y-dimension of the input neuron layer
  neuron_params_vec[0].push_back(BinaryModel->AddInputNeuronGroup(neuron_params_0_1));

  // Set-up background noise neurons; these ensure no 'dead' neurons following plasticity by guarenteeing a random input to every neuron
  patterned_poisson_input_spiking_neuron_parameters_struct* back_input_neuron_params = new patterned_poisson_input_spiking_neuron_parameters_struct();
  back_input_neuron_params->group_shape[0] = x_dim;    // x-dimension of the input neuron layer
  back_input_neuron_params->group_shape[1] = y_dim;   // y-dimension of the input neuron layer
  int back_input_layer_ID = BinaryModel->AddInputNeuronGroup(back_input_neuron_params);

  int total_number_of_input_neurons = (neuron_params_0_0->group_shape[0]*neuron_params_0_0->group_shape[1] 
    + neuron_params_0_1->group_shape[0]*neuron_params_0_1->group_shape[1] 
    + back_input_neuron_params->group_shape[0]*back_input_neuron_params->group_shape[1]);

  // SETTING UP NEURON GROUPS
  lif_spiking_neuron_parameters_struct * excitatory_population_params = new lif_spiking_neuron_parameters_struct();
  excitatory_population_params->group_shape[0] = x_dim;
  excitatory_population_params->group_shape[1] = y_dim;
  excitatory_population_params->resting_potential_v0 = -0.06f;
  excitatory_population_params->absolute_refractory_period = 0.002f;
  excitatory_population_params->threshold_for_action_potential_spike = -0.05f;
  excitatory_population_params->somatic_capacitance_Cm = 200.0*pow(10, -12);
  excitatory_population_params->somatic_leakage_conductance_g0 = 10.0*pow(10, -9);

  std::cout << "New layer ID is " << neuron_params_vec[0][0] << "\n";
  std::cout << "New layer ID is " << neuron_params_vec[0][1] << "\n";

  // Iteratively create all the additional layers of excitatory neurons, storing their IDs in a 2D vector-of-a-vector
  // Note the input neurons are the 0th layer, specified earlier
  for (int ii = 1; ii < 4; ii++){
    // As neuron_params_vec is a vector-of-a-vector without a defined size, need to add an element to the base vector
    neuron_params_vec.push_back(std::vector<int>());
    for (int jj = 0; jj < 2; jj++){ // Iterate through the left and right hand sides of each layer
      //Add an element to the inner vector, and assign the desired value
      neuron_params_vec[ii].push_back(BinaryModel->AddNeuronGroup(excitatory_population_params));
      std::cout << "New layer ID is " << neuron_params_vec[ii][jj] << "\n";
    }
  }


  // SETTING UP SYNAPSES

  std::cout << "\n\n.......\nBuilding feed-forward connectivity...\n.......\n\n";

  // FEED-FORWARD CONNECTIONS
  // Create vector-of-a-vector-of-a-vector to store the synapse structure for feed-forward connections
  // Note the type is actually specified as the synapses parameter structure during the vector initialization
  // The first dimension corresponds to the layer, the second to the source layer side (left or right), and the receiving layer side (left or right)
  std::vector<std::vector<std::vector<conductance_spiking_synapse_parameters_struct*>>> ff_synapse_params_vec;

  //Create the feed-forward synapses for from the input neurons to the first layer, for the left hand side source, with left hand projections
  ff_synapse_params_vec.push_back(std::vector<std::vector<conductance_spiking_synapse_parameters_struct*>>());
  ff_synapse_params_vec[0].push_back(std::vector<conductance_spiking_synapse_parameters_struct*>());
  ff_synapse_params_vec[0][0].push_back(new conductance_spiking_synapse_parameters_struct());
  ff_synapse_params_vec[0][0][0]->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  ff_synapse_params_vec[0][0][0]->delay_range[0] = 10.0*timestep;
  ff_synapse_params_vec[0][0][0]->delay_range[1] = 10.0*timestep; //NB that as the delays will be set later from loaded files, these values are arbitrary, albeit required by Spike
  ff_synapse_params_vec[0][0][0]->decay_term_tau_g = 0.0017f;  // Seconds (Conductance Parameter)
  ff_synapse_params_vec[0][0][0]->reversal_potential_Vhat = 0.0*pow(10.0, -3);
  ff_synapse_params_vec[0][0][0]->connectivity_type = CONNECTIVITY_TYPE_PAIRWISE;

  // Create the right hand projections of the above
  ff_synapse_params_vec[0][0].push_back(new conductance_spiking_synapse_parameters_struct());
  std::memcpy(ff_synapse_params_vec[0][0][1], ff_synapse_params_vec[0][0][0], sizeof(* ff_synapse_params_vec[0][0][0]));

  //Create the above two equivalents, but for the right hand side of the source
  ff_synapse_params_vec[0].push_back(std::vector<conductance_spiking_synapse_parameters_struct*>());
  ff_synapse_params_vec[0][1].push_back(new conductance_spiking_synapse_parameters_struct());
  std::memcpy(ff_synapse_params_vec[0][1][0], ff_synapse_params_vec[0][0][0], sizeof(* ff_synapse_params_vec[0][0][0])); //Left hand projections
  ff_synapse_params_vec[0][1].push_back(new conductance_spiking_synapse_parameters_struct());
  std::memcpy(ff_synapse_params_vec[0][1][1], ff_synapse_params_vec[0][0][0], sizeof(* ff_synapse_params_vec[0][0][0])); //Right hand projections

  // Iteratively allocate the additional feed-forward synapses
  for (int ii = 1; ii < 3; ii++){ //Note the first element (input layer) has already been asigned above
    ff_synapse_params_vec.push_back(std::vector<std::vector<conductance_spiking_synapse_parameters_struct*>>()); //Add element for next layer
    for (int jj = 0; jj < 2; jj++){ // Iterate through the left and right hand sides of each layer sending the connections
      ff_synapse_params_vec[ii].push_back(std::vector<conductance_spiking_synapse_parameters_struct*>());
      for (int kk = 0; kk < 2; kk++){ // Iterate through the left and right hand sides of each layer receiving the connections
        ff_synapse_params_vec[ii][jj].push_back(new conductance_spiking_synapse_parameters_struct());
        std::memcpy(ff_synapse_params_vec[ii][jj][kk], ff_synapse_params_vec[0][0][0], sizeof(* ff_synapse_params_vec[0][0][0])); //Note all the ff synapses share the same base properties
      }
    }
  }

  // Iteratively load connectivity data; note this loop has to be run separately as higher layers are called, which would otherwise not yet have been defined
  for (int ii = 0; ii < 3; ii++){ // Iterate through the layers
    for (int jj = 0; jj < 2; jj++){ // Iterate through the left and right hand sides of each layer sending the connections
      for (int kk = 0; kk < 2; kk++){ // Iterate through the left and right hand sides of each layer sending the connections
          //std::cout << "\n\nCurrent projecting group ID is " << neuron_params_vec[ii][jj] << ", sending to group ID " << neuron_params_vec[ii+1][kk] << "\n";
          //std::cout << "File being called is " << ("Connectivity_Data_ff_" + std::to_string(ii) + std::to_string(jj) + std::to_string(kk) + ".syn") << "\n";
          connect_from_python(neuron_params_vec[ii][jj],
          neuron_params_vec[ii+1][kk],
          ff_synapse_params_vec[ii][jj][kk],
          ("Connectivity_Data_ff_" + std::to_string(ii) + std::to_string(jj) + std::to_string(kk) + ".syn"),
          BinaryModel);
          //std::cout << "The number of synapses is " << ff_synapse_params_vec[ii][jj][kk]->pairwise_connect_presynaptic.size() << "\n";
      }
    }
  }

  // //Check all connectivity data has been assigned to parameter structures as expected by printing to screen
  // for (int ii = 0; ii < 4; ii++){ // Iterate through the layers
  //   for (int jj = 0; jj < 2; jj++){ // Iterate through the left and right hand sides of each layer sending the connections
  //     for (int kk = 0; kk < 2; kk++){ // Iterate through the left and right hand sides of each layer sending the connections
  //         for (int ll = 100; ll < 105; ll++){
  //           printf("Pre ID %d, post ID %d, weight %f, delay %f\n", ff_synapse_params_vec[ii][jj][kk]->pairwise_connect_presynaptic[ll],
  //             ff_synapse_params_vec[ii][jj][kk]->pairwise_connect_postsynaptic[ll],
  //             ff_synapse_params_vec[ii][jj][kk]->pairwise_connect_weight[ll],
  //             ff_synapse_params_vec[ii][jj][kk]->pairwise_connect_delay[ll]);
  //       }
  //     }
  //   }
  // }

  std::cout << "\n\n.......\nBuilding lateral connectivity...\n.......\n\n";




  // LATERAL CONNECTIONS
  

  // *** inbuilt Spike method - limited to one synapse and therefore one delay per conneciton ***


  conductance_spiking_synapse_parameters_struct * lateral_excitation_ipsilateral_parameters = new conductance_spiking_synapse_parameters_struct();
  lateral_excitation_ipsilateral_parameters->weight_range[0] = lower_weight_limit; 
  lateral_excitation_ipsilateral_parameters->weight_range[1] = upper_weight_limit;
  lateral_excitation_ipsilateral_parameters->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  lateral_excitation_ipsilateral_parameters->delay_range[0] = 10.0*timestep; //Delays range from 1 to 2 ms for inhibitory connectivity
  lateral_excitation_ipsilateral_parameters->delay_range[1] = 100.0*timestep;
  lateral_excitation_ipsilateral_parameters->decay_term_tau_g = 0.0017f;  // Seconds (Conductance Parameter)
  lateral_excitation_ipsilateral_parameters->reversal_potential_Vhat = 0.0*pow(10.0, -3);
  lateral_excitation_ipsilateral_parameters->connectivity_type = CONNECTIVITY_TYPE_RANDOM;
  lateral_excitation_ipsilateral_parameters->random_connectivity_probability = competitive_connection_prob;

  conductance_spiking_synapse_parameters_struct * lateral_excitation_contralateral_parameters = new conductance_spiking_synapse_parameters_struct();
  std::memcpy(lateral_excitation_contralateral_parameters, lateral_excitation_ipsilateral_parameters, sizeof(* lateral_excitation_ipsilateral_parameters));
  lateral_excitation_contralateral_parameters->random_connectivity_probability = (1 - competitive_connection_prob);

  // Iteratively create the inhibitory synapses
  for (int ii = 0; ii < 3; ii++){ // Iterate through the layers
    for (int jj = 0; jj < 2; jj++){ // Iterate through the left and right hand sides of each layer sending the connections
      for (int kk = 0; kk < 2; kk++){ // Iterate through the left and right hand sides of each layer sending the connections
          if (jj==kk){
            BinaryModel->AddSynapseGroup(neuron_params_vec[ii+1][jj], neuron_params_vec[ii+1][kk], lateral_excitation_ipsilateral_parameters);
            //std::cout << "\n\nCurrent projecting group ID is " << neuron_params_vec[ii+1][jj] << ", sending to group ID " << neuron_params_vec[ii+1][kk] << ", with ipsilateral probability\n";
          }
          else{
            BinaryModel->AddSynapseGroup(neuron_params_vec[ii+1][jj], neuron_params_vec[ii+1][kk], lateral_excitation_contralateral_parameters);  
            //std::cout << "\n\nCurrent projecting group ID is " << neuron_params_vec[ii+1][jj] << ", sending to group ID " << neuron_params_vec[ii+1][kk] << ", with contralateral probability\n";
          }
      }
    }
  }




  // *** connectivity file loading method: prevents the simulation from running ***


  // // Create vector-of-a-vector-of-a-vector to store the synapse structure for lateral connections
  // // The first dimension corresponds to the layer, the second to the source side (left or right), and the receiving side (left or right)
  // std::vector<std::vector<std::vector<conductance_spiking_synapse_parameters_struct*>>> lat_synapse_params_vec;

  // //Create the lateral synapses from the first layer, for the left hand source, with left hand projections
  // lat_synapse_params_vec.push_back(std::vector<std::vector<conductance_spiking_synapse_parameters_struct*>>());
  // lat_synapse_params_vec[0].push_back(std::vector<conductance_spiking_synapse_parameters_struct*>());
  // lat_synapse_params_vec[0][0].push_back(new conductance_spiking_synapse_parameters_struct());
  // lat_synapse_params_vec[0][0][0]->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  // lat_synapse_params_vec[0][0][0]->delay_range[0] = 10.0*timestep;
  // lat_synapse_params_vec[0][0][0]->delay_range[1] = 10.0*timestep; //NB that as the delays will be set later from loaded files, these values are arbitrary, albeit required by Spike
  // lat_synapse_params_vec[0][0][0]->decay_term_tau_g = 0.0017f;  // Seconds (Conductance Parameter)
  // lat_synapse_params_vec[0][0][0]->reversal_potential_Vhat = 0.0*pow(10.0, -3);
  // lat_synapse_params_vec[0][0][0]->connectivity_type = CONNECTIVITY_TYPE_PAIRWISE;

  // // Create the right hand projections of the above
  // lat_synapse_params_vec[0][0].push_back(new conductance_spiking_synapse_parameters_struct());
  // std::memcpy(lat_synapse_params_vec[0][0][1], lat_synapse_params_vec[0][0][0], sizeof(* lat_synapse_params_vec[0][0][0]));

  // //Create the above two equivalents, but for the right hand  source
  // lat_synapse_params_vec[0].push_back(std::vector<conductance_spiking_synapse_parameters_struct*>());
  // lat_synapse_params_vec[0][1].push_back(new conductance_spiking_synapse_parameters_struct());
  // std::memcpy(lat_synapse_params_vec[0][1][0], lat_synapse_params_vec[0][0][0], sizeof(* lat_synapse_params_vec[0][0][0])); //Left hand projections
  // lat_synapse_params_vec[0][1].push_back(new conductance_spiking_synapse_parameters_struct());
  // std::memcpy(lat_synapse_params_vec[0][1][1], lat_synapse_params_vec[0][0][0], sizeof(* lat_synapse_params_vec[0][0][0])); //Right hand projections

  // // Iteratively allocate the additional lateral synapses
  // for (int ii = 1; ii < 3; ii++){ //Note the first element (the first layer) has already been assigned above
  //   lat_synapse_params_vec.push_back(std::vector<std::vector<conductance_spiking_synapse_parameters_struct*>>()); //Add element for next layer
  //   for (int jj = 0; jj < 2; jj++){ // Iterate through the left and right hand sides of each layer sending the connections
  //     lat_synapse_params_vec[ii].push_back(std::vector<conductance_spiking_synapse_parameters_struct*>());
  //     for (int kk = 0; kk < 2; kk++){ // Iterate through the left and right hand sides of each layer receiving the connections
  //       lat_synapse_params_vec[ii][jj].push_back(new conductance_spiking_synapse_parameters_struct());
  //       //std::cout << "Additional lateral layer synapse is " << ii << jj << kk << "\n";
  //       std::memcpy(lat_synapse_params_vec[ii][jj][kk], lat_synapse_params_vec[0][0][0], sizeof(* lat_synapse_params_vec[0][0][0]));
  //     }
  //   }
  // }

  // // Iteratively load connectivity data
  // for (int ii = 0; ii < 2; ii++){ // Iterate through the layers
  //   for (int jj = 0; jj < 2; jj++){ // Iterate through the left and right hand sides of each layer sending the connections
  //     for (int kk = 0; kk < 2; kk++){ // Iterate through the left and right hand sides of each layer sending the connections
  //         std::cout << "\n\nCurrent projecting group ID is " << neuron_params_vec[ii+1][jj] << ", sending to group ID " << neuron_params_vec[ii+1][kk] << "\n";
  //         std::cout << "File being called is " << ("Connectivity_Data_lat_" + std::to_string(ii) + std::to_string(jj) + std::to_string(kk) + ".syn") << "\n";
  //         connect_from_python(neuron_params_vec[ii+1][jj], //Note the input layer is skipped due to the addition of 1
  //         neuron_params_vec[ii+1][kk],
  //         lat_synapse_params_vec[ii][jj][kk], //lat_synapse_params_vec[ii][jj][kk], //note synapse vector is indexed with ii, not ii+1
  //         ("Connectivity_Data_lat_" + std::to_string(ii) + std::to_string(jj) + std::to_string(kk) + ".syn"),
  //         BinaryModel);
  //         std::cout << "The number of synapses is " << lat_synapse_params_vec[ii][jj][kk]->pairwise_connect_presynaptic.size() << "\n";
  //     }
  //   }
  // }

  // //Check all connectivity data has been assigned ot parameter structures as expected by printing to screen
  // for (int ii = 0; ii < 3; ii++){ // Iterate through the layers
  //   for (int jj = 0; jj < 2; jj++){ // Iterate through the left and right hand sides of each layer sending the connections
  //     for (int kk = 0; kk < 2; kk++){ // Iterate through the left and right hand sides of each layer sending the connections
  //         for (int ll = 100; ll < 105; ll++){
  //           printf("Pre ID %d, post ID %d, weight %f, delay %f\n", lat_synapse_params_vec[ii][jj][kk]->pairwise_connect_presynaptic[ll],
  //             lat_synapse_params_vec[ii][jj][kk]->pairwise_connect_postsynaptic[ll],
  //             lat_synapse_params_vec[ii][jj][kk]->pairwise_connect_weight[ll],
  //             lat_synapse_params_vec[ii][jj][kk]->pairwise_connect_delay[ll]);
  //       }
  //     }
  //   }
  // }

  std::cout << "\n\n.......\nBuilding inhibitory and background connectivity...\n.......\n\n";

  // INHIBITORY CONNECTIONS
  conductance_spiking_synapse_parameters_struct * non_Dale_inhibition_ipsilateral_parameters = new conductance_spiking_synapse_parameters_struct();
  non_Dale_inhibition_ipsilateral_parameters->weight_range[0] = lower_weight_limit * exc_inh_weight_ratio; 
  non_Dale_inhibition_ipsilateral_parameters->weight_range[1] = upper_weight_limit * exc_inh_weight_ratio;
  // *** NB the inhibitory weights must be positive if the reversal potential is set appropriately (e.g. negative 80mv)
  non_Dale_inhibition_ipsilateral_parameters->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  non_Dale_inhibition_ipsilateral_parameters->delay_range[0] = 10.0*timestep; //Delays range from 1 to 2 ms for inhibitory connectivity
  non_Dale_inhibition_ipsilateral_parameters->delay_range[1] = 10.0*timestep;
  non_Dale_inhibition_ipsilateral_parameters->decay_term_tau_g = 0.0017f;  // Seconds (Conductance Parameter)
  non_Dale_inhibition_ipsilateral_parameters->reversal_potential_Vhat = -80.0*pow(10.0, -3);
  non_Dale_inhibition_ipsilateral_parameters->connectivity_type = CONNECTIVITY_TYPE_RANDOM;
  non_Dale_inhibition_ipsilateral_parameters->random_connectivity_probability = (1 - competitive_connection_prob);

  conductance_spiking_synapse_parameters_struct * non_Dale_inhibition_contralateral_parameters = new conductance_spiking_synapse_parameters_struct();
  std::memcpy(non_Dale_inhibition_contralateral_parameters, non_Dale_inhibition_ipsilateral_parameters, sizeof(* non_Dale_inhibition_ipsilateral_parameters));
  non_Dale_inhibition_contralateral_parameters->random_connectivity_probability = competitive_connection_prob;

  // Iteratively create the inhibitory synapses
  for (int ii = 0; ii < 3; ii++){ // Iterate through the layers
    for (int jj = 0; jj < 2; jj++){ // Iterate through the left and right hand sides of each layer sending the connections
      for (int kk = 0; kk < 2; kk++){ // Iterate through the left and right hand sides of each layer sending the connections
          if (jj==kk){
            BinaryModel->AddSynapseGroup(neuron_params_vec[ii+1][jj], neuron_params_vec[ii+1][kk], non_Dale_inhibition_ipsilateral_parameters);
            //std::cout << "\n\nCurrent projecting group ID is " << neuron_params_vec[ii+1][jj] << ", sending to group ID " << neuron_params_vec[ii+1][kk] << ", with ipsilateral probability\n";
          }
          else{
            BinaryModel->AddSynapseGroup(neuron_params_vec[ii+1][jj], neuron_params_vec[ii+1][kk], non_Dale_inhibition_contralateral_parameters);  
            //std::cout << "\n\nCurrent projecting group ID is " << neuron_params_vec[ii+1][jj] << ", sending to group ID " << neuron_params_vec[ii+1][kk] << ", with contralateral probability\n";
          }
      }
    }
  }


  // BACKGROUND CONNECTIONS
  //Create synapses for the background noise input to all other neurons
  conductance_spiking_synapse_parameters_struct * back_input_to_all = new conductance_spiking_synapse_parameters_struct();
  back_input_to_all->weight_range[0] = lower_weight_limit;   // Create uniform distributions of weights between the upper and lower bound
  back_input_to_all->weight_range[1] = upper_weight_limit; //NB the weight range is simply the initialization
  back_input_to_all->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  back_input_to_all->delay_range[0] = 10.0*timestep;
  back_input_to_all->delay_range[1] = 100.0*timestep;
  back_input_to_all->decay_term_tau_g = 0.0017f;  // Seconds (Conductance Parameter)
  back_input_to_all->reversal_potential_Vhat = 0.0*pow(10.0, -3);
  back_input_to_all->connectivity_type = CONNECTIVITY_TYPE_ALL_TO_ALL;

  // Iteratively create the background input synapses
  for (int ii = 0; ii < 3; ii++){ // Iterate through the layers
    for (int jj = 0; jj < 2; jj++){ // Iterate through the left and right hand sides of each layer
      //Note the first (input) layers are skipped, and each side of each layer only needs to receive one connection
      BinaryModel->AddSynapseGroup(back_input_layer_ID, neuron_params_vec[ii+1][jj], back_input_to_all);
      //std::cout << "\n\nBackground input neurons are sending to group ID " << neuron_params_vec[ii+1][jj] << "\n";
    }
  }

  std::cout << "\n\n.......\nAll synapses created...\n.......\n\n";

  /*
      ADD INPUT STIMULI TO THE PATTERNED POISSON NEURONS CLASS
  */

  std::cout << "\n\n.......\nAssigning stimuli to the network...\n.......\n\n";

  //Initialize array for input firing rates; note that althought it is a 2D input in the model, this is represented in Spike as a 1D array
  int total_input_size = (x_dim * y_dim * num_images * 2);
  std::vector<float> input_rates(total_input_size, 1.0); //Initializes an array of ones

  //Set the first and last section of input_rates to values such that they will be the 'on' stimuli (note the zero's and one's are inverted later)
  for (int ii = 0; ii < (x_dim * y_dim); ++ii){
    input_rates[ii] = 0.0f;
  }
  for (int jj = (x_dim * y_dim * 3); jj < (x_dim * y_dim * 4); ++jj){
    input_rates[jj] = 0.0f;
  }

  //Uncomment the following section to test that firing rates for each stimulus have maintained their 2D structure by printing to screen
  for (int ii = 0; ii < num_images; ++ii){
    std::cout << "\n\n\n\n*** Stimulus " << (ii+1) << "***\n\n";
    //Iterate through each row
    for (int jj = 0; jj < y_dim*2; ++jj){
      //Iterate through each column in a row
      for (int kk = 0; kk < x_dim; ++kk){
        std::cout << input_rates[(2 * ii * x_dim * y_dim) + jj*y_dim + kk];
      }
      std::cout << "\n";
    }
  }
  

  //Invert firing rate values (i.e. 0's and 1's) so that stimuli are the active neurons, and multiply by baseline firing rate
  for (int ii = 0; ii < total_input_size; ++ii){
    input_rates[ii] = ((input_rates[ii] - 1.10) * -1) * input_firing_rate; //Results in a stimuli firing rate that is 1.10*baseline, and a background firing rate that is 0.1*baseline
  }


  /*** Assign firing rates to stimuli ***/

  //Initialize an array of integers to hold the stimulus ID values
  int stimuli_array[num_images];
  //Initialize a temporary array for holding stimulus firing rates
  float temp_stimulus_array[total_number_of_input_neurons]; //

  //Iterate through each image/stimulus
  for (int ii = 0; ii < num_images; ++ii){
    //Iterate through each image's firing rates and assign to a temporary array
    for (int jj = 0; jj < (2*x_dim*y_dim); jj++){
      temp_stimulus_array[jj] = input_rates[(2 * ii * x_dim * y_dim) + jj];
    }

    // Add the firing rate of the background neurons that input to all others
    for (int kk = (2*x_dim*y_dim); kk < total_number_of_input_neurons; kk++){
      temp_stimulus_array[kk] = background_firing_rate;
    }
    stimuli_array[ii] = patterned_poisson_input_neurons->add_stimulus(temp_stimulus_array, total_number_of_input_neurons);
  }


  BinaryModel->finalise_model();
  float simtime = 0.2f; //This should be long enough to allow any recursive signalling to finish propagating

  std::cout << "\n\n.......\nModel finalized and ready for simulating...\n.......\n\n";

  
  //    RUN THE SIMULATION WITH TRAINING
  
  BinaryModel->spiking_synapses->save_weights_as_binary("./", "Initial_Sandbox_Network");

  // Loop through a certain number of epoch's of presentation
  for (int ii = 0; ii < training_epochs; ++ii) {
    // Within each epoch, loop through each stimulus 
    //*** Eventually this order should probably be randomized ***
    for (int jj = 0; jj < num_images; ++jj){
      BinaryModel->reset_state(); //Re-set the activity of the network, but not e.g. weights and connectivity
      patterned_poisson_input_neurons->select_stimulus(stimuli_array[jj]);
      BinaryModel->run(simtime, 0); //the second argument is a boolean determining if STDP is on or off
    }
    //Save a snapshot of the model's current weights to enable looking for convergence 
    BinaryModel->spiking_synapses->save_weights_as_binary("./", "Epoch" + std::to_string(ii) + "Sandbox_Network");

  }


  spike_monitor_main->reset_state(); //Dumps all recorded spikes
  spike_monitor_input->reset_state();
  BinaryModel->reset_time(); //Resets the internal clock to 0

  /*
      RUN THE SIMULATION AFTER TRAINING WITH FIRST STIMULUS
  */

  // Loop through a certain number of epoch's of presentation
  for (int ii = 0; ii < display_epochs; ++ii) {

    BinaryModel->reset_state(); //Re-set the activity of the network, but not e.g. weights and connectivity
    patterned_poisson_input_neurons->select_stimulus(stimuli_array[0]);
    BinaryModel->run(simtime, 0); //the second argument is a boolean determining if STDP is on or off

  }

  spike_monitor_main->save_spikes_as_binary("./", "output_spikes_posttraining_stim1");
  spike_monitor_main->save_spikes_as_txt("./", "output_spikes_posttraining_stim1");
  
  spike_monitor_input->save_spikes_as_binary("./", "input_Poisson_stim1");


  spike_monitor_main->reset_state(); //Dumps all recorded spikes
  spike_monitor_input->reset_state();
  BinaryModel->reset_time(); //Resets the internal clock to 0

  /*
      RUN THE SIMULATION AFTER TRAINING WITH SECOND STIMULUS
  */

  // Loop through a certain number of epoch's of presentation
  for (int ii = 0; ii < display_epochs; ++ii) {

    BinaryModel->reset_state(); //Re-set the activity of the network, but not e.g. weights and connectivity
    patterned_poisson_input_neurons->select_stimulus(stimuli_array[1]);
    BinaryModel->run(simtime, 0); 

  }

  spike_monitor_main->save_spikes_as_binary("./", "output_spikes_posttraining_stim2");
  spike_monitor_main->save_spikes_as_txt("./", "output_spikes_posttraining_stim2");

  spike_monitor_input->save_spikes_as_binary("./", "input_Poisson_stim2");
  BinaryModel->spiking_synapses->save_weights_as_binary("./", "Sandbox_Network");
  BinaryModel->spiking_synapses->save_connectivity_as_binary("./", "Sandbox_Network");


  return 0;
}