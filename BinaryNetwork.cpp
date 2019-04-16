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

  // Initialize model parameters
  int x_dim = 5;
  int y_dim = 5;
  int num_images = 2; 
  int input_firing_rate = 20; //approximate firing rate of input stimuli; note multiplier used later to generate actual stimuli
  int background_firing_rate = 100; //approximate firing rate of noisy neurons feeding into all layers, and preventing dead neurons
  // Set up the simulator with a timestep at which the neuron, synapse and STDP properties will be calculated 
  float timestep = 0.0001;  // In seconds
  int training_epochs = 5; // Number of epochs to have STDP active
  int display_epochs = 5; // Number of epochs where the each stimulus is presented with STDP inactive
  float competitive_connection_prob = 0.75; // Probability parameter that controls how the two competing halves of the network are connected
  float lower_weight_limit = 0.1;
  float upper_weight_limit = 0.3;
  float exc_inh_weight_ratio = 3.0; //parameter that determines how much stronger inhibitory synapses are than excitatory synapses
  BinaryModel->SetTimestep(timestep);

  // Choose an input neuron type
  PatternedPoissonInputSpikingNeurons* patterned_poisson_input_neurons = new PatternedPoissonInputSpikingNeurons();
  // Choose your neuron type
  LIFSpikingNeurons* lif_spiking_neurons = new LIFSpikingNeurons();
  // Choose your synapse type
  ConductanceSpikingSynapses * conductance_spiking_synapses = new ConductanceSpikingSynapses();

  // Allocate your chosen components to the simulator
  BinaryModel->input_spiking_neurons = patterned_poisson_input_neurons;
  BinaryModel->spiking_neurons = lif_spiking_neurons;
  BinaryModel->spiking_synapses = conductance_spiking_synapses;

  // *** Allocate chosen plasticity rule
  custom_stdp_plasticity_parameters_struct * Excit_STDP_PARAMS = new custom_stdp_plasticity_parameters_struct;
  Excit_STDP_PARAMS->a_plus = 0.2f;
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


  /*
      SETUP PROPERTIES AND CREATE NETWORK:
    Note: 
    All Neuron, Synapse and STDP types have associated parameters structures.
    These structures are defined in the header file for that class and allow us to set properties.
  */

  // SETTING UP INPUT NEURONS
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
  // Creating an LIF parameter structure for an excitatory neuron population and an inhibitory
  // 1 x 100 Layer
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
  for (int ii = 1; ii < 4; ii++){ //Note the first element (input layer) has already been asigned above
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
          std::cout << "\n\nCurrent projecting group ID is " << neuron_params_vec[ii][jj] << ", sending to group ID " << neuron_params_vec[ii+1][kk] << "\n";
          std::cout << "File being called is " << ("Connectivity_Data_ff_" + std::to_string(ii) + std::to_string(jj) + std::to_string(kk) + ".syn") << "\n";
          connect_from_python(neuron_params_vec[ii][jj],
          neuron_params_vec[ii+1][kk],
          ff_synapse_params_vec[ii][jj][kk],
          ("Connectivity_Data_ff_" + std::to_string(ii) + std::to_string(jj) + std::to_string(kk) + ".syn"),
          BinaryModel);
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

  std::cout << "\n\n.......\nNow doing lateral connectivity...\n.......\n\n";

  // LATERAL CONNECTIONS
  // Create vector-of-a-vector-of-a-vector to store the synapse structure for lateral connections
  // The first dimension corresponds to the layer, the second to the source side (left or right), and the receiving side (left or right)
  std::vector<std::vector<std::vector<conductance_spiking_synapse_parameters_struct*>>> lat_synapse_params_vec;

  //Create the lateral synapses from the first layer, for the left hand source, with left hand projections
  lat_synapse_params_vec.push_back(std::vector<std::vector<conductance_spiking_synapse_parameters_struct*>>());
  lat_synapse_params_vec[0].push_back(std::vector<conductance_spiking_synapse_parameters_struct*>());
  lat_synapse_params_vec[0][0].push_back(new conductance_spiking_synapse_parameters_struct());
  lat_synapse_params_vec[0][0][0]->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  lat_synapse_params_vec[0][0][0]->delay_range[0] = 10.0*timestep;
  lat_synapse_params_vec[0][0][0]->delay_range[1] = 10.0*timestep; //NB that as the delays will be set later from loaded files, these values are arbitrary, albeit required by Spike
  lat_synapse_params_vec[0][0][0]->decay_term_tau_g = 0.0017f;  // Seconds (Conductance Parameter)
  lat_synapse_params_vec[0][0][0]->reversal_potential_Vhat = 0.0*pow(10.0, -3);
  lat_synapse_params_vec[0][0][0]->connectivity_type = CONNECTIVITY_TYPE_PAIRWISE;

  // Create the right hand projections of the above
  lat_synapse_params_vec[0][0].push_back(new conductance_spiking_synapse_parameters_struct());
  std::memcpy(lat_synapse_params_vec[0][0][1], lat_synapse_params_vec[0][0][0], sizeof(* lat_synapse_params_vec[0][0][0]));

  //Create the above two equivalents, but for the right hand  source
  lat_synapse_params_vec[0].push_back(std::vector<conductance_spiking_synapse_parameters_struct*>());
  lat_synapse_params_vec[0][1].push_back(new conductance_spiking_synapse_parameters_struct());
  std::memcpy(lat_synapse_params_vec[0][1][0], lat_synapse_params_vec[0][0][0], sizeof(* lat_synapse_params_vec[0][0][0])); //Left hand projections
  lat_synapse_params_vec[0][1].push_back(new conductance_spiking_synapse_parameters_struct());
  std::memcpy(lat_synapse_params_vec[0][1][1], lat_synapse_params_vec[0][0][0], sizeof(* lat_synapse_params_vec[0][0][0])); //Right hand projections

  // Iteratively allocate the additional feed-forward synapses
  for (int ii = 1; ii < 3; ii++){ //Note the first element (input layer) has already been asigned above, and the input layer has no lateral connections, so there are fewer ii iterations
    lat_synapse_params_vec.push_back(std::vector<std::vector<conductance_spiking_synapse_parameters_struct*>>()); //Add element for next layer
    for (int jj = 0; jj < 2; jj++){ // Iterate through the left and right hand sides of each layer sending the connections
      lat_synapse_params_vec[ii].push_back(std::vector<conductance_spiking_synapse_parameters_struct*>());
      for (int kk = 0; kk < 2; kk++){ // Iterate through the left and right hand sides of each layer receiving the connections
        lat_synapse_params_vec[ii][jj].push_back(new conductance_spiking_synapse_parameters_struct());
        std::cout << "Additional lateral layer synapse is " << ii << jj << kk << "\n";
        std::memcpy(lat_synapse_params_vec[ii][jj][kk], lat_synapse_params_vec[0][0][0], sizeof(* lat_synapse_params_vec[0][0][0])); //Note all the ff synapses share the same base properties
      }
    }
  }

  // Iteratively load connectivity data
  for (int ii = 0; ii < 3; ii++){ // Iterate through the layers
    for (int jj = 0; jj < 2; jj++){ // Iterate through the left and right hand sides of each layer sending the connections
      for (int kk = 0; kk < 2; kk++){ // Iterate through the left and right hand sides of each layer sending the connections
          std::cout << "\n\nCurrent projecting group ID is " << neuron_params_vec[ii+1][jj] << ", sending to group ID " << neuron_params_vec[ii+1][kk] << "\n";
          std::cout << "File being called is " << ("Connectivity_Data_lat_" + std::to_string(ii) + std::to_string(jj) + std::to_string(kk) + ".syn") << "\n";
          connect_from_python(neuron_params_vec[ii+1][jj], //Note the input layer is skipped due to the addition of 1
          neuron_params_vec[ii+1][kk],
          lat_synapse_params_vec[ii][jj][kk], //lat_synapse_params_vec[ii][jj][kk], //note synapse vector is indexed with ii, not ii+1
          ("Connectivity_Data_lat_" + std::to_string(ii) + std::to_string(jj) + std::to_string(kk) + ".syn"),
          BinaryModel);
      }
    }
  }

  //Check all connectivity data has been assigned ot parameter structures as expected by printing to screen
  for (int ii = 0; ii < 3; ii++){ // Iterate through the layers
    for (int jj = 0; jj < 2; jj++){ // Iterate through the left and right hand sides of each layer sending the connections
      for (int kk = 0; kk < 2; kk++){ // Iterate through the left and right hand sides of each layer sending the connections
          for (int ll = 100; ll < 105; ll++){
            printf("Pre ID %d, post ID %d, weight %f, delay %f\n", lat_synapse_params_vec[ii][jj][kk]->pairwise_connect_presynaptic[ll],
              lat_synapse_params_vec[ii][jj][kk]->pairwise_connect_postsynaptic[ll],
              lat_synapse_params_vec[ii][jj][kk]->pairwise_connect_weight[ll],
              lat_synapse_params_vec[ii][jj][kk]->pairwise_connect_delay[ll]);
        }
      }
    }
  }

  // //Iterate through each layer's excitatory connections between layers

  // connect_from_python(first_excitatory_left_neuron_layer_ID,
  // first_excitatory_left_neuron_layer_ID,
  // excitatory_to_excitatory_ipsilateral_parameters,
  // "Connectivity_Data_Excit_Excit_Ipsi.syn",
  // BinaryModel);

  // // Check successful loading by printing to screen:

  // //Test that the data has been loaded as expected
  // std::cout << "\nCalling the pre-synapse member after file-closing, the number of synapses is " << input_to_ipsilateral_parameters->pairwise_connect_presynaptic.size() << "\n";
  // std::cout << "\nSample of pre-synapse IDs: \n";
  // for (int ii = 1000; ii < 1005; ++ii){
  //     std::cout << input_to_ipsilateral_parameters->pairwise_connect_presynaptic[ii] << " ";
  // }

  // //Test that the data has been loaded as expected
  // std::cout << "\nCalling the post-synapse member after file-closing, the number of synapses is " << input_to_ipsilateral_parameters->pairwise_connect_postsynaptic.size() << "\n";
  // std::cout << "\nSample of post-synapse IDs: \n";
  // for (int ii = 1000; ii < 1005; ++ii){
  //     std::cout << input_to_ipsilateral_parameters->pairwise_connect_postsynaptic[ii] << " ";
  // }

  // //Test that the data has been loaded as expected
  // std::cout << "\nCalling the weights member after file-closing, the number of synapses is " << input_to_ipsilateral_parameters->pairwise_connect_weight.size() << "\n";
  // std::cout << "\nSample of weights: \n";
  // for (int ii = 1000; ii < 1005; ++ii){
  //     std::cout << input_to_ipsilateral_parameters->pairwise_connect_weight[ii] << " ";
  // }

  // //Test that the data has been loaded as expected
  // std::cout << "\nCalling the delays member after file-closing, the number of synapses is " << input_to_ipsilateral_parameters->pairwise_connect_delay.size() << "\n";
  // std::cout << "\nSample of delays: \n";
  // for (int ii = 1000; ii < 1005; ++ii){
  //     std::cout << input_to_ipsilateral_parameters->pairwise_connect_delay[ii] << " ";
  // }

  // //Test that the data has been loaded as expected
  // std::cout << "\nCalling the pre-synapse member after file-closing, the number of synapses is " << input_to_ipsilateral_parameters->pairwise_connect_presynaptic.size() << "\n";
  // std::cout << "\nSample of pre-synapse IDs: \n";
  // for (int ii = 1000; ii < 1005; ++ii){
  //     std::cout << excitatory_to_excitatory_ipsilateral_parameters->pairwise_connect_presynaptic[ii] << " ";
  // }

  // //Test that the data has been loaded as expected
  // std::cout << "\nCalling the post-synapse member after file-closing, the number of synapses is " << input_to_ipsilateral_parameters->pairwise_connect_postsynaptic.size() << "\n";
  // std::cout << "\nSample of post-synapse IDs: \n";
  // for (int ii = 1000; ii < 1005; ++ii){
  //     std::cout << excitatory_to_excitatory_ipsilateral_parameters->pairwise_connect_postsynaptic[ii] << " ";
  // }

  // //Test that the data has been loaded as expected
  // std::cout << "\nCalling the weights member after file-closing, the number of synapses is " << input_to_ipsilateral_parameters->pairwise_connect_weight.size() << "\n";
  // std::cout << "\nSample of weights: \n";
  // for (int ii = 1000; ii < 1005; ++ii){
  //     std::cout << excitatory_to_excitatory_ipsilateral_parameters->pairwise_connect_weight[ii] << " ";
  // }

  // //Test that the data has been loaded as expected
  // std::cout << "\nCalling the delays member after file-closing, the number of synapses is " << input_to_ipsilateral_parameters->pairwise_connect_delay.size() << "\n";
  // std::cout << "\nSample of delays: \n";
  // for (int ii = 1000; ii < 1005; ++ii){
  //     std::cout << excitatory_to_excitatory_ipsilateral_parameters->pairwise_connect_delay[ii] << " ";
  // }

  // std::cout << "\nThe weight-scaling constant for the input-to-excitatory synapses is " << input_to_ipsilateral_parameters->weight_scaling_constant;
  // std::cout << "\nThe weight-scaling constant for the excitatory-to-excitatory synapses is " << excitatory_to_excitatory_ipsilateral_parameters->weight_scaling_constant;


  // //Create synapses for the background noise input to all other neurons
  // conductance_spiking_synapse_parameters_struct * back_input_to_all = new conductance_spiking_synapse_parameters_struct();
  // back_input_to_all->weight_range[0] = lower_weight_limit;   // Create uniform distributions of weights between the upper and lower bound
  // back_input_to_all->weight_range[1] = upper_weight_limit; //NB the weight range is simply the initialization
  // back_input_to_all->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  // back_input_to_all->delay_range[0] = 10.0*timestep;
  // back_input_to_all->delay_range[1] = 100.0*timestep;
  // back_input_to_all->decay_term_tau_g = 0.0017f;  // Seconds (Conductance Parameter)
  // back_input_to_all->reversal_potential_Vhat = 0.0*pow(10.0, -3);
  // // std::memcpy(back_input_to_all, input_to_ipsilateral_parameters, sizeof(* input_to_ipsilateral_parameters));
  // back_input_to_all->connectivity_type = CONNECTIVITY_TYPE_ALL_TO_ALL;


  // // CREATING SYNAPSES
  // // When creating synapses, the ids of the presynaptic and postsynaptic populations are all that are required
  // // Note: Input neuron populations cannot be post-synaptic on any synapse

  // // *** FIRST LAYER SYNAPSES ***
  // // Input (from input neurons)
  // // BinaryModel->AddSynapseGroup(input_left_layer_ID, first_excitatory_left_neuron_layer_ID, input_to_ipsilateral_parameters);
  // // // BinaryModel->AddSynapseGroup(input_left_layer_ID, first_excitatory_right_neuron_layer_ID, input_to_contralateral_parameters);
  // // // BinaryModel->AddSynapseGroup(input_right_layer_ID, first_excitatory_left_neuron_layer_ID, input_to_contralateral_parameters);
  // // BinaryModel->AddSynapseGroup(input_right_layer_ID, first_excitatory_right_neuron_layer_ID, input_to_ipsilateral_parameters);

  // // *** BACKGROUND NOISE SYNAPSES
  // //Create synapses between the background neurons and all others
  // BinaryModel->AddSynapseGroup(back_input_layer_ID, first_excitatory_left_neuron_layer_ID, back_input_to_all);
  // BinaryModel->AddSynapseGroup(back_input_layer_ID, first_excitatory_right_neuron_layer_ID, back_input_to_all);


  // /*
  //     ADD INPUT STIMULI TO THE PATTERNED POISSON NEURONS CLASS
  // */

  // //Initialize array for input firing rates; note that althought it is a 2D input in the model, this is represented in Spike as a 1D array, n*n long
  // int total_input_size = (x_dim * y_dim * num_images * 2);
  // std::vector<float> input_rates(total_input_size, 1.0);

  // for (int ii = 0; ii < (x_dim * y_dim); ++ii){
  //   input_rates[ii] = 0.0f;
  // }
  // for (int jj = (x_dim * y_dim * 3); jj < (x_dim * y_dim * 4); ++jj){
  //   input_rates[jj] = 0.0f;
  // }

  // //Uncomment the following section to test that firing rates for each stimulus have maintained their 2D structure
  
  // //Test that the firing rates have maintained their correct x-y structure by printing to screen
  // for (int ii = 0; ii < num_images; ++ii){
  //   std::cout << "\n\n\n\n*** Stimulus " << (ii+1) << "***\n\n";
  //   //Iterate through each row
  //   for (int jj = 0; jj < y_dim*2; ++jj){
  //     //Iterate through each column in a row
  //     for (int kk = 0; kk < x_dim; ++kk){
  //       std::cout << input_rates[(2 * ii * x_dim * y_dim) + jj*y_dim + kk];
  //     }
  //     std::cout << "\n";
  //   }
  // }
  
  


  // //Invert firing rate values (i.e. 0's and 1's) so that stimuli are the active neurons, and multiply by baseline firing rate
  // for (int ii = 0; ii < total_input_size; ++ii){
  //   input_rates[ii] = ((input_rates[ii] - 1.10) * -1) * input_firing_rate; //Results in a stimuli firing rate that is 1.10*baseline, and a background firing rate that is 0.1*baseline
  // }


  // /*** Assign firing rates to stimuli ***/

  // //Initialize an array of integers to hold the stimulus ID values
  // int stimuli_array[num_images];
  // //Initialize a temporary array for holding stimulus firing rates
  // float temp_stimulus_array[total_number_of_input_neurons]; //

  // //Iterate through each image
  // for (int ii = 0; ii < num_images; ++ii){
  //   //Iterate through each image's firing rates and assign to a temporary array
  //   for (int jj = 0; jj < (2*x_dim*y_dim); jj++){
  //     temp_stimulus_array[jj] = input_rates[(2 * ii * x_dim * y_dim) + jj];
  //   }

  //   // Add the firing rate of the background neurons that input to all others
  //   for (int kk = (2*x_dim*y_dim); kk < total_number_of_input_neurons; kk++){
  //     temp_stimulus_array[kk] = background_firing_rate;
  //   }
  //   stimuli_array[ii] = patterned_poisson_input_neurons->add_stimulus(temp_stimulus_array, total_number_of_input_neurons);
  // }


  // BinaryModel->finalise_model();
  // float simtime = 0.2f; //This should be long enough to allow any recursive signalling to finish propagating

  
  // //    RUN THE SIMULATION WITH TRAINING
  
  // BinaryModel->spiking_synapses->save_weights_as_binary("./", "Initial_Sandbox_Network");

  // // Loop through a certain number of epoch's of presentation
  // for (int ii = 0; ii < training_epochs; ++ii) {
  //   // Within each epoch, loop through each stimulus 
  //   //*** Eventually this order should probably be randomized ***
  //   for (int jj = 0; jj < num_images; ++jj){
  //     BinaryModel->reset_state(); //Re-set the activity of the network, but not e.g. weights and connectivity
  //     patterned_poisson_input_neurons->select_stimulus(stimuli_array[jj]);
  //     BinaryModel->run(simtime, 1); //the second argument ensures STDP is on or off
  //   }
  //   //Save a snapshot of the model's current weights to enable looking for convergence 
  //   BinaryModel->spiking_synapses->save_weights_as_binary("./", "Epoch" + std::to_string(ii) + "Sandbox_Network");

  // }


  // spike_monitor_main->reset_state(); //Dumps all recorded spikes
  // spike_monitor_input->reset_state();
  // BinaryModel->reset_time(); //Resets the internal clock to 0

  // /*
  //     RUN THE SIMULATION AFTER TRAINING WITH FIRST STIMULUS
  // */

  // // Loop through a certain number of epoch's of presentation
  // for (int ii = 0; ii < display_epochs; ++ii) {

  //   BinaryModel->reset_state(); //Re-set the activity of the network, but not e.g. weights and connectivity
  //   patterned_poisson_input_neurons->select_stimulus(stimuli_array[0]);
  //   BinaryModel->run(simtime, 0); //the second argument determines if STDP is on or off

  // }

  // spike_monitor_main->save_spikes_as_binary("./", "output_spikes_posttraining_stim1");
  // spike_monitor_main->save_spikes_as_txt("./", "output_spikes_posttraining_stim1");
  
  // spike_monitor_input->save_spikes_as_binary("./", "input_Poisson_stim1");


  // spike_monitor_main->reset_state(); //Dumps all recorded spikes
  // spike_monitor_input->reset_state();
  // BinaryModel->reset_time(); //Resets the internal clock to 0

  // /*
  //     RUN THE SIMULATION AFTER TRAINING WITH SECOND STIMULUS
  // */

  // // Loop through a certain number of epoch's of presentation
  // for (int ii = 0; ii < display_epochs; ++ii) {

  //   BinaryModel->reset_state(); //Re-set the activity of the network, but not e.g. weights and connectivity
  //   patterned_poisson_input_neurons->select_stimulus(stimuli_array[1]);
  //   BinaryModel->run(simtime, 0); //the second argument determines if STDP is on or off

  // }

  // spike_monitor_main->save_spikes_as_binary("./", "output_spikes_posttraining_stim2");
  // spike_monitor_main->save_spikes_as_txt("./", "output_spikes_posttraining_stim2");

  // spike_monitor_input->save_spikes_as_binary("./", "input_Poisson_stim2");
  // BinaryModel->spiking_synapses->save_weights_as_binary("./", "Sandbox_Network");
  // BinaryModel->spiking_synapses->save_connectivity_as_binary("./", "Sandbox_Network");


  return 0;
}