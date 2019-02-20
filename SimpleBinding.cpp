
#include "Spike/Spike.hpp"
#include <array>
#include <iostream>

// The function which will autorun when the executable is created
int main (int argc, char *argv[]){

  /*
      CHOOSE THE COMPONENTS OF YOUR SIMULATION
  */

  // Create an instance of the Model
  SpikingModel* ExampleModel = new SpikingModel();
  /* Explanation of above notation:
    ExampleModel is intiliazed as a pointer to an object of class SpikingModel
    The 'new' operator is essentially the C++ equivalent of 'malloc' allocates memory for the un-named object, and returns the pointer to this object,
    or if it is an array, the first element. The memory allocation performed by new is with 'dynamic storage duration', such that the lifetime of the 
    object isn't limited to the scope in which it was created. This is also known as allocating memory to the 'heap' (as opposed to the stack)
    and as such memory *de*-allocation is critical in order to prevent a memory leak/'garbage' building up
  */

  // Initialize model parameters
  int x_dim = 32;
  int y_dim = 32;
  int num_images = 2; 
  int baseline_firing_rate = 10;
  // Set up the simulator with a timestep at which the neuron, synapse and STDP properties will be calculated 
  float timestep = 0.0001;  // In seconds
  ExampleModel->SetTimestep(timestep);
  /* Explanation of above notation:
    The arrow notation references the member (SetTimeStep) of the object that is pointed to by the pointer, ExampleModel
    SetTimeStep, a member of the SpikingModel class, takes an argument, in this case timestep.
  */


  // Choose an input neuron type
  PatternedPoissonInputSpikingNeurons* patterned_poisson_input_neurons = new PatternedPoissonInputSpikingNeurons();

  // Choose your neuron type
  LIFSpikingNeurons* lif_spiking_neurons = new LIFSpikingNeurons();

  // Choose your synapse type
  ConductanceSpikingSynapses * conductance_spiking_synapses = new ConductanceSpikingSynapses();

  // Allocate your chosen components to the simulator
  ExampleModel->input_spiking_neurons = patterned_poisson_input_neurons;
  ExampleModel->spiking_neurons = lif_spiking_neurons;
  ExampleModel->spiking_synapses = conductance_spiking_synapses;

  // *** Allocate chosen plasticity rule
  // weightdependent_stdp_plasticity_parameters_struct * WDSTDP_PARAMS = new weightdependent_stdp_plasticity_parameters_struct;
  // WDSTDP_PARAMS->a_plus = 1.0;
  // WDSTDP_PARAMS->a_minus = 1.0;
  // WDSTDP_PARAMS->tau_plus = 0.02;
  // WDSTDP_PARAMS->tau_minus = 0.02;
  // WDSTDP_PARAMS->lambda = 1.0f*powf(10.0, -2);
  // WDSTDP_PARAMS->alpha = 2.02;
  // WDSTDP_PARAMS->w_max = 0.6; //w_max determines the maximum value the weight variable can take during learning

  // WeightDependentSTDPPlasticity * weightdependent_stdp = new WeightDependentSTDPPlasticity((SpikingSynapses *) conductance_spiking_synapses, (SpikingNeurons *) lif_spiking_neurons, (SpikingNeurons *) patterned_poisson_input_neurons, (stdp_plasticity_parameters_struct *) WDSTDP_PARAMS);  
  
  // ExampleModel->AddPlasticityRule(weightdependent_stdp);

  /*
      ADD ANY ACTIVITY MONITORS
  */
  SpikingActivityMonitor* spike_monitor_main = new SpikingActivityMonitor(lif_spiking_neurons);
  ExampleModel->AddActivityMonitor(spike_monitor_main);

  // Add activity monitor for poisson input neurons
  SpikingActivityMonitor* spike_monitor_input = new SpikingActivityMonitor(patterned_poisson_input_neurons);
  ExampleModel->AddActivityMonitor(spike_monitor_input);


  /*
      SETUP PROPERTIES AND CREATE NETWORK:
    Note: 
    All Neuron, Synapse and STDP types have associated parameters structures.
    These structures are defined in the header file for that class and allow us to set properties.
  */

  // SETTING UP INPUT NEURONS
  // Creating an input neuron parameter structure
  patterned_poisson_input_spiking_neuron_parameters_struct* input_neuron_params = new patterned_poisson_input_spiking_neuron_parameters_struct();
  // Setting the dimensions of the input neuron layer
  input_neuron_params->group_shape[0] = x_dim;    // x-dimension of the input neuron layer
  input_neuron_params->group_shape[1] = y_dim;   // y-dimension of the input neuron layer
  // Create a group of input neurons. This function returns the ID of the input neuron group
  int input_layer_ID = ExampleModel->AddInputNeuronGroup(input_neuron_params);

  // SETTING UP NEURON GROUPS
  // Creating an LIF parameter structure for an excitatory neuron population and an inhibitory
  // 1 x 100 Layer
  lif_spiking_neuron_parameters_struct * excitatory_population_params = new lif_spiking_neuron_parameters_struct();
  excitatory_population_params->group_shape[0] = x_dim;
  excitatory_population_params->group_shape[1] = y_dim;
  excitatory_population_params->resting_potential_v0 = -0.074f;
  excitatory_population_params->threshold_for_action_potential_spike = -0.053f;
  excitatory_population_params->somatic_capacitance_Cm = 500.0*pow(10, -12);
  excitatory_population_params->somatic_leakage_conductance_g0 = 25.0*pow(10, -9);

  lif_spiking_neuron_parameters_struct * inhibitory_population_params = new lif_spiking_neuron_parameters_struct();
  inhibitory_population_params->group_shape[0] = 12;
  inhibitory_population_params->group_shape[1] = 12;
  inhibitory_population_params->resting_potential_v0 = -0.082f;
  inhibitory_population_params->threshold_for_action_potential_spike = -0.053f;
  inhibitory_population_params->somatic_capacitance_Cm = 214.0*pow(10, -12);
  inhibitory_population_params->somatic_leakage_conductance_g0 = 18.0*pow(10, -9);

  // Create populations of excitatory and inhibitory neurons
  int first_excitatory_neuron_layer_ID = ExampleModel->AddNeuronGroup(excitatory_population_params);
  int first_inhibitory_neuron_layer_ID = ExampleModel->AddNeuronGroup(inhibitory_population_params);

  // *** Create additional layers - NB that is shares the same properties as the first layer
  // int second_excitatory_neuron_layer_ID = ExampleModel->AddNeuronGroup(excitatory_population_params);
  // int second_inhibitory_neuron_layer_ID = ExampleModel->AddNeuronGroup(inhibitory_population_params);
  // int third_excitatory_neuron_layer_ID = ExampleModel->AddNeuronGroup(excitatory_population_params);
  // int third_inhibitory_neuron_layer_ID = ExampleModel->AddNeuronGroup(inhibitory_population_params);


  // SETTING UP SYNAPSES
  // Creating a synapses parameter structure for connections from the input neurons to the excitatory neurons
  conductance_spiking_synapse_parameters_struct* input_to_excitatory_parameters = new conductance_spiking_synapse_parameters_struct();
  input_to_excitatory_parameters->weight_range[0] = 0.1f;   // Create uniform distributions of weights between the upper and lower bound
  input_to_excitatory_parameters->weight_range[1] = 0.2f; //NB the weight range is simply the initialization
  input_to_excitatory_parameters->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  input_to_excitatory_parameters->delay_range[0] = 10.0*timestep;
  input_to_excitatory_parameters->delay_range[1] = 100.0*timestep;
  input_to_excitatory_parameters->decay_term_tau_g = 0.005f;  // Seconds (Conductance Parameter)
  // input_to_excitatory_parameters->connectivity_type = CONNECTIVITY_TYPE_ONE_TO_ONE;
  input_to_excitatory_parameters->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
  input_to_excitatory_parameters->gaussian_synapses_standard_deviation = 5.0; //connects neurons with specified Gaussian SD
  input_to_excitatory_parameters->max_number_of_connections_per_pair = 5;
  input_to_excitatory_parameters->gaussian_synapses_per_postsynaptic_neuron = 10;

  // *** Add plasticity to input synapses
  // input_to_excitatory_parameters->plasticity_vec.push_back(weightdependent_stdp);


  // Creating a set of synapse parameters for connections from the excitatory neurons to the inhibitory neurons *within a layer*
  conductance_spiking_synapse_parameters_struct * excitatory_to_inhibitory_parameters = new conductance_spiking_synapse_parameters_struct();
  excitatory_to_inhibitory_parameters->weight_range[0] = 0.1f;
  excitatory_to_inhibitory_parameters->weight_range[1] = 0.2f;
  excitatory_to_inhibitory_parameters->weight_scaling_constant = inhibitory_population_params->somatic_leakage_conductance_g0;
  excitatory_to_inhibitory_parameters->delay_range[0] = 10.0*timestep; //Delays range from 1 to 2 ms for inhibitory connectivity
  excitatory_to_inhibitory_parameters->delay_range[1] = 20.0*timestep;
  excitatory_to_inhibitory_parameters->decay_term_tau_g = 0.005f;  // Seconds (Conductance Parameter)
  excitatory_to_inhibitory_parameters->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
  excitatory_to_inhibitory_parameters->gaussian_synapses_standard_deviation = 5.0; //connects neurons with specified Gaussian SD
  excitatory_to_inhibitory_parameters->max_number_of_connections_per_pair = 5;
  excitatory_to_inhibitory_parameters->gaussian_synapses_per_postsynaptic_neuron = 10;

  // Creating a set of synapse parameters from the inhibitory neurons to the excitatory neurons *within a layer*
  conductance_spiking_synapse_parameters_struct * inhibitory_to_excitatory_parameters = new conductance_spiking_synapse_parameters_struct();
  inhibitory_to_excitatory_parameters->weight_range[0] = -0.4f;
  inhibitory_to_excitatory_parameters->weight_range[1] = -0.2f;
  inhibitory_to_excitatory_parameters->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  inhibitory_to_excitatory_parameters->delay_range[0] = 10.0*timestep; //Delays range from 1 to 2 ms for inhibitory connectivity
  inhibitory_to_excitatory_parameters->delay_range[1] = 20.0*timestep;
  inhibitory_to_excitatory_parameters->decay_term_tau_g = 0.005f;  // Seconds (Conductance Parameter)
  inhibitory_to_excitatory_parameters->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
  inhibitory_to_excitatory_parameters->gaussian_synapses_standard_deviation = 1.0; //connects neurons with specified Gaussian SD
  inhibitory_to_excitatory_parameters->max_number_of_connections_per_pair = 5;
  inhibitory_to_excitatory_parameters->gaussian_synapses_per_postsynaptic_neuron = 10;
  
  // Creating a set of synapse parameters for connections from the excitatory neurons back to the excitatory neurons *within a layer*
  conductance_spiking_synapse_parameters_struct * excitatory_to_excitatory_parameters = new conductance_spiking_synapse_parameters_struct();
  excitatory_to_excitatory_parameters->weight_range[0] = 0.01f;
  excitatory_to_excitatory_parameters->weight_range[1] = 0.05f;
  excitatory_to_excitatory_parameters->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  excitatory_to_excitatory_parameters->delay_range[0] = 10.0*timestep; //Delays range from 1 to 10 ms for excitatory connectivity
  excitatory_to_excitatory_parameters->delay_range[1] = 100.0*timestep;
  excitatory_to_excitatory_parameters->decay_term_tau_g = 0.005f;  // Seconds (Conductance Parameter)
  excitatory_to_excitatory_parameters->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
  excitatory_to_excitatory_parameters->gaussian_synapses_standard_deviation = 1.0; //connects neurons with specified Gaussian SD
  excitatory_to_excitatory_parameters->max_number_of_connections_per_pair = 5;
  excitatory_to_excitatory_parameters->gaussian_synapses_per_postsynaptic_neuron = 10;

  // Creating a set of synapse parameters for connections from the excitatory neurons *in a lower layer to the layer above*
  // conductance_spiking_synapse_parameters_struct * lower_to_upper_parameters = new conductance_spiking_synapse_parameters_struct();
  // lower_to_upper_parameters->weight_range[0] = 0.0f;
  // lower_to_upper_parameters->weight_range[1] = 0.5f;
  // lower_to_upper_parameters->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  // lower_to_upper_parameters->delay_range[0] = 10.0*timestep; //Delays range from 1 to 10 ms for excitatory connectivity
  // lower_to_upper_parameters->delay_range[1] = 100.0*timestep;
  // lower_to_upper_parameters->decay_term_tau_g = 0.005f;  // Seconds (Conductance Parameter)
  // lower_to_upper_parameters->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
  // lower_to_upper_parameters->gaussian_synapses_standard_deviation = 1.0; //connects neurons with specified Gaussian SD
  // lower_to_upper_parameters->max_number_of_connections_per_pair = 5;
  // lower_to_upper_parameters->gaussian_synapses_per_postsynaptic_neuron = 50;


  // *** Add plasticity to excitatory to excitatory synapses (w/in layers), and the excitatory connections projecting up layers
  // excitatory_to_excitatory_parameters->plasticity_vec.push_back(weightdependent_stdp);
  // lower_to_upper_parameters->plasticity_vec.push_back(weightdependent_stdp);


  // CREATING SYNAPSES
  // When creating synapses, the ids of the presynaptic and postsynaptic populations are all that are required
  // Note: Input neuron populations cannot be post-synaptic on any synapse
  ExampleModel->AddSynapseGroup(input_layer_ID, first_excitatory_neuron_layer_ID, input_to_excitatory_parameters);
  ExampleModel->AddSynapseGroup(first_excitatory_neuron_layer_ID, first_inhibitory_neuron_layer_ID, excitatory_to_inhibitory_parameters);
  ExampleModel->AddSynapseGroup(first_inhibitory_neuron_layer_ID, first_excitatory_neuron_layer_ID, inhibitory_to_excitatory_parameters);
  ExampleModel->AddSynapseGroup(first_excitatory_neuron_layer_ID, first_excitatory_neuron_layer_ID, excitatory_to_excitatory_parameters);

  // *** Add synapses relevant to the additional layers
  // ExampleModel->AddSynapseGroup(first_excitatory_neuron_layer_ID, second_excitatory_neuron_layer_ID, lower_to_upper_parameters);
  // ExampleModel->AddSynapseGroup(second_excitatory_neuron_layer_ID, second_inhibitory_neuron_layer_ID, excitatory_to_inhibitory_parameters);
  // ExampleModel->AddSynapseGroup(second_inhibitory_neuron_layer_ID, second_excitatory_neuron_layer_ID, inhibitory_to_excitatory_parameters);
  // ExampleModel->AddSynapseGroup(second_excitatory_neuron_layer_ID, second_excitatory_neuron_layer_ID, excitatory_to_excitatory_parameters);

  // ExampleModel->AddSynapseGroup(second_excitatory_neuron_layer_ID, third_excitatory_neuron_layer_ID, lower_to_upper_parameters);
  // ExampleModel->AddSynapseGroup(third_excitatory_neuron_layer_ID, third_inhibitory_neuron_layer_ID, excitatory_to_inhibitory_parameters);
  // ExampleModel->AddSynapseGroup(third_inhibitory_neuron_layer_ID, third_excitatory_neuron_layer_ID, inhibitory_to_excitatory_parameters);
  // ExampleModel->AddSynapseGroup(third_excitatory_neuron_layer_ID, third_excitatory_neuron_layer_ID, excitatory_to_excitatory_parameters);


  /*
      ADD INPUT STIMULI TO THE PATTERNED POISSON NEURONS CLASS
  */

  //Initialize array for input firing rates; note that althought it is a 2D input in the model, this is represented in Spike as a 1D array, n*n long
  int total_input_size = (x_dim * y_dim * num_images);
  std::vector<float> input_rates(total_input_size);


  //Load binary file containing firing rates
  ifstream firing_rates_file;
  firing_rates_file.open("sub_file.gbo", ios::binary);
  if (firing_rates_file.is_open()){ //checks the binary fire successfully opened
    std::cout << "Firing rates file opened, extracting rates...\n";

    firing_rates_file.seekg(0, std::ios::end);
    int num_elements = firing_rates_file.tellg() / sizeof(float); //tellg will return the total number of bytes (indicated by the final read position in the file, which was obtained by seekg and ::end in the line above)
    assert(num_elements == total_input_size); //Check the size of the number of firing rates in the input file is as expected
    firing_rates_file.seekg(0, std::ios::beg);

    firing_rates_file.read(reinterpret_cast<char*>(&input_rates[0]), num_elements*sizeof(float)); //char types are a single byte in C++, so the char pointer is a way of...
    // enforcing C++ to read the file byte-by-byte; this does not however change the type of e..g input_rates from whatever it was initialized to before

  }
  else{
    std::cout << "Error: Issue opening file containing input firing rates\n";
  }

  firing_rates_file.close(); //closes file when input firing rates have been successfully copied


  //Uncomment the following section to test that firing rates for each stimulus have maintained their 2D structure
  
  //Test that the firing rates have maintained their correct x-y structure by printing to screen
  for (int ii = 0; ii < num_images; ++ii){
    std::cout << "\n\n\n\n*** Stimulus " << (ii+1) << "***\n\n";
    //Iterate through each row
    for (int jj = 0; jj < y_dim; ++jj){
      //Iterate through each column in a row
      for (int kk = 0; kk < x_dim; ++kk){
        std::cout << input_rates[(ii * x_dim * y_dim) + jj*32 + kk];
      }
      std::cout << "\n";
    }
  }
  
  


  //Invert firing rate values (i.e. 0's and 1's) so that stimuli are the active neurons, and multiply by baseline firing rate
  for (int ii = 0; ii < total_input_size; ++ii){
    input_rates[ii] = ((input_rates[ii] - 1.25) * -1) * baseline_firing_rate; //Results in a stimuli firing rate that is 1.25*baseline, and a background firing rate that is 0.25*baseline
  }


  /*** Assign firing rates to stimuli ***/

  //Initialize an array of integers to hold the stimulus ID values
  int stimuli_array[num_images];
  //Initialize a temporary array for holding stimulus firing rates
  float temp_stimulus_array[x_dim*y_dim]; //

  //Iterate through each image
  for (int ii = 0; ii < num_images; ++ii){
    //Iterate through each image's firing rates and assign to a temporary array
    for (int jj = 0; jj < (x_dim*y_dim); jj++){
      temp_stimulus_array[jj] = input_rates[(ii * x_dim * y_dim) + jj];
    }
    stimuli_array[ii] = patterned_poisson_input_neurons->add_stimulus(temp_stimulus_array, x_dim*y_dim);
  }


  ExampleModel->finalise_model();
  float simtime = 1.0f; //This should be long enough to allow any recursive signalling to finish propagating

  /*
      RUN THE SIMULATION BEFORE TRAINING
  */

  // // Loop through a certain number of epoch's of presentation
  // for (int ii = 0; ii < 20; ++ii) {
  //   // Within each epoch, loop through each stimulus 
  //   //*** Eventually this order should probably be randomized ***
  //   for (int jj = 0; jj < num_images; ++jj){
  //     ExampleModel->reset_state(); //Re-set the activity of the network, but not e.g. weights and connectivity
  //     patterned_poisson_input_neurons->select_stimulus(stimuli_array[jj]);
  //     ExampleModel->run(simtime, 0); //the second argument ensures STDP is on or off
  //   }

  // }

  //spike_monitor_main->save_spikes_as_txt("./", "output_spikes_pretraining");

  /*
      RUN THE SIMULATION WITH TRAINING
  */

  // Loop through a certain number of epoch's of presentation
  // for (int ii = 0; ii < 10; ++ii) {
  //   // Within each epoch, loop through each stimulus 
  //   //*** Eventually this order should probably be randomized ***
  //   for (int jj = 0; jj < num_images; ++jj){
  //     ExampleModel->reset_state(); //Re-set the activity of the network, but not e.g. weights and connectivity
  //     patterned_poisson_input_neurons->select_stimulus(stimuli_array[jj]);
  //     ExampleModel->run(simtime, 1); //the second argument ensures STDP is on or off
  //   }

  // }


  // spike_monitor_main->reset_state(); //Dumps all recorded spikes
  // // spike_monitor_input->reset_state();
  // ExampleModel->reset_time(); //Resets the internal clock to 0

  /*
      RUN THE SIMULATION AFTER TRAINING WITH FIRST STIMULUS
  */

  // Loop through a certain number of epoch's of presentation
  for (int ii = 0; ii < 2; ++ii) {

    ExampleModel->reset_state(); //Re-set the activity of the network, but not e.g. weights and connectivity
    patterned_poisson_input_neurons->select_stimulus(stimuli_array[0]);
    ExampleModel->run(simtime, 0); //the second argument determines if STDP is on or off

  }

  spike_monitor_main->save_spikes_as_binary("./", "output_spikes_posttraining_stim1");
  spike_monitor_main->save_spikes_as_txt("./", "output_spikes_posttraining_stim1");
  
  spike_monitor_input->save_spikes_as_binary("./", "input_Poisson_stim1");


  spike_monitor_main->reset_state(); //Dumps all recorded spikes
  spike_monitor_input->reset_state();
  ExampleModel->reset_time(); //Resets the internal clock to 0

  /*
      RUN THE SIMULATION AFTER TRAINING WITH SECOND STIMULUS
  */

  // Loop through a certain number of epoch's of presentation
  for (int ii = 0; ii < 2; ++ii) {

    ExampleModel->reset_state(); //Re-set the activity of the network, but not e.g. weights and connectivity
    patterned_poisson_input_neurons->select_stimulus(stimuli_array[1]);
    ExampleModel->run(simtime, 0); //the second argument determines if STDP is on or off

  }

  spike_monitor_main->save_spikes_as_binary("./", "output_spikes_posttraining_stim2");
  spike_monitor_main->save_spikes_as_txt("./", "output_spikes_posttraining_stim2");

  spike_monitor_input->save_spikes_as_binary("./", "input_Poisson_stim2");
  ExampleModel->spiking_synapses->save_weights_as_binary("./", "Sandbox_Network");
  ExampleModel->spiking_synapses->save_connectivity_as_binary("./", "Sandbox_Network");


  return 0;
}