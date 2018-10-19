
#include "Spike/Spike.hpp"

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

  // Set up the simulator with a timestep at which the neuron, synapse and STDP properties will be calculated 
  float timestep = 0.0001;  // In seconds
  ExampleModel->SetTimestep(timestep);
  /* Explanation of above notation:
    The arrow notation references the member (SetTimeStep) of the object that is pointed to by the pointer, ExampleModel
    SetTimeStep, a member of the SpikingModel class, takes an argument, in this case timestep.
  */


  // Choose an input neuron type
  GeneratorInputSpikingNeurons* generator_input_neurons = new GeneratorInputSpikingNeurons();

  // Choose your neuron type
  LIFSpikingNeurons* lif_spiking_neurons = new LIFSpikingNeurons();

  // Choose your synapse type
  ConductanceSpikingSynapses * conductance_spiking_synapses = new ConductanceSpikingSynapses();

  // Allocate your chosen components to the simulator
  ExampleModel->input_spiking_neurons = generator_input_neurons;
  ExampleModel->spiking_neurons = lif_spiking_neurons;
  ExampleModel->spiking_synapses = conductance_spiking_synapses;

  // *** Allocate chosen plasticity rule
  weightdependent_stdp_plasticity_parameters_struct * WDSTDP_PARAMS = new weightdependent_stdp_plasticity_parameters_struct;
  WDSTDP_PARAMS->a_plus = 1.0;
  WDSTDP_PARAMS->a_minus = 1.0;
  WDSTDP_PARAMS->tau_plus = 0.02;
  WDSTDP_PARAMS->tau_minus = 0.02;
  WDSTDP_PARAMS->lambda = 1.0f*powf(10.0, -2);
  WDSTDP_PARAMS->alpha = 2.02;
  WDSTDP_PARAMS->w_max = 0.3*powf(10.0, -3);

  WeightDependentSTDPPlasticity * weightdependent_stdp = new WeightDependentSTDPPlasticity((SpikingSynapses *) conductance_spiking_synapses, (SpikingNeurons *) lif_spiking_neurons, (SpikingNeurons *) generator_input_neurons, (stdp_plasticity_parameters_struct *) WDSTDP_PARAMS);  
  
  ExampleModel->AddPlasticityRule(weightdependent_stdp);

  /*
      ADD ANY ACTIVITY MONITORS OR PLASTICITY RULES YOU WISH FOR 
  */
  SpikingActivityMonitor* spike_monitor_main = new SpikingActivityMonitor(lif_spiking_neurons);
  ExampleModel->AddActivityMonitor(spike_monitor_main);

  // Add activity monitor for generator neurons
  SpikingActivityMonitor* spike_monitor_input = new SpikingActivityMonitor(generator_input_neurons);
  ExampleModel->AddActivityMonitor(spike_monitor_input);


  /*
      SETUP PROPERTIES AND CREATE NETWORK:
    
    Note: 
    All Neuron, Synapse and STDP types have associated parameters structures.
    These structures are defined in the header file for that class and allow us to set properties.
  */

  // SETTING UP INPUT NEURONS
  // Creating an input neuron parameter structure
  generator_input_spiking_neuron_parameters_struct* input_neuron_params = new generator_input_spiking_neuron_parameters_struct();
  // Setting the dimensions of the input neuron layer
  input_neuron_params->group_shape[0] = 1;    // x-dimension of the input neuron layer
  input_neuron_params->group_shape[1] = 10;   // y-dimension of the input neuron layer
  // Create a group of input neurons. This function returns the ID of the input neuron group
  int input_layer_ID = ExampleModel->AddInputNeuronGroup(input_neuron_params);

  // SETTING UP NEURON GROUPS
  // Creating an LIF parameter structure for an excitatory neuron population and an inhibitory
  // 1 x 100 Layer
  lif_spiking_neuron_parameters_struct * excitatory_population_params = new lif_spiking_neuron_parameters_struct();
  excitatory_population_params->group_shape[0] = 1;
  excitatory_population_params->group_shape[1] = 80;
  excitatory_population_params->resting_potential_v0 = -0.074f;
  excitatory_population_params->threshold_for_action_potential_spike = -0.053f;
  excitatory_population_params->somatic_capacitance_Cm = 500.0*pow(10, -12);
  excitatory_population_params->somatic_leakage_conductance_g0 = 25.0*pow(10, -9);

  lif_spiking_neuron_parameters_struct * inhibitory_population_params = new lif_spiking_neuron_parameters_struct();
  inhibitory_population_params->group_shape[0] = 1;
  inhibitory_population_params->group_shape[1] = 20;
  inhibitory_population_params->resting_potential_v0 = -0.082f;
  inhibitory_population_params->threshold_for_action_potential_spike = -0.053f;
  inhibitory_population_params->somatic_capacitance_Cm = 214.0*pow(10, -12);
  inhibitory_population_params->somatic_leakage_conductance_g0 = 18.0*pow(10, -9);

  // Create populations of excitatory and inhibitory neurons
  int excitatory_neuron_layer_ID = ExampleModel->AddNeuronGroup(excitatory_population_params);
  int inhibitory_neuron_layer_ID = ExampleModel->AddNeuronGroup(inhibitory_population_params);


  // SETTING UP SYNAPSES
  // Creating a synapses parameter structure for connections from the input neurons to the excitatory neurons
  conductance_spiking_synapse_parameters_struct* input_to_excitatory_parameters = new conductance_spiking_synapse_parameters_struct();
  input_to_excitatory_parameters->weight_range[0] = 0.5f;   // Create uniform distributions of weights [0.5, 10.0]
  input_to_excitatory_parameters->weight_range[1] = 10.0f;
  input_to_excitatory_parameters->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  input_to_excitatory_parameters->delay_range[0] = 10.0*timestep;   //Delays range from 1 to 10 ms for excitatory connectivity
  input_to_excitatory_parameters->delay_range[1] = 100.0*timestep;

  // *** Add plasticity to input synapses
  input_to_excitatory_parameters->plasticity_vec.push_back(weightdependent_stdp);

  // The connectivity types for synapses include:
    // CONNECTIVITY_TYPE_ALL_TO_ALL
    // CONNECTIVITY_TYPE_ONE_TO_ONE
    // CONNECTIVITY_TYPE_RANDOM
    // CONNECTIVITY_TYPE_PAIRWISE
  input_to_excitatory_parameters->connectivity_type = CONNECTIVITY_TYPE_RANDOM;
  input_to_excitatory_parameters->random_connectivity_probability = 0.8;
  //input_to_excitatory_parameters->plasticity_vec.push_back(STDP_RULE);


  // Creating a set of synapse parameters for connections from the excitatory neurons to the inhibitory neurons
  conductance_spiking_synapse_parameters_struct * excitatory_to_inhibitory_parameters = new conductance_spiking_synapse_parameters_struct();
  excitatory_to_inhibitory_parameters->weight_range[0] = 10.0f;
  excitatory_to_inhibitory_parameters->weight_range[1] = 10.0f;
  excitatory_to_inhibitory_parameters->weight_scaling_constant = inhibitory_population_params->somatic_leakage_conductance_g0;
  excitatory_to_inhibitory_parameters->delay_range[0] = 10.0*timestep; //Delays range from 1 to 2 ms for inhibitory connectivity
  excitatory_to_inhibitory_parameters->delay_range[1] = 20.0*timestep;
  excitatory_to_inhibitory_parameters->connectivity_type = CONNECTIVITY_TYPE_RANDOM;
  excitatory_to_inhibitory_parameters->random_connectivity_probability = 0.5; //connects neurons with 10% probability 

  // Creating a set of synapse parameters from the inhibitory neurons to the excitatory neurons
  conductance_spiking_synapse_parameters_struct * inhibitory_to_excitatory_parameters = new conductance_spiking_synapse_parameters_struct();
  inhibitory_to_excitatory_parameters->weight_range[0] = -5.0f;
  inhibitory_to_excitatory_parameters->weight_range[1] = -2.5f;
  inhibitory_to_excitatory_parameters->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  inhibitory_to_excitatory_parameters->delay_range[0] = 10.0*timestep; //Delays range from 1 to 2 ms for inhibitory connectivity
  inhibitory_to_excitatory_parameters->delay_range[1] = 20.0*timestep;
  inhibitory_to_excitatory_parameters->connectivity_type = CONNECTIVITY_TYPE_ALL_TO_ALL;
  
  // Creating a set of synapse parameters for connections from the excitatory neurons back to the excitatory neurons
  conductance_spiking_synapse_parameters_struct * excitatory_to_excitatory_parameters = new conductance_spiking_synapse_parameters_struct();
  excitatory_to_excitatory_parameters->weight_range[0] = 10.0f;
  excitatory_to_excitatory_parameters->weight_range[1] = 10.0f;
  excitatory_to_excitatory_parameters->weight_scaling_constant = inhibitory_population_params->somatic_leakage_conductance_g0;
  excitatory_to_excitatory_parameters->delay_range[0] = 10.0*timestep; //Delays range from 1 to 10 ms for excitatory connectivity
  excitatory_to_excitatory_parameters->delay_range[1] = 100.0*timestep;
  excitatory_to_excitatory_parameters->connectivity_type = CONNECTIVITY_TYPE_RANDOM;
  excitatory_to_excitatory_parameters->random_connectivity_probability = 0.8; //connects neurons with 10% probability

  // *** Add plasticity to excitatory to excitatory synapses
  excitatory_to_excitatory_parameters->plasticity_vec.push_back(weightdependent_stdp);

  // CREATING SYNAPSES
  // When creating synapses, the ids of the presynaptic and postsynaptic populations are all that are required
  // Note: Input neuron populations cannot be post-synaptic on any synapse
  ExampleModel->AddSynapseGroup(input_layer_ID, excitatory_neuron_layer_ID, input_to_excitatory_parameters);
  ExampleModel->AddSynapseGroup(excitatory_neuron_layer_ID, inhibitory_neuron_layer_ID, excitatory_to_inhibitory_parameters);
  ExampleModel->AddSynapseGroup(inhibitory_neuron_layer_ID, excitatory_neuron_layer_ID, inhibitory_to_excitatory_parameters);
  ExampleModel->AddSynapseGroup(excitatory_neuron_layer_ID, excitatory_neuron_layer_ID, excitatory_to_excitatory_parameters);


  /*
      ADD INPUT STIMULI TO THE GENERATOR NEURONS CLASS
  */
  // First stimulus is the 'ascending' pattern; pattern takes place over 10 ms.
  int s1_num_spikes = 10;
  int s1_neuron_ids[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  float s1_spike_times[10] = {0.001f, 0.002f, 0.003f, 0.004f, 0.005f, 0.006f, 0.007f, 0.008f, 0.009f, 0.01f};
  // Adding this stimulus to the input neurons
  int first_stimulus = generator_input_neurons->add_stimulus(s1_num_spikes, s1_neuron_ids, s1_spike_times);
  // Creating a second stimulus (descending pattern)
  int s2_num_spikes = 10;
  int s2_neuron_ids[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  float s2_spike_times[10] = {0.01f, 0.009f, 0.008f, 0.007f, 0.006f, 0.005f, 0.004f, 0.003f, 0.002f, 0.001f};
  int second_stimulus = generator_input_neurons->add_stimulus(s2_num_spikes, s2_neuron_ids, s2_spike_times);
  

  /*
      RUN THE SIMULATION
  */

  // The only argument to run is the number of seconds
  ExampleModel->finalise_model();
  float simtime = 0.03f; //This should be long enough to allow any recursive signalling to finish propagating

  // Very first 'loop' is executed outside, such that stimulus_onset_adjustment is incremented appropriately, without repeating a logical operator
  generator_input_neurons->select_stimulus(first_stimulus);
  ExampleModel->run(simtime);

  // Loop through a series of stimulus presentations, alternating the chosen stimulus in each example; note due to the above, loop begins at 1
  for (int ii = 1; ii < 10; ++ii) {
    ExampleModel->reset_state(); //Re-set the activity of the network, but not e.g. weights and connectivity
    //On even loops, present stimulus 1
    if (ii%2 == 0) {
      generator_input_neurons->stimulus_onset_adjustment += simtime;
      generator_input_neurons->select_stimulus(first_stimulus);
      ExampleModel->run(simtime);
    }
    else {
      generator_input_neurons->stimulus_onset_adjustment += simtime; //simtime ensures the stimulus begins at the start of the new simulation run, rather than e.g. not running because the spike times for each neuron are in the past
      generator_input_neurons->select_stimulus(second_stimulus);
      ExampleModel->run(simtime);
    }

  }

  spike_monitor_main->save_spikes_as_binary("./", "main_spikes");
  spike_monitor_input->save_spikes_as_binary("./", "input_spikes"); //Save the input neurons spiking activity

  //ExampleModel->spiking_synapses->save_connectivity_as_txt("./");

  return 0;
}