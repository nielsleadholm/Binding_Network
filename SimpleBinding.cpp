
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
  weightdependent_stdp_plasticity_parameters_struct * WDSTDP_PARAMS = new weightdependent_stdp_plasticity_parameters_struct;
  WDSTDP_PARAMS->a_plus = 1.0;
  WDSTDP_PARAMS->a_minus = 1.0;
  WDSTDP_PARAMS->tau_plus = 0.02;
  WDSTDP_PARAMS->tau_minus = 0.02;
  WDSTDP_PARAMS->lambda = 1.0f*powf(10.0, -2);
  WDSTDP_PARAMS->alpha = 2.02;
  WDSTDP_PARAMS->w_max = 0.3*powf(10.0, -3);

  WeightDependentSTDPPlasticity * weightdependent_stdp = new WeightDependentSTDPPlasticity((SpikingSynapses *) conductance_spiking_synapses, (SpikingNeurons *) lif_spiking_neurons, (SpikingNeurons *) patterned_poisson_input_neurons, (stdp_plasticity_parameters_struct *) WDSTDP_PARAMS);  
  
  ExampleModel->AddPlasticityRule(weightdependent_stdp);

  /*
      ADD ANY ACTIVITY MONITORS OR PLASTICITY RULES YOU WISH FOR 
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
  input_neuron_params->group_shape[0] = 1;    // x-dimension of the input neuron layer
  input_neuron_params->group_shape[1] = 27;   // y-dimension of the input neuron layer
  // Create a group of input neurons. This function returns the ID of the input neuron group
  int input_layer_ID = ExampleModel->AddInputNeuronGroup(input_neuron_params);

  // SETTING UP NEURON GROUPS
  // Creating an LIF parameter structure for an excitatory neuron population and an inhibitory
  // 1 x 100 Layer
  lif_spiking_neuron_parameters_struct * excitatory_population_params = new lif_spiking_neuron_parameters_struct();
  excitatory_population_params->group_shape[0] = 1;
  excitatory_population_params->group_shape[1] = 27;
  excitatory_population_params->resting_potential_v0 = -0.074f;
  excitatory_population_params->threshold_for_action_potential_spike = -0.053f;
  excitatory_population_params->somatic_capacitance_Cm = 500.0*pow(10, -12);
  excitatory_population_params->somatic_leakage_conductance_g0 = 25.0*pow(10, -9);

  lif_spiking_neuron_parameters_struct * inhibitory_population_params = new lif_spiking_neuron_parameters_struct();
  inhibitory_population_params->group_shape[0] = 1;
  inhibitory_population_params->group_shape[1] = 6;
  inhibitory_population_params->resting_potential_v0 = -0.082f;
  inhibitory_population_params->threshold_for_action_potential_spike = -0.053f;
  inhibitory_population_params->somatic_capacitance_Cm = 214.0*pow(10, -12);
  inhibitory_population_params->somatic_leakage_conductance_g0 = 18.0*pow(10, -9);

  // Create populations of excitatory and inhibitory neurons
  int first_excitatory_neuron_layer_ID = ExampleModel->AddNeuronGroup(excitatory_population_params);
  int first_inhibitory_neuron_layer_ID = ExampleModel->AddNeuronGroup(inhibitory_population_params);

  // *** Create additional layers - NB that is shares the same properties as the first layer
  int second_excitatory_neuron_layer_ID = ExampleModel->AddNeuronGroup(excitatory_population_params);
  int second_inhibitory_neuron_layer_ID = ExampleModel->AddNeuronGroup(inhibitory_population_params);
  int third_excitatory_neuron_layer_ID = ExampleModel->AddNeuronGroup(excitatory_population_params);
  int third_inhibitory_neuron_layer_ID = ExampleModel->AddNeuronGroup(inhibitory_population_params);

  // SETTING UP SYNAPSES
  // Creating a synapses parameter structure for connections from the input neurons to the excitatory neurons
  conductance_spiking_synapse_parameters_struct* input_to_excitatory_parameters = new conductance_spiking_synapse_parameters_struct();
  input_to_excitatory_parameters->weight_range[0] = 0.5f;   // Create uniform distributions of weights [0.5, 10.0]
  input_to_excitatory_parameters->weight_range[1] = 10.0f;
  input_to_excitatory_parameters->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  input_to_excitatory_parameters->delay_range[0] = 10.0*timestep;   //Delays range from 1 to 10 ms for excitatory connectivity
  input_to_excitatory_parameters->delay_range[1] = 100.0*timestep;
  input_to_excitatory_parameters->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
  input_to_excitatory_parameters->gaussian_synapses_standard_deviation = 2.0; //connects neurons with specified Gaussian SD
  input_to_excitatory_parameters->max_number_of_connections_per_pair = 5;
  input_to_excitatory_parameters->gaussian_synapses_per_postsynaptic_neuron = 4;

  // *** Add plasticity to input synapses
  input_to_excitatory_parameters->plasticity_vec.push_back(weightdependent_stdp);


  // Creating a set of synapse parameters for connections from the excitatory neurons to the inhibitory neurons *within a layer*
  conductance_spiking_synapse_parameters_struct * excitatory_to_inhibitory_parameters = new conductance_spiking_synapse_parameters_struct();
  excitatory_to_inhibitory_parameters->weight_range[0] = 10.0f;
  excitatory_to_inhibitory_parameters->weight_range[1] = 10.0f;
  excitatory_to_inhibitory_parameters->weight_scaling_constant = inhibitory_population_params->somatic_leakage_conductance_g0;
  excitatory_to_inhibitory_parameters->delay_range[0] = 10.0*timestep; //Delays range from 1 to 2 ms for inhibitory connectivity
  excitatory_to_inhibitory_parameters->delay_range[1] = 20.0*timestep;
  excitatory_to_inhibitory_parameters->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
  excitatory_to_inhibitory_parameters->gaussian_synapses_standard_deviation = 2.0; //connects neurons with specified Gaussian SD
  excitatory_to_inhibitory_parameters->max_number_of_connections_per_pair = 5;
  excitatory_to_inhibitory_parameters->gaussian_synapses_per_postsynaptic_neuron = 4;

  // Creating a set of synapse parameters from the inhibitory neurons to the excitatory neurons *within a layer*
  conductance_spiking_synapse_parameters_struct * inhibitory_to_excitatory_parameters = new conductance_spiking_synapse_parameters_struct();
  inhibitory_to_excitatory_parameters->weight_range[0] = -5.0f;
  inhibitory_to_excitatory_parameters->weight_range[1] = -2.5f;
  inhibitory_to_excitatory_parameters->weight_scaling_constant = excitatory_population_params->somatic_leakage_conductance_g0;
  inhibitory_to_excitatory_parameters->delay_range[0] = 10.0*timestep; //Delays range from 1 to 2 ms for inhibitory connectivity
  inhibitory_to_excitatory_parameters->delay_range[1] = 20.0*timestep;
  inhibitory_to_excitatory_parameters->connectivity_type = CONNECTIVITY_TYPE_ALL_TO_ALL;
  
  // Creating a set of synapse parameters for connections from the excitatory neurons back to the excitatory neurons *within a layer*
  conductance_spiking_synapse_parameters_struct * excitatory_to_excitatory_parameters = new conductance_spiking_synapse_parameters_struct();
  excitatory_to_excitatory_parameters->weight_range[0] = 10.0f;
  excitatory_to_excitatory_parameters->weight_range[1] = 10.0f;
  excitatory_to_excitatory_parameters->weight_scaling_constant = inhibitory_population_params->somatic_leakage_conductance_g0;
  excitatory_to_excitatory_parameters->delay_range[0] = 10.0*timestep; //Delays range from 1 to 10 ms for excitatory connectivity
  excitatory_to_excitatory_parameters->delay_range[1] = 100.0*timestep;
  excitatory_to_excitatory_parameters->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
  excitatory_to_excitatory_parameters->gaussian_synapses_standard_deviation = 2.0; //connects neurons with specified Gaussian SD
  excitatory_to_excitatory_parameters->max_number_of_connections_per_pair = 5;
  excitatory_to_excitatory_parameters->gaussian_synapses_per_postsynaptic_neuron = 4;

  // Creating a set of synapse parameters for connections from the excitatory neurons *in a lower layer to the layer above*
  conductance_spiking_synapse_parameters_struct * lower_to_upper_parameters = new conductance_spiking_synapse_parameters_struct();
  lower_to_upper_parameters->weight_range[0] = 10.0f;
  lower_to_upper_parameters->weight_range[1] = 10.0f;
  lower_to_upper_parameters->weight_scaling_constant = inhibitory_population_params->somatic_leakage_conductance_g0;
  lower_to_upper_parameters->delay_range[0] = 10.0*timestep; //Delays range from 1 to 10 ms for excitatory connectivity
  lower_to_upper_parameters->delay_range[1] = 100.0*timestep;
  lower_to_upper_parameters->connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
  lower_to_upper_parameters->gaussian_synapses_standard_deviation = 2.0; //connects neurons with specified Gaussian SD
  lower_to_upper_parameters->max_number_of_connections_per_pair = 5;
  lower_to_upper_parameters->gaussian_synapses_per_postsynaptic_neuron = 4;


  // *** Add plasticity to excitatory to excitatory synapses (w/in layers), and the excitatory connections projecting up layers
  excitatory_to_excitatory_parameters->plasticity_vec.push_back(weightdependent_stdp);
  lower_to_upper_parameters->plasticity_vec.push_back(weightdependent_stdp);


  // CREATING SYNAPSES
  // When creating synapses, the ids of the presynaptic and postsynaptic populations are all that are required
  // Note: Input neuron populations cannot be post-synaptic on any synapse
  ExampleModel->AddSynapseGroup(input_layer_ID, first_excitatory_neuron_layer_ID, input_to_excitatory_parameters);
  ExampleModel->AddSynapseGroup(first_excitatory_neuron_layer_ID, first_inhibitory_neuron_layer_ID, excitatory_to_inhibitory_parameters);
  ExampleModel->AddSynapseGroup(first_inhibitory_neuron_layer_ID, first_excitatory_neuron_layer_ID, inhibitory_to_excitatory_parameters);
  ExampleModel->AddSynapseGroup(first_excitatory_neuron_layer_ID, first_excitatory_neuron_layer_ID, excitatory_to_excitatory_parameters);

  // *** Add synapses relevant to the second and third layers
  ExampleModel->AddSynapseGroup(first_excitatory_neuron_layer_ID, second_excitatory_neuron_layer_ID, lower_to_upper_parameters);
  ExampleModel->AddSynapseGroup(second_excitatory_neuron_layer_ID, second_inhibitory_neuron_layer_ID, excitatory_to_inhibitory_parameters);
  ExampleModel->AddSynapseGroup(second_inhibitory_neuron_layer_ID, second_excitatory_neuron_layer_ID, inhibitory_to_excitatory_parameters);
  ExampleModel->AddSynapseGroup(second_excitatory_neuron_layer_ID, second_excitatory_neuron_layer_ID, excitatory_to_excitatory_parameters);

  ExampleModel->AddSynapseGroup(second_excitatory_neuron_layer_ID, third_excitatory_neuron_layer_ID, lower_to_upper_parameters);
  ExampleModel->AddSynapseGroup(third_excitatory_neuron_layer_ID, third_inhibitory_neuron_layer_ID, excitatory_to_inhibitory_parameters);
  ExampleModel->AddSynapseGroup(third_inhibitory_neuron_layer_ID, third_excitatory_neuron_layer_ID, inhibitory_to_excitatory_parameters);
  ExampleModel->AddSynapseGroup(third_excitatory_neuron_layer_ID, third_excitatory_neuron_layer_ID, excitatory_to_excitatory_parameters);


  /*
      ADD INPUT STIMULI TO THE PATTERNED POISSON NEURONS CLASS
  */
  // First stimulus is the 'ascending' pattern; pattern takes place over 10 ms.
  float s1_poisson_rates[27] = {10.0f, 15.0f, 10.0f, 15.0f, 10.0f, 15.0f, 10.0f, 15.0f, 10.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f};
  // Adding this stimulus to the input neurons
  int first_stimulus = patterned_poisson_input_neurons->add_stimulus(s1_poisson_rates, 27);
  // Creating a second stimulus (descending pattern)
 float s2_poisson_rates[27] = {2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 10.0f, 10.0f, 10.0f, 10.0f, 10.0f, 15.0f, 15.0f, 15.0f, 15.0f};
  int second_stimulus = patterned_poisson_input_neurons->add_stimulus(s2_poisson_rates, 27);
  

  /*
      RUN THE SIMULATION
  */

  // The only argument to run is the number of seconds
  ExampleModel->finalise_model();
  float simtime = 0.05f; //This should be long enough to allow any recursive signalling to finish propagating

  // Very first 'loop' is executed outside, such that stimulus_onset_adjustment is incremented appropriately, without repeating a logical operator
  patterned_poisson_input_neurons->select_stimulus(first_stimulus);
  ExampleModel->run(simtime);

  // Loop through a series of stimulus presentations, alternating the chosen stimulus in each example; note due to the above, loop begins at 1
  for (int ii = 1; ii < 10; ++ii) {
    ExampleModel->reset_state(); //Re-set the activity of the network, but not e.g. weights and connectivity
    //On even loops, present stimulus 1
    if (ii%2 == 0) {
      //NB stimulus adjusment is not relevant for Poisson input: patterned_poisson_input_neurons->stimulus_onset_adjustment += simtime;
      patterned_poisson_input_neurons->select_stimulus(first_stimulus);
      ExampleModel->run(simtime);
    }
    else {
      //patterned_poisson_input_neurons->stimulus_onset_adjustment += simtime; //simtime ensures the stimulus begins at the start of the new simulation run, rather than e.g. not running because the spike times for each neuron are in the past
      patterned_poisson_input_neurons->select_stimulus(second_stimulus);
      ExampleModel->run(simtime);
    }

  }

  spike_monitor_main->save_spikes_as_binary("./", "main_spikes");
  spike_monitor_input->save_spikes_as_binary("./", "input_spikes"); //Save the input neurons spiking activity

  //ExampleModel->spiking_synapses->save_connectivity_as_txt("./");

  return 0;
}