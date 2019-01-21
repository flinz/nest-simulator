/*
 *  ht_neuron.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef HT_NEURON_COMPTE_H
#define HT_NEURON_COMPTE_H

#include "archiving_node.h"
#include <vector>
#include <string>
#include "stringdatum.h"

#ifdef HAVE_GSL_1_11

#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/* BeginDocumentation
   Name: ht_neuron_compte - Neuron model modified from ht_neuron
  for Seeholzer et al (2019). It is inspired by the model by
   Compte (2000), albeit with linear synapses and no voltage dependence
   of the NMDA synapses.

   Description:
   This model neuron implements the neuron model described in [1], and was based
   of the implementation of the ht_neuron.

   The most important properties are:

   - Integrate-and-fire with fixed threshold, fixed refractory period, and hard reset
   - AMPA, NMDA, GABA conductance-based synapses with
     exponential time course

   Parameters:
   V_m  -  membrane potential
   g_L, E_L - conductances and reversal potentials leak current

   {AMPA,NMDA,GABA}_{E_rev,g_peak,Tau_1}
   - reversal potentials, peak conductances and time constants for synapses
     (Tau_1: decay time)

   receptor_types - dictionary mapping synapse names to ports on neuron model
   recordables - list of recordable quantities.

   Author: Alexander Seeholzer

   Sends: SpikeEvent

   Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

   FirstVersion: 2015

   References:
   [1] A. Seeholzer, M. Deger, and W. Gerstner, “Stability of working memory in continuous attractor
  networks under the control of short-term plasticity,” bioRxiv, p. 424515, Sep. 2018.

   SeeAlso: ht_neuron
*/

namespace nest
{
/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
 */
extern "C" int ht_neuron_compte_dynamics( double, const double*, double*, void* );

class ht_neuron_compte : public Archiving_Node
{
public:
  ht_neuron_compte();
  ht_neuron_compte( const ht_neuron_compte& );
  ~ht_neuron_compte();

  /**
   * Import sets of overloaded virtual functions.
   * We need to explicitly include sets of overloaded
   * virtual functions into the current scope.
   * According to the SUN C++ FAQ, this is the correct
   * way of doing things, although all other compilers
   * happily live without.
   */

  using Node::handle;
  using Node::handles_test_event;

  port send_test_event( Node&, rport, synindex, bool );

  void handle( SpikeEvent& e );
  void handle( CurrentEvent& e );
  void handle( DataLoggingRequest& );

  port handles_test_event( SpikeEvent&, rport );
  port handles_test_event( CurrentEvent&, rport );
  port handles_test_event( DataLoggingRequest&, rport );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );

private:
  /**
   * Synapse types to connect to
   * @note Excluded upper and lower bounds are defined as INF_, SUP_.
   *       Excluding port 0 avoids accidental connections.
   */
  enum SynapseTypes
  {
    INF_SPIKE_RECEPTOR = 0,
    AMPA,
    NMDA,
    GABA,
    SUP_SPIKE_RECEPTOR
  };

  void init_state_( const Node& proto );
  void init_buffers_();
  void calibrate();

  void update( Time const&, const long_t, const long_t );

  // END Boilerplate function declarations ----------------------------

  // Friends --------------------------------------------------------

  // make dynamics function quasi-member
  friend int ht_neuron_compte_dynamics( double, const double*, double*, void* );

  // ----------------------------------------------------------------

  /**
   * Independent parameters of the model.
   */
  struct Parameters_
  {
    // Leaks
    double_t E_L;
    double_t g_L;
    double_t C_m;

    double_t V_reset_; //!< Reset Potential in mV
    double_t t_ref_;   //!< Refractory period in ms

    Parameters_();

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
    void set( const DictionaryDatum& ); //!< Set values from dicitonary

    // Parameters for synapse of type AMPA, GABA, GABA_B and NMDA
    double_t AMPA_g_peak;
    double_t AMPA_Tau_1; // ms
    double_t AMPA_E_rev; // mV

    double_t NMDA_g_peak;
    double_t NMDA_Tau_1; // ms
    double_t NMDA_E_rev; // mV

    double_t GABA_g_peak;
    double_t GABA_Tau_1; // ms
    double_t GABA_E_rev; // mV
  };

  // ----------------------------------------------------------------

  /**
   * State variables of the model.
   */
public:
  struct State_
  {

    // y_ = [V, Theta, Synapses]
    enum StateVecElems_
    {
      VM = 0,
      G_AMPA,
      G_NMDA,
      G_GABA,
      STATE_VEC_SIZE
    };

    double_t y_[ STATE_VEC_SIZE ]; //!< neuron state, must be C-array for GSL solver
    int_t r_;                      //!< number of refractory steps remaining

    State_();
    State_( const Parameters_& p );
    State_( const State_& s );
    ~State_();

    State_& operator=( const State_& s );

    void get( DictionaryDatum& ) const;
    void set( const DictionaryDatum&, const Parameters_& );
  };

private:
  // These friend declarations must be precisely here.
  friend class RecordablesMap< ht_neuron_compte >;
  friend class UniversalDataLogger< ht_neuron_compte >;


  // ----------------------------------------------------------------

  /**
   * Buffers of the model.
   */
  struct Buffers_
  {
    Buffers_( ht_neuron_compte& );
    Buffers_( const Buffers_&, ht_neuron_compte& );

    UniversalDataLogger< ht_neuron_compte > logger_;

    /** buffers and sums up incoming spikes/currents */
    std::vector< RingBuffer > spike_inputs_;
    RingBuffer currents_;

    /** GSL ODE stuff */
    gsl_odeiv_step* s_;    //!< stepping function
    gsl_odeiv_control* c_; //!< adaptive stepsize control function
    gsl_odeiv_evolve* e_;  //!< evolution function
    gsl_odeiv_system sys_; //!< struct describing system

    // IntergrationStep_ should be reset with the neuron on ResetNetwork,
    // but remain unchanged during calibration. Since it is initialized with
    // step_, and the resolution cannot change after nodes have been created,
    // it is safe to place both here.
    double_t step_;          //!< step size in ms
    double IntegrationStep_; //!< current integration time step, updated by GSL

    /**
     * Input current injected by CurrentEvent.
     * This variable is used to transport the current applied into the
     * _dynamics function computing the derivative of the state vector.
     * It must be a part of Buffers_, since it is initialized once before
     * the first simulation, but not modified before later Simulate calls.
     */
    double_t I_stim_;
  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
    //! size of conductance steps for arriving spikes
    std::vector< double_t > cond_steps_;

    //! Duration of potassium current.
    // int_t    PotassiumRefractoryCounts_;
    int_t RefractoryCounts_;
  };


  // readout functions, can use template for vector elements
  template < State_::StateVecElems_ elem >
  double_t
  get_y_elem_() const
  {
    return S_.y_[ elem ];
  }

  static RecordablesMap< ht_neuron_compte > recordablesMap_;

  Parameters_ P_;
  State_ S_;
  Variables_ V_;
  Buffers_ B_;
};

inline port
ht_neuron_compte::send_test_event( Node& target, rport receptor_type, synindex, bool )
{
  SpikeEvent e;
  e.set_sender( *this );

  return target.handles_test_event( e, receptor_type );
}

inline port
ht_neuron_compte::handles_test_event( SpikeEvent&, rport receptor_type )
{
  assert( B_.spike_inputs_.size() == 3 );

  if ( !( INF_SPIKE_RECEPTOR < receptor_type && receptor_type < SUP_SPIKE_RECEPTOR ) )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
    return 0;
  }
  else
  {
    return receptor_type - 1;
  }
}

inline port
ht_neuron_compte::handles_test_event( CurrentEvent&, rport receptor_type )
{

  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline port
ht_neuron_compte::handles_test_event( DataLoggingRequest& dlr, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}
}

#endif // HAVE_GSL
#endif // HT_NEURON_COMPTE_H
