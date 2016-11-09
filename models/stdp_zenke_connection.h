/*
 *  stdp_triplet_connection.h
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

#ifndef STDP_ZENKE_CONNECTION_H
#define STDP_ZENKE_CONNECTION_H

/* BeginDocumentation
  Name: stdp_zenke_synapse - Synapse type with spike-timing dependent
                               plasticity (triplets).

  Description:
    stdp_zenke_synapse is a connection with spike time dependent
    plasticity accounting for spike triplet effects (as defined in [1]).

  STDP examples:
    pair-based   Aplus_triplet = Aminus_triplet = 0.0
    triplet      Aplus_triplet = Aminus_triplet = 1.0

  Parameters:
    tau_plus           double - time constant of short presynaptic trace
                              - (tau_plus of [1])
    tau_plus_triplet   double - time constant of long presynaptic trace
                              - (tau_x of [1])
    Aplus              double - weight of pair potentiation rule
                              - (A_plus_2 of [1])
    Aplus_triplet      double - weight of triplet potentiation rule
                              - (A_plus_3 of [1])
    Aminus             double - weight of pair depression rule
                                (A_minus_2 of [1])
    Aminus_triplet     double - weight of triplet depression rule
                              - (A_minus_3 of [1])
    Wmax               double - maximum allowed weight

  States:
    Kplus              double: pre-synaptic trace (r_1 of [1])
    Kplus_triplet      double: triplet pre-synaptic trace (r_2 of [1])

  Transmits: SpikeEvent

  References:
    [1] J.-P. Pfister & W. Gerstner (2006) Triplets of Spikes in a Model
        of Spike Timing-Dependent Plasticity.  The Journal of Neuroscience
        26(38):9673-9682; doi:10.1523/JNEUROSCI.1425-06.2006

  Notes:
    - Presynaptic traces r_1 and r_2 of [1] are stored in the connection as
      Kplus and Kplus_triplet and decay with time-constants tau_plus and
      tau_plus_triplet, respectively.
    - Postsynaptic traces o_1 and o_2 of [1] are acquired from the post-synaptic
      neuron states Kminus_ and triplet_Kminus_ which decay on time-constants
      tau_minus and tau_minus_triplet, respectively. These two time-constants
      can be set as properties of the postsynaptic neuron.
    - This version implements the 'all-to-all' spike interaction of [1]. The
      'nearest-spike' interaction of [1] can currently not be implemented
      without changing the postsynaptic archiving-node (clip the traces to a
      maximum of 1).

  FirstVersion: Nov 2007
  Author: Abigail Morrison, Eilif Muller, Alexander Seeholzer, Teo Stocco
  Adapted by: Philipp Weidel
  SeeAlso: stdp_zenke_synapse_hpc, synapsedict, stdp_synapse, static_synapse
*/

// C-header for math.h since copysign() is in C99 but not C++98
#include <math.h>
#include "connection.h"

namespace nest
{
// connections are templates of target identifier type
// (used for pointer / target index addressing)
// derived from generic connection template
template < typename targetidentifierT >
class STDPZenkeConnection : public Connection< targetidentifierT >
{

public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  STDPZenkeConnection();

  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  STDPZenkeConnection( const STDPZenkeConnection& );

  /**
   * Default Destructor.
   */
  ~STDPZenkeConnection()
  {
  }

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places
  // these functions are used. Since ConnectionBase depends on the template
  // parameter, they are not automatically found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_delay;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param t_lastspike Point in time of last spike sent.
   * \param cp common properties of all synapses (empty).
   */
  void send( Event& e,
    thread t,
    double t_lastspike,
    const CommonSynapseProperties& cp );

  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport )
    {
      return invalid_port_;
    }
  };

  /*
   * This function calls check_connection on the sender and checks if the
   * receiver accepts the event type and receptor type requested by the sender.
   * Node::check_connection() will either confirm the receiver port by returning
   * true or false if the connection should be ignored.
   * We have to override the base class' implementation, since for STDP
   * connections we have to call register_stdp_connection on the target neuron
   * to inform the Archiver to collect spikes for this connection.
   *
   * \param s The source node
   * \param r The target node
   * \param receptor_type The ID of the requested receptor type
   * \param t_lastspike last spike emitted by presynaptic neuron
   */
  void
  check_connection( Node& s,
    Node& t,
    rport receptor_type,
    double t_lastspike,
    const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;

    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    t.register_stdp_connection( t_lastspike - get_delay() );
  }

  void
  set_weight( double w )
  {
    weight_ = w;
  }

private:

  // Parameters
  double A_;
  double P_;
  double beta_;
  double delta_;
  double Wmax_;
  double dt_fast_;
  double dt_slow_;  

  // States
  double weight_;
  double w_tilde_;
  double w_tilde_last_t_;
  double C_;
  double z_j_plus_;
  double z_i_ht_;
  
  // Timescales
  double tau_plus_;
  double tau_slow_;
  double tau_ht_;
  double tau_hom_;
  
};

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param t The thread on which this connection is stored.
 * \param t_lastspike Time point of last spike emitted
 * \param cp Common properties object, containing the stdp parameters.
 */
template < typename targetidentifierT >
inline void
STDPZenkeConnection< targetidentifierT >::send( Event& e,
  thread t,
  double t_lastspike,
  const CommonSynapseProperties& )
{

  double t_spike = e.get_stamp().get_ms();
  double dendritic_delay = get_delay();
  Node* target = get_target( t );

  // get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque< histentry >::iterator start;
  std::deque< histentry >::iterator finish;
  
  // integration of synapse state starts from the last spike received
  long t_last_postspike = t_lastspike + dendritic_delay;

  target->get_history(
    t_lastspike - dendritic_delay, t_spike - dendritic_delay, &start, &finish );

  // facilitation due to post-synaptic spikes since last pre-synaptic spike
  while ( start != finish )
  {
    // post-synaptic spike is delayed by dendritic_delay so that
    // it is effectively late by that much at the synapse.
    long delta = Time( Time::ms( start->t_ + dendritic_delay - t_last_postspike ) ).get_steps();

    long_t delta_done = 0;
    while ( delta_done < delta )
    {
      // how many steps are left to be processed?
      long_t delta_this = delta - delta_done;

      // state variable integration, only for existing contacts
      if ( delta_this > 0 )
      {

        // increment the step counter of the loop over delta
        delta_done += delta_this;
      }
    }

    t_last_postspike = start->t_ + dendritic_delay;
    ++start;
  }

  long remaining_delta = Time( Time::ms( t_spike - t_last_postspike ) ).get_steps();

  // integrate remaining
  
  e.set_receiver( *target );
  e.set_weight( weight_ );
  e.set_delay( get_delay_steps() );
  e.set_rport( get_rport() );
  e();
}

// Defaults come from reference [1] data fitting and table 3.
template < typename targetidentifierT >
STDPZenkeConnection< targetidentifierT >::STDPZenkeConnection()
  : ConnectionBase()
  , A_( 1.0 )
  , P_( 1.0 )
  , beta_( 1.0 )
  , delta_( 1.0 )
  , Wmax_( 100.0 )
  , dt_slow_( 1200.0 )
  , dt_fast_( 0.1 )

  , weight_( 1.0 )
  , w_tilde_( 1.0 )
  , w_tilde_last_t_( 0. )
  , C_( 1.0 )
  , z_j_plus_( 1.0 )
  , z_i_ht_( 1.0 )
  
  , tau_plus_( 100.0 )
  , tau_slow_( 100.0 )
  , tau_ht_( 100.0 )
  , tau_hom_( 100.0 )
{
}

template < typename targetidentifierT >
STDPZenkeConnection< targetidentifierT >::STDPZenkeConnection(
  const STDPZenkeConnection< targetidentifierT >& rhs )
  : ConnectionBase( rhs )
  , A_( rhs.A_ )
  , P_( rhs.P_ )
  , beta_( rhs.beta_ )
  , delta_( rhs.delta_ )
  , Wmax_( rhs.Wmax_ )
  , dt_slow_( rhs.dt_slow_ )
  , dt_fast_( rhs.dt_fast_ )
  , weight_( rhs.weight_ )
  , w_tilde_( rhs.w_tilde_ )
  , w_tilde_last_t_( rhs.w_tilde_last_t_ )
  , C_( rhs.C_ )
  , z_j_plus_( rhs.z_j_plus_ )
  , z_i_ht_( rhs.z_i_ht_ )
  , tau_plus_( rhs.tau_plus_ )
  , tau_slow_( rhs.tau_slow_ )
  , tau_ht_( rhs.tau_ht_ )
  , tau_hom_( rhs.tau_hom_ )
{
}

template < typename targetidentifierT >
void
STDPZenkeConnection< targetidentifierT >::get_status(
  DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, "A", A_ );
  def< double >( d, "P", P_ );
  def< double >( d, "beta", beta_ );
  def< double >( d, "delta", delta_ );
  def< double >( d, "Wmax", Wmax_ );
  def< double >( d, "dt_slow_", dt_slow_ );
  def< double >( d, "dt_fast_", dt_fast_ );
  def< double >( d, names::weight, weight_ );
  def< double >( d, "w_tilde", w_tilde_ );
  def< double >( d, "w_tilde_last_t", w_tilde_last_t_ );
  def< double >( d, "C", C_ );
  def< double >( d, "z_j_plus", z_j_plus_ );
  def< double >( d, "z_i_ht", z_i_ht_ );
  def< double >( d, "tau_plus", tau_plus_ );
  def< double >( d, "tau_slow", tau_slow_ );
  def< double >( d, "tau_ht", tau_ht_ );
  def< double >( d, "tau_hom", tau_hom_ );
}

template < typename targetidentifierT >
void
STDPZenkeConnection< targetidentifierT >::set_status(
  const DictionaryDatum& d,
  ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, "A", A_ );
  updateValue< double >( d, "P", P_ );
  updateValue< double >( d, "beta", beta_ );
  updateValue< double >( d, "delta", delta_ );
  updateValue< double >( d, "Wmax", Wmax_ );
  updateValue< double >( d, "dt_slow_", dt_slow_ );
  updateValue< double >( d, "dt_fast_", dt_fast_ );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, "w_tilde", w_tilde_ );
  updateValue< double >( d, "w_tilde_last_t", w_tilde_last_t_ );
  updateValue< double >( d, "C", C_ );
  updateValue< double >( d, "z_j_plus", z_j_plus_ );
  updateValue< double >( d, "z_i_ht", z_i_ht_ );
  updateValue< double >( d, "tau_plus", tau_plus_ );
  updateValue< double >( d, "tau_slow", tau_slow_ );
  updateValue< double >( d, "tau_ht", tau_ht_ );
  updateValue< double >( d, "tau_hom", tau_hom_ );

  // check if weight_ and Wmax_ has the same sign
  if ( not( ( ( weight_ >= 0 ) - ( weight_ < 0 ) )
         == ( ( Wmax_ >= 0 ) - ( Wmax_ < 0 ) ) ) )
  {
    throw BadProperty( "Weight and Wmax must have same sign." );
  }

  // check other states
}

} // of namespace nest

#endif // of #ifnupdateValue STDP_TRIPLET_CONNECTION_H
