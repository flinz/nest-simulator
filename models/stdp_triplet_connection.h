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

#ifndef STDP_TRIPLET_CONNECTION_H
#define STDP_TRIPLET_CONNECTION_H

/* BeginDocumentation
  Name: stdp_triplet_synapse - Synapse type for spike-timing dependent 
  plasticity accounting for spike triplets as described in [1].

  Description:
   stdp_triplet_synapse is a connector to create synapses with spike time 
   dependent plasticity accounting for spike triplets (as defined in [1]). 
   
   Here, a multiplicative weight dependence is added (in contrast to [1]) 
   to depression resulting in a stable weight distribution.
  
  STDP examples:
   pair-based   Aplus_triplet_ = Aminus_triplet = 0.0
   triplet      Aplus_triplet_ = Aminus_triplet = 1.0
 
  Parameters:
   tau_plus          time constant of STDP window, potentiation
                     (tau_minus defined in post-synaptic neuron)
   tau_x             time constant of triplet potentiation
   tau_y             time constant of triplet depression
   Aplus             weight of pair potentiation rule
   Aminus            weight of pair depression rule
   Aplus_triplet_    weight of triplet potentiation rule
   Aminus_triplet    weight of triplet depression rule

  Transmits: SpikeEvent

  References:
   [1] J.-P. Pfister & W. Gerstner (2006) Triplets of Spikes in a Model 
   of Spike Timing-Dependent Plasticity.  The Journal of Neuroscience 
   26(38):9673-9682; doi:10.1523/JNEUROSCI.1425-06.2006

  FirstVersion: Nov 2007
  Author: Moritz Helias, Abigail Morrison, Eilif Muller, Alex Seeholzer, Teo Stocco
  SeeAlso: synapsedict, stdp_synapse, tsodyks_synapse, static_synapse
*/

#include <cmath>
#include "connection.h"

namespace nest
{
// connections are templates of target identifier type 
// (used for pointer / target index addressing)
// derived from generic connection template
template < typename targetidentifierT >
class STDPTripletConnection : public Connection< targetidentifierT >
{ 

public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  STDPTripletConnection();

  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  STDPTripletConnection(const STDPTripletConnection &);

  /**
   * Default Destructor.
   */
  ~STDPTripletConnection() {}

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
  void send( Event& e, thread t, double_t t_lastspike, const CommonSynapseProperties& cp );

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
   * This function calls check_connection on the sender and checks if the receiver
   * accepts the event type and receptor type requested by the sender.
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
    double_t t_lastspike,
    const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;

    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    t.register_stdp_connection( t_lastspike - get_delay() );
  }

  void
  set_weight( double_t w )
  {
    weight_ = w;
  } 

 private:
  inline double_t // TBD
  facilitate_( double_t w, double_t kplus, double_t ky )
  {
    return w + kplus * ( Aplus_ + Aplus_triplet_ * ky );
  }

  inline double_t
  depress_( double_t w, double_t kminus, double_t Kplus_triplet_ )
  {
    double new_w = w - kminus * ( Aminus_ + Aminus_triplet_ * Kplus_triplet_ );
    return new_w > 0.0 ? new_w : 0.0; // TBD
  } // max weight

  // data members of each connection
  double_t weight_;
  double_t tau_plus_;
  double_t tau_plus_triplet;
  double_t Aplus_;
  double_t Aminus_;
  double_t Aplus_triplet__;
  double_t Aminus_triplet_;
  double_t Kplus_;
  double_t Kplus_triplet__;
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
STDPTripletConnection< targetidentifierT >::send( Event& e,
  thread t,
  double_t t_lastspike,
  const CommonSynapseProperties& )
{
  // synapse STDP depressing/facilitation dynamics
  
  double_t t_spike = e.get_stamp().get_ms();
  // t_lastspike_ = 0 initially

  // use accessor functions (inherited from Connection< >) to obtain delay and target
  Node* target = get_target( t );
  double_t dendritic_delay = get_delay(); 
    
  //get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque<histentry>::iterator start;
  std::deque<histentry>::iterator finish;    
  target->get_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay, &start, &finish);
  //facilitation due to post-synaptic spikes since last pre-synaptic spike
  double_t minus_dt;
  double_t ky;
  while (start != finish)
  {      
    // post-synaptic spike is delayed by dendritic_delay so that
    // it is effectively late by that much at the synapse.
    minus_dt = t_lastspike - (start->t_ + dendritic_delay);
    // subtract 1.0 yields the triplet_Kminus value just prior to
    // the post synaptic spike, implementing the t-epsilon in 
    // Pfister et al, 2006
    ky = start->triplet_Kminus_ - 1.0;
    ++start;
	  
	if ( minus_dt == 0 )
	{
      continue;
	}
	  
    weight_ = facilitate_( weight_, Kplus_ * std::exp( minus_dt / tau_plus_ ), ky );
  }

  //depression due to new pre-synaptic spike
  Kplus_triplet__ *= std::exp( ( t_lastspike - t_spike ) / tau_plus_triplet );

  // dendritic delay means we must look back in time by that amount
  // for determining the K value, because the K value must propagate
  // out to the synapse
  weight_ = depress_( weight_, target->get_K_value( t_spike - dendritic_delay ), Kplus_triplet__ );

  Kplus_triplet__ += 1.0;
  Kplus_ = Kplus_ * std::exp( ( t_lastspike - t_spike ) / tau_plus_ ) + 1.0;
	
  e.set_receiver( *target );
  e.set_weight( weight_ );
  // use accessor functions (inherited from Connection< >) to obtain delay in steps and rport
  e.set_delay( get_delay_steps() );
  e.set_rport( get_rport() );
  e();
}

// Defaults come from reference [1] data fitting and table 3.
template < typename targetidentifierT >
STDPTripletConnection< targetidentifierT >::STDPTripletConnection()
  : ConnectionBase()
  , weight_( 1.0 )
  , tau_plus_( 16.8 )
  , tau_plus_triplet( 101.0 )
  , Aplus_( 5e-10 )
  , Aminus_( 7e-3 )
  , Aplus_triplet__( 6.2e-3 )
  , Aminus_triplet_( 2.3e-4 )
  , Kplus_( 0.0 )
  , Kplus_triplet__( 0.0 )
{
}

template < typename targetidentifierT >
STDPTripletConnection< targetidentifierT >::STDPTripletConnection(
  const STDPTripletConnection< targetidentifierT >& rhs )
  : ConnectionBase(rhs)
  , weight_( rhs.weight_ )
  , tau_plus_( rhs.tau_plus_ )
  , tau_plus_triplet( rhs.tau_plus_triplet )
  , Aplus_( rhs.Aplus_ )
  , Aminus_( rhs.Aminus_ )
  , Aplus_triplet__( rhs.Aplus_triplet__ )
  , Aminus_triplet_( rhs.Aminus_triplet_)
  , Kplus_( rhs.Kplus_ )
  , Kplus_triplet__( rhs.Kplus_triplet__ )
{
}

template < typename targetidentifierT >
void
STDPTripletConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double_t >( d, names::weight, weight_ );
  def< double_t >( d, "tau_plus", tau_plus_ );
  def< double_t >( d, "tau_x", tau_plus_triplet );
  def< double_t >( d, "Aplus", Aplus_ );
  def< double_t >( d, "Aminus", Aminus_ );
  def< double_t >( d, "Aplus_triplet_", Aplus_triplet__ );
  def< double_t >( d, "Aminus_triplet", Aminus_triplet_ );
  def< double_t >( d, "Kplus", Kplus_ );
  def<double_t>( d, "Kplus_triplet_", Kplus_triplet__ );
} // TBD names ?

template < typename targetidentifierT >
void
STDPTripletConnection< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double_t >( d, names::weight, weight_ );
  updateValue< double_t >( d, "tau_plus", tau_plus_ );
  updateValue< double_t >( d, "tau_x", tau_plus_triplet );
  updateValue< double_t >( d, "Aplus", Aplus_ );
  updateValue< double_t >( d, "Aminus", Aminus_ );
  updateValue< double_t >( d, "Aplus_triplet_", Aplus_triplet__ );
  updateValue< double_t >( d, "Aminus_triplet", Aminus_triplet_ );
  updateValue< double_t >( d, "Kplus", Kplus_ );
  updateValue<double_t>( d, "Kplus_triplet_", Kplus_triplet__ );
}

} // of namespace nest

#endif // of #ifndef STDP_TRIPLET_CONNECTION_H
