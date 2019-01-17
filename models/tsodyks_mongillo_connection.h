/*
 *  tsodyks2_connection.h
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

#ifndef TSODYKS_MONGILLO_CONNECTION_H
#define TSODYKS_MONGILLO_CONNECTION_H


/* BeginDocumentation
  Name: tsodyks_mongillo_synapse - Synapse type with short term plasticity.

  Description:
   This synapse model implements synaptic short-term depression and short-term facilitation
   according to [1].

   Parameters:
     The following parameters can be set in the status dictionary:
     U          double - probability of release increment (U1) [0,1], default=0.5
     u          double - Maximum probability of release (U_se) [0,1], default=0.5
     x          double - current scaling factor of the weight, default=U
     tau_rec    double - time constant for depression in ms, default=800 ms
     tau_fac    double - time constant for facilitation in ms, default=0 (off)

  Remarks:
    It gives the same results as Tsodyks2 synapse [2], but implements
    the equations of [1] with a baseline decay to U for the facilitation variable.
    The main difference is that at a spike, the values BEFORE the spike-triggered update
    are transmitted (x_trans, u_trans), whereas in [2] the updated values are transmitted.

  References:
   [1] G. Mongillo, O. Barak, and M. Tsodyks, “Synaptic theory of working memory.,”
       Science, vol. 319, no. 5869, pp. 1543–1546, Mar. 2008.
   [2] Tsodyks, M. V., & Markram, H. (1997). The neural code between neocortical pyramidal neurons
       depends on neurotransmitter release probability. PNAS, 94(2), 719-23.

  Transmits: SpikeEvent

  FirstVersion: 2017
  Author: Marc-Oliver Gewaltig, based on tsodyks2_synapse by Marc-Oliver Gewaltig
  SeeAlso: tsodyks_synapse, tsodyks2_synapse, synapsedict, stdp_synapse, static_synapse
*/


/**
 * Class representing a synapse with Tsodyks short term plasticity, based on the ODE
 * A suitable Connector containing these connections can be obtained from the template
 * GenericConnector.
 */
#include "connection.h"
#include <cmath>

namespace nest
{

template < typename targetidentifierT >
class TsodyksMongilloConnection : public Connection< targetidentifierT >
{
public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  TsodyksMongilloConnection();

  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  TsodyksMongilloConnection( const TsodyksMongilloConnection& );

  /**
   * Default Destructor.
   */
  ~TsodyksMongilloConnection()
  {
  }

  // Explicitly declare all methods inherited from the dependent base ConnectionBase.
  // This avoids explicit name prefixes in all places these functions are used.
  // Since ConnectionBase depends on the template parameter, they are not automatically
  // found in the base class.
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
   * \param cp Common properties to all synapses (empty).
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


  void
  check_connection( Node& s, Node& t, rport receptor_type, double_t, const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );
  }

  void
  set_weight( double_t w )
  {
    weight_ = w;
  }


private:
  double_t weight_;
  double_t U_;       //!< unit increment of a facilitating synapse
  double_t u_;       //!< dynamic value of probability of release
  double_t x_;       //!< current fraction of the synaptic weight
  double_t tau_rec_; //!< [ms] time constant for recovery
  double_t tau_fac_; //!< [ms] time constant for facilitation
};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
template < typename targetidentifierT >
inline void
TsodyksMongilloConnection< targetidentifierT >::send( Event& e,
  thread t,
  double_t t_lastspike,
  const CommonSynapseProperties& )
{
  Node* target = get_target( t );

  double_t h = e.get_stamp().get_ms() - t_lastspike;
  double_t x_decay = std::exp( -h / tau_rec_ );
  double_t u_decay = ( tau_fac_ < 1.0e-10 ) ? 0.0 : std::exp( -h / tau_fac_ );

  // Changes from Tsodyks2Connection: we transmit the values before the spike-triggered increase,
  // and only decay the state variables.
  
  // These are the values before the per-spike jump, which get transmitted
  double_t u_trans = u_ * u_decay + U_ * (1. - u_decay); // first decay
  double_t x_trans = x_ * x_decay + 1. * (1. - x_decay); // first decay

  // increment the values after the spike which are propagated to next timestep
  u_ = u_trans + U_*(1.0-u_trans);
  x_ = x_trans - u_trans * x_trans;

  // We use the current values for the spike number n.
  e.set_receiver( *target );
  e.set_weight( u_trans * x_trans * weight_ );
  // send the spike to the target
  e.set_delay( get_delay_steps() );
  e.set_rport( get_rport() );
  e();
}

template < typename targetidentifierT >
TsodyksMongilloConnection< targetidentifierT >::TsodyksMongilloConnection()
  : ConnectionBase()
  , weight_( 1.0 )
  , U_( 0.5 )
  , u_( U_ )
  , x_( U_ )
  , tau_rec_( 800.0 )
  , tau_fac_( 0.0 )
{
}

template < typename targetidentifierT >
TsodyksMongilloConnection< targetidentifierT >::TsodyksMongilloConnection( const TsodyksMongilloConnection& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , U_( rhs.U_ )
  , u_( rhs.u_ )
  , x_( rhs.x_ )
  , tau_rec_( rhs.tau_rec_ )
  , tau_fac_( rhs.tau_fac_ )
{
}


template < typename targetidentifierT >
void
TsodyksMongilloConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double_t >( d, names::weight, weight_ );

  def< double_t >( d, names::dU, U_ );
  def< double_t >( d, names::u, u_ );
  def< double_t >( d, names::tau_rec, tau_rec_ );
  def< double_t >( d, names::tau_fac, tau_fac_ );
  def< double_t >( d, names::x, x_ );
  def< long_t >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
TsodyksMongilloConnection< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double_t >( d, names::weight, weight_ );

  updateValue< double_t >( d, names::dU, U_ );
  if ( U_ > 1.0 || U_ < 0.0 )
    throw BadProperty( "U must be in [0,1]." );

  updateValue< double_t >( d, names::u, u_ );
  if ( u_ > 1.0 || u_ < 0.0 )
    throw BadProperty( "u must be in [0,1]." );

  updateValue< double_t >( d, names::tau_rec, tau_rec_ );
  if ( tau_rec_ <= 0.0 )
    throw BadProperty( "tau_rec must be > 0." );

  updateValue< double_t >( d, names::tau_fac, tau_fac_ );
  if ( tau_fac_ < 0.0 )
    throw BadProperty( "tau_fac must be >= 0." );

  updateValue< double_t >( d, names::x, x_ );
}

} // namespace

#endif // TSODYKS_MONGILLO_CONNECTION_H
