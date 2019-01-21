/*
 *  ht_neuron_compte.cpp
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

#include "ht_neuron_compte.h"

#ifdef HAVE_GSL_1_11

#include "universal_data_logger_impl.h"

#include <cmath>


namespace nest
{

RecordablesMap< ht_neuron_compte > ht_neuron_compte::recordablesMap_;

template <>
void
RecordablesMap< ht_neuron_compte >::create()
{
  insert_( "V_m", &ht_neuron_compte::get_y_elem_< ht_neuron_compte::State_::VM > );
  insert_( Name( "g_AMPA" ), &ht_neuron_compte::get_y_elem_< ht_neuron_compte::State_::G_AMPA > );
  insert_( Name( "g_NMDA" ), &ht_neuron_compte::get_y_elem_< ht_neuron_compte::State_::G_NMDA > );
  insert_( Name( "g_GABA" ), &ht_neuron_compte::get_y_elem_< ht_neuron_compte::State_::G_GABA > );
}

/* ----------------------------------------------------------------
 * Iteration function
 * ---------------------------------------------------------------- */

extern "C" inline int
ht_neuron_compte_dynamics( double, const double y[], double f[], void* pnode )
{

  // shorthand
  typedef nest::ht_neuron_compte::State_ S;

  // get access to node so we can almost work as in a member class
  assert( pnode );
  nest::ht_neuron_compte& node = *( reinterpret_cast< nest::ht_neuron_compte* >( pnode ) );

  // easier access to membrane potential
  const double_t& V = y[ S::VM ];

  // Synaptic channels
  double_t I_syn = 0;

  // Calculate sum of all synaptic channels.
  // Sign convention: For each current, write I = - g * ( V - E )
  //    then dV/dt ~ Sum(I) / C_m

  // Compte 2000, AMPA has maximal conductance and exponentially decaying gating variable which is
  // increased on spike arrival by 1
  I_syn += -y[ S::G_AMPA ] * node.P_.AMPA_g_peak * ( V - node.P_.AMPA_E_rev );

  // Compte 2000, DG_NMDA=s, G_NMDA=x=spike receiver
  I_syn += -y[ S::G_NMDA ] * node.P_.NMDA_g_peak * ( V - node.P_.NMDA_E_rev );
  I_syn += -y[ S::G_GABA ] * node.P_.GABA_g_peak * ( V - node.P_.GABA_E_rev );

  // leak currents
  const double_t I_L = -node.P_.g_L * ( V - node.P_.E_L );

  // decays
  f[ S::VM ] = ( I_L + I_syn + node.B_.I_stim_ ) / node.P_.C_m;
  f[ S::G_AMPA ] = -y[ S::G_AMPA ] / node.P_.AMPA_Tau_1;
  f[ S::G_NMDA ] = -y[ S::G_NMDA ] / node.P_.NMDA_Tau_1;
  f[ S::G_GABA ] = -y[ S::G_GABA ] / node.P_.GABA_Tau_1;

  return GSL_SUCCESS;
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::ht_neuron_compte::Parameters_::Parameters_()
  : E_L( -60.0 )
  , // 30 mV
  g_L( 1.0 )
  , // 0.2
  C_m( 10.0 )
  , // ms
  AMPA_g_peak( 3.1 )
  , // nS
  AMPA_Tau_1( 2.0 )
  , // ms
  AMPA_E_rev( 0.0 )
  , // mV
  NMDA_g_peak( 0.075 )
  , NMDA_Tau_1( 80.0 )
  , // ms
  NMDA_E_rev( 0.0 )
  , // mV
  GABA_g_peak( 0.33 )
  , GABA_Tau_1( 10.0 )
  , // mss
  GABA_E_rev( -70.0 )
  , // mV
  V_reset_( -50.0 )
  ,             // mV
  t_ref_( 2.0 ) // ms
{
}

nest::ht_neuron_compte::State_::State_()
  : r_( 0 )
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = 0;
  }
}

nest::ht_neuron_compte::State_::State_( const Parameters_& p )
  : r_( 0 )
{
  y_[ VM ] = p.E_L;
  for ( size_t i = 2; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = 0.0;
  }
}

nest::ht_neuron_compte::State_::State_( const State_& s )
  : r_( s.r_ )
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = s.y_[ i ];
  }
}

nest::ht_neuron_compte::State_& nest::ht_neuron_compte::State_::operator=( const State_& s )
{
  if ( this == &s )
  {
    return *this;
  }
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = s.y_[ i ];
  }
  r_ = s.r_;
  return *this;
}

nest::ht_neuron_compte::State_::~State_()
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::ht_neuron_compte::Parameters_::get( DictionaryDatum& d ) const
{
  def< double_t >( d, "E_L", E_L );
  def< double_t >( d, "g_L", g_L );
  def< double_t >( d, "C_m", C_m );
  def< double_t >( d, "AMPA_g_peak", AMPA_g_peak );
  def< double_t >( d, "AMPA_Tau_1", AMPA_Tau_1 );
  def< double_t >( d, "AMPA_E_rev", AMPA_E_rev );
  def< double_t >( d, "NMDA_g_peak", NMDA_g_peak );
  def< double_t >( d, "NMDA_Tau_1", NMDA_Tau_1 );
  def< double_t >( d, "NMDA_E_rev", NMDA_E_rev );
  def< double_t >( d, "GABA_g_peak", GABA_g_peak );
  def< double_t >( d, "GABA_Tau_1", GABA_Tau_1 );
  def< double_t >( d, "GABA_E_rev", GABA_E_rev );
  def< double_t >( d, "V_reset", V_reset_ );
  def< double_t >( d, "t_ref", t_ref_ );
}

void
nest::ht_neuron_compte::Parameters_::set( const DictionaryDatum& d )
{
  updateValue< double_t >( d, "E_L", E_L );
  updateValue< double_t >( d, "g_L", g_L );
  updateValue< double_t >( d, "C_m", C_m );
  updateValue< double_t >( d, "AMPA_g_peak", AMPA_g_peak );
  updateValue< double_t >( d, "AMPA_Tau_1", AMPA_Tau_1 );
  updateValue< double_t >( d, "AMPA_E_rev", AMPA_E_rev );
  updateValue< double_t >( d, "NMDA_g_peak", NMDA_g_peak );
  updateValue< double_t >( d, "NMDA_Tau_1", NMDA_Tau_1 );
  updateValue< double_t >( d, "NMDA_E_rev", NMDA_E_rev );
  updateValue< double_t >( d, "GABA_g_peak", GABA_g_peak );
  updateValue< double_t >( d, "GABA_Tau_1", GABA_Tau_1 );
  updateValue< double_t >( d, "GABA_E_rev", GABA_E_rev );
  updateValue< double_t >( d, "V_reset", V_reset_ );
  updateValue< double_t >( d, "t_ref", t_ref_ );
}

void
nest::ht_neuron_compte::State_::get( DictionaryDatum& d ) const
{
  def< double_t >( d, "V_m", y_[ VM ] ); // Membrane potential
}

void
nest::ht_neuron_compte::State_::set( const DictionaryDatum& d, const Parameters_& )
{
  updateValue< double_t >( d, "V_m", y_[ VM ] );
}

nest::ht_neuron_compte::Buffers_::Buffers_( ht_neuron_compte& n )
  : logger_( n )
  , spike_inputs_( std::vector< RingBuffer >( SUP_SPIKE_RECEPTOR - 1 ) )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

nest::ht_neuron_compte::Buffers_::Buffers_( const Buffers_&, ht_neuron_compte& n )
  : logger_( n )
  , spike_inputs_( std::vector< RingBuffer >( SUP_SPIKE_RECEPTOR - 1 ) )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest::ht_neuron_compte::ht_neuron_compte()
  : Archiving_Node()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  recordablesMap_.create();
}

nest::ht_neuron_compte::ht_neuron_compte( const ht_neuron_compte& n )
  : Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

nest::ht_neuron_compte::~ht_neuron_compte()
{
  // GSL structs may not be initialized, so we need to protect destruction.
  if ( B_.e_ )
  {
    gsl_odeiv_evolve_free( B_.e_ );
  }
  if ( B_.c_ )
  {
    gsl_odeiv_control_free( B_.c_ );
  }
  if ( B_.s_ )
  {
    gsl_odeiv_step_free( B_.s_ );
  }
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
nest::ht_neuron_compte::init_state_( const Node& proto )
{
  const ht_neuron_compte& pr = downcast< ht_neuron_compte >( proto );
  S_ = pr.S_;
}

void
nest::ht_neuron_compte::init_buffers_()
{
  // Reset spike buffers.
  for ( std::vector< RingBuffer >::iterator it = B_.spike_inputs_.begin();
        it != B_.spike_inputs_.end();
        ++it )
  {
    it->clear(); // include resize
  }

  B_.currents_.clear(); // include resize

  B_.logger_.reset();

  Archiving_Node::clear_history();

  B_.step_ = Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

  if ( B_.s_ == 0 )
  {
    B_.s_ = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_step_reset( B_.s_ );
  }

  if ( B_.c_ == 0 )
  {
    B_.c_ = gsl_odeiv_control_y_new( 1e-3, 0.0 );
  }
  else
  {
    gsl_odeiv_control_init( B_.c_, 1e-3, 0.0, 1.0, 0.0 );
  }

  if ( B_.e_ == 0 )
  {
    B_.e_ = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_evolve_reset( B_.e_ );
  }

  B_.sys_.function = ht_neuron_compte_dynamics;
  B_.sys_.jacobian = 0;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );

  B_.I_stim_ = 0.0;
}

void
nest::ht_neuron_compte::calibrate()
{
  B_.logger_.init(); // ensures initialization in case mm connected after Simulate

  V_.cond_steps_.resize( SUP_SPIKE_RECEPTOR - 1 );

  V_.cond_steps_[ AMPA - 1 ] =
    1.; // Compte 2000, AMPA gating variable increases by 1 at spike arrival

  V_.cond_steps_[ NMDA - 1 ] =
    1.; // Compte 2000, NMDA gating variable x is increases by 1 at spike arrival

  V_.cond_steps_[ GABA - 1 ] =
    1.; // Compte 2000, GABA gating variable x is increases by 1 at spike arrival

  V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref_ ) ).get_steps();
  assert( V_.RefractoryCounts_ >= 0 );
}

void
nest::ht_neuron_compte::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d );
  Archiving_Node::get_status( d );

  DictionaryDatum receptor_type = new Dictionary();

  ( *receptor_type )[ "AMPA" ] = AMPA;
  ( *receptor_type )[ "NMDA" ] = NMDA;
  ( *receptor_type )[ "GABA" ] = GABA;

  ( *d )[ "receptor_types" ] = receptor_type;
  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

void
nest::ht_neuron_compte::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors
  ptmp.set( d );         // throws if BadProperty
  State_ stmp = S_;      // temporary copy in case of errors
  stmp.set( d, ptmp );   // throws if BadProperty

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  Archiving_Node::set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void
ht_neuron_compte::update( Time const& origin, const long_t from, const long_t to )
{
  assert( to >= 0 && ( delay ) from < Scheduler::get_min_delay() );
  assert( from < to );

  for ( long_t lag = from; lag < to; ++lag )
  {
    double tt = 0.0; // it's all relative!

    // adaptive step integration
    while ( tt < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply( B_.e_,
        B_.c_,
        B_.s_,
        &B_.sys_,             // system of ODE
        &tt,                  // from t...
        B_.step_,             // ...to t=t+h
        &B_.IntegrationStep_, // integration window (written on!)
        S_.y_ );              // neuron state

      if ( status != GSL_SUCCESS )
      {
        throw GSLSolverFailure( get_name(), status );
      }
    }

    // Add new spikes to node state array
    for ( size_t i = 0; i < B_.spike_inputs_.size(); ++i )
    {
      // state variable indices are 0=VM, G_AMPA, G_NMDA, G_GABA
      // spike inputs are 0=G_AMPA, G_NMDA, G_GABA
      // therefore we skip the first state variable by setting i+1
      S_.y_[ i + 1 ] += V_.cond_steps_[ i ] * B_.spike_inputs_[ i ].get_value( lag );
    }

    // A spike is generated when the membrane potential (V) exceeds
    // the threshold (Theta).
    if ( S_.r_ )
    { // neuron is absolute refractory
      --S_.r_;
      S_.y_[ State_::VM ] = P_.V_reset_;
    }
    else
      // HARD CODED THRESHOLD
      if ( S_.y_[ State_::VM ] >= -50. )
    {
      // Set Vm to the reset potential
      S_.y_[ State_::VM ] = P_.V_reset_;

      // Set neuron refractory
      S_.r_ = V_.RefractoryCounts_;

      set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

      SpikeEvent se;
      network()->send( *this, se, lag );
    }

    // set new input current
    B_.I_stim_ = B_.currents_.get_value( lag );

    B_.logger_.record_data( origin.get_steps() + lag );
  }
}

void
nest::ht_neuron_compte::handle( SpikeEvent& e )
{
  assert( e.get_delay() > 0 );
  assert( e.get_rport() < static_cast< int_t >( B_.spike_inputs_.size() ) );

  B_.spike_inputs_[ e.get_rport() ].add_value(
    e.get_rel_delivery_steps( network()->get_slice_origin() ),
    e.get_weight() * e.get_multiplicity() );
}

void
nest::ht_neuron_compte::handle( CurrentEvent& e )
{
  assert( e.get_delay() > 0 );

  const double_t I = e.get_current();
  const double_t w = e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value( e.get_rel_delivery_steps( network()->get_slice_origin() ), w * I );
}

void
nest::ht_neuron_compte::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}
}

#endif // HAVE_GSL
