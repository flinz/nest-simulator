/*
 *  simulation_manager.cpp
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

#include "simulation_manager.h"

// C includes:
#include <sys/time.h>

// C++ includes:
#include <vector>

// Includes from libnestutil:
#include "compose.hpp"

// Includes from nestkernel:
#include "kernel_manager.h"
#include "sibling_container.h"

// Includes from sli:
#include "dictutils.h"
#include "psignal.h"

nest::SimulationManager::SimulationManager()
  : simulating_( false )
  , clock_( Time::tic( 0L ) )
  , slice_( 0L )
  , to_do_( 0L )
  , to_do_total_( 0L )
  , from_step_( 0L )
  , to_step_( 0L ) // consistent with to_do_ == 0
  , t_real_( 0L )
  , terminate_( false )
  , simulated_( false )
  , print_time_( false )
  , use_wfr_( true )
  , wfr_comm_interval_( 1.0 )
  , wfr_tol_( 0.0001 )
  , wfr_max_iterations_( 15 )
  , wfr_interpolation_order_( 3 )
{
}

void
nest::SimulationManager::initialize()
{
  // set resolution, ensure clock is calibrated to new resolution
  Time::reset_resolution();
  clock_.calibrate();

  simulated_ = false;
}

void
nest::SimulationManager::finalize()
{
  nest::Time::reset_to_defaults();

  clock_.set_to_zero(); // ensures consistent state
  to_do_ = 0;
  slice_ = 0;
  from_step_ = 0;
  to_step_ = 0; // consistent with to_do_ = 0
}

void
nest::SimulationManager::set_status( const DictionaryDatum& d )
{
  // Create an instance of time converter here to capture the current
  // representation of time objects: TICS_PER_MS and TICS_PER_STEP
  // will be stored in time_converter.
  // This object can then be used to convert times in steps
  // (e.g. Connection::delay_) or tics to the new representation.
  // We pass this object to ConnectionManager::calibrate to update
  // all time objects in the connection system to the new representation.
  // MH 08-04-14
  TimeConverter time_converter;

  double time;
  if ( updateValue< double >( d, "time", time ) )
  {
    if ( time != 0.0 )
      throw BadProperty( "The simulation time can only be set to 0.0." );

    if ( clock_ > TimeZero )
    {
      // reset only if time has passed
      LOG( M_WARNING,
        "SimulationManager::set_status",
        "Simulation time reset to t=0.0. Resetting the simulation time is not "
        "fully supported in NEST at present. Some spikes may be lost, and "
        "stimulating devices may behave unexpectedly. PLEASE REVIEW YOUR "
        "SIMULATION OUTPUT CAREFULLY!" );

      clock_ = Time::step( 0 );
      from_step_ = 0;
      slice_ = 0;
      // clear all old spikes
      kernel().event_delivery_manager.configure_spike_buffers();
    }
  }

  updateValue< bool >( d, "print_time", print_time_ );

  // tics_per_ms and resolution must come after local_num_thread /
  // total_num_threads because they might reset the network and the time
  // representation
  double tics_per_ms = 0.0;
  bool tics_per_ms_updated =
    updateValue< double >( d, "tics_per_ms", tics_per_ms );
  double resd = 0.0;
  bool res_updated = updateValue< double >( d, "resolution", resd );

  if ( tics_per_ms_updated || res_updated )
  {
    if ( kernel().node_manager.size() > 1 ) // root always exists
    {
      LOG( M_ERROR,
        "SimulationManager::set_status",
        "Cannot change time representation after nodes have been created. "
        "Please call ResetKernel first." );
      throw KernelException();
    }
    else if ( has_been_simulated() ) // someone may have simulated empty network
    {
      LOG( M_ERROR,
        "SimulationManager::set_status",
        "Cannot change time representation after the network has been "
        "simulated. Please call ResetKernel first." );
      throw KernelException();
    }
    else if ( kernel().connection_manager.get_num_connections() != 0 )
    {
      LOG( M_ERROR,
        "SimulationManager::set_status",
        "Cannot change time representation after connections have been "
        "created. Please call ResetKernel first." );
      throw KernelException();
    }
    else if ( res_updated && tics_per_ms_updated ) // only allow TICS_PER_MS to
                                                   // be changed together with
                                                   // resolution
    {
      if ( resd < 1.0 / tics_per_ms )
      {
        LOG( M_ERROR,
          "SimulationManager::set_status",
          "Resolution must be greater than or equal to one tic. Value "
          "unchanged." );
        throw KernelException();
      }
      else
      {
        nest::Time::set_resolution( tics_per_ms, resd );
        // adjust to new resolution
        clock_.calibrate();
        // adjust delays in the connection system to new resolution
        kernel().connection_manager.calibrate( time_converter );
        kernel().model_manager.calibrate( time_converter );
        LOG( M_INFO,
          "SimulationManager::set_status",
          "tics per ms and resolution changed." );

        // make sure that wfr communication interval is always greater or equal
        // to resolution if no wfr is used explicitly set wfr_comm_interval
        // to resolution because communication in every step is needed
        if ( wfr_comm_interval_ < Time::get_resolution().get_ms()
          || not use_wfr_ )
        {
          wfr_comm_interval_ = Time::get_resolution().get_ms();
        }
      }
    }
    else if ( res_updated ) // only resolution changed
    {
      if ( resd < Time::get_ms_per_tic() )
      {
        LOG( M_ERROR,
          "SimulationManager::set_status",
          "Resolution must be greater than or equal to one tic. Value "
          "unchanged." );
        throw KernelException();
      }
      else
      {
        Time::set_resolution( resd );
        clock_.calibrate(); // adjust to new resolution
        // adjust delays in the connection system to new resolution
        kernel().connection_manager.calibrate( time_converter );
        kernel().model_manager.calibrate( time_converter );
        LOG( M_INFO,
          "SimulationManager::set_status",
          "Temporal resolution changed." );

        // make sure that wfr communication interval is always greater or equal
        // to resolution if no wfr is used explicitly set wfr_comm_interval
        // to resolution because communication in every step is needed
        if ( wfr_comm_interval_ < Time::get_resolution().get_ms()
          || not use_wfr_ )
        {
          wfr_comm_interval_ = Time::get_resolution().get_ms();
        }
      }
    }
    else
    {
      LOG( M_ERROR,
        "SimulationManager::set_status",
        "change of tics_per_step requires simultaneous specification of "
        "resolution." );
      throw KernelException();
    }
  }

  // The decision whether the waveform relaxation is used
  // must be set before nodes are created.
  // Important: wfr_comm_interval_ may change depending on use_wfr_
  bool wfr;
  if ( updateValue< bool >( d, "use_wfr", wfr ) )
  {
    if ( kernel().node_manager.size() > 1 )
    {
      LOG( M_ERROR,
        "SimulationManager::set_status",
        "Cannot enable/disable usage of waveform relaxation after nodes have "
        "been created. Please call ResetKernel first." );
      throw KernelException();
    }
    else
    {
      use_wfr_ = wfr;
      // if no wfr is used explicitly set wfr_comm_interval to resolution
      // because communication in every step is needed
      if ( not use_wfr_ )
      {
        wfr_comm_interval_ = Time::get_resolution().get_ms();
      }
    }
  }

  // wfr_comm_interval_ can only be changed if use_wfr_ is true and before
  // connections are created. If use_wfr_ is false wfr_comm_interval_ is set to
  // the resolution whenever the resolution changes.
  double wfr_interval;
  if ( updateValue< double >( d, "wfr_comm_interval", wfr_interval ) )
  {
    if ( not use_wfr_ )
    {
      LOG( M_ERROR,
        "SimulationManager::set_status",
        "Cannot set waveform communication interval when usage of waveform "
        "relaxation is disabled. Set use_wfr to true first." );
      throw KernelException();
    }
    else if ( kernel().connection_manager.get_num_connections() != 0 )
    {
      LOG( M_ERROR,
        "SimulationManager::set_status",
        "Cannot change waveform communication interval after connections have "
        "been created. Please call ResetKernel first." );
      throw KernelException();
    }
    else if ( wfr_interval < Time::get_resolution().get_ms() )
    {
      LOG( M_ERROR,
        "SimulationManager::set_status",
        "Communication interval of the waveform relaxation must be greater or "
        "equal to the resolution of the simulation." );
      throw KernelException();
    }
    else
    {
      LOG( M_INFO,
        "SimulationManager::set_status",
        "Waveform communication interval changed successfully. " );
      wfr_comm_interval_ = wfr_interval;
    }
  }

  // set the convergence tolerance for the waveform relaxation method
  double tol;
  if ( updateValue< double >( d, "wfr_tol", tol ) )
  {
    if ( tol < 0.0 )
      LOG( M_ERROR,
        "SimulationManager::set_status",
        "Tolerance must be zero or positive" );
    else
      wfr_tol_ = tol;
  }

  // set the maximal number of iterations for the waveform relaxation method
  long max_iter;
  if ( updateValue< long >( d, "wfr_max_iterations", max_iter ) )
  {
    if ( max_iter <= 0 )
      LOG( M_ERROR,
        "SimulationManager::set_status",
        "Maximal number of iterations  for the waveform relaxation must be "
        "positive. To disable waveform relaxation set use_wfr instead." );
    else
      wfr_max_iterations_ = max_iter;
  }

  // set the interpolation order for the waveform relaxation method
  long interp_order;
  if ( updateValue< long >( d, "wfr_interpolation_order", interp_order ) )
  {
    if ( ( interp_order < 0 ) || ( interp_order == 2 ) || ( interp_order > 3 ) )
      LOG( M_ERROR,
        "SimulationManager::set_status",
        "Interpolation order must be 0, 1, or 3." );
    else
      wfr_interpolation_order_ = interp_order;
  }
}

void
nest::SimulationManager::get_status( DictionaryDatum& d )
{
  def< double >( d, "ms_per_tic", Time::get_ms_per_tic() );
  def< double >( d, "tics_per_ms", Time::get_tics_per_ms() );
  def< long >( d, "tics_per_step", Time::get_tics_per_step() );
  def< double >( d, "resolution", Time::get_resolution().get_ms() );

  def< double >( d, "T_min", Time::min().get_ms() );
  def< double >( d, "T_max", Time::max().get_ms() );

  def< double >( d, "time", get_time().get_ms() );
  def< long >( d, "to_do", to_do_ );
  def< bool >( d, "print_time", print_time_ );

  def< bool >( d, "use_wfr", use_wfr_ );
  def< double >( d, "wfr_comm_interval", wfr_comm_interval_ );
  def< double >( d, "wfr_tol", wfr_tol_ );
  def< long >( d, "wfr_max_iterations", wfr_max_iterations_ );
  def< long >( d, "wfr_interpolation_order", wfr_interpolation_order_ );
}

void
nest::SimulationManager::simulate( Time const& t )
{
  assert( kernel().is_initialized() );

  t_real_ = 0;
  t_slice_begin_ = timeval();
  t_slice_end_ = timeval();

  if ( t == Time::ms( 0.0 ) )
    return;

  if ( t < Time::step( 1 ) )
  {
    LOG( M_ERROR,
      "SimulationManager::simulate",
      String::compose( "Simulation time must be >= %1 ms (one time step).",
           Time::get_resolution().get_ms() ) );
    throw KernelException();
  }

  if ( t.is_finite() )
  {
    Time time1 = clock_ + t;
    if ( !time1.is_finite() )
    {
      std::string msg = String::compose(
        "A clock overflow will occur after %1 of %2 ms. Please reset network "
        "clock first!",
        ( Time::max() - clock_ ).get_ms(),
        t.get_ms() );
      LOG( M_ERROR, "SimulationManager::simulate", msg );
      throw KernelException();
    }
  }
  else
  {
    std::string msg = String::compose(
      "The requested simulation time exceeds the largest time NEST can handle "
      "(T_max = %1 ms). Please use a shorter time!",
      Time::max().get_ms() );
    LOG( M_ERROR, "SimulationManager::simulate", msg );
    throw KernelException();
  }

  to_do_ += t.get_steps();
  to_do_total_ = to_do_;

  const size_t num_active_nodes = prepare_simulation_();

  // from_step_ is not touched here.  If we are at the beginning
  // of a simulation, it has been reset properly elsewhere.  If
  // a simulation was ended and is now continued, from_step_ will
  // have the proper value.  to_step_ is set as in advance_time().

  delay end_sim = from_step_ + to_do_;
  if ( kernel().connection_manager.get_min_delay() < end_sim )
    to_step_ =
      kernel()
        .connection_manager.get_min_delay(); // update to end of time slice
  else
    to_step_ = end_sim; // update to end of simulation time

  // Warn about possible inconsistencies, see #504.
  // This test cannot come any earlier, because we first need to compute
  // min_delay_
  // above.
  if ( t.get_steps() % kernel().connection_manager.get_min_delay() != 0 )
    LOG( M_WARNING,
      "SimulationManager::simulate",
      "The requested simulation time is not an integer multiple of the minimal "
      "delay in the network. This may result in inconsistent results under the "
      "following conditions: (i) A network contains more than one source of "
      "randomness, e.g., two different poisson_generators, and (ii) Simulate "
      "is called repeatedly with simulation times that are not multiples of "
      "the minimal delay." );

  resume_( num_active_nodes );

  finalize_simulation_();
}

void
nest::SimulationManager::resume_( size_t num_active_nodes )
{
  assert( kernel().is_initialized() );

  std::ostringstream os;
  double t_sim = to_do_ * Time::get_resolution().get_ms();

  os << "Number of local nodes: " << num_active_nodes << std::endl;
  os << "Simulaton time (ms): " << t_sim;

#ifdef _OPENMP
  os << std::endl
     << "Number of OpenMP threads: " << kernel().vp_manager.get_num_threads();
#else
  os << std::endl
     << "Not using OpenMP";
#endif

#ifdef HAVE_MPI
  os << std::endl
     << "Number of MPI processes: " << kernel().mpi_manager.get_num_processes();
#else
  os << std::endl
     << "Not using MPI";
#endif

  LOG( M_INFO, "SimulationManager::resume", os.str() );


  terminate_ = false;

  if ( to_do_ == 0 )
    return;

  if ( print_time_ )
  {
    // TODO: Remove direct output
    std::cout << std::endl;
    print_progress_();
  }

  simulating_ = true;
  simulated_ = true;

  update_();

  simulating_ = false;

  if ( print_time_ )
    std::cout << std::endl;

  kernel().mpi_manager.synchronize();

  if ( terminate_ )
  {
    LOG( M_ERROR,
      "SimulationManager::resume",
      "Exiting on error or user signal." );
    LOG( M_ERROR,
      "SimulationManager::resume",
      "SimulationManager: Use 'ResumeSimulation' to resume." );

    if ( SLIsignalflag != 0 )
    {
      SystemSignal signal( SLIsignalflag );
      SLIsignalflag = 0;
      throw signal;
    }
    else
      throw SimulationError();
  }

  LOG( M_INFO, "SimulationManager::resume", "Simulation finished." );
}

size_t
nest::SimulationManager::prepare_simulation_()
{
  assert( to_do_ != 0 ); // This is checked in simulate()

  // Reset profiling timers and counters within event_delivery_manager
  kernel().event_delivery_manager.reset_timers_counters();

  // find shortest and longest delay across all MPI processes
  // this call sets the member variables
  kernel().connection_manager.update_delay_extrema_();
  kernel().event_delivery_manager.init_moduli();

  // Check for synchronicity of global rngs over processes.
  // We need to do this ahead of any simulation in case random numbers
  // have been consumed on the SLI level.
  if ( kernel().mpi_manager.get_num_processes() > 1 )
  {
    if ( !kernel().mpi_manager.grng_synchrony(
           kernel().rng_manager.get_grng()->ulrand( 100000 ) ) )
    {
      LOG( M_ERROR,
        "SimulationManager::simulate",
        "Global Random Number Generators are not synchronized prior to "
        "simulation." );
      throw KernelException();
    }
  }

  // if at the beginning of a simulation, set up spike buffers
  if ( !simulated_ )
    kernel().event_delivery_manager.configure_spike_buffers();

  kernel().node_manager.ensure_valid_thread_local_ids();
  const size_t num_active_nodes = kernel().node_manager.prepare_nodes();

  kernel().model_manager.create_secondary_events_prototypes();

  // we have to do enter_runtime after prepre_nodes, since we use
  // calibrate to map the ports of MUSIC devices, which has to be done
  // before enter_runtime
  if ( !simulated_ ) // only enter the runtime mode once
  {
    double tick = Time::get_resolution().get_ms()
      * kernel().connection_manager.get_min_delay();
    kernel().music_manager.enter_runtime( tick );
  }

  return num_active_nodes;
}

bool
nest::SimulationManager::wfr_update_( Node* n )
{
  return ( n->wfr_update( clock_, from_step_, to_step_ ) );
}

void
nest::SimulationManager::update_()
{
  // to store done values of the different threads
  std::vector< bool > done;
  bool done_all = true;
  delay old_to_step;

  std::vector< lockPTR< WrappedThreadException > > exceptions_raised(
    kernel().vp_manager.get_num_threads() );
// parallel section begins
#pragma omp parallel
  {
    const int thrd = kernel().vp_manager.get_thread_id();

    do
    {
      if ( print_time_ )
        gettimeofday( &t_slice_begin_, NULL );

      if ( kernel().sp_manager.is_structural_plasticity_enabled()
        && ( clock_.get_steps() + from_step_ )
            % kernel().sp_manager.get_structural_plasticity_update_interval()
          == 0 )
      {
        for ( std::vector< Node* >::const_iterator i =
                kernel().node_manager.get_nodes_on_thread( thrd ).begin();
              i != kernel().node_manager.get_nodes_on_thread( thrd ).end();
              ++i )
        {
          ( *i )->update_synaptic_elements(
            Time( Time::step( clock_.get_steps() + from_step_ ) ).get_ms() );
        }
#pragma omp barrier
#pragma omp single
        {
          kernel().sp_manager.update_structural_plasticity();
        }
        // Remove 10% of the vacant elements
        for ( std::vector< Node* >::const_iterator i =
                kernel().node_manager.get_nodes_on_thread( thrd ).begin();
              i != kernel().node_manager.get_nodes_on_thread( thrd ).end();
              ++i )
        {
          ( *i )->decay_synaptic_elements_vacant();
        }
      }


      if ( from_step_ == 0 ) // deliver only at beginning of slice
      {
        kernel().event_delivery_manager.deliver_events( thrd );
#ifdef HAVE_MUSIC
// advance the time of music by one step (min_delay * h) must
// be done after deliver_events_() since it calls
// music_event_out_proxy::handle(), which hands the spikes over to
// MUSIC *before* MUSIC time is advanced

// wait until all threads are done -> synchronize
#pragma omp barrier
// the following block is executed by the master thread only
// the other threads are enforced to wait at the end of the block
#pragma omp master
        {
          // advance the time of music by one step (min_delay * h) must
          // be done after deliver_events_() since it calls
          // music_event_out_proxy::handle(), which hands the spikes over to
          // MUSIC *before* MUSIC time is advanced
          if ( slice_ > 0 )
            kernel().music_manager.advance_music_time();

          // the following could be made thread-safe
          kernel().music_manager.update_music_event_handlers(
            clock_, from_step_, to_step_ );
        }
// end of master section, all threads have to synchronize at this point
#pragma omp barrier
#endif
      }

      // preliminary update of nodes that use waveform relaxtion
      if ( kernel().node_manager.any_node_uses_wfr() )
      {
#pragma omp single
        {
          // if the end of the simulation is in the middle
          // of a min_delay_ step, we need to make a complete
          // step in the wfr_update and only do
          // the partial step in the final update
          // needs to be done in omp single since to_step_ is a scheduler
          // variable
          old_to_step = to_step_;
          if ( to_step_ < kernel().connection_manager.get_min_delay() )
            to_step_ = kernel().connection_manager.get_min_delay();
        }

        bool max_iterations_reached = true;
        const std::vector< Node* >& thread_local_wfr_nodes =
          kernel().node_manager.get_wfr_nodes_on_thread( thrd );
        for ( long n = 0; n < wfr_max_iterations_; ++n )
        {
          bool done_p = true;

          // this loop may be empty for those threads
          // that do not have any nodes requiring wfr_update
          for ( std::vector< Node* >::const_iterator i =
                  thread_local_wfr_nodes.begin();
                i != thread_local_wfr_nodes.end();
                ++i )
            done_p = wfr_update_( *i ) && done_p;

// add done value of thread p to done vector
#pragma omp critical
          done.push_back( done_p );
// parallel section ends, wait until all threads are done -> synchronize
#pragma omp barrier

// the following block is executed by a single thread
// the other threads wait at the end of the block
#pragma omp single
          {
            // set done_all
            for ( size_t i = 0; i < done.size(); i++ )
              done_all = done[ i ] && done_all;

            // gather SecondaryEvents (e.g. GapJunctionEvents)
            kernel().event_delivery_manager.gather_events( done_all );

            // reset done and done_all
            //(needs to be in the single threaded part)
            done_all = true;
            done.clear();
          }

          // deliver SecondaryEvents generated during wfr_update
          // returns the done value over all threads
          done_p = kernel().event_delivery_manager.deliver_events( thrd );

          if ( done_p )
          {
            max_iterations_reached = false;
            break;
          }
        } // of for (wfr_max_iterations) ...

#pragma omp single
        {
          to_step_ = old_to_step;
          if ( max_iterations_reached )
          {
            std::string msg = String::compose(
              "Maximum number of iterations reached at interval %1-%2 ms",
              clock_.get_ms(),
              clock_.get_ms() + to_step_ * Time::get_resolution().get_ms() );
            LOG( M_WARNING, "SimulationManager::wfr_update", msg );
          }
        }

      } // of if(any_node_uses_wfr)
      // end of preliminary update

      const std::vector< Node* >& thread_local_nodes =
        kernel().node_manager.get_nodes_on_thread( thrd );
      for (
        std::vector< Node* >::const_iterator node = thread_local_nodes.begin();
        node != thread_local_nodes.end();
        ++node )
      {
        // We update in a parallel region. Therefore, we need to catch
        // exceptions here and then handle them after the parallel region.
        try
        {
          if ( not( *node )->is_frozen() )
            ( *node )->update( clock_, from_step_, to_step_ );
        }
        catch ( std::exception& e )
        {
          // so throw the exception after parallel region
          exceptions_raised.at( thrd ) = lockPTR< WrappedThreadException >(
            new WrappedThreadException( e ) );
          terminate_ = true;
        }
      }

// parallel section ends, wait until all threads are done -> synchronize
#pragma omp barrier

// the following block is executed by the master thread only
// the other threads are enforced to wait at the end of the block
#pragma omp master
      {
        // gather only at end of slice
        if ( to_step_ == kernel().connection_manager.get_min_delay() )
          kernel().event_delivery_manager.gather_events( true );

        advance_time_();

        if ( SLIsignalflag != 0 )
        {
          LOG( M_INFO,
            "SimulationManager::update",
            "Simulation exiting on user signal." );
          terminate_ = true;
        }

        if ( print_time_ )
        {
          gettimeofday( &t_slice_end_, NULL );
          print_progress_();
        }
      }
// end of master section, all threads have to synchronize at this point
#pragma omp barrier

    } while ( ( to_do_ != 0 ) && ( !terminate_ ) );

    // End of the slice, we update the number of synaptic element
    for ( std::vector< Node* >::const_iterator i =
            kernel().node_manager.get_nodes_on_thread( thrd ).begin();
          i != kernel().node_manager.get_nodes_on_thread( thrd ).end();
          ++i )
    {
      ( *i )->update_synaptic_elements(
        Time( Time::step( clock_.get_steps() + to_step_ ) ).get_ms() );
    }

  } // end of #pragma parallel omp

  // check if any exceptions have been raised
  for ( index thrd = 0; thrd < kernel().vp_manager.get_num_threads(); ++thrd )
    if ( exceptions_raised.at( thrd ).valid() )
      throw WrappedThreadException( *( exceptions_raised.at( thrd ) ) );
}

void
nest::SimulationManager::finalize_simulation_()
{
  if ( not simulated_ )
    return;

  // Check for synchronicity of global rngs over processes
  // TODO: This seems double up, there is such a test at end of simulate()
  if ( kernel().mpi_manager.get_num_processes() > 1 )
    if ( !kernel().mpi_manager.grng_synchrony(
           kernel().rng_manager.get_grng()->ulrand( 100000 ) ) )
    {
      LOG( M_ERROR,
        "SimulationManager::simulate",
        "Global Random Number Generators are not synchronized after "
        "simulation." );
      throw KernelException();
    }

  kernel().node_manager.finalize_nodes();
}

void
nest::SimulationManager::reset_network()
{
  if ( not has_been_simulated() )
    return; // nothing to do


  kernel().event_delivery_manager.clear_pending_spikes();

  kernel().node_manager.reset_nodes_state();

  // ConnectionManager doesn't support resetting dynamic synapses yet
  LOG( M_WARNING,
    "SimulationManager::ResetNetwork",
    "Synapses with internal dynamics (facilitation, STDP) are not reset.\n"
    "This will be implemented in a future version of NEST." );
}

void
nest::SimulationManager::advance_time_()
{
  // time now advanced time by the duration of the previous step
  to_do_ -= to_step_ - from_step_;

  // advance clock, update modulos, slice counter only if slice completed
  if ( ( delay ) to_step_ == kernel().connection_manager.get_min_delay() )
  {
    clock_ += Time::step( kernel().connection_manager.get_min_delay() );
    ++slice_;
    kernel().event_delivery_manager.update_moduli();
    from_step_ = 0;
  }
  else
    from_step_ = to_step_;

  long end_sim = from_step_ + to_do_;

  if ( kernel().connection_manager.get_min_delay() < ( delay ) end_sim )
    // update to end of time slice
    to_step_ = kernel().connection_manager.get_min_delay();
  else
    to_step_ = end_sim; // update to end of simulation time

  assert( to_step_ - from_step_
    <= ( long ) kernel().connection_manager.get_min_delay() );
}

void
nest::SimulationManager::print_progress_()
{
  double rt_factor = 0.0;

  if ( t_slice_end_.tv_sec != 0 )
  {
    // usec
    long t_real_s = ( t_slice_end_.tv_sec - t_slice_begin_.tv_sec ) * 1e6;
    // usec
    t_real_ += t_real_s + ( t_slice_end_.tv_usec - t_slice_begin_.tv_usec );
    // ms
    double t_real_acc = ( t_real_ ) / 1000.;
    double t_sim_acc =
      ( to_do_total_ - to_do_ ) * Time::get_resolution().get_ms();
    rt_factor = t_sim_acc / t_real_acc;
  }

  int percentage = ( 100 - int( float( to_do_ ) / to_do_total_ * 100 ) );

  std::cout << "\r" << std::setw( 3 ) << std::right << percentage << " %: "
            << "network time: " << std::fixed << std::setprecision( 1 )
            << clock_.get_ms() << " ms, "
            << "realtime factor: " << std::setprecision( 4 ) << rt_factor
            << std::resetiosflags( std::ios_base::floatfield );
  std::flush( std::cout );
}

nest::Time const
nest::SimulationManager::get_previous_slice_origin() const
{
  return clock_ - Time::step( kernel().connection_manager.get_min_delay() );
}
