#include <aspect/global.h>
#include "composition_time_step.h"

namespace aspect
{
  namespace TimeStepping
  {
    template <int dim>
    double
    CompositionTimeStep<dim>::execute()
    {
      double min_local_composition_timestep = std::numeric_limits<double>::max();

      const QIterated<dim> quadrature_formula (QTrapezoid<1>(),
                                               this->get_parameters().stokes_velocity_degree);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values |
                               update_gradients |
                               update_quadrature_points);

      const unsigned int n_q_points = quadrature_formula.size();

      double current_time_step;
      if (this->get_timestep_number() > 0) current_time_step = this->get_timestep();
      else current_time_step = 1e30; // or return ?

      MaterialModel::MaterialModelInputs<dim> in(n_q_points,
                                                 this->introspection().n_compositional_fields);

      std::vector<double> reactions(n_q_points, numbers::signaling_nan<double>()); // Why?
      const unsigned int porosity_index = this->introspection().compositional_index_for_name("porosity");

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            in.reinit(fe_values,
                      cell,
                      this->introspection(),
                      this->get_solution());
            
            // TODO iterate over all compositions or choose some specific?

            // reactions vector contains accummulated reactions over time step (integral, not rate; see helper_functions.cc)
            fe_values[this->introspection().extractors.compositional_fields[porosity_index]].get_function_values(this->get_reaction_vector(),
            reactions);

            for (unsigned int q=0; q<n_q_points; ++q)
              {                
                // convert accummulated porosity change to rate:
                double reaction_rate = std::max(std::abs(reactions[q]),1e-30)/current_time_step;
                //std::cout << reaction_rate << " ";
                // if necessary, reduce time step to get maximum porosity change:
                min_local_composition_timestep = std::min(min_local_composition_timestep,
                                                             max_change/reaction_rate);
                //min_local_composition_timestep = 1e2*year_in_seconds; // TEST
              }
          }

      const double min_composition_timestep = Utilities::MPI::min (min_local_composition_timestep, this->get_mpi_communicator());

      AssertThrow (min_composition_timestep > 0,
                   ExcMessage("The time step length for the each time step needs to be positive, "
                              "but the computed step length was: " + std::to_string(min_composition_timestep) + ". "
                              "Please check."));

      std::cout << "Composition time step calculated:" << min_composition_timestep/year_in_seconds << " yrs. \n";
      return min_composition_timestep;
    }

    template <int dim>
    void
    CompositionTimeStep<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Time stepping");
      {
        prm.enter_subsection("Composition time step");

        prm.declare_entry("Max porosity change", "0.01",
                          Patterns::Double (0.),
                          "Maximum absolute change of porosity over a time step.");

        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    CompositionTimeStep<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Time stepping");
      {
        prm.enter_subsection("Composition time step");

        max_change = prm.get_double("Max porosity change");
        
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace TimeStepping
  {
    ASPECT_REGISTER_TIME_STEPPING_MODEL(CompositionTimeStep,
                                        "composition time step",
                                        "TO CHANGE This model computes the composition time step "
                                        "from reaction rates. Requires operator splitting")
  }
}