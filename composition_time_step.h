#ifndef _aspect_time_stepping_composition_time_step_h
#define _aspect_time_stepping_composition_time_step_h

#include <aspect/time_stepping/interface.h>


namespace aspect
{
  namespace TimeStepping
  {
    using namespace dealii;

    /**
     * Compute the maximum time step based on the current solution
     * of the compositional field and
     * return this as the time step.
     *
     * @ingroup TimeStepping
     */
    template <int dim>
    class CompositionTimeStep : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        CompositionTimeStep () = default;


        /**
         * @copydoc aspect::TimeStepping::Interface<dim>::execute()
         */
        double
        execute() override;
      /**
       * Declare the parameters this class takes through input files.
       */
      static void
      declare_parameters(ParameterHandler &prm);

      /**
       * Read the parameters this class declares from the parameter file.
       */
      void
      parse_parameters(ParameterHandler &prm) override;
      /**
       * @}
       */
      
      private:
        double max_change;

    };
  }
}

#endif